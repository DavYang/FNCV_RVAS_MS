#!/usr/bin/env python3
"""Sample common SNPs from the AoU ACAF splitMT for background variant set.

Reads the pre-built ACAF threshold split MatrixTable directly, filters to
EUR ancestry samples, and samples a target number of variants proportionally
across autosomes (chr1-22). Each chromosome is processed independently with
partition-pruned reads and written as a separate Hail MT to GCS.

Usage:
    python 00_sample_common_SNPs_mt.py
"""

import math
import os
import sys
import json
import time
import hail as hl
import hailtop.fs as hfs
from datetime import datetime
from utils import load_config, setup_logger

logger = setup_logger("SNP_sampling_mt")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
# GRCh38 autosome sizes
CHR_SIZES = {
    'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559,
    'chr4': 190214555, 'chr5': 181538259, 'chr6': 170805979,
    'chr7': 159345973, 'chr8': 145138636, 'chr9': 138394717,
    'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
    'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189,
    'chr16': 90338345, 'chr17': 83257441, 'chr18': 80373285,
    'chr19': 58617616, 'chr20': 64444167, 'chr21': 46709983,
    'chr22': 50818468,
}

# Estimated ACAF variant density: ~99M variants / 2.88 Gbp autosomal genome
# Calibrated from chr21 test: 1,787,412 variants / 46,709,983 bp = 0.03826
ACAF_VARIANTS_PER_BP = 0.03826

# Overshoot buffer for Bernoulli sampling (sample slightly more, trim after)
SAMPLING_OVERSHOOT = 1.08

# Maximum chunk size in base pairs (~50 Mbp, similar to chr21 which succeeded)
CHUNK_SIZE_BP = 50_000_000


def _cleanup_gcs_path(path: str) -> None:
    """Recursively remove a GCS path, logging but not raising on failure."""
    try:
        hfs.rmtree(path)
        logger.info(f"Cleaned up: {path}")
    except Exception as e:
        logger.warning(f"Failed to clean up {path}: {e}")

# ---------------------------------------------------------------------------
# Global caches
# ---------------------------------------------------------------------------
_eur_ids_cache = None
_eur_set_cache = None


def _load_eur_sample_ids(config: dict) -> list:
    """Load EUR ancestry sample IDs from the ancestry TSV (cached)."""
    global _eur_ids_cache
    if _eur_ids_cache is None:
        logger.info("Loading EUR ancestry sample IDs ...")
        ancestry_ht = hl.import_table(
            config['inputs']['ancestry_pred'],
            impute=True,
            types={'research_id': hl.tstr, 'ancestry_pred': hl.tstr},
        )
        eur_ht = ancestry_ht.filter(ancestry_ht.ancestry_pred == 'eur')
        _eur_ids_cache = eur_ht.aggregate(hl.agg.collect(eur_ht.research_id))
        logger.info(f"Cached {len(_eur_ids_cache)} EUR sample IDs")
    return _eur_ids_cache


def _get_eur_set(config: dict) -> hl.expr.SetExpression:
    """Return a Hail literal Set of EUR sample IDs (cached)."""
    global _eur_set_cache
    if _eur_set_cache is None:
        ids = _load_eur_sample_ids(config)
        _eur_set_cache = hl.literal(set(ids))
        logger.info("Created Hail Set for EUR filtering")
    return _eur_set_cache


# ---------------------------------------------------------------------------
# Chromosome targets
# ---------------------------------------------------------------------------
def calculate_chromosome_targets(target_snps: int) -> dict:
    """Distribute target SNPs proportionally across autosomes (chr1-22).

    Args:
        target_snps: Total number of SNPs to sample genome-wide.

    Returns:
        dict mapping chromosome name to per-chromosome target count.
    """
    total_size = sum(CHR_SIZES.values())
    chr_targets = {}
    for chrom, size in CHR_SIZES.items():
        chr_targets[chrom] = max(1, int(target_snps * (size / total_size)))

    # Adjust rounding so the total matches target_snps exactly
    current_total = sum(chr_targets.values())
    remaining = target_snps - current_total
    if remaining != 0:
        sorted_chroms = sorted(chr_targets, key=lambda c: chr_targets[c], reverse=True)
        for i, chrom in enumerate(sorted_chroms):
            if i >= abs(remaining):
                break
            chr_targets[chrom] += 1 if remaining > 0 else -1

    logger.info(f"Per-chromosome targets (autosomes): {chr_targets}")
    logger.info(f"Total target: {sum(chr_targets.values()):,}")
    return chr_targets


# ---------------------------------------------------------------------------
# Per-chromosome processing
# ---------------------------------------------------------------------------
def _make_chunk_intervals(chrom: str) -> list:
    """Split a chromosome into sub-intervals of CHUNK_SIZE_BP.

    Args:
        chrom: Chromosome name, e.g. 'chr1'.

    Returns:
        List of (start, end) tuples in 1-based coordinates.
    """
    chr_size = CHR_SIZES[chrom]
    n_chunks = math.ceil(chr_size / CHUNK_SIZE_BP)
    intervals = []
    for i in range(n_chunks):
        start = i * CHUNK_SIZE_BP + 1
        end = min((i + 1) * CHUNK_SIZE_BP, chr_size)
        intervals.append((start, end))
    return intervals


def process_chromosome(
    chrom: str,
    chr_target: int,
    mt_path: str,
    output_dir: str,
    config: dict,
) -> dict:
    """Sample variants for a single chromosome from the ACAF splitMT.

    Large chromosomes are split into ~50 Mbp chunks (chr21-sized) to avoid
    OOM. Each chunk is processed independently as a small Spark job, then
    the tiny sampled chunks are combined.

    Pipeline per chunk:
      1. read_matrix_table  (lazy)
      2. filter_intervals   (partition pruning to chunk)
      3. filter_cols EUR    (lazy predicate)
      4. select_entries GT  (drop unused fields)
      5. sample_rows        (lazy, fraction from estimated density)
      6. write chunk MT     (small Spark job)

    After all chunks:
      7. union_rows tiny chunk MTs -> write final MT
      8. trim if overshot

    Args:
        chrom: Chromosome name, e.g. 'chr21'.
        chr_target: Number of SNPs to sample for this chromosome.
        mt_path: GCS path to the ACAF splitMT.
        output_dir: GCS base output directory for this run.
        config: Loaded config dict.

    Returns:
        Summary dict for this chromosome.
    """
    seed = config['sampling']['random_seed']
    eur_set = _get_eur_set(config)

    logger.info(f"\n{'='*50}")
    logger.info(f"Processing {chrom} (target: {chr_target:,} SNPs)")
    logger.info(f"{'='*50}")
    chr_start = time.time()

    # Compute sampling fraction from estimated ACAF density
    chr_size = CHR_SIZES[chrom]
    estimated_pool = int(chr_size * ACAF_VARIANTS_PER_BP)
    fraction = min(1.0, (chr_target / estimated_pool) * SAMPLING_OVERSHOOT)

    # Split chromosome into manageable chunks
    chunks = _make_chunk_intervals(chrom)
    logger.info(f"{chrom}: {len(chunks)} chunks, estimated pool ~{estimated_pool:,}, "
                f"sampling fraction {fraction:.6f}")

    chr_output_dir = f"{output_dir}/{chrom}"
    chunks_dir = f"{chr_output_dir}/_chunks"
    chunk_paths = []

    try:
        for ci, (start, end) in enumerate(chunks):
            chunk_label = f"{chrom}_chunk{ci}"
            chunk_mt_path = f"{chunks_dir}/{chunk_label}.mt"

            logger.info(f"  {chunk_label}: {start:,}-{end:,} "
                        f"(~{(end - start) / 1e6:.0f} Mbp)")
            chunk_start = time.time()

            # Lazy read + partition-pruned interval filter for this chunk
            mt = hl.read_matrix_table(mt_path)
            interval = hl.parse_locus_interval(
                f"{chrom}:{start}-{end}", reference_genome='GRCh38'
            )
            mt = hl.filter_intervals(mt, [interval])

            # Filter to EUR samples, drop unused row/entry fields
            mt = mt.filter_cols(eur_set.contains(mt.s))
            mt = mt.select_entries('GT')
            mt = mt.select_rows()

            # Sample (lazy) and write this chunk
            if fraction < 1.0:
                mt_chunk = mt.sample_rows(fraction, seed=seed + ci)
            else:
                mt_chunk = mt

            mt_chunk.write(chunk_mt_path, overwrite=True)
            chunk_time = time.time() - chunk_start
            logger.info(f"  {chunk_label}: wrote in {chunk_time:.1f}s")
            chunk_paths.append(chunk_mt_path)

    except Exception as e:
        logger.error(f"{chrom}: Chunk processing failed at chunk {ci}: {e}")
        _cleanup_gcs_path(chunks_dir)
        raise

    # Combine all tiny chunk MTs into the final chromosome MT
    logger.info(f"{chrom}: Combining {len(chunk_paths)} chunks ...")
    combine_start = time.time()

    mt_combined = hl.read_matrix_table(chunk_paths[0])
    for cp in chunk_paths[1:]:
        mt_combined = mt_combined.union_rows(hl.read_matrix_table(cp))

    mt_output_path = f"{chr_output_dir}/{chrom}_background_snps.mt"
    mt_combined.write(mt_output_path, overwrite=True)
    combine_time = time.time() - combine_start
    logger.info(f"{chrom}: Combined in {combine_time:.1f}s")

    # Clean up temporary chunk MTs now that final MT is written
    _cleanup_gcs_path(chunks_dir)

    # Count the final MT (cheap -- tiny)
    mt_final = hl.read_matrix_table(mt_output_path)
    final_count = mt_final.count_rows()
    logger.info(f"{chrom}: {final_count:,} SNPs (target: {chr_target:,})")

    # Trim if overshot
    if final_count > chr_target:
        logger.info(f"{chrom}: Trimming {final_count:,} -> {chr_target:,}")
        mt_trimmed = mt_final.head(chr_target)
        mt_trimmed.write(mt_output_path, overwrite=True)
        final_count = chr_target

    chr_time = time.time() - chr_start

    # Write JSON summary
    summary = {
        'chromosome': chrom,
        'target_snps': chr_target,
        'estimated_pool': estimated_pool,
        'sampling_fraction': round(fraction, 6),
        'n_chunks': len(chunks),
        'final_sampled_variants': final_count,
        'output_path': mt_output_path,
        'combine_time_seconds': round(combine_time, 1),
        'total_time_seconds': round(chr_time, 1),
        'timestamp': datetime.now().isoformat(),
    }
    summary_path = f"{chr_output_dir}/{chrom}_sampling_summary.json"
    with hl.hadoop_open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)

    logger.info(f"{chrom}: Done -- {final_count:,} SNPs in {chr_time:.1f}s")
    return summary


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    config = load_config("config/config.json")

    # Initialize Hail
    hl.init(
        log='/tmp/hail_snp_sampling_mt.log',
        spark_conf={
            'spark.driver.memory': '12g',
            'spark.executor.memory': '8g',
            'spark.network.timeout': '600s',
            'spark.executor.heartbeatInterval': '120s',
            'spark.driver.maxResultSize': '4g',
        },
    )
    hl.default_reference('GRCh38')

    # Resolve workspace bucket
    workspace_bucket = os.environ.get('WORKSPACE_BUCKET')
    if not workspace_bucket:
        logger.error("WORKSPACE_BUCKET environment variable not set")
        sys.exit(1)
    if not workspace_bucket.startswith('gs://'):
        workspace_bucket = f"gs://{workspace_bucket}"

    logger.info(f"Workspace bucket: {workspace_bucket}")

    # splitMT path
    mt_path = config['inputs']['wgs_matrix_table']
    logger.info(f"Source MT: {mt_path}")

    # Output directory
    target_snps = config['sampling']['target_total_snps']
    if target_snps >= 1_000_000:
        snp_label = f"{target_snps // 1_000_000}M"
    elif target_snps >= 1_000:
        snp_label = f"{target_snps // 1_000}K"
    else:
        snp_label = str(target_snps)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    base_filename = f"{snp_label}_background_snps"
    output_dir = (f"{workspace_bucket}/results/FNCV_RVAS_MS/"
                  f"{base_filename}_{timestamp}")
    logger.info(f"Output directory: {output_dir}")

    # Test mode vs production
    test_mode = config['params'].get('test_mode', False)
    test_chromosome = config['params'].get('test_chromosome', None)

    chr_targets = calculate_chromosome_targets(target_snps)

    try:
        if test_mode and test_chromosome:
            # ----- Test mode: single chromosome -----
            logger.info(f"=== TEST MODE: {test_chromosome} only ===")
            if test_chromosome not in chr_targets:
                logger.error(f"{test_chromosome} not in autosome targets")
                sys.exit(1)

            chr_target = chr_targets[test_chromosome]
            logger.info(f"Target for {test_chromosome}: {chr_target:,}")

            summary = process_chromosome(
                test_chromosome, chr_target, mt_path, output_dir, config,
            )

            if summary:
                test_summary = {
                    'test_mode': True,
                    'test_chromosome': test_chromosome,
                    'target_snps': chr_target,
                    'final_sampled_variants': summary['final_sampled_variants'],
                    'output_directory': output_dir,
                    'timestamp': datetime.now().isoformat(),
                }
                path = f"{output_dir}/{base_filename}_test_summary.json"
                with hl.hadoop_open(path, 'w') as f:
                    json.dump(test_summary, f, indent=2)

                logger.info(f"=== TEST COMPLETE: {test_chromosome} ===")
                logger.info(f"Sampled {summary['final_sampled_variants']:,} SNPs")
                logger.info(f"Results: {output_dir}")
            else:
                logger.error(f"TEST FAILED: no variants for {test_chromosome}")

        else:
            # ----- Production mode: all autosomes -----
            logger.info("=== PRODUCTION MODE: chr1-22 ===")
            all_summaries = []
            total_sampled = 0
            failed_chroms = []

            # Load progress for resume support
            progress_path = f"{output_dir}/_progress.json"
            completed_chroms = set()
            try:
                with hl.hadoop_open(progress_path, 'r') as f:
                    progress = json.load(f)
                    completed_chroms = set(progress.get('completed', []))
                if completed_chroms:
                    logger.info(f"Resuming: {len(completed_chroms)} chromosomes "
                                f"already completed: {sorted(completed_chroms)}")
            except Exception:
                pass

            for chrom, chr_target in chr_targets.items():
                if chrom in completed_chroms:
                    logger.info(f"{chrom}: Already completed, skipping")
                    continue

                try:
                    summary = process_chromosome(
                        chrom, chr_target, mt_path, output_dir, config,
                    )
                    if summary:
                        all_summaries.append(summary)
                        total_sampled += summary['final_sampled_variants']
                        completed_chroms.add(chrom)

                        # Update progress file after each chromosome
                        progress_data = {
                            'completed': sorted(completed_chroms),
                            'last_updated': datetime.now().isoformat(),
                        }
                        with hl.hadoop_open(progress_path, 'w') as f:
                            json.dump(progress_data, f, indent=2)

                except Exception as e:
                    logger.error(f"{chrom}: FAILED - {e}")
                    failed_chroms.append(chrom)
                    continue

            overall = {
                'target_snps': target_snps,
                'chromosome_targets': chr_targets,
                'total_sampled_variants': total_sampled,
                'chromosomes_processed': len(all_summaries),
                'chromosomes_failed': failed_chroms,
                'output_directory': output_dir,
                'timestamp': datetime.now().isoformat(),
            }
            path = f"{output_dir}/{base_filename}_overall_summary.json"
            with hl.hadoop_open(path, 'w') as f:
                json.dump(overall, f, indent=2)

            logger.info("=== PRODUCTION COMPLETE ===")
            logger.info(f"Chromosomes completed: {len(all_summaries)}")
            if failed_chroms:
                logger.warning(f"Chromosomes failed: {failed_chroms}")
            logger.info(f"Total SNPs: {total_sampled:,}")
            logger.info(f"Results: {output_dir}")

    except Exception as e:
        logger.error(f"Pipeline failed: {e}")
        raise


if __name__ == "__main__":
    main()
