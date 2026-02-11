#!/usr/bin/env python3
"""Build background SNP PLINK files for Regenie null model (Phase 1).

Two-pass pipeline to sample ~500K common SNPs from the AoU ACAF splitMT
and export per-chromosome PLINK files for EUR ancestry samples.

Pass 1 (cheap): Sample loci from the rows Table (no entry data).
Pass 2 (targeted): Extract genotypes at sampled loci, export PLINK.

Usage:
    python 01_build_background_snps_plink.py
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

logger = setup_logger("background_snps_plink")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
AUTOSOMES = [f"chr{i}" for i in range(1, 23)]

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

# Estimated ACAF variant density (calibrated from chr21 test run)
ACAF_VARIANTS_PER_BP = 0.03826

# Overshoot buffer for Bernoulli sampling (trim to exact target after)
SAMPLING_OVERSHOOT = 1.10

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


def _cleanup_gcs_path(path: str) -> None:
    """Recursively remove a GCS path, logging but not raising on failure."""
    try:
        hfs.rmtree(path)
        logger.info(f"Cleaned up: {path}")
    except Exception as e:
        logger.warning(f"Failed to clean up {path}: {e}")


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
    for chrom in AUTOSOMES:
        size = CHR_SIZES[chrom]
        chr_targets[chrom] = max(1, int(target_snps * (size / total_size)))

    current_total = sum(chr_targets.values())
    remaining = target_snps - current_total
    if remaining != 0:
        sorted_chroms = sorted(chr_targets, key=lambda c: chr_targets[c], reverse=True)
        for i, chrom in enumerate(sorted_chroms):
            if i >= abs(remaining):
                break
            chr_targets[chrom] += 1 if remaining > 0 else -1

    logger.info(f"Total target: {sum(chr_targets.values()):,} across {len(chr_targets)} autosomes")
    return chr_targets


# ---------------------------------------------------------------------------
# Pass 1: Sample loci from rows Table (no entry data)
# ---------------------------------------------------------------------------
def pass1_sample_loci(
    mt_path: str,
    output_ht_path: str,
    chr_targets: dict,
    config: dict,
) -> dict:
    """Sample loci proportionally across autosomes using rows Table only.

    This pass never touches entry data. It reads the row keys of the ACAF
    splitMT, applies Bernoulli sampling per chromosome, and writes a small
    Hail Table of selected (locus, alleles) pairs.

    Args:
        mt_path: GCS path to ACAF splitMT.
        output_ht_path: GCS path to write the sampled loci Table.
        chr_targets: Dict of chromosome -> target SNP count.
        config: Loaded config dict.

    Returns:
        Summary dict with per-chromosome counts.
    """
    seed = config['sampling']['random_seed']
    avoid_hla = config['sampling'].get('avoid_hla', False)
    hla_region = config['sampling'].get('hla_region', {})

    logger.info("=" * 60)
    logger.info("PASS 1: Sampling loci from rows Table (no entry data)")
    logger.info("=" * 60)
    pass1_start = time.time()

    # Read only the rows table -- zero entry data
    logger.info(f"Reading rows Table from {mt_path} ...")
    mt = hl.read_matrix_table(mt_path)
    ht_rows = mt.rows().select()
    logger.info("Rows Table loaded (no entry data accessed)")

    per_chr_tables = []
    per_chr_counts = {}

    for chrom in AUTOSOMES:
        chr_target = chr_targets.get(chrom, 0)
        if chr_target == 0:
            continue

        chr_start = time.time()

        # Filter to this chromosome using interval pruning
        interval = hl.parse_locus_interval(
            f"{chrom}:1-{CHR_SIZES[chrom]}", reference_genome='GRCh38'
        )
        ht_chr = hl.filter_intervals(ht_rows, [interval])

        # Exclude HLA region on chr6 if configured
        if avoid_hla and chrom == 'chr6' and hla_region:
            hla_chrom = f"chr{hla_region['chrom']}"
            hla_start = hla_region['start']
            hla_end = hla_region['end']
            hla_interval = hl.parse_locus_interval(
                f"{hla_chrom}:{hla_start}-{hla_end}", reference_genome='GRCh38'
            )
            ht_chr = ht_chr.filter(
                ~hl.is_defined(hla_interval.contains(ht_chr.locus))
                | ~hla_interval.contains(ht_chr.locus)
            )
            logger.info(f"  {chrom}: Excluded HLA region "
                        f"{hla_chrom}:{hla_start}-{hla_end}")

        # Compute sampling fraction from estimated density
        estimated_pool = int(CHR_SIZES[chrom] * ACAF_VARIANTS_PER_BP)
        fraction = min(1.0, (chr_target / estimated_pool) * SAMPLING_OVERSHOOT)

        # Bernoulli sampling using rand_unif on the rows Table
        chr_seed = seed + AUTOSOMES.index(chrom)
        ht_sampled = ht_chr.filter(hl.rand_unif(0.0, 1.0, seed=chr_seed) < fraction)

        per_chr_tables.append(ht_sampled)
        chr_time = time.time() - chr_start
        logger.info(f"  {chrom}: target={chr_target:,}, "
                    f"est_pool={estimated_pool:,}, "
                    f"fraction={fraction:.6f} ({chr_time:.1f}s)")

    # Union all per-chromosome sampled Tables
    logger.info("Combining per-chromosome sampled loci ...")
    ht_combined = per_chr_tables[0]
    for ht in per_chr_tables[1:]:
        ht_combined = ht_combined.union(ht)

    # Write the combined sampled loci Table
    logger.info(f"Writing sampled loci Table to {output_ht_path} ...")
    ht_combined.write(output_ht_path, overwrite=True)

    # Read back and count per chromosome
    ht_final = hl.read_table(output_ht_path)
    total_count = ht_final.count()
    logger.info(f"Pass 1 complete: {total_count:,} loci sampled")

    # Get per-chromosome breakdown
    chr_counts = ht_final.group_by(
        contig=ht_final.locus.contig
    ).aggregate(n=hl.agg.count())
    chr_counts_dict = {row.contig: row.n for row in chr_counts.collect()}

    for chrom in AUTOSOMES:
        n = chr_counts_dict.get(chrom, 0)
        target = chr_targets.get(chrom, 0)
        per_chr_counts[chrom] = {'sampled': n, 'target': target}
        logger.info(f"  {chrom}: {n:,} sampled (target: {target:,})")

    # Trim if total overshoots global target
    global_target = sum(chr_targets.values())
    if total_count > global_target:
        logger.info(f"Trimming {total_count:,} -> {global_target:,} loci")
        ht_trimmed = ht_final.head(global_target)
        ht_trimmed.write(output_ht_path, overwrite=True)
        total_count = global_target

    pass1_time = time.time() - pass1_start
    summary = {
        'pass': 1,
        'total_loci_sampled': total_count,
        'global_target': global_target,
        'per_chromosome': per_chr_counts,
        'output_ht_path': output_ht_path,
        'time_seconds': round(pass1_time, 1),
        'timestamp': datetime.now().isoformat(),
    }

    logger.info(f"Pass 1 done: {total_count:,} loci in {pass1_time:.1f}s")
    return summary


# ---------------------------------------------------------------------------
# Pass 2: Extract genotypes at sampled loci, export PLINK per chromosome
# ---------------------------------------------------------------------------
def pass2_export_plink(
    mt_path: str,
    loci_ht_path: str,
    output_dir: str,
    chr_targets: dict,
    config: dict,
) -> dict:
    """Extract EUR genotypes at sampled loci and export per-chromosome PLINK.

    For each chromosome:
      1. filter_intervals for partition pruning
      2. semi_join_rows against sampled loci (tiny Table)
      3. filter_cols to EUR samples
      4. select_entries('GT'), select_rows(), select_cols()
      5. export_plink with FID=IID=research_id

    Args:
        mt_path: GCS path to ACAF splitMT.
        loci_ht_path: GCS path to sampled loci Hail Table from Pass 1.
        output_dir: GCS base directory for PLINK output.
        chr_targets: Dict of chromosome -> target SNP count.
        config: Loaded config dict.

    Returns:
        Summary dict with per-chromosome variant/sample counts.
    """
    eur_set = _get_eur_set(config)

    logger.info("=" * 60)
    logger.info("PASS 2: Extracting genotypes and exporting PLINK")
    logger.info("=" * 60)
    pass2_start = time.time()

    # Load the sampled loci Table
    logger.info(f"Reading sampled loci from {loci_ht_path} ...")
    ht_loci = hl.read_table(loci_ht_path)

    plink_dir = f"{output_dir}/null_model_data"
    per_chr_results = {}
    failed_chroms = []

    # Load progress for resume support
    progress_path = f"{output_dir}/_pass2_progress.json"
    completed_chroms = set()
    try:
        with hl.hadoop_open(progress_path, 'r') as f:
            progress = json.load(f)
            completed_chroms = set(progress.get('completed', []))
        if completed_chroms:
            logger.info(f"Resuming Pass 2: {len(completed_chroms)} chromosomes "
                        f"already completed: {sorted(completed_chroms)}")
    except Exception:
        pass

    for chrom in AUTOSOMES:
        if chrom not in chr_targets:
            continue
        if chrom in completed_chroms:
            logger.info(f"  {chrom}: Already completed, skipping")
            continue

        try:
            chr_start = time.time()
            logger.info(f"\n  Processing {chrom} ...")

            # Filter sampled loci to this chromosome
            ht_chr_loci = ht_loci.filter(ht_loci.locus.contig == chrom)

            # Read MT with partition pruning for this chromosome
            mt = hl.read_matrix_table(mt_path)
            interval = hl.parse_locus_interval(
                f"{chrom}:1-{CHR_SIZES[chrom]}", reference_genome='GRCh38'
            )
            mt = hl.filter_intervals(mt, [interval])

            # Restrict to sampled loci only (the key operation)
            mt = mt.semi_join_rows(ht_chr_loci)

            # Filter to EUR samples
            mt = mt.filter_cols(eur_set.contains(mt.s))

            # Drop all unused fields to minimize data
            mt = mt.select_entries('GT')
            mt = mt.select_rows()
            mt = mt.select_cols()

            # Export PLINK with FID = IID = research_id (mt.s)
            plink_prefix = f"{plink_dir}/{chrom}_background"
            logger.info(f"  {chrom}: Exporting PLINK to {plink_prefix} ...")
            hl.export_plink(
                mt,
                plink_prefix,
                fam_id=mt.s,
                ind_id=mt.s,
            )

            chr_time = time.time() - chr_start
            per_chr_results[chrom] = {
                'plink_prefix': plink_prefix,
                'time_seconds': round(chr_time, 1),
                'status': 'success',
            }
            completed_chroms.add(chrom)
            logger.info(f"  {chrom}: PLINK exported in {chr_time:.1f}s")

            # Update progress after each chromosome
            progress_data = {
                'completed': sorted(completed_chroms),
                'last_updated': datetime.now().isoformat(),
            }
            with hl.hadoop_open(progress_path, 'w') as f:
                json.dump(progress_data, f, indent=2)

        except Exception as e:
            logger.error(f"  {chrom}: FAILED - {e}")
            failed_chroms.append(chrom)
            continue

    pass2_time = time.time() - pass2_start
    summary = {
        'pass': 2,
        'chromosomes_completed': len(completed_chroms),
        'chromosomes_failed': failed_chroms,
        'plink_output_dir': plink_dir,
        'per_chromosome': per_chr_results,
        'time_seconds': round(pass2_time, 1),
        'timestamp': datetime.now().isoformat(),
    }

    logger.info(f"\nPass 2 done: {len(completed_chroms)} chromosomes in {pass2_time:.1f}s")
    if failed_chroms:
        logger.warning(f"Failed chromosomes: {failed_chroms}")

    return summary


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    config = load_config("config/config.json")

    # Initialize Hail with hardened Spark config
    hl.init(
        log='/tmp/hail_background_snps_plink.log',
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

    # Source MT path
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

    # Paths for intermediate and final outputs
    loci_ht_path = f"{output_dir}/sampled_loci.ht"

    # Test mode vs production
    test_mode = config['params'].get('test_mode', False)
    test_chromosome = config['params'].get('test_chromosome', None)

    # Calculate per-chromosome targets
    chr_targets = calculate_chromosome_targets(target_snps)

    # In test mode, restrict to a single chromosome
    if test_mode and test_chromosome:
        logger.info(f"=== TEST MODE: {test_chromosome} only ===")
        if test_chromosome not in chr_targets:
            logger.error(f"{test_chromosome} not in autosome targets")
            sys.exit(1)
        chr_targets = {test_chromosome: chr_targets[test_chromosome]}
        logger.info(f"Target for {test_chromosome}: {chr_targets[test_chromosome]:,}")
    else:
        logger.info("=== PRODUCTION MODE: chr1-22 ===")

    try:
        # --- Pass 1: Sample loci (rows only, no entry data) ---
        if hfs.exists(f"{loci_ht_path}/_SUCCESS") or hfs.exists(f"{loci_ht_path}/metadata.json.gz"):
            logger.info(f"Sampled loci Table already exists at {loci_ht_path}, skipping Pass 1")
            pass1_summary = {'skipped': True, 'output_ht_path': loci_ht_path}
        else:
            pass1_summary = pass1_sample_loci(
                mt_path, loci_ht_path, chr_targets, config,
            )

            # Write Pass 1 summary
            p1_summary_path = f"{output_dir}/pass1_summary.json"
            with hl.hadoop_open(p1_summary_path, 'w') as f:
                json.dump(pass1_summary, f, indent=2, default=str)
            logger.info(f"Pass 1 summary written to {p1_summary_path}")

        # --- Pass 2: Extract genotypes + export PLINK ---
        pass2_summary = pass2_export_plink(
            mt_path, loci_ht_path, output_dir, chr_targets, config,
        )

        # Write Pass 2 summary
        p2_summary_path = f"{output_dir}/pass2_summary.json"
        with hl.hadoop_open(p2_summary_path, 'w') as f:
            json.dump(pass2_summary, f, indent=2, default=str)
        logger.info(f"Pass 2 summary written to {p2_summary_path}")

        # Write overall summary
        overall = {
            'target_snps': target_snps,
            'test_mode': test_mode,
            'test_chromosome': test_chromosome if test_mode else None,
            'output_directory': output_dir,
            'plink_dir': f"{output_dir}/null_model_data",
            'loci_table': loci_ht_path,
            'pass1': pass1_summary,
            'pass2': pass2_summary,
            'timestamp': datetime.now().isoformat(),
        }
        overall_path = f"{output_dir}/{base_filename}_overall_summary.json"
        with hl.hadoop_open(overall_path, 'w') as f:
            json.dump(overall, f, indent=2, default=str)

        logger.info("=" * 60)
        if test_mode:
            logger.info(f"=== TEST COMPLETE: {test_chromosome} ===")
        else:
            logger.info("=== PRODUCTION COMPLETE ===")
        logger.info(f"PLINK files: {output_dir}/null_model_data/")
        logger.info(f"Sampled loci Table: {loci_ht_path}")
        logger.info(f"Overall summary: {overall_path}")
        logger.info("=" * 60)

    except Exception as e:
        logger.error(f"Pipeline failed: {e}")
        raise


if __name__ == "__main__":
    main()
