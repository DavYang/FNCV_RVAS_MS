#!/usr/bin/env python3
"""Build background SNP PLINK files for Regenie Step 1 null model.

Samples ~500K common SNPs from the AoU ACAF splitMT and exports
per-chromosome PLINK files for EUR ancestry samples.

Pass 1 (cheap): Sample loci from the rows Table (no entry data).
Pass 2 (targeted): Extract genotypes at sampled loci, checkpoint,
                    then export_plink per chromosome.

Usage:
    python 01_build_background_snps.py
"""

import json
import os
import sys
import time
import traceback

import hail as hl
import hailtop.fs as hfs
from datetime import datetime

from utils import load_config, setup_logger

logger = setup_logger("background_snps")


def _fmt_elapsed(seconds: float) -> str:
    """Format elapsed seconds as human-readable string."""
    if seconds < 60:
        return f"{seconds:.0f}s"
    elif seconds < 3600:
        return f"{seconds / 60:.1f}m"
    else:
        h = int(seconds // 3600)
        m = int((seconds % 3600) // 60)
        return f"{h}h {m}m"


def _log_progress(
    current: int, total: int, label: str, start_time: float
) -> None:
    """Log a progress line with ETA."""
    elapsed = time.time() - start_time
    pct = (current / total * 100) if total > 0 else 0
    if current > 0:
        eta = (elapsed / current) * (total - current)
        eta_str = _fmt_elapsed(eta)
    else:
        eta_str = "?"
    logger.info(
        f"  [{current}/{total}] ({pct:.0f}%) {label} "
        f"| elapsed: {_fmt_elapsed(elapsed)} | ETA: {eta_str}"
    )


def _log_memory() -> None:
    """Log current system memory usage if /proc/meminfo is available."""
    try:
        with open('/proc/meminfo', 'r') as f:
            lines = {l.split(':')[0]: l.split(':')[1].strip()
                     for l in f.readlines() if ':' in l}
        total = lines.get('MemTotal', '?')
        avail = lines.get('MemAvailable', '?')
        logger.info(f"  [memory] total={total}, available={avail}")
    except Exception:
        pass

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

ACAF_VARIANTS_PER_BP = 0.03826
SAMPLING_OVERSHOOT = 1.02

# ---------------------------------------------------------------------------
# EUR sample Table (cached across chromosomes)
# ---------------------------------------------------------------------------
_eur_ht_cache = None


def _get_eur_samples_ht(config: dict, tmp_dir: str) -> hl.Table:
    """Return a checkpointed Hail Table of EUR sample IDs for semi_join_cols.

    Uses a Table-based approach instead of hl.literal(set(...)) to avoid
    embedding ~234K sample IDs into the Spark DAG, which crashes the JVM.
    Checkpointed once and reused across per-chromosome iterations.

    Args:
        config: Loaded config dict.
        tmp_dir: GCS path for temporary checkpoint files.

    Returns:
        Hail Table keyed by 's' (sample ID string).
    """
    global _eur_ht_cache
    if _eur_ht_cache is None:
        logger.info("Loading EUR ancestry sample IDs as Table ...")
        ancestry_ht = hl.import_table(
            config['inputs']['ancestry_pred'],
            impute=True,
            types={'research_id': hl.tstr, 'ancestry_pred': hl.tstr},
        )
        eur_ht = ancestry_ht.filter(ancestry_ht.ancestry_pred == 'eur')
        eur_ht = eur_ht.select(s=eur_ht.research_id).key_by('s')
        eur_ht_path = f"{tmp_dir}/eur_samples.ht"
        eur_ht = eur_ht.checkpoint(eur_ht_path, overwrite=True)
        n_eur = eur_ht.count()
        logger.info(f"EUR samples: {n_eur:,} (checkpointed to {eur_ht_path})")
        _eur_ht_cache = eur_ht
    return _eur_ht_cache


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
        Dict mapping chromosome name to per-chromosome target count.
    """
    total_size = sum(CHR_SIZES.values())
    chr_targets = {}
    for chrom in AUTOSOMES:
        size = CHR_SIZES[chrom]
        chr_targets[chrom] = max(1, int(target_snps * (size / total_size)))

    current_total = sum(chr_targets.values())
    remaining = target_snps - current_total
    if remaining != 0:
        sorted_chroms = sorted(
            chr_targets, key=lambda c: chr_targets[c], reverse=True
        )
        for i, chrom in enumerate(sorted_chroms):
            if i >= abs(remaining):
                break
            chr_targets[chrom] += 1 if remaining > 0 else -1

    logger.info(
        f"Total target: {sum(chr_targets.values()):,} "
        f"across {len(chr_targets)} autosomes"
    )
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
    _log_memory()

    logger.info(f"Reading rows Table from {mt_path} ...")
    mt = hl.read_matrix_table(mt_path)
    ht_rows = mt.rows().select()
    logger.info("Rows Table loaded (no entry data accessed)")

    chroms_to_process = [c for c in AUTOSOMES if chr_targets.get(c, 0) > 0]
    total_chroms = len(chroms_to_process)
    logger.info(f"Chromosomes to sample: {total_chroms}")

    per_chr_counts = {}
    per_chr_ht_paths = []
    completed_count = 0

    for chrom in chroms_to_process:
        chr_target = chr_targets[chrom]
        chr_ht_path = f"{output_ht_path}_parts/{chrom}.ht"
        per_chr_ht_paths.append(chr_ht_path)

        # Skip if this chromosome was already written (resume support)
        if (hfs.exists(f"{chr_ht_path}/_SUCCESS")
                or hfs.exists(f"{chr_ht_path}/metadata.json.gz")):
            ht_chr_final = hl.read_table(chr_ht_path)
            n = ht_chr_final.count()
            per_chr_counts[chrom] = {'sampled': n, 'target': chr_target}
            completed_count += 1
            logger.info(
                f"  {chrom}: {n:,} already sampled "
                f"(target: {chr_target:,}), skipping "
                f"[{completed_count}/{total_chroms}]"
            )
            continue

        chr_start = time.time()
        logger.info(f"  {chrom}: Starting sampling (target: {chr_target:,}) ...")

        try:
            interval = hl.parse_locus_interval(
                f"{chrom}:1-{CHR_SIZES[chrom]}", reference_genome='GRCh38'
            )
            ht_chr = hl.filter_intervals(ht_rows, [interval])

            if avoid_hla and chrom == 'chr6' and hla_region:
                hla_chrom = f"chr{hla_region['chrom']}"
                hla_start = hla_region['start']
                hla_end = hla_region['end']
                hla_interval = hl.parse_locus_interval(
                    f"{hla_chrom}:{hla_start}-{hla_end}",
                    reference_genome='GRCh38',
                )
                ht_chr = ht_chr.filter(~hla_interval.contains(ht_chr.locus))
                logger.info(
                    f"  {chrom}: Excluded HLA region "
                    f"{hla_chrom}:{hla_start}-{hla_end}"
                )

            estimated_pool = int(CHR_SIZES[chrom] * ACAF_VARIANTS_PER_BP)
            fraction = min(
                1.0, (chr_target / estimated_pool) * SAMPLING_OVERSHOOT
            )

            chr_seed = seed + AUTOSOMES.index(chrom)
            ht_sampled = ht_chr.filter(
                hl.rand_unif(0.0, 1.0, seed=chr_seed) < fraction
            )

            logger.info(
                f"  {chrom}: Writing sampled loci "
                f"(fraction={fraction:.6f}, est_pool={estimated_pool:,}) ..."
            )
            ht_sampled.write(chr_ht_path, overwrite=True)

            logger.info(f"  {chrom}: Counting sampled loci ...")
            ht_chr_final = hl.read_table(chr_ht_path)
            n = ht_chr_final.count()
            per_chr_counts[chrom] = {'sampled': n, 'target': chr_target}

            completed_count += 1
            chr_time = time.time() - chr_start
            _log_progress(
                completed_count, total_chroms,
                f"{chrom}: {n:,} sampled (target: {chr_target:,}) "
                f"in {_fmt_elapsed(chr_time)}",
                pass1_start,
            )
            _log_memory()

        except Exception as e:
            logger.error(
                f"  {chrom}: FAILED during sampling - {e}\n"
                f"{traceback.format_exc()}"
            )
            raise

    # Union all per-chromosome Tables into the final combined Table
    logger.info("Combining per-chromosome sampled loci into single Table ...")
    combine_start = time.time()
    chr_tables = [hl.read_table(p) for p in per_chr_ht_paths]
    ht_combined = chr_tables[0]
    for ht in chr_tables[1:]:
        ht_combined = ht_combined.union(ht)
    ht_combined.write(output_ht_path, overwrite=True)
    logger.info(
        f"Combined Table written in {_fmt_elapsed(time.time() - combine_start)}"
    )

    total_count = sum(c['sampled'] for c in per_chr_counts.values())

    # Summary table
    logger.info("")
    logger.info("Pass 1 Summary:")
    logger.info(f"  {'Chrom':<8} {'Sampled':>10} {'Target':>10} {'Diff':>8}")
    logger.info(f"  {'-'*8} {'-'*10} {'-'*10} {'-'*8}")
    for chrom in AUTOSOMES:
        if chrom in per_chr_counts:
            c = per_chr_counts[chrom]
            diff = c['sampled'] - c['target']
            sign = '+' if diff >= 0 else ''
            logger.info(
                f"  {chrom:<8} {c['sampled']:>10,} {c['target']:>10,} "
                f"{sign}{diff:>7,}"
            )
    logger.info(f"  {'TOTAL':<8} {total_count:>10,} {sum(chr_targets.values()):>10,}")
    logger.info("")

    pass1_time = time.time() - pass1_start
    global_target = sum(chr_targets.values())
    summary = {
        'pass': 1,
        'total_loci_sampled': total_count,
        'global_target': global_target,
        'per_chromosome': per_chr_counts,
        'output_ht_path': output_ht_path,
        'time_seconds': round(pass1_time, 1),
        'timestamp': datetime.now().isoformat(),
    }

    logger.info(
        f"Pass 1 COMPLETE: {total_count:,} loci sampled "
        f"in {_fmt_elapsed(pass1_time)}"
    )
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
      2. semi_join_rows against sampled loci
      3. semi_join_cols to EUR samples (Table-based, no hl.literal)
      4. select_entries('GT'), drop row/col fields
      5. checkpoint to break DAG
      6. export_plink with FID=IID=research_id

    Args:
        mt_path: GCS path to ACAF splitMT.
        loci_ht_path: GCS path to sampled loci Hail Table from Pass 1.
        output_dir: GCS base directory for output.
        chr_targets: Dict of chromosome -> target SNP count.
        config: Loaded config dict.

    Returns:
        Summary dict with per-chromosome results.
    """
    logger.info("Loading EUR samples Table ...")
    eur_load_start = time.time()
    eur_samples_ht = _get_eur_samples_ht(config, f"{output_dir}/tmp")
    logger.info(f"EUR samples loaded in {_fmt_elapsed(time.time() - eur_load_start)}")

    logger.info("=" * 60)
    logger.info("PASS 2: Extracting genotypes and exporting PLINK")
    logger.info("=" * 60)
    pass2_start = time.time()
    _log_memory()

    logger.info(f"Reading sampled loci from {loci_ht_path} ...")
    ht_loci = hl.read_table(loci_ht_path)

    plink_dir = f"{output_dir}/null_model_data"
    per_chr_results = {}
    failed_chroms = []

    # Resume support
    progress_path = f"{output_dir}/_pass2_progress.json"
    completed_chroms = set()
    try:
        with hl.hadoop_open(progress_path, 'r') as f:
            progress = json.load(f)
            completed_chroms = set(progress.get('completed', []))
        if completed_chroms:
            logger.info(
                f"Resuming Pass 2: {len(completed_chroms)} chromosomes "
                f"already done: {sorted(completed_chroms)}"
            )
    except Exception:
        pass

    chroms_to_process = [c for c in AUTOSOMES if c in chr_targets]
    total_chroms = len(chroms_to_process)
    export_count = 0
    logger.info(f"Chromosomes to export: {total_chroms}")

    for chrom in chroms_to_process:
        if chrom in completed_chroms:
            export_count += 1
            logger.info(
                f"  {chrom}: Already completed, skipping "
                f"[{export_count}/{total_chroms}]"
            )
            continue

        checkpoint_path = f"{output_dir}/tmp/{chrom}_checkpoint.mt"
        try:
            chr_start = time.time()
            logger.info("")
            logger.info(f"  {chrom}: === Starting PLINK export ===")
            _log_memory()

            logger.info(f"  {chrom}: Step 1/5 - Filtering loci and MT by interval ...")
            step_start = time.time()
            chr_interval = hl.parse_locus_interval(
                f"{chrom}:1-{CHR_SIZES[chrom]}",
                reference_genome='GRCh38',
            )
            ht_chr_loci = hl.filter_intervals(ht_loci, [chr_interval])
            mt = hl.read_matrix_table(mt_path)
            mt = hl.filter_intervals(mt, [chr_interval])
            logger.info(f"  {chrom}: Step 1/5 done ({_fmt_elapsed(time.time() - step_start)})")

            logger.info(f"  {chrom}: Step 2/5 - Joining rows (sampled loci) and cols (EUR) ...")
            step_start = time.time()
            mt = mt.semi_join_rows(ht_chr_loci)
            mt = mt.semi_join_cols(eur_samples_ht)
            mt = mt.select_entries('GT')
            mt = mt.select_rows()
            mt = mt.select_cols()
            mt = mt.naive_coalesce(max(1, mt.n_partitions() // 50))
            logger.info(f"  {chrom}: Step 2/5 done ({_fmt_elapsed(time.time() - step_start)})")

            logger.info(f"  {chrom}: Step 3/5 - Checkpointing filtered MT (breaks DAG) ...")
            step_start = time.time()
            mt = mt.checkpoint(checkpoint_path, overwrite=True)
            n_variants = mt.count_rows()
            n_samples = mt.count_cols()
            logger.info(
                f"  {chrom}: Step 3/5 done - {n_variants:,} variants x "
                f"{n_samples:,} samples ({_fmt_elapsed(time.time() - step_start)})"
            )

            logger.info(f"  {chrom}: Step 4/5 - Exporting PLINK ...")
            step_start = time.time()
            plink_prefix = f"{plink_dir}/{chrom}_background"
            hl.export_plink(
                mt,
                plink_prefix,
                fam_id=mt.s,
                ind_id=mt.s,
            )
            logger.info(f"  {chrom}: Step 4/5 done ({_fmt_elapsed(time.time() - step_start)})")

            logger.info(f"  {chrom}: Step 5/5 - Cleaning up checkpoint ...")
            _cleanup_gcs_path(checkpoint_path)

            chr_time = time.time() - chr_start
            per_chr_results[chrom] = {
                'plink_prefix': plink_prefix,
                'n_variants': n_variants,
                'n_samples': n_samples,
                'time_seconds': round(chr_time, 1),
                'status': 'success',
            }
            completed_chroms.add(chrom)
            export_count += 1

            _log_progress(
                export_count, total_chroms,
                f"{chrom}: {n_variants:,} variants exported "
                f"in {_fmt_elapsed(chr_time)}",
                pass2_start,
            )
            _log_memory()

            # Persist progress
            progress_data = {
                'completed': sorted(completed_chroms),
                'last_updated': datetime.now().isoformat(),
            }
            with hl.hadoop_open(progress_path, 'w') as f:
                json.dump(progress_data, f, indent=2)

        except Exception as e:
            logger.error(
                f"  {chrom}: FAILED - {e}\n"
                f"{traceback.format_exc()}"
            )
            _cleanup_gcs_path(checkpoint_path)
            failed_chroms.append(chrom)
            continue

    # Summary table
    pass2_time = time.time() - pass2_start
    logger.info("")
    logger.info("Pass 2 Summary:")
    logger.info(f"  {'Chrom':<8} {'Variants':>10} {'Samples':>10} {'Time':>10} {'Status':>10}")
    logger.info(f"  {'-'*8} {'-'*10} {'-'*10} {'-'*10} {'-'*10}")
    for chrom in chroms_to_process:
        if chrom in per_chr_results:
            r = per_chr_results[chrom]
            logger.info(
                f"  {chrom:<8} {r.get('n_variants', '?'):>10,} "
                f"{r.get('n_samples', '?'):>10,} "
                f"{_fmt_elapsed(r['time_seconds']):>10} "
                f"{'OK':>10}"
            )
        elif chrom in completed_chroms:
            logger.info(f"  {chrom:<8} {'(resumed)':>10} {'':>10} {'':>10} {'OK':>10}")
        elif chrom in failed_chroms:
            logger.info(f"  {chrom:<8} {'':>10} {'':>10} {'':>10} {'FAILED':>10}")
    logger.info("")

    summary = {
        'pass': 2,
        'chromosomes_completed': len(completed_chroms),
        'chromosomes_failed': failed_chroms,
        'plink_output_dir': plink_dir,
        'per_chromosome': per_chr_results,
        'time_seconds': round(pass2_time, 1),
        'timestamp': datetime.now().isoformat(),
    }

    logger.info(
        f"Pass 2 COMPLETE: {len(completed_chroms)}/{total_chroms} chromosomes "
        f"in {_fmt_elapsed(pass2_time)}"
    )
    if failed_chroms:
        logger.warning(f"FAILED chromosomes ({len(failed_chroms)}): {failed_chroms}")

    return summary


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    pipeline_start = time.time()
    config = load_config("config/config.json")

    logger.info("=" * 60)
    logger.info("PIPELINE START")
    logger.info("=" * 60)
    logger.info(f"Timestamp: {datetime.now().isoformat()}")
    logger.info(f"PID: {os.getpid()}")
    _log_memory()

    hl.init(log='/tmp/hail_background_snps.log')
    hl.default_reference('GRCh38')

    workspace_bucket = os.environ.get('WORKSPACE_BUCKET')
    if not workspace_bucket:
        logger.error("WORKSPACE_BUCKET environment variable not set")
        sys.exit(1)
    if not workspace_bucket.startswith('gs://'):
        workspace_bucket = f"gs://{workspace_bucket}"

    mt_path = config['inputs']['wgs_matrix_table']
    target_snps = config['sampling']['target_total_snps']
    test_mode = config['params'].get('test_mode', False)
    test_chromosome = config['params'].get('test_chromosome', None)

    # Log configuration
    logger.info("Configuration:")
    logger.info(f"  Workspace bucket: {workspace_bucket}")
    logger.info(f"  Source MT: {mt_path}")
    logger.info(f"  Target SNPs: {target_snps:,}")
    logger.info(f"  Avoid HLA: {config['sampling'].get('avoid_hla', False)}")
    logger.info(f"  Random seed: {config['sampling']['random_seed']}")
    logger.info(f"  Test mode: {test_mode}")
    if test_mode:
        logger.info(f"  Test chromosome: {test_chromosome}")

    if target_snps >= 1_000_000:
        snp_label = f"{target_snps // 1_000_000}M"
    elif target_snps >= 1_000:
        snp_label = f"{target_snps // 1_000}K"
    else:
        snp_label = str(target_snps)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = (
        f"{workspace_bucket}/results/FNCV_RVAS_MS/"
        f"{snp_label}_background_snps_{timestamp}"
    )
    logger.info(f"  Output directory: {output_dir}")
    logger.info("")

    chr_targets = calculate_chromosome_targets(target_snps)

    if test_mode and test_chromosome:
        logger.info(f"=== TEST MODE: {test_chromosome} only ===")
        if test_chromosome not in chr_targets:
            logger.error(f"{test_chromosome} not in autosome targets")
            sys.exit(1)
        chr_targets = {test_chromosome: chr_targets[test_chromosome]}
        logger.info(
            f"Target for {test_chromosome}: "
            f"{chr_targets[test_chromosome]:,}"
        )
    else:
        logger.info("=== PRODUCTION MODE: chr1-22 ===")

    loci_ht_path = f"{output_dir}/sampled_loci.ht"

    try:
        # --- Pass 1: Sample loci (rows only) ---
        if (hfs.exists(f"{loci_ht_path}/_SUCCESS")
                or hfs.exists(f"{loci_ht_path}/metadata.json.gz")):
            logger.info(
                f"Sampled loci Table exists at {loci_ht_path}, "
                "skipping Pass 1"
            )
            pass1_summary = {
                'skipped': True,
                'output_ht_path': loci_ht_path,
            }
        else:
            pass1_summary = pass1_sample_loci(
                mt_path, loci_ht_path, chr_targets, config,
            )

        p1_path = f"{output_dir}/pass1_summary.json"
        with hl.hadoop_open(p1_path, 'w') as f:
            json.dump(pass1_summary, f, indent=2, default=str)
        logger.info(f"Pass 1 summary written to: {p1_path}")

        # --- Pass 2: Extract genotypes + export PLINK ---
        pass2_summary = pass2_export_plink(
            mt_path, loci_ht_path, output_dir, chr_targets, config,
        )

        try:
            p2_path = f"{output_dir}/pass2_summary.json"
            with hl.hadoop_open(p2_path, 'w') as f:
                json.dump(pass2_summary, f, indent=2, default=str)

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
            overall_path = f"{output_dir}/{snp_label}_background_snps_summary.json"
            with hl.hadoop_open(overall_path, 'w') as f:
                json.dump(overall, f, indent=2, default=str)
        except Exception as summary_err:
            logger.warning(
                f"Failed to write summary JSONs: {summary_err}\n"
                f"{traceback.format_exc()}"
            )

        pipeline_time = time.time() - pipeline_start
        logger.info("")
        logger.info("=" * 60)
        if test_mode:
            logger.info(f"=== TEST COMPLETE: {test_chromosome} ===")
        else:
            logger.info("=== PRODUCTION COMPLETE ===")
        logger.info(f"Total pipeline time: {_fmt_elapsed(pipeline_time)}")
        logger.info(f"PLINK files: {output_dir}/null_model_data/")
        logger.info(f"Sampled loci: {loci_ht_path}")
        logger.info(f"Summaries: {output_dir}/")
        n_failed = len(pass2_summary.get('chromosomes_failed', []))
        n_done = pass2_summary.get('chromosomes_completed', 0)
        logger.info(f"Chromosomes: {n_done} succeeded, {n_failed} failed")
        logger.info("=" * 60)

    except Exception as e:
        pipeline_time = time.time() - pipeline_start
        logger.error(
            f"Pipeline FAILED after {_fmt_elapsed(pipeline_time)}: {e}\n"
            f"{traceback.format_exc()}"
        )
        raise


if __name__ == "__main__":
    main()
