#!/usr/bin/env python3
"""Build background SNP PLINK files for Regenie Step 1 null model.

Processes a SINGLE chromosome per invocation: samples loci from the AoU
ACAF splitMT rows Table, then extracts EUR genotypes and exports PLINK.
The bash wrapper calls this script once per chromosome (chr1-22), giving
each invocation a fresh Hail/JVM session to avoid accumulated state crashes.

Pipeline per chromosome:
  1. count_qc_pool()      - Count real post-filter variant pool (rows-only, no GT).
                            Uses pre-computed info.AF and variant_qc.call_rate.
                            Result cached to pool_count.txt to avoid recomputation.
  2. Retry loop:
     a. sample_loci()     - Bernoulli sample from QC-passing rows at fraction
                            computed from real pool size * oversample_factor.
     b. extract_and_count() - Extract EUR GTs, apply safety call-rate filter,
                            checkpoint, count. No PLINK export on failed attempts.
     c. If count >= target: export_plink_final() once, break.
        Else: boost fraction, discard checkpoint, retry.

Usage:
    python 01_build_background_snps.py --chrom chr21 --target 11700 \
        --output-dir gs://bucket/results/FNCV_RVAS_MS/500K_background_snps \
        [--force-rerun]
"""

import argparse
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

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
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

DEFAULT_MIN_CALL_RATE = 0.95
DEFAULT_HAIL_MIN_MAF = 0.01
DEFAULT_HAIL_MIN_VARIANT_CALL_RATE = 0.99
DEFAULT_RESAMPLE_MAX_ATTEMPTS = 5
DEFAULT_OVERSAMPLE_FACTOR = 2.5


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
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


def _cleanup_gcs_path(path: str) -> None:
    """Recursively remove a GCS path, logging but not raising on failure."""
    try:
        hfs.rmtree(path)
        logger.info(f"  Cleaned up: {path}")
    except Exception as e:
        logger.warning(f"  Failed to clean up {path}: {e}")


def _get_eur_samples_ht(config: dict, shared_dir: str) -> hl.Table:
    """Return a Hail Table of EUR sample IDs for semi_join_cols.

    Uses a Table-based approach instead of hl.literal(set(...)) to avoid
    embedding ~234K sample IDs into the Spark DAG, which crashes the JVM.
    The Table is written to a shared path so it persists across chromosome
    invocations without re-computing.

    Args:
        config: Loaded config dict.
        shared_dir: GCS path for shared intermediate files.

    Returns:
        Hail Table keyed by 's' (sample ID string).
    """
    eur_ht_path = f"{shared_dir}/eur_samples.ht"

    if (hfs.exists(f"{eur_ht_path}/_SUCCESS")
            or hfs.exists(f"{eur_ht_path}/metadata.json.gz")):
        logger.info(f"  EUR samples Table exists at {eur_ht_path}")
        eur_ht = hl.read_table(eur_ht_path)
        n_eur = eur_ht.count()
        logger.info(f"  EUR samples: {n_eur:,}")
        return eur_ht

    logger.info("  Creating EUR samples Table from ancestry TSV ...")
    ancestry_ht = hl.import_table(
        config['inputs']['ancestry_pred'],
        impute=True,
        types={'research_id': hl.tstr, 'ancestry_pred': hl.tstr},
    )
    eur_ht = ancestry_ht.filter(ancestry_ht.ancestry_pred == 'eur')
    eur_ht = eur_ht.select(s=eur_ht.research_id).key_by('s')
    eur_ht.write(eur_ht_path, overwrite=True)
    eur_ht = hl.read_table(eur_ht_path)
    n_eur = eur_ht.count()
    logger.info(f"  EUR samples: {n_eur:,} (written to {eur_ht_path})")
    return eur_ht


# ---------------------------------------------------------------------------
# Step 0: Count the real QC-passing variant pool for this chromosome
# ---------------------------------------------------------------------------
def count_qc_pool(
    mt_path: str,
    chrom: str,
    pool_count_path: str,
    config: dict,
    force: bool = False,
) -> int:
    """Count QC-passing variants for this chromosome using rows Table only.

    Applies the same pre-filters used in sample_loci() (global info.AF proxy
    for MAF, pre-computed variant_qc.call_rate, HLA exclusion on chr6) but
    counts rather than samples. No entry (genotype) data is read.

    The count is cached to pool_count_path so that retries within the same
    invocation do not re-scan; also safe to reuse across re-invocations unless
    force=True.

    Args:
        mt_path: GCS path to ACAF splitMT.
        chrom: Chromosome name, e.g. 'chr21'.
        pool_count_path: GCS path to a plain-text file caching the count.
        config: Loaded config dict.
        force: If True, ignore cached count and recompute.

    Returns:
        Integer count of QC-passing variants in the pool.
    """
    # Check cache first (skip if force=True)
    if not force and hfs.exists(pool_count_path):
        try:
            with hl.hadoop_open(pool_count_path, 'r') as f:
                cached = int(f.read().strip())
            logger.info(
                f"  Step 0: Loaded cached pool count: {cached:,} "
                f"(from {pool_count_path})"
            )
            return cached
        except Exception as e:
            logger.warning(f"  Step 0: Could not read cached pool count ({e}); recomputing")

    step_start = time.time()
    hail_min_maf = config['sampling'].get('hail_min_maf', DEFAULT_HAIL_MIN_MAF)
    hail_min_cr = config['sampling'].get(
        'hail_min_variant_call_rate', DEFAULT_HAIL_MIN_VARIANT_CALL_RATE
    )
    avoid_hla = config['sampling'].get('avoid_hla', False)
    hla_region = config['sampling'].get('hla_region', {})

    logger.info(
        f"  Step 0: Counting QC-passing pool for {chrom} "
        f"(info.AF > {hail_min_maf}, call_rate >= {hail_min_cr}) ..."
    )
    _log_memory()

    t0 = time.time()
    mt = hl.read_matrix_table(mt_path)
    ht_rows = mt.rows()

    chr_size = CHR_SIZES[chrom]
    interval = hl.parse_locus_interval(
        f"{chrom}:1-{chr_size}", reference_genome='GRCh38'
    )
    ht_chr = hl.filter_intervals(ht_rows, [interval])
    ht_chr = ht_chr.filter(
        (ht_chr.info.AF[0] > hail_min_maf) &
        (ht_chr.info.AF[0] < (1.0 - hail_min_maf)) &
        (ht_chr.variant_qc.call_rate >= hail_min_cr)
    )

    if avoid_hla and chrom == 'chr6' and hla_region:
        hla_chrom = f"chr{hla_region['chrom']}"
        hla_start = hla_region['start']
        hla_end = hla_region['end']
        hla_interval = hl.parse_locus_interval(
            f"{hla_chrom}:{hla_start}-{hla_end}",
            reference_genome='GRCh38',
        )
        ht_chr = ht_chr.filter(~hla_interval.contains(ht_chr.locus))
        logger.info(f"  Excluded HLA region {hla_chrom}:{hla_start}-{hla_end}")

    n_pool = ht_chr.count()
    logger.info(
        f"  Step 0 DONE: pool={n_pool:,} variants "
        f"({_fmt_elapsed(time.time() - step_start)})"
    )

    # Cache to GCS
    try:
        with hl.hadoop_open(pool_count_path, 'w') as f:
            f.write(str(n_pool))
        logger.info(f"  Pool count cached to {pool_count_path}")
    except Exception as e:
        logger.warning(f"  Could not cache pool count ({e}); will recompute on re-run")

    _log_memory()
    return n_pool


# ---------------------------------------------------------------------------
# Step 1: Sample loci for this chromosome (rows Table only, no entry data)
# ---------------------------------------------------------------------------
def sample_loci(
    mt_path: str,
    chrom: str,
    target: int,
    output_ht_path: str,
    config: dict,
    fraction: float,
    seed_offset: int = 0,
) -> dict:
    """Sample loci on a single chromosome using Bernoulli sampling.

    Reads only the rows Table of the ACAF splitMT (no entry data touched).
    Pre-filters by pre-computed info.AF (global AF > hail_min_maf as a proxy
    for EUR MAF) and variant_qc.call_rate. HLA exclusion is applied on chr6
    if configured. The retry loop lives in main(); this function performs a
    single sampling attempt at the given fraction.

    Args:
        mt_path: GCS path to ACAF splitMT.
        chrom: Chromosome name, e.g. 'chr21'.
        target: Target number of SNPs (used for logging only).
        output_ht_path: GCS path to write the sampled loci Table.
        config: Loaded config dict.
        fraction: Bernoulli sampling fraction (computed in main() from real pool).
        seed_offset: Added to the base chromosome seed to vary draws on retries.

    Returns:
        Dict with 'n_sampled', 'target', 'fraction', 'time_seconds'.
    """
    step_start = time.time()
    seed = config['sampling']['random_seed']
    avoid_hla = config['sampling'].get('avoid_hla', False)
    hla_region = config['sampling'].get('hla_region', {})
    hail_min_maf = config['sampling'].get('hail_min_maf', DEFAULT_HAIL_MIN_MAF)
    hail_min_cr = config['sampling'].get(
        'hail_min_variant_call_rate', DEFAULT_HAIL_MIN_VARIANT_CALL_RATE
    )

    logger.info(f"  Step 1: Sampling loci (target={target:,}, fraction={fraction:.6f})")
    logger.info(
        f"  Pre-filters: info.AF > {hail_min_maf} (global AF MAF proxy), "
        f"variant_qc.call_rate >= {hail_min_cr} (rows-only, no GT scan)"
    )
    _log_memory()

    logger.info(f"  Reading rows Table from {mt_path} ...")
    t0 = time.time()
    mt = hl.read_matrix_table(mt_path)
    ht_rows = mt.rows()

    chr_size = CHR_SIZES[chrom]
    interval = hl.parse_locus_interval(
        f"{chrom}:1-{chr_size}", reference_genome='GRCh38'
    )
    ht_chr = hl.filter_intervals(ht_rows, [interval])

    # Pre-filter by pre-computed annotations (rows-only, zero genotype scan)
    # info.AF is global (pan-ancestry) -- used as an approximation of EUR MAF.
    # Variants with global AF > hail_min_maf (0.1) are very likely to also have
    # MAF > 0.01-0.05 in EUR, which is the downstream PLINK --maf target.
    ht_chr = ht_chr.filter(
        (ht_chr.info.AF[0] > hail_min_maf) &
        (ht_chr.info.AF[0] < (1.0 - hail_min_maf)) &
        (ht_chr.variant_qc.call_rate >= hail_min_cr)
    )
    logger.info(f"  Rows Table read and pre-filter applied ({_fmt_elapsed(time.time() - t0)})")

    if avoid_hla and chrom == 'chr6' and hla_region:
        hla_chrom = f"chr{hla_region['chrom']}"
        hla_start = hla_region['start']
        hla_end = hla_region['end']
        hla_interval = hl.parse_locus_interval(
            f"{hla_chrom}:{hla_start}-{hla_end}",
            reference_genome='GRCh38',
        )
        ht_chr = ht_chr.filter(~hla_interval.contains(ht_chr.locus))
        logger.info(f"  Excluded HLA region {hla_chrom}:{hla_start}-{hla_end}")

    autosomes = [f"chr{i}" for i in range(1, 23)]
    chr_seed = seed + autosomes.index(chrom) + seed_offset

    logger.info(f"  Sampling: seed={chr_seed}")
    ht_sampled = ht_chr.filter(
        hl.rand_unif(0.0, 1.0, seed=chr_seed) < fraction
    )
    ht_sampled = ht_sampled.select()  # drop all row fields except key
    ht_sampled.write(output_ht_path, overwrite=True)
    ht_sampled = hl.read_table(output_ht_path)
    n_sampled = ht_sampled.count()
    logger.info(f"  Sampled {n_sampled:,} loci")

    elapsed = time.time() - step_start
    logger.info(
        f"  Step 1 DONE: {n_sampled:,} loci sampled "
        f"(target={target:,}) in {_fmt_elapsed(elapsed)}"
    )
    _log_memory()

    return {
        'n_sampled': n_sampled,
        'target': target,
        'fraction': round(fraction, 6),
        'time_seconds': round(elapsed, 1),
    }


# ---------------------------------------------------------------------------
# Step 2a: Extract EUR genotypes, apply call-rate filter, checkpoint + count
# ---------------------------------------------------------------------------
def extract_and_count(
    mt_path: str,
    chrom: str,
    loci_ht_path: str,
    eur_samples_ht: hl.Table,
    checkpoint_path: str,
    config: dict,
) -> tuple:
    """Filter MT to sampled loci + EUR samples, apply call-rate filter,
    checkpoint, and count QC-passing variants. Does NOT export PLINK.

    Separating count from export lets the retry loop skip the expensive
    export on failed attempts, paying that cost only once when target is met.

    Args:
        mt_path: GCS path to ACAF splitMT.
        chrom: Chromosome name, e.g. 'chr21'.
        loci_ht_path: GCS path to sampled loci Table for this chromosome.
        eur_samples_ht: Hail Table of EUR sample IDs keyed by 's'.
        checkpoint_path: GCS path for the intermediate checkpoint MT.
        config: Loaded config dict.

    Returns:
        Tuple of (mt, n_passing, n_samples, elapsed_seconds) where mt is the
        checkpointed MatrixTable ready for export.
    """
    step_start = time.time()
    min_call_rate = config['sampling'].get('min_call_rate', DEFAULT_MIN_CALL_RATE)

    logger.info(f"  Step 2a: Extracting EUR genotypes (call_rate >= {min_call_rate})")
    _log_memory()

    logger.info(f"  2a-interval: Filtering MT to {chrom} interval ...")
    t0 = time.time()
    chr_size = CHR_SIZES[chrom]
    chr_interval = hl.parse_locus_interval(
        f"{chrom}:1-{chr_size}", reference_genome='GRCh38'
    )
    ht_loci = hl.read_table(loci_ht_path)
    mt = hl.read_matrix_table(mt_path)
    mt = hl.filter_intervals(mt, [chr_interval])
    logger.info(f"  2a-interval done ({_fmt_elapsed(time.time() - t0)})")

    logger.info(f"  2a-join: Joining to sampled loci + EUR samples ...")
    t0 = time.time()
    mt = mt.semi_join_rows(ht_loci)
    mt = mt.semi_join_cols(eur_samples_ht)
    mt = mt.select_entries('GT')
    mt = mt.select_rows()
    mt = mt.select_cols()
    logger.info(f"  2a-join done ({_fmt_elapsed(time.time() - t0)})")

    logger.info(f"  2a-filter: Applying call_rate >= {min_call_rate} on EUR GTs ...")
    t0 = time.time()
    mt = mt.filter_rows(
        hl.agg.fraction(hl.is_defined(mt.GT)) >= min_call_rate
    )
    logger.info(f"  2a-filter done ({_fmt_elapsed(time.time() - t0)})")

    logger.info(f"  2a-coalesce: Repartitioning to 200 partitions ...")
    t0 = time.time()
    mt = mt.naive_coalesce(200)
    logger.info(f"  2a-coalesce done ({_fmt_elapsed(time.time() - t0)})")

    logger.info(f"  2a-checkpoint: Checkpointing (breaks DAG for count) ...")
    t0 = time.time()
    mt = mt.checkpoint(checkpoint_path, overwrite=True)
    n_passing = mt.count_rows()
    n_samples = mt.count_cols()
    logger.info(
        f"  2a-checkpoint done: {n_passing:,} variants x {n_samples:,} samples "
        f"({_fmt_elapsed(time.time() - t0)})"
    )

    elapsed = time.time() - step_start
    logger.info(
        f"  Step 2a DONE: {n_passing:,} variants pass EUR QC "
        f"in {_fmt_elapsed(elapsed)}"
    )
    _log_memory()

    return mt, n_passing, n_samples, round(elapsed, 1)


# ---------------------------------------------------------------------------
# Step 2b: Export PLINK from already-checkpointed MatrixTable
# ---------------------------------------------------------------------------
def export_plink_final(
    mt: hl.MatrixTable,
    plink_prefix: str,
    checkpoint_path: str,
) -> dict:
    """Export PLINK files from a checkpointed MatrixTable and clean up.

    Called only after extract_and_count() confirms the variant count meets the
    target, so the expensive checkpoint+export is paid exactly once.

    Args:
        mt: Checkpointed MatrixTable returned by extract_and_count().
        plink_prefix: Output prefix for PLINK files (.bed/.bim/.fam).
        checkpoint_path: GCS path of the checkpoint MT to clean up after export.

    Returns:
        Dict with 'plink_prefix', 'time_seconds'.
    """
    step_start = time.time()

    logger.info(f"  Step 2b: Exporting PLINK to {plink_prefix} ...")
    _log_memory()

    t0 = time.time()
    hl.export_plink(mt, plink_prefix, fam_id=mt.s, ind_id=mt.s)
    logger.info(f"  2b-export done ({_fmt_elapsed(time.time() - t0)})")

    logger.info(f"  2b-cleanup: Removing checkpoint {checkpoint_path} ...")
    _cleanup_gcs_path(checkpoint_path)

    elapsed = time.time() - step_start
    logger.info(f"  Step 2b DONE: PLINK exported in {_fmt_elapsed(elapsed)}")
    _log_memory()

    return {
        'plink_prefix': plink_prefix,
        'time_seconds': round(elapsed, 1),
    }


# ---------------------------------------------------------------------------
# Main: process a single chromosome end-to-end
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Process one chromosome: sample loci + export PLINK"
    )
    parser.add_argument(
        '--chrom', required=True,
        help='Chromosome to process, e.g. chr21'
    )
    parser.add_argument(
        '--target', required=True, type=int,
        help='Number of SNPs to sample on this chromosome'
    )
    parser.add_argument(
        '--output-dir', required=True,
        help='GCS base output directory (shared across chromosomes)'
    )
    parser.add_argument(
        '--config', default='config/config.json',
        help='Path to config JSON file (default: config/config.json)'
    )
    parser.add_argument(
        '--force-rerun', action='store_true', default=False,
        help='Ignore existing summary.json and pool_count.txt; overwrite all outputs'
    )
    return parser.parse_args()


def main() -> None:
    """Run the single-chromosome pipeline."""
    args = parse_args()
    chrom = args.chrom
    target = args.target
    output_dir = args.output_dir
    config_path = args.config
    force_rerun = args.force_rerun

    run_start = time.time()

    logger.info("=" * 60)
    logger.info(f"SINGLE-CHROMOSOME PIPELINE: {chrom}")
    logger.info("=" * 60)
    logger.info(f"  Timestamp   : {datetime.now().isoformat()}")
    logger.info(f"  PID         : {os.getpid()}")
    logger.info(f"  Chrom       : {chrom}")
    logger.info(f"  Target      : {target:,}")
    logger.info(f"  Output dir  : {output_dir}")
    logger.info(f"  Config      : {config_path}")
    logger.info(f"  Force rerun : {force_rerun}")
    _log_memory()

    if chrom not in CHR_SIZES:
        logger.error(f"Unknown chromosome: {chrom}")
        sys.exit(1)

    config = load_config(config_path)
    mt_path = config['inputs']['wgs_matrix_table']

    chrom_dir = f"{output_dir}/{chrom}"
    shared_dir = f"{output_dir}/shared"
    tmp_dir = f"{output_dir}/tmp"
    summary_path = f"{chrom_dir}/summary.json"
    pool_count_path = f"{chrom_dir}/pool_count.txt"

    logger.info(f"  Source MT   : {mt_path}")
    logger.info("")

    # Initialize Hail (fresh session for this chromosome)
    # Must init before hfs.exists() since hailtop.fs needs Hadoop/GCS access
    hail_log = f"/tmp/hail_{chrom}.log"
    logger.info(f"  Initializing Hail (log: {hail_log}) ...")
    hl.init(log=hail_log)
    hl.default_reference('GRCh38')
    logger.info("  Hail initialized")

    # Resume: skip if summary.json status=success AND not force_rerun
    if not force_rerun and hfs.exists(summary_path):
        try:
            with hl.hadoop_open(summary_path, 'r') as f:
                prior = json.load(f)
            if prior.get('status') == 'success':
                logger.info(f"  summary.json exists with status=success at {summary_path}")
                logger.info(f"  Skipping {chrom} (already completed). Use --force-rerun to override.")
                hl.stop()
                return
            else:
                logger.warning(
                    f"  summary.json exists but status={prior.get('status')!r} — "
                    f"re-running {chrom}"
                )
        except Exception as read_err:
            logger.warning(
                f"  Could not read prior summary.json ({read_err}); re-running {chrom}"
            )

    summary = {
        'chrom': chrom,
        'target': target,
        'status': 'running',
        'start_time': datetime.now().isoformat(),
        'force_rerun': force_rerun,
    }

    try:
        # --- EUR samples (shared, written once, reused) ---
        logger.info("  Loading EUR samples Table ...")
        eur_samples_ht = _get_eur_samples_ht(config, shared_dir)

        # --- Step 0: Count real QC-passing pool for this chromosome ---
        # Uses pre-computed info.AF (global AF > hail_min_maf as EUR MAF proxy)
        # and variant_qc.call_rate. Rows-only scan, no genotype data read.
        # Cached to pool_count.txt; recomputed only if force_rerun=True.
        oversample_factor = config['sampling'].get(
            'oversample_factor', DEFAULT_OVERSAMPLE_FACTOR
        )
        n_pool = count_qc_pool(
            mt_path, chrom, pool_count_path, config, force=force_rerun
        )
        if n_pool == 0:
            raise RuntimeError(
                f"{chrom}: QC-passing pool is empty after info.AF and call_rate filters. "
                f"Check hail_min_maf and hail_min_variant_call_rate in config."
            )

        logger.info(
            f"  Pool: {n_pool:,} variants | target: {target:,} | "
            f"oversample_factor: {oversample_factor}"
        )

        # --- Steps 1 + 2: Sample loci, extract GTs, check count, export ---
        # Fraction is computed from the REAL pool size, not an estimated constant.
        # export_plink_final() is called only when the target is met, so the
        # expensive checkpoint+export cost is paid at most once per chromosome.
        loci_ht_path = f"{chrom_dir}/sampled_loci.ht"
        plink_prefix = f"{chrom_dir}/{chrom}_background"
        checkpoint_path = f"{tmp_dir}/{chrom}_checkpoint.mt"

        max_attempts = config['sampling'].get(
            'resample_max_attempts', DEFAULT_RESAMPLE_MAX_ATTEMPTS
        )

        fraction = min(1.0, (target * oversample_factor) / n_pool)

        loci_summary = None
        export_summary = None
        n_passing = 0
        n_samples = 0
        attempt = 0

        for attempt in range(1, max_attempts + 1):
            logger.info(
                f"  --- Attempt {attempt}/{max_attempts} "
                f"(fraction={fraction:.6f}) ---"
            )
            seed_offset = (attempt - 1) * 100

            # Step 1: sample QC-passing loci (rows-only, no GT scan)
            loci_summary = sample_loci(
                mt_path, chrom, target, loci_ht_path, config,
                fraction=fraction,
                seed_offset=seed_offset,
            )

            # Step 2a: extract EUR GTs, apply safety call-rate filter,
            # checkpoint, and count — no PLINK export yet
            mt, n_passing, n_samples, extract_elapsed = extract_and_count(
                mt_path, chrom, loci_ht_path,
                eur_samples_ht, checkpoint_path, config,
            )

            logger.info(
                f"  Attempt {attempt}: {n_passing:,} variants passed EUR QC "
                f"(target={target:,})"
            )

            if n_passing >= target:
                logger.info(
                    f"  Target met ({n_passing:,} >= {target:,}) on attempt {attempt}"
                )
                # Step 2b: export PLINK only now that target is confirmed
                export_summary = export_plink_final(mt, plink_prefix, checkpoint_path)
                break

            # Target not met: discard checkpoint, boost fraction, retry
            logger.warning(
                f"  Attempt {attempt}: {n_passing:,} < target {target:,}."
            )
            _cleanup_gcs_path(checkpoint_path)

            if attempt < max_attempts:
                boost = (target / n_passing) * 1.1 if n_passing > 0 else 2.0
                fraction = min(1.0, fraction * boost)
                logger.warning(
                    f"  Boosting fraction to {fraction:.6f} for attempt {attempt + 1}"
                )

        if n_passing < target:
            logger.warning(
                f"  {chrom}: reached {n_passing:,} variants after "
                f"{max_attempts} attempt(s); below target {target:,}. "
                f"Proceeding with available variants."
            )
            if export_summary is None:
                # Final attempt did not export; export whatever we have
                export_summary = export_plink_final(mt, plink_prefix, checkpoint_path)

        loci_summary['attempts'] = attempt
        export_summary['attempts'] = attempt
        export_summary['n_variants'] = n_passing
        export_summary['n_samples'] = n_samples

        summary['pool_size'] = n_pool
        summary['oversample_factor'] = oversample_factor
        summary['loci'] = loci_summary
        summary['plink'] = export_summary
        summary['status'] = 'success'

    except Exception as e:
        summary['status'] = 'failed'
        summary['error'] = str(e)
        summary['traceback'] = traceback.format_exc()
        logger.error(
            f"  {chrom} FAILED: {e}\n{traceback.format_exc()}"
        )
    finally:
        run_elapsed = time.time() - run_start
        summary['end_time'] = datetime.now().isoformat()
        summary['time_seconds'] = round(run_elapsed, 1)

        # Write summary JSON
        try:
            with hl.hadoop_open(summary_path, 'w') as f:
                json.dump(summary, f, indent=2, default=str)
            logger.info(f"  Summary written to {summary_path}")
        except Exception as write_err:
            logger.warning(f"  Failed to write summary: {write_err}")

        logger.info("")
        logger.info("=" * 60)
        logger.info(
            f"  {chrom}: {summary['status'].upper()} "
            f"in {_fmt_elapsed(run_elapsed)}"
        )
        logger.info("=" * 60)

        # Stop Hail / JVM
        try:
            hl.stop()
            logger.info("  Hail stopped")
        except Exception:
            pass

    if summary['status'] == 'failed':
        sys.exit(1)


if __name__ == "__main__":
    main()
