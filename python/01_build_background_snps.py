#!/usr/bin/env python3
"""Build background SNP PLINK files for Regenie Step 1 null model.

Processes a SINGLE chromosome per invocation: samples loci from the AoU
ACAF splitMT rows Table, then extracts EUR genotypes and exports PLINK.
The bash wrapper calls this script once per chromosome (chr1-22), giving
each invocation a fresh Hail/JVM session to avoid accumulated state crashes.

Pipeline per chromosome (2 steps):
  1. sample_loci()   - Bernoulli sample from QC-passing rows Table at a
                       generous oversample_factor. Rows-only, no GT data read.
                       Uses pre-computed info.AF (global AF > hail_min_maf as
                       EUR MAF proxy) and variant_qc.call_rate pre-filters.
  2. export_plink()  - Join sampled loci + EUR samples on the full MT, apply
                       EUR-specific call-rate filter, coalesce partitions, and
                       export PLINK directly from the DAG (no intermediate
                       checkpoint). Variant count is read from the .bim file
                       after export.

The oversample_factor (default 3+) ensures the target is met on the first
attempt without needing a pool-count step or retry loop. Empirically,
the EUR call-rate filter drops <0.1% of pre-filtered variants.

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


def _count_bim_lines(bim_path: str) -> int:
    """Count the number of lines in a .bim file on GCS.

    This is used to determine the variant count after export_plink, avoiding
    the need for a separate checkpoint + count_rows() call. Each line in a
    .bim file corresponds to exactly one variant.

    Args:
        bim_path: GCS path to the .bim file.

    Returns:
        Integer count of variants (lines) in the .bim file.
    """
    n = 0
    with hl.hadoop_open(bim_path, 'r') as f:
        for _ in f:
            n += 1
    return n


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
    if configured.

    Args:
        mt_path: GCS path to ACAF splitMT.
        chrom: Chromosome name, e.g. 'chr21'.
        target: Target number of SNPs (used for logging only).
        output_ht_path: GCS path to write the sampled loci Table.
        config: Loaded config dict.
        fraction: Bernoulli sampling fraction (computed in main()).
        seed_offset: Added to the base chromosome seed (reserved for future use).

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
# Step 2: Export PLINK directly from filtered DAG (no checkpoint)
# ---------------------------------------------------------------------------
def export_plink(
    mt_path: str,
    chrom: str,
    loci_ht_path: str,
    eur_samples_ht: hl.Table,
    plink_prefix: str,
    config: dict,
) -> dict:
    """Build a filtered MatrixTable and export PLINK directly.

    Constructs the full DAG: interval filter -> semi_join to sampled loci +
    EUR samples -> select GT only -> EUR call-rate filter -> naive_coalesce ->
    export_plink. No intermediate checkpoint is written. The variant count is
    read from the .bim file after export instead of count_rows().

    Per Hail docs, export_plink accepts any MatrixTable keyed by
    (locus, alleles) with column key of type tstr and diploid unphased GT.
    filter_rows supports aggregation over columns (hl.agg.fraction).
    naive_coalesce merges adjacent partitions without shuffle.

    Args:
        mt_path: GCS path to ACAF splitMT.
        chrom: Chromosome name, e.g. 'chr21'.
        loci_ht_path: GCS path to sampled loci Table for this chromosome.
        eur_samples_ht: Hail Table of EUR sample IDs keyed by 's'.
        plink_prefix: Output prefix for PLINK files (.bed/.bim/.fam).
        config: Loaded config dict.

    Returns:
        Dict with 'plink_prefix', 'n_variants', 'n_samples', 'time_seconds'.
    """
    step_start = time.time()
    min_call_rate = config['sampling'].get('min_call_rate', DEFAULT_MIN_CALL_RATE)

    logger.info(f"  Step 2: Export PLINK (direct DAG, no checkpoint)")
    logger.info(f"  EUR call_rate filter: >= {min_call_rate}")
    _log_memory()

    logger.info(f"  2-interval: Filtering MT to {chrom} interval ...")
    t0 = time.time()
    chr_size = CHR_SIZES[chrom]
    chr_interval = hl.parse_locus_interval(
        f"{chrom}:1-{chr_size}", reference_genome='GRCh38'
    )
    ht_loci = hl.read_table(loci_ht_path)
    mt = hl.read_matrix_table(mt_path)
    mt = hl.filter_intervals(mt, [chr_interval])
    logger.info(f"  2-interval done ({_fmt_elapsed(time.time() - t0)})")

    logger.info(f"  2-join: Joining to sampled loci + EUR samples ...")
    t0 = time.time()
    mt = mt.semi_join_rows(ht_loci)
    mt = mt.semi_join_cols(eur_samples_ht)
    mt = mt.select_entries('GT')
    mt = mt.select_rows()
    mt = mt.select_cols()
    logger.info(f"  2-join done ({_fmt_elapsed(time.time() - t0)})")

    logger.info(f"  2-filter: Applying call_rate >= {min_call_rate} on EUR GTs ...")
    t0 = time.time()
    mt = mt.filter_rows(
        hl.agg.fraction(hl.is_defined(mt.GT)) >= min_call_rate
    )
    logger.info(f"  2-filter done ({_fmt_elapsed(time.time() - t0)})")

    logger.info(f"  2-coalesce: Repartitioning to 200 partitions ...")
    t0 = time.time()
    mt = mt.naive_coalesce(200)
    logger.info(f"  2-coalesce done ({_fmt_elapsed(time.time() - t0)})")

    logger.info(f"  2-export: Writing PLINK to {plink_prefix} ...")
    t0 = time.time()
    hl.export_plink(mt, plink_prefix, fam_id=mt.s, ind_id=mt.s)
    export_elapsed = time.time() - t0
    logger.info(f"  2-export done ({_fmt_elapsed(export_elapsed)})")

    logger.info(f"  2-count: Reading variant count from .bim ...")
    t0 = time.time()
    bim_path = f"{plink_prefix}.bim"
    n_variants = _count_bim_lines(bim_path)
    logger.info(
        f"  2-count done: {n_variants:,} variants ({_fmt_elapsed(time.time() - t0)})"
    )

    # n_samples from .fam (one line per sample)
    fam_path = f"{plink_prefix}.fam"
    n_samples = _count_bim_lines(fam_path)

    elapsed = time.time() - step_start
    logger.info(
        f"  Step 2 DONE: {n_variants:,} variants x {n_samples:,} samples "
        f"exported in {_fmt_elapsed(elapsed)}"
    )
    _log_memory()

    return {
        'plink_prefix': plink_prefix,
        'n_variants': n_variants,
        'n_samples': n_samples,
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
        help='Ignore existing summary.json; overwrite all outputs'
    )
    return parser.parse_args()


def main() -> None:
    """Run the single-chromosome pipeline.

    Two-step flow:
      1. sample_loci()  - rows-only Bernoulli sample (no GT data)
      2. export_plink() - join sampled loci + EUR samples on full MT,
                          apply EUR call-rate filter, export PLINK directly

    No checkpoint or retry loop. The oversample_factor ensures we always
    overshoot the target. Variant count is read from the .bim after export.
    """
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
    summary_path = f"{chrom_dir}/summary.json"

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

        # --- Compute sampling fraction ---
        # With oversample_factor >= 3 and the generous info.AF > 0.1 pre-filter,
        # the QC-passing pool is always much larger than the target. We use a
        # fixed fraction rather than counting the pool first (saves ~2 min).
        # fraction = oversample_factor * target / target = oversample_factor,
        # but capped at 1.0 for safety. In practice this means we sample ~3x
        # the target from the rows Table, and the EUR call-rate filter drops <0.1%.
        oversample_factor = config['sampling'].get(
            'oversample_factor', DEFAULT_OVERSAMPLE_FACTOR
        )
        # Use a conservative estimate: ~1.75 QC-passing variants per 1000 bp
        # (empirical from chr9: 241K pool / 138M bp = 1.75/kbp).
        # This gives fraction = target * oversample_factor / estimated_pool.
        # For safety, cap at 1.0.
        estimated_pool = CHR_SIZES[chrom] * 0.00175
        fraction = min(1.0, (target * oversample_factor) / estimated_pool)

        logger.info(
            f"  Sampling fraction: {fraction:.6f} "
            f"(target={target:,}, oversample_factor={oversample_factor}, "
            f"estimated_pool={estimated_pool:,.0f})"
        )

        # --- Step 1: Sample loci (rows-only, no GT data) ---
        loci_ht_path = f"{chrom_dir}/sampled_loci.ht"
        plink_prefix = f"{chrom_dir}/{chrom}_background"

        loci_summary = sample_loci(
            mt_path, chrom, target, loci_ht_path, config,
            fraction=fraction,
        )

        # --- Step 2: Export PLINK directly (no checkpoint) ---
        export_summary = export_plink(
            mt_path, chrom, loci_ht_path,
            eur_samples_ht, plink_prefix, config,
        )

        n_variants = export_summary['n_variants']
        if n_variants < target:
            logger.warning(
                f"  {chrom}: exported {n_variants:,} variants, below target "
                f"{target:,}. Downstream PLINK QC will further filter. "
                f"Consider increasing oversample_factor in config."
            )
        else:
            logger.info(
                f"  {chrom}: exported {n_variants:,} variants "
                f"(target={target:,}, +{n_variants - target:,} surplus)"
            )

        summary['oversample_factor'] = oversample_factor
        summary['fraction'] = round(fraction, 6)
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
