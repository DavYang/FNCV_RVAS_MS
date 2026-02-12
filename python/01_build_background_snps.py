#!/usr/bin/env python3
"""Build background SNP PLINK files for Regenie Step 1 null model.

Processes a SINGLE chromosome per invocation: samples loci from the AoU
ACAF splitMT rows Table, then extracts EUR genotypes and exports PLINK.
The bash wrapper calls this script once per chromosome (chr1-22), giving
each invocation a fresh Hail/JVM session to avoid accumulated state crashes.

Usage:
    python 01_build_background_snps.py --chrom chr21 --target 11700 \
        --output-dir gs://bucket/results/FNCV_RVAS_MS/500K_background_snps
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

ACAF_VARIANTS_PER_BP = 0.03826
SAMPLING_OVERSHOOT = 1.02


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
# Step 1: Sample loci for this chromosome (rows Table only, no entry data)
# ---------------------------------------------------------------------------
def sample_loci(
    mt_path: str,
    chrom: str,
    target: int,
    output_ht_path: str,
    config: dict,
) -> dict:
    """Sample loci on a single chromosome using Bernoulli sampling.

    Reads only the rows Table of the ACAF splitMT (no entry data touched).
    Applies HLA exclusion on chr6 if configured.

    Args:
        mt_path: GCS path to ACAF splitMT.
        chrom: Chromosome name, e.g. 'chr21'.
        target: Target number of SNPs to sample on this chromosome.
        output_ht_path: GCS path to write the sampled loci Table.
        config: Loaded config dict.

    Returns:
        Dict with 'n_sampled', 'target', 'time_seconds'.
    """
    step_start = time.time()
    seed = config['sampling']['random_seed']
    avoid_hla = config['sampling'].get('avoid_hla', False)
    hla_region = config['sampling'].get('hla_region', {})

    logger.info(f"  Step 1: Sampling loci (target={target:,})")
    _log_memory()

    logger.info(f"  Reading rows Table from {mt_path} ...")
    mt = hl.read_matrix_table(mt_path)
    ht_rows = mt.rows().select()

    chr_size = CHR_SIZES[chrom]
    interval = hl.parse_locus_interval(
        f"{chrom}:1-{chr_size}", reference_genome='GRCh38'
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
        logger.info(f"  Excluded HLA region {hla_chrom}:{hla_start}-{hla_end}")

    estimated_pool = int(chr_size * ACAF_VARIANTS_PER_BP)
    fraction = min(1.0, (target / estimated_pool) * SAMPLING_OVERSHOOT)

    autosomes = [f"chr{i}" for i in range(1, 23)]
    chr_seed = seed + autosomes.index(chrom)
    ht_sampled = ht_chr.filter(
        hl.rand_unif(0.0, 1.0, seed=chr_seed) < fraction
    )

    logger.info(
        f"  Sampling fraction={fraction:.6f}, est_pool={estimated_pool:,}, "
        f"seed={chr_seed}"
    )
    ht_sampled.write(output_ht_path, overwrite=True)

    ht_sampled = hl.read_table(output_ht_path)
    n_sampled = ht_sampled.count()

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
# Step 2: Extract genotypes at sampled loci, export PLINK
# ---------------------------------------------------------------------------
def export_plink(
    mt_path: str,
    chrom: str,
    loci_ht_path: str,
    plink_prefix: str,
    eur_samples_ht: hl.Table,
    tmp_dir: str,
) -> dict:
    """Filter MT to sampled loci + EUR samples, checkpoint, export PLINK.

    Args:
        mt_path: GCS path to ACAF splitMT.
        chrom: Chromosome name, e.g. 'chr21'.
        loci_ht_path: GCS path to sampled loci Table for this chromosome.
        plink_prefix: Output prefix for PLINK files (.bed/.bim/.fam).
        eur_samples_ht: Hail Table of EUR sample IDs keyed by 's'.
        tmp_dir: GCS path for checkpoint files.

    Returns:
        Dict with 'n_variants', 'n_samples', 'time_seconds'.
    """
    step_start = time.time()
    checkpoint_path = f"{tmp_dir}/{chrom}_checkpoint.mt"

    logger.info(f"  Step 2: Extracting genotypes and exporting PLINK")
    _log_memory()

    logger.info(f"  2a: Filtering MT by interval ...")
    t0 = time.time()
    chr_size = CHR_SIZES[chrom]
    chr_interval = hl.parse_locus_interval(
        f"{chrom}:1-{chr_size}", reference_genome='GRCh38'
    )
    ht_loci = hl.read_table(loci_ht_path)
    mt = hl.read_matrix_table(mt_path)
    mt = hl.filter_intervals(mt, [chr_interval])
    logger.info(f"  2a done ({_fmt_elapsed(time.time() - t0)})")

    logger.info(f"  2b: Joining rows (sampled loci) and cols (EUR) ...")
    t0 = time.time()
    mt = mt.semi_join_rows(ht_loci)
    mt = mt.semi_join_cols(eur_samples_ht)
    mt = mt.select_entries('GT')
    mt = mt.select_rows()
    mt = mt.select_cols()
    logger.info(f"  2b done ({_fmt_elapsed(time.time() - t0)})")

    logger.info(f"  2c: Checkpointing filtered MT (breaks DAG) ...")
    t0 = time.time()
    mt = mt.checkpoint(checkpoint_path, overwrite=True)
    n_variants = mt.count_rows()
    n_samples = mt.count_cols()
    logger.info(
        f"  2c done: {n_variants:,} variants x {n_samples:,} samples "
        f"({_fmt_elapsed(time.time() - t0)})"
    )

    logger.info(f"  2d: Exporting PLINK to {plink_prefix} ...")
    t0 = time.time()
    hl.export_plink(mt, plink_prefix, fam_id=mt.s, ind_id=mt.s)
    logger.info(f"  2d done ({_fmt_elapsed(time.time() - t0)})")

    logger.info(f"  2e: Cleaning up checkpoint ...")
    _cleanup_gcs_path(checkpoint_path)

    elapsed = time.time() - step_start
    logger.info(
        f"  Step 2 DONE: PLINK exported in {_fmt_elapsed(elapsed)}"
    )
    _log_memory()

    return {
        'n_variants': n_variants,
        'n_samples': n_samples,
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
    return parser.parse_args()


def main() -> None:
    """Run the single-chromosome pipeline."""
    args = parse_args()
    chrom = args.chrom
    target = args.target
    output_dir = args.output_dir
    config_path = args.config

    run_start = time.time()

    logger.info("=" * 60)
    logger.info(f"SINGLE-CHROMOSOME PIPELINE: {chrom}")
    logger.info("=" * 60)
    logger.info(f"  Timestamp : {datetime.now().isoformat()}")
    logger.info(f"  PID       : {os.getpid()}")
    logger.info(f"  Chrom     : {chrom}")
    logger.info(f"  Target    : {target:,}")
    logger.info(f"  Output dir: {output_dir}")
    logger.info(f"  Config    : {config_path}")
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

    # Resume: skip if summary.json already exists for this chromosome
    if hfs.exists(summary_path):
        logger.info(f"  summary.json already exists at {summary_path}")
        logger.info(f"  Skipping {chrom} (already completed)")
        return

    logger.info(f"  Source MT : {mt_path}")
    logger.info("")

    # Initialize Hail (fresh session for this chromosome)
    hail_log = f"/tmp/hail_{chrom}.log"
    logger.info(f"  Initializing Hail (log: {hail_log}) ...")
    hl.init(log=hail_log)
    hl.default_reference('GRCh38')
    logger.info("  Hail initialized")

    summary = {
        'chrom': chrom,
        'target': target,
        'status': 'running',
        'start_time': datetime.now().isoformat(),
    }

    try:
        # --- EUR samples (shared, written once, reused) ---
        logger.info("  Loading EUR samples Table ...")
        eur_samples_ht = _get_eur_samples_ht(config, shared_dir)

        # --- Step 1: Sample loci ---
        loci_ht_path = f"{chrom_dir}/sampled_loci.ht"
        if (hfs.exists(f"{loci_ht_path}/_SUCCESS")
                or hfs.exists(f"{loci_ht_path}/metadata.json.gz")):
            logger.info(f"  Sampled loci Table exists, skipping Step 1")
            ht_loci = hl.read_table(loci_ht_path)
            n_sampled = ht_loci.count()
            loci_summary = {
                'n_sampled': n_sampled,
                'target': target,
                'skipped': True,
            }
        else:
            loci_summary = sample_loci(
                mt_path, chrom, target, loci_ht_path, config
            )

        summary['loci'] = loci_summary

        # --- Step 2: Export PLINK ---
        plink_prefix = f"{chrom_dir}/{chrom}_background"
        bed_path = f"{plink_prefix}.bed"

        if hfs.exists(bed_path):
            logger.info(f"  PLINK .bed exists at {bed_path}, skipping Step 2")
            plink_summary = {'skipped': True, 'plink_prefix': plink_prefix}
        else:
            plink_summary = export_plink(
                mt_path, chrom, loci_ht_path,
                plink_prefix, eur_samples_ht, tmp_dir,
            )

        summary['plink'] = plink_summary
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
