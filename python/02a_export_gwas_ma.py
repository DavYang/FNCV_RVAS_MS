#!/usr/bin/env python3
"""Export per-chromosome GCTA-COJO .ma summary statistics from AllxAll GWAS HT.

GCTA requires full-chromosome .ma files (all SNPs, not just per-window) to
correctly estimate phenotypic variance. This script reads the AllxAll GWAS
Hail Table, resolves field names, and writes one .ma file per autosome.

.ma format (tab-separated, with header):
    SNP  A1  A2  freq  b  se  p  N

    SNP  = variant ID (chrN:pos:ref:alt, matching Hail export_plink .bim)
    A1   = effect allele (alt)
    A2   = other allele (ref)
    freq = frequency of A1
    b    = effect size (log OR for case-control)
    se   = standard error
    p    = p-value
    N    = sample size

Usage:
    python python/02a_export_gwas_ma.py --config config/config.json
    python python/02a_export_gwas_ma.py --config config/config.json --force
    python python/02a_export_gwas_ma.py --config config/config.json --chrom chr21
"""

import argparse
import json
import logging
import os
import sys
import time
from datetime import datetime
from typing import Dict, List

import hail as hl
import pandas as pd

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

def setup_logger(name: str) -> logging.Logger:
    """Set up a logger writing to stdout with timestamps."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)],
    )
    return logging.getLogger(name)


logger = setup_logger("export_gwas_ma")


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

def load_config(config_path: str) -> dict:
    """Load JSON config file."""
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Config file not found: {config_path}")
    with open(config_path, "r") as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_FIELD_CANDIDATES: Dict[str, List[str]] = {
    'beta':   ['beta', 'effect_size', 'b'],
    'se':     ['standard_error', 'se', 'stderr'],
    'pval':   ['p_value', 'p_value_EUR', 'pval', 'p'],
    'af':     ['AF', 'AF_EUR', 'allele_frequency', 'freq'],
    'n':      ['n_complete_samples', 'N', 'n_samples', 'n'],
}

CHR_SIZES: Dict[str, int] = {
    'chr1': 248956422,  'chr2': 242193529,  'chr3': 198295559,
    'chr4': 190214555,  'chr5': 181538259,  'chr6': 170805979,
    'chr7': 159345973,  'chr8': 145138636,  'chr9': 138394717,
    'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
    'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189,
    'chr16': 90338345,  'chr17': 83257441,  'chr18': 80373285,
    'chr19': 58617616,  'chr20': 64444167,  'chr21': 46709983,
    'chr22': 50818468,
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _fmt_elapsed(seconds: float) -> str:
    """Format elapsed seconds as human-readable string."""
    if seconds < 60:
        return f"{seconds:.1f}s"
    m, s = divmod(int(seconds), 60)
    if m < 60:
        return f"{m}m {s:02d}s"
    h, m = divmod(m, 60)
    return f"{h}h {m:02d}m {s:02d}s"


def _resolve_field(ht: hl.Table, canonical: str, candidates: List[str]) -> str:
    """Return the first matching field name from candidates present in ht.row."""
    row_fields = set(ht.row)
    for name in candidates:
        if name in row_fields:
            logger.info(f"  Field '{canonical}' resolved to '{name}'")
            return name
    raise ValueError(
        f"Could not find field '{canonical}' in GWAS HT. "
        f"Tried: {candidates}. Available fields: {sorted(row_fields)}"
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Export per-chromosome GCTA-COJO .ma files from GWAS HT"
    )
    parser.add_argument(
        '--config', default='config/config.json',
        help='Path to config.json (default: config/config.json)'
    )
    parser.add_argument(
        '--output-dir', default=None,
        help='Override output directory for .ma files'
    )
    parser.add_argument(
        '--chrom', default=None,
        help='Process a single chromosome (e.g. chr21). Default: all autosomes.'
    )
    parser.add_argument(
        '--force', action='store_true',
        help='Overwrite existing .ma files'
    )
    return parser.parse_args()


def main() -> None:
    """Export per-chromosome .ma files from GWAS HT."""
    args = parse_args()

    run_start = time.time()
    logger.info("=" * 60)
    logger.info("02a: Export per-chromosome GWAS .ma files for GCTA-COJO")
    logger.info("=" * 60)
    logger.info(f"  Timestamp : {datetime.now().isoformat()}")
    logger.info(f"  Config    : {args.config}")

    config = load_config(args.config)

    gwas_path = config['inputs']['phenotype_gwas']
    output_dir = args.output_dir or os.path.join(
        config['outputs'].get('locus_def_dir', 'results/2-locus_definition'),
        'ma'
    )
    os.makedirs(output_dir, exist_ok=True)

    logger.info(f"  GWAS HT   : {gwas_path}")
    logger.info(f"  Output dir: {output_dir}")
    logger.info(f"  Force     : {args.force}")
    logger.info(f"  Chrom     : {args.chrom or 'all autosomes'}")
    logger.info("")

    # Determine chromosomes to process
    if args.chrom:
        chrom = args.chrom if args.chrom.startswith('chr') else f"chr{args.chrom}"
        if chrom not in CHR_SIZES:
            logger.error(f"Unknown chromosome: {chrom}")
            sys.exit(1)
        chroms = [chrom]
    else:
        chroms = [f"chr{i}" for i in range(1, 23)]

    # Initialize Hail
    hail_log = '/tmp/hail_export_ma.log'
    logger.info(f"  Initializing Hail (log: {hail_log}) ...")
    hl.init(log=hail_log)
    hl.default_reference('GRCh38')
    logger.info("  Hail initialized")
    logger.info("")

    try:
        # Load GWAS HT
        logger.info(f"Loading GWAS HT from {gwas_path} ...")
        t0 = time.time()
        gwas_ht = hl.read_table(gwas_path)
        logger.info(f"  GWAS HT loaded ({_fmt_elapsed(time.time() - t0)})")

        # Resolve field names
        logger.info("Resolving GWAS field names ...")
        field_map = {
            canonical: _resolve_field(gwas_ht, canonical, candidates)
            for canonical, candidates in _FIELD_CANDIDATES.items()
        }
        logger.info("")

        # Process each chromosome
        n_exported = 0
        n_skipped = 0

        for chrom in chroms:
            ma_path = os.path.join(output_dir, f"{chrom}.ma")

            # Resume check
            if os.path.exists(ma_path) and not args.force:
                n_lines = sum(1 for _ in open(ma_path)) - 1
                logger.info(f"[{chrom}] SKIP (exists) -- {n_lines:,} variants in {ma_path}")
                n_skipped += 1
                continue

            logger.info(f"[{chrom}] Exporting .ma file ...")
            t0 = time.time()

            # Filter to chromosome
            chr_size = CHR_SIZES[chrom]
            interval = hl.parse_locus_interval(
                f"{chrom}:1-{chr_size}", reference_genome='GRCh38'
            )
            chr_ht = hl.filter_intervals(gwas_ht, [interval])

            # Annotate with .ma columns
            chr_ht = chr_ht.annotate(
                SNP=(
                    hl.str(chr_ht.locus.contig) + ':' +
                    hl.str(chr_ht.locus.position) + ':' +
                    chr_ht.alleles[0] + ':' + chr_ht.alleles[1]
                ),
                A1=chr_ht.alleles[1],
                A2=chr_ht.alleles[0],
                freq=hl.float64(chr_ht[field_map['af']]),
                b=hl.float64(chr_ht[field_map['beta']]),
                se=hl.float64(chr_ht[field_map['se']]),
                p=hl.float64(chr_ht[field_map['pval']]),
                N=hl.int32(chr_ht[field_map['n']]),
            )

            # Select only .ma columns and export to pandas
            df = chr_ht.select('SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N').to_pandas()

            if df.empty:
                logger.warning(f"[{chrom}] No variants found -- skipping")
                continue

            # Filter out rows with missing values
            n_before = len(df)
            df = df.dropna(subset=['SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N'])
            n_dropped = n_before - len(df)
            if n_dropped > 0:
                logger.info(f"[{chrom}] Dropped {n_dropped:,} rows with missing values")

            # Write .ma file
            df.to_csv(ma_path, sep='\t', index=False)
            elapsed = time.time() - t0
            logger.info(
                f"[{chrom}] Exported {len(df):,} variants to {ma_path} "
                f"({_fmt_elapsed(elapsed)})"
            )
            n_exported += 1

        logger.info("")
        logger.info("=" * 60)
        logger.info(f"  Exported: {n_exported}  Skipped: {n_skipped}")
        logger.info("=" * 60)

    finally:
        try:
            hl.stop()
            logger.info("  Hail stopped")
        except Exception:
            pass

    elapsed = time.time() - run_start
    logger.info(f"02a complete in {_fmt_elapsed(elapsed)}")


if __name__ == "__main__":
    main()
