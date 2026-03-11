#!/usr/bin/env python3
"""Export per-chromosome summary statistics (.ma format) from the AllxAll GWAS HT.

Reads the AllxAll EUR MS GWAS Hail Table, resolves field names, and writes one
.ma file per autosome. These files are consumed by 02c_export_top_gwas_snps.py
which filters to SNPs at p < gwas_p_threshold (default 5e-6) and produces
top_gwas_snps.tsv for cross-referencing against published EUR MS GWAS Catalog
loci in the windowed locus definition step (02e_define_loci_catalog.py).

.ma format (tab-separated, with header):
    SNP  A1  A2  freq  b  se  p  N

    SNP  = variant ID (chrN:pos:ref:alt)
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
from typing import Dict, List, Optional

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
    'beta':   ['BETA', 'beta', 'effect_size', 'b'],
    'se':     ['SE', 'standard_error', 'se', 'stderr'],
    'pval':   ['Pvalue', 'p_value', 'p_value_EUR', 'pval', 'p'],
    'af':     ['AF_Allele2', 'AF', 'AF_EUR', 'allele_frequency', 'freq'],
    'n':      ['n_complete_samples', 'N', 'n_samples', 'n'],
    'ac':     ['AC_Allele2', 'AC', 'allele_count'],
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


def _resolve_field(
    ht: hl.Table, canonical: str, candidates: List[str], required: bool = True,
) -> Optional[str]:
    """Return the first matching field name from candidates present in ht.row.

    Args:
        ht: Hail Table to search.
        canonical: Canonical field name for logging.
        candidates: Ordered list of candidate field names to try.
        required: If True, raise ValueError when no match found.

    Returns:
        Matched field name, or None if not required and no match found.
    """
    row_fields = set(ht.row)
    for name in candidates:
        if name in row_fields:
            logger.info(f"  Field '{canonical}' resolved to '{name}'")
            return name
    if required:
        raise ValueError(
            f"Could not find field '{canonical}' in GWAS HT. "
            f"Tried: {candidates}. Available fields: {sorted(row_fields)}"
        )
    logger.info(f"  Field '{canonical}' not found (optional) -- tried: {candidates}")
    return None


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
    logger.info("02a: Export per-chromosome GWAS .ma files")
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

        # Resolve field names (n and ac are optional; N derived from AC/AF if needed)
        logger.info("Resolving GWAS field names ...")
        field_map: Dict[str, Optional[str]] = {}
        for canonical, candidates in _FIELD_CANDIDATES.items():
            required = canonical not in ('n', 'ac')
            field_map[canonical] = _resolve_field(
                gwas_ht, canonical, candidates, required=required,
            )

        # Determine N strategy
        if field_map['n'] is not None:
            n_strategy = 'direct'
            logger.info(f"  N strategy: direct from field '{field_map['n']}'")
        elif field_map['ac'] is not None and field_map['af'] is not None:
            n_strategy = 'derived'
            logger.info(
                f"  N strategy: derived from {field_map['ac']} / (2 * {field_map['af']})"
            )
        else:
            raise ValueError(
                "Cannot determine sample size (N). GWAS HT has neither a direct N "
                "field nor AC + AF fields to derive it. "
                f"Available fields: {sorted(set(gwas_ht.row))}"
            )
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

            # Derive N expression based on strategy
            if n_strategy == 'direct':
                n_expr = hl.int32(chr_ht[field_map['n']])
            else:
                # N = AC_Allele2 / (2 * AF_Allele2) -- diploid sample count
                af_val = hl.float64(chr_ht[field_map['af']])
                ac_val = hl.float64(chr_ht[field_map['ac']])
                n_expr = hl.or_missing(
                    af_val > 0,
                    hl.int32(hl.floor(ac_val / (2.0 * af_val) + 0.5)),
                )

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
                N=n_expr,
            )

            # Drop keys so locus/alleles don't leak into output, then select .ma cols
            ma_cols = ['SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N']
            df = chr_ht.key_by().select(*ma_cols).to_pandas()
            df = df[ma_cols]

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
    logger.info("Next step: bash bash/02c_export_top_gwas_snps.sh")


if __name__ == "__main__":
    main()
