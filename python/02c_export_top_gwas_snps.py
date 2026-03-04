#!/usr/bin/env python3
"""Export top GWAS SNPs above a configurable p-value threshold.

Reads per-chromosome .ma files (produced by 02a_export_gwas_ma.py), filters to
SNPs with p < threshold, and writes a single consolidated TSV sorted by p-value.
No Hail or Spark required -- pure pandas.

Output columns:
    chrom, pos, ref, alt, SNP, A1, A2, freq, beta, se, p, N

Usage:
    python python/02c_export_top_gwas_snps.py --config config/config.json
    python python/02c_export_top_gwas_snps.py --config config/config.json --p-threshold 5e-8
    python python/02c_export_top_gwas_snps.py --config config/config.json --p-threshold 1e-5 --force
"""

import argparse
import logging
import os
import sys
import time
from datetime import datetime
from typing import List, Optional

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


logger = setup_logger("export_top_snps")


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


def parse_snp_id(snp_series: pd.Series) -> pd.DataFrame:
    """Parse chr:pos:ref:alt SNP IDs into component columns.

    Args:
        snp_series: pandas Series of SNP IDs in chr:pos:ref:alt format.

    Returns:
        DataFrame with columns: chrom, pos, ref, alt.
    """
    parts = snp_series.str.split(":", n=3, expand=True)
    parts.columns = ["chrom", "pos", "ref", "alt"]
    parts["pos"] = pd.to_numeric(parts["pos"], errors="coerce").astype("Int64")
    return parts


def load_config(config_path: str) -> dict:
    """Load JSON config file."""
    import json
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Config file not found: {config_path}")
    with open(config_path, "r") as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Export top GWAS SNPs above a p-value threshold from .ma files"
    )
    parser.add_argument(
        "--config", default="config/config.json",
        help="Path to config.json (default: config/config.json)"
    )
    parser.add_argument(
        "--p-threshold", type=float, default=None,
        help="P-value threshold (default: use gwas_p_threshold from config)"
    )
    parser.add_argument(
        "--ma-dir", default=None,
        help="Override directory containing per-chromosome .ma files"
    )
    parser.add_argument(
        "--output", default=None,
        help="Output TSV path (default: <locus_def_dir>/top_gwas_snps.tsv)"
    )
    parser.add_argument(
        "--force", action="store_true",
        help="Overwrite existing output file"
    )
    return parser.parse_args()


def main() -> None:
    """Read per-chromosome .ma files, filter by p-value, write consolidated TSV."""
    args = parse_args()
    run_start = time.time()

    logger.info("=" * 60)
    logger.info("02c: Export top GWAS SNPs")
    logger.info("=" * 60)
    logger.info(f"  Timestamp : {datetime.now().isoformat()}")
    logger.info(f"  Config    : {args.config}")

    config = load_config(args.config)

    # Resolve parameters
    p_threshold: float = args.p_threshold or config["params"].get("gwas_p_threshold", 5e-6)
    ma_dir: str = args.ma_dir or config["outputs"].get("gwas_ma_dir", "results/2-locus_definition/ma")
    locus_def_dir: str = config["outputs"].get("locus_def_dir", "results/2-locus_definition")
    output_path: str = args.output or os.path.join(locus_def_dir, "top_gwas_snps.tsv")

    logger.info(f"  P threshold: {p_threshold:.1e}")
    logger.info(f"  .ma dir    : {ma_dir}")
    logger.info(f"  Output     : {output_path}")
    logger.info(f"  Force      : {args.force}")
    logger.info("")

    # Check for existing output
    if os.path.exists(output_path) and not args.force:
        logger.error(f"Output already exists: {output_path}. Use --force to overwrite.")
        sys.exit(1)

    # Validate .ma directory
    if not os.path.isdir(ma_dir):
        logger.error(f".ma directory not found: {ma_dir}")
        logger.error("  Run 02a_export_gwas_ma first.")
        sys.exit(1)

    # Collect per-chromosome results
    chroms: List[str] = [f"chr{i}" for i in range(1, 23)]
    dfs: List[pd.DataFrame] = []
    total_variants: int = 0
    missing_chroms: List[str] = []

    for chrom in chroms:
        ma_path = os.path.join(ma_dir, f"{chrom}.ma")
        if not os.path.exists(ma_path):
            logger.warning(f"  [{chrom}] .ma file not found: {ma_path}")
            missing_chroms.append(chrom)
            continue

        df = pd.read_csv(ma_path, sep="\t")
        n_total = len(df)
        total_variants += n_total

        # Filter by p-value threshold
        df = df[df["p"] < p_threshold].copy()
        n_sig = len(df)

        if n_sig > 0:
            logger.info(f"  [{chrom}] {n_sig:,} / {n_total:,} SNPs pass p < {p_threshold:.1e}")
            dfs.append(df)
        else:
            logger.info(f"  [{chrom}] 0 / {n_total:,} SNPs pass threshold")

    if missing_chroms:
        logger.warning(f"  Missing .ma files for: {', '.join(missing_chroms)}")

    if not dfs:
        logger.warning("No SNPs found below the threshold across all chromosomes.")
        logger.warning("Writing empty output file with header only.")
        header_cols = [
            "chrom", "pos", "ref", "alt", "SNP", "A1", "A2",
            "freq", "beta", "se", "p", "N",
        ]
        pd.DataFrame(columns=header_cols).to_csv(output_path, sep="\t", index=False)
        logger.info(f"  Empty output written to {output_path}")
        return

    # Concatenate all chromosomes
    combined = pd.concat(dfs, ignore_index=True)

    # Parse SNP ID into component columns
    snp_parts = parse_snp_id(combined["SNP"])
    combined = pd.concat([snp_parts, combined], axis=1)

    # Rename .ma columns for clarity
    combined = combined.rename(columns={"b": "beta"})

    # Sort by p-value ascending
    combined = combined.sort_values("p", ascending=True).reset_index(drop=True)

    # Select and order output columns
    out_cols = ["chrom", "pos", "ref", "alt", "SNP", "A1", "A2", "freq", "beta", "se", "p", "N"]
    combined = combined[out_cols]

    # Write output
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    combined.to_csv(output_path, sep="\t", index=False)

    elapsed = time.time() - run_start
    logger.info("")
    logger.info("=" * 60)
    logger.info(f"  Total variants scanned : {total_variants:,}")
    logger.info(f"  SNPs at p < {p_threshold:.1e}   : {len(combined):,}")
    logger.info(f"  Chromosomes with hits  : {len(dfs)}")
    logger.info(f"  Output written to      : {output_path}")
    logger.info(f"  Completed in           : {_fmt_elapsed(elapsed)}")
    logger.info("=" * 60)

    # Print top 20 hits to log for quick inspection
    logger.info("")
    logger.info("Top 20 SNPs by p-value:")
    logger.info("-" * 100)
    top = combined.head(20)
    for _, row in top.iterrows():
        logger.info(
            f"  {row['chrom']}:{row['pos']}\t"
            f"p={row['p']:.2e}\t"
            f"beta={row['beta']:.4f}\t"
            f"se={row['se']:.4f}\t"
            f"freq={row['freq']:.4f}\t"
            f"{row['SNP']}"
        )
    logger.info("-" * 100)


if __name__ == "__main__":
    main()
