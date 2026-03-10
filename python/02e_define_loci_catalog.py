#!/usr/bin/env python3
"""Define RVAS target loci by cross-referencing AoU top GWAS SNPs against EUR MS GWAS Catalog loci.

Loads the filtered GWAS Catalog EUR-only MS loci (GRCh38, produced by
02d_parse_gwas_catalog.py), positionally cross-references against AoU
top GWAS SNPs (p < gwas_p_threshold), and generates fixed-width ±250kb windows
around validated SNPs. Overlapping windows are merged, the MHC region is
excluded, and results are exported as target_loci.bed.

All input/output paths are read from config/config.json:
    outputs.locus_def_dir              -- base output directory
    outputs.loci_bed                   -- path for target_loci.bed
    params.gwas_catalog_hg38_path      -- local path to gwas_catalog_ms_eur_hg38.tsv
    params.gwas_catalog_gcs_path       -- GCS URI to download catalog from if not local
    params.catalog_validation_window_bp -- ±bp for positional matching (default 250000)
    params.locus_flank_bp              -- ±bp window around each validated SNP (default 250000)
    params.mhc_interval                -- MHC region to exclude (default chr6:25000000-35000000)

Matching logic:
    For each AoU top SNP (GRCh38 position), the SNP is "validated" if at least
    one EUR MS catalog SNP lies within ±catalog_validation_window_bp on the same
    chromosome. Unmatched SNPs are written to novel_snps.tsv and excluded from
    locus generation.

Output files:
    <locus_def_dir>/gwas_catalog_validated_snps.tsv
    <locus_def_dir>/novel_snps.tsv
    <locus_def_dir>/target_loci.bed

Usage (run from project root):
    python python/02e_define_loci_catalog.py --config config/config.json
    python python/02e_define_loci_catalog.py --config config/config.json --force
"""

import argparse
import json
import logging
import os
import sys
import time
from datetime import datetime
from typing import Dict, List, Optional, Tuple

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


logger = setup_logger("define_loci_catalog")


def _fmt_elapsed(seconds: float) -> str:
    """Format elapsed seconds as human-readable string."""
    if seconds < 60:
        return f"{seconds:.1f}s"
    m, s = divmod(int(seconds), 60)
    if m < 60:
        return f"{m}m {s:02d}s"
    h, m = divmod(m, 60)
    return f"{h}h {m:02d}m {s:02d}s"


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

def load_config(config_path: str) -> dict:
    """Load JSON config file.

    Args:
        config_path: Path to config.json.

    Returns:
        Parsed config dict.
    """
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Config file not found: {config_path}")
    with open(config_path, "r") as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# SNP ID parsing
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------

def load_top_snps(top_snps_path: str) -> pd.DataFrame:
    """Load AoU top GWAS SNPs TSV.

    Args:
        top_snps_path: Path to top_gwas_snps.tsv.

    Returns:
        DataFrame with columns: chrom, pos, ref, alt, SNP, p, beta, se, freq.
    """
    logger.info(f"Loading top GWAS SNPs: {top_snps_path}")
    df = pd.read_csv(top_snps_path, sep="\t", dtype={"chrom": str})
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce").astype("Int64")
    df["p"] = pd.to_numeric(df["p"], errors="coerce")
    df["chrom"] = df["chrom"].astype(str).str.strip()
    if not df["chrom"].str.startswith("chr").all():
        df["chrom"] = "chr" + df["chrom"].str.lstrip("chr")
    logger.info(f"  Loaded {len(df):,} top GWAS SNPs")
    return df


def load_catalog(catalog_path: str) -> pd.DataFrame:
    """Load the filtered GWAS Catalog EUR MS loci TSV (GRCh38).

    Args:
        catalog_path: Path to gwas_catalog_ms_eur_hg38.tsv.

    Returns:
        DataFrame with columns: chrom, pos_hg38, rsid, p_value, etc.
    """
    logger.info(f"Loading GWAS Catalog EUR MS loci: {catalog_path}")
    df = pd.read_csv(catalog_path, sep="\t", dtype={"chrom": str})
    df["pos_hg38"] = pd.to_numeric(df["pos_hg38"], errors="coerce").astype("Int64")
    df["p_value"] = pd.to_numeric(df["p_value"], errors="coerce")
    df["chrom"] = df["chrom"].astype(str).str.strip()
    if not df["chrom"].str.startswith("chr").all():
        df["chrom"] = "chr" + df["chrom"].str.lstrip("chr")
    logger.info(f"  Loaded {len(df):,} catalog EUR MS loci")
    return df


# ---------------------------------------------------------------------------
# Cross-reference: positional matching
# ---------------------------------------------------------------------------

def cross_reference(
    top_snps: pd.DataFrame,
    catalog: pd.DataFrame,
    window_bp: int,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Cross-reference AoU top SNPs against catalog loci by position.

    For each AoU top SNP, find the nearest catalog SNP on the same chromosome
    within ±window_bp. SNPs with a hit are validated; the rest are unmatched.

    Args:
        top_snps: AoU top GWAS SNPs with chrom, pos columns (GRCh38).
        catalog: GWAS Catalog EUR MS loci with chrom, pos_hg38 columns (GRCh38).
        window_bp: Half-window size in base pairs for positional matching.

    Returns:
        Tuple of (validated_df, unmatched_df). validated_df has extra columns:
        catalog_rsid, catalog_pos_hg38, catalog_p_value, catalog_pubmed_id,
        distance_bp.
    """
    logger.info(
        f"Cross-referencing {len(top_snps):,} AoU SNPs vs "
        f"{len(catalog):,} catalog loci (window ±{window_bp:,} bp)..."
    )

    validated_rows = []
    unmatched_rows = []

    # Group catalog by chromosome for fast lookup
    catalog_by_chrom: Dict[str, pd.DataFrame] = {}
    for chrom, grp in catalog.groupby("chrom"):
        catalog_by_chrom[str(chrom)] = grp.reset_index(drop=True)

    for _, snp in top_snps.iterrows():
        chrom = str(snp["chrom"])
        pos = int(snp["pos"])

        if chrom not in catalog_by_chrom:
            unmatched_rows.append(snp.to_dict())
            continue

        cat_chrom = catalog_by_chrom[chrom]
        in_window = cat_chrom[
            (cat_chrom["pos_hg38"] >= pos - window_bp) &
            (cat_chrom["pos_hg38"] <= pos + window_bp)
        ]

        if len(in_window) == 0:
            unmatched_rows.append(snp.to_dict())
            continue

        # Pick nearest catalog hit
        in_window = in_window.copy()
        in_window["_dist"] = (in_window["pos_hg38"] - pos).abs()
        nearest = in_window.loc[in_window["_dist"].idxmin()]

        row = snp.to_dict()
        row["catalog_rsid"] = nearest.get("rsid", None)
        row["catalog_pos_hg38"] = int(nearest["pos_hg38"])
        row["catalog_p_value"] = nearest.get("p_value", None)
        row["catalog_pubmed_id"] = nearest.get("pubmed_id", None)
        row["distance_bp"] = int(nearest["_dist"])
        validated_rows.append(row)

    validated_df = pd.DataFrame(validated_rows) if validated_rows else pd.DataFrame()
    unmatched_df = pd.DataFrame(unmatched_rows) if unmatched_rows else pd.DataFrame()

    logger.info(
        f"  Validated (catalog hit within ±{window_bp:,} bp): {len(validated_df):,}"
    )
    logger.info(f"  Unmatched (no catalog hit): {len(unmatched_df):,}")

    return validated_df, unmatched_df


# ---------------------------------------------------------------------------
# Window generation and merging
# ---------------------------------------------------------------------------

def generate_windows(validated: pd.DataFrame, flank_bp: int) -> pd.DataFrame:
    """Generate ±flank_bp windows around each validated AoU SNP position.

    Args:
        validated: Validated SNPs dataframe with chrom and pos columns (GRCh38).
        flank_bp: One-sided flank in base pairs.

    Returns:
        DataFrame with columns: chrom, start (0-based), end (1-based), snp_id, p.
    """
    windows = validated[["chrom", "pos", "SNP", "p"]].copy()
    windows["start"] = (windows["pos"] - flank_bp).clip(lower=1) - 1  # 0-based
    windows["end"] = windows["pos"] + flank_bp
    windows = windows.rename(columns={"SNP": "snp_id"})
    return windows[["chrom", "start", "end", "snp_id", "p"]].copy()


def merge_windows(windows: pd.DataFrame) -> pd.DataFrame:
    """Merge overlapping windows, assigning locus_id from the lowest-p SNP.

    Args:
        windows: DataFrame with chrom, start (0-based), end, snp_id, p.

    Returns:
        Merged DataFrame with chrom, start, end, lead_snp columns.
    """
    logger.info(f"Merging {len(windows):,} windows...")

    chrom_order = [f"chr{i}" for i in range(1, 23)]
    windows["chrom"] = pd.Categorical(windows["chrom"], categories=chrom_order, ordered=True)
    windows = windows.sort_values(["chrom", "start"]).reset_index(drop=True)

    merged_rows = []

    for chrom, grp in windows.groupby("chrom", observed=True):
        grp = grp.sort_values("start").reset_index(drop=True)
        cur_start = int(grp.loc[0, "start"])
        cur_end = int(grp.loc[0, "end"])
        cur_snps = [(float(grp.loc[0, "p"]), str(grp.loc[0, "snp_id"]))]

        for i in range(1, len(grp)):
            row_start = int(grp.loc[i, "start"])
            row_end = int(grp.loc[i, "end"])
            row_p = float(grp.loc[i, "p"])
            row_snp = str(grp.loc[i, "snp_id"])

            if row_start <= cur_end:
                cur_end = max(cur_end, row_end)
                cur_snps.append((row_p, row_snp))
            else:
                lead_snp = min(cur_snps, key=lambda x: x[0])[1]
                merged_rows.append({
                    "chrom": str(chrom),
                    "start": cur_start,
                    "end": cur_end,
                    "lead_snp": lead_snp,
                })
                cur_start = row_start
                cur_end = row_end
                cur_snps = [(row_p, row_snp)]

        lead_snp = min(cur_snps, key=lambda x: x[0])[1]
        merged_rows.append({
            "chrom": str(chrom),
            "start": cur_start,
            "end": cur_end,
            "lead_snp": lead_snp,
        })

    merged = pd.DataFrame(merged_rows)
    logger.info(f"  {len(windows):,} windows -> {len(merged):,} merged loci")
    return merged


# ---------------------------------------------------------------------------
# MHC exclusion
# ---------------------------------------------------------------------------

def exclude_mhc(
    merged: pd.DataFrame, mhc_interval: str
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Exclude windows overlapping the MHC region.

    Args:
        merged: Merged windows DataFrame with chrom, start, end.
        mhc_interval: MHC region string e.g. "chr6:25000000-35000000".

    Returns:
        Tuple of (clean_df, mhc_df) where mhc_df contains excluded windows.
    """
    chrom_part, range_part = mhc_interval.split(":")
    mhc_chrom = chrom_part.strip()
    mhc_start, mhc_end = [int(x) for x in range_part.split("-")]

    mhc_mask = (
        (merged["chrom"] == mhc_chrom) &
        (merged["start"] < mhc_end) &
        (merged["end"] > mhc_start)
    )

    mhc_windows = merged[mhc_mask].copy()
    clean = merged[~mhc_mask].copy()

    if len(mhc_windows) > 0:
        logger.info(
            f"Excluded {len(mhc_windows):,} window(s) overlapping MHC ({mhc_interval})"
        )
    else:
        logger.info("No windows overlap MHC region.")

    return clean, mhc_windows


# ---------------------------------------------------------------------------
# Assign locus IDs
# ---------------------------------------------------------------------------

def assign_locus_ids(merged: pd.DataFrame) -> pd.DataFrame:
    """Assign locus_id based on lead SNP chrom:pos.

    Args:
        merged: Merged windows with chrom, start, end, lead_snp columns.

    Returns:
        DataFrame with locus_id column added.
    """
    def _make_locus_id(row: pd.Series) -> str:
        parts = str(row["lead_snp"]).split(":")
        if len(parts) >= 2:
            return f"{parts[0]}_{parts[1]}"
        return f"{row['chrom']}_{row['start']}"

    merged = merged.copy()
    merged["locus_id"] = merged.apply(_make_locus_id, axis=1)
    return merged


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Define RVAS loci by cross-referencing AoU GWAS SNPs vs "
            "GWAS Catalog EUR MS loci. All paths resolved from config/config.json."
        )
    )
    parser.add_argument(
        "--config", default="config/config.json",
        help="Path to config.json (default: config/config.json)"
    )
    parser.add_argument(
        "--catalog", default=None,
        help="Override path to gwas_catalog_ms_eur_hg38.tsv (bypasses config)"
    )
    parser.add_argument(
        "--top-snps", default=None,
        help="Override path to top_gwas_snps.tsv (bypasses config default)"
    )
    parser.add_argument(
        "--force", action="store_true",
        help="Overwrite existing output files"
    )
    return parser.parse_args()


def main() -> None:
    """Cross-reference AoU top SNPs vs GWAS Catalog, generate RVAS loci BED."""
    args = parse_args()
    run_start = time.time()

    logger.info("=" * 60)
    logger.info("02e: Define Loci via GWAS Catalog Cross-Reference")
    logger.info("=" * 60)
    logger.info(f"  Timestamp : {datetime.now().isoformat()}")
    logger.info(f"  Config    : {args.config}")

    config = load_config(args.config)
    params = config.get("params", {})
    outputs = config.get("outputs", {})

    # Resolve paths from config
    locus_def_dir = outputs.get("locus_def_dir", "results/2-locus_definition")
    os.makedirs(locus_def_dir, exist_ok=True)

    top_snps_path = args.top_snps or os.path.join(locus_def_dir, "top_gwas_snps.tsv")

    catalog_path = args.catalog or params.get(
        "gwas_catalog_hg38_path",
        os.path.join(locus_def_dir, "gwas_catalog_ms_eur_hg38.tsv"),
    )

    validation_window = int(params.get("catalog_validation_window_bp", 250_000))
    flank_bp = int(params.get("locus_flank_bp", 250_000))
    mhc_interval = params.get("mhc_interval", "chr6:25000000-35000000")
    gcs_path = params.get("gwas_catalog_gcs_path", "")

    validated_out = os.path.join(locus_def_dir, "gwas_catalog_validated_snps.tsv")
    novel_out = os.path.join(locus_def_dir, "novel_snps.tsv")
    bed_out = outputs.get("loci_bed", os.path.join(locus_def_dir, "target_loci.bed"))

    logger.info(f"  Top SNPs  : {top_snps_path}")
    logger.info(f"  Catalog   : {catalog_path}")
    logger.info(f"  Val window: +/-{validation_window:,} bp")
    logger.info(f"  Flank     : +/-{flank_bp:,} bp")
    logger.info(f"  MHC       : {mhc_interval}")
    logger.info(f"  BED out   : {bed_out}")
    logger.info("")

    # Pre-flight
    if not os.path.exists(top_snps_path):
        logger.error(f"Top SNPs file not found: {top_snps_path}")
        logger.error("Run: bash bash/02c_export_top_gwas_snps.sh")
        sys.exit(1)

    if not os.path.exists(catalog_path):
        logger.error(f"GWAS Catalog file not found: {catalog_path}")
        if gcs_path:
            logger.error(
                f"Download from GCS with:\n"
                f"  gsutil cp {gcs_path} {catalog_path}"
            )
        else:
            logger.error(
                "Run bash/02d_parse_gwas_catalog.sh on HPC, then upload to GCS "
                "and set params.gwas_catalog_gcs_path in config.json."
            )
        sys.exit(1)

    if os.path.exists(bed_out) and not args.force:
        logger.info(f"Output already exists (use --force to overwrite): {bed_out}")
        sys.exit(0)

    # Load data
    top_snps = load_top_snps(top_snps_path)
    catalog = load_catalog(catalog_path)

    # Cross-reference
    validated, unmatched = cross_reference(top_snps, catalog, validation_window)

    if len(validated) == 0:
        logger.error(
            "No AoU top SNPs validated against catalog. "
            "Check that both files use the same chromosome naming (chr-prefix) "
            "and that positions are GRCh38."
        )
        sys.exit(1)

    # Write validated and novel SNP tables
    validated.to_csv(validated_out, sep="\t", index=False)
    logger.info(f"Written: {validated_out}  ({len(validated):,} validated SNPs)")

    unmatched.to_csv(novel_out, sep="\t", index=False)
    logger.info(f"Written: {novel_out}  ({len(unmatched):,} unmatched SNPs)")

    # Generate windows
    windows = generate_windows(validated, flank_bp)

    # Merge overlapping windows
    merged = merge_windows(windows)

    # Exclude MHC
    clean, mhc_excluded = exclude_mhc(merged, mhc_interval)

    if len(mhc_excluded) > 0:
        mhc_out = os.path.join(locus_def_dir, "mhc_excluded_windows.tsv")
        mhc_excluded.to_csv(mhc_out, sep="\t", index=False)
        logger.info(f"Written MHC-excluded windows: {mhc_out}")

    # Assign locus IDs
    clean = assign_locus_ids(clean)

    # Write BED (0-based, half-open: chrom, start, end, locus_id)
    chrom_order = [f"chr{i}" for i in range(1, 23)]
    bed = clean[["chrom", "start", "end", "locus_id"]].copy()
    bed["chrom"] = pd.Categorical(bed["chrom"], categories=chrom_order, ordered=True)
    bed = bed.sort_values(["chrom", "start"]).reset_index(drop=True)
    bed["chrom"] = bed["chrom"].astype(str)
    bed.to_csv(bed_out, sep="\t", index=False, header=False)
    logger.info(f"Written: {bed_out}  ({len(bed):,} loci)")

    # Summary
    logger.info("")
    logger.info("Summary:")
    logger.info(f"  AoU top SNPs loaded          : {len(top_snps):,}")
    logger.info(f"  Validated (catalog support)  : {len(validated):,}")
    logger.info(f"  Unmatched (novel)            : {len(unmatched):,}")
    logger.info(f"  Windows pre-merge            : {len(windows):,}")
    logger.info(f"  Loci after merge             : {len(merged):,}")
    logger.info(f"  MHC excluded                 : {len(mhc_excluded):,}")
    logger.info(f"  Final RVAS loci              : {len(bed):,}")
    logger.info("")
    logger.info(f"Total runtime: {_fmt_elapsed(time.time() - run_start)}")
    logger.info("Done.")


if __name__ == "__main__":
    main()
