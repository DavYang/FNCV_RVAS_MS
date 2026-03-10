#!/usr/bin/env python3
"""Parse the GWAS Catalog full associations TSV to extract EUR-only MS loci.

Filters to Multiple Sclerosis associations from European-only GWAS studies
(multi-ancestry studies excluded entirely), lifts positions from GRCh37 to
GRCh38 using pyliftover, and writes a clean TSV to the locus definition output
directory for use by 02e_define_loci_catalog.py on the AoU workbench.

All input/output paths are read from config/config.json:
    inputs.gwas_catalog_raw  -- full GWAS Catalog associations TSV (GRCh37)
    inputs.liftover_chain    -- hg19ToHg38.over.chain.gz (downloaded if absent)
    outputs.locus_def_dir    -- output directory
    params.gwas_catalog_p_threshold  -- p-value threshold (default 5e-8)
    params.gwas_catalog_gcs_path     -- GCS destination for upload

Filtering criteria:
    - DISEASE/TRAIT contains "multiple sclerosis" (case-insensitive)
    - P-VALUE < gwas_catalog_p_threshold
    - CHR_ID in 1-22 (autosomes only)
    - CHR_POS non-null
    - EUR-only study: INITIAL SAMPLE SIZE contains "European ancestry" AND
      does NOT mention African, Asian, Japanese, Hispanic, admixed, Latino,
      South Asian, East Asian, multi, etc. Excludes multi-ancestry GWAS entirely.

Output columns:
    chrom, pos_hg37, pos_hg38, rsid, trait, p_value, pubmed_id, sample_size

Output file:
    <locus_def_dir>/gwas_catalog_ms_eur_hg38.tsv

Usage (run from project root):
    python python/02d_parse_gwas_catalog.py --config config/config.json
    python python/02d_parse_gwas_catalog.py --config config/config.json --force
"""

import argparse
import json
import logging
import os
import sys
import time
from datetime import datetime
from typing import List

import pandas as pd


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

CHAIN_URL = (
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
)
CHUNK_SIZE = 50_000

# Presence of any of these keywords in INITIAL SAMPLE SIZE means the study
# is multi-ancestry and must be excluded entirely.
NON_EUR_KEYWORDS: List[str] = [
    "african",
    "asian",
    "japanese",
    "hispanic",
    "admixed",
    "latino",
    "latina",
    "south asian",
    "east asian",
    "multi",
    "chinese",
    "korean",
    "indigenous",
    "amerindian",
    "native american",
]

AUTOSOMES = {str(i) for i in range(1, 23)}


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


logger = setup_logger("parse_gwas_catalog")


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
# Chain file
# ---------------------------------------------------------------------------

def ensure_chain_file(chain_path: str) -> None:
    """Download hg19ToHg38 chain file from UCSC if not already present.

    Args:
        chain_path: Destination path for the chain file.
    """
    if os.path.exists(chain_path):
        logger.info(f"Chain file found: {chain_path}")
        return
    logger.info("Chain file not found. Downloading from UCSC...")
    import urllib.request
    chain_dir = os.path.dirname(chain_path)
    if chain_dir:
        os.makedirs(chain_dir, exist_ok=True)
    urllib.request.urlretrieve(CHAIN_URL, chain_path)
    logger.info(f"Downloaded chain file: {chain_path}")


# ---------------------------------------------------------------------------
# Filtering helpers
# ---------------------------------------------------------------------------

def is_eur_only(sample_size: str) -> bool:
    """Return True if INITIAL SAMPLE SIZE describes a EUR-only study.

    Requires "European ancestry" to be present and rejects any row whose
    sample size mentions a non-EUR ancestry group, excluding multi-ancestry
    GWAS studies entirely.

    Args:
        sample_size: Free-text INITIAL SAMPLE SIZE field from GWAS Catalog.

    Returns:
        True if the study is EUR-only.
    """
    if not isinstance(sample_size, str) or not sample_size.strip():
        return False
    s_lower = sample_size.lower()
    if "european ancestry" not in s_lower and "european" not in s_lower:
        return False
    for kw in NON_EUR_KEYWORDS:
        if kw in s_lower:
            return False
    return True


def filter_catalog(df: pd.DataFrame, p_threshold: float) -> pd.DataFrame:
    """Apply all filtering criteria to a catalog chunk.

    Args:
        df: Raw GWAS Catalog dataframe chunk (dtype=str).
        p_threshold: Maximum p-value to retain.

    Returns:
        Filtered dataframe.
    """
    n_start = len(df)

    # MS trait filter (case-insensitive partial match)
    mask_trait = df["DISEASE/TRAIT"].str.contains(
        "multiple sclerosis", case=False, na=False
    )
    df = df[mask_trait].copy()
    logger.debug(f"  After MS trait filter: {len(df):,} / {n_start:,}")

    if df.empty:
        return df

    # P-value filter
    df["P-VALUE"] = pd.to_numeric(df["P-VALUE"], errors="coerce")
    df = df[df["P-VALUE"] < p_threshold].copy()
    logger.debug(f"  After p < {p_threshold}: {len(df):,}")

    if df.empty:
        return df

    # Autosome filter
    df["CHR_ID"] = df["CHR_ID"].astype(str).str.strip()
    df = df[df["CHR_ID"].isin(AUTOSOMES)].copy()
    logger.debug(f"  After autosome filter: {len(df):,}")

    if df.empty:
        return df

    # Non-null CHR_POS
    df["CHR_POS"] = pd.to_numeric(df["CHR_POS"], errors="coerce")
    df = df[df["CHR_POS"].notna()].copy()
    df["CHR_POS"] = df["CHR_POS"].astype(int)
    logger.debug(f"  After non-null CHR_POS: {len(df):,}")

    if df.empty:
        return df

    # EUR-only study filter
    eur_mask = df["INITIAL SAMPLE SIZE"].apply(is_eur_only)
    df = df[eur_mask].copy()
    logger.debug(f"  After EUR-only filter: {len(df):,}")

    return df


# ---------------------------------------------------------------------------
# Liftover
# ---------------------------------------------------------------------------

def liftover_positions(df: pd.DataFrame, chain_path: str) -> pd.DataFrame:
    """Lift CHR_POS from GRCh37 to GRCh38 using pyliftover.

    Rows that fail liftover are dropped and counted. Positions are 1-based
    in both input (CHR_POS) and output (pos_hg38).

    Args:
        df: Filtered catalog dataframe with CHR_ID and CHR_POS (GRCh37, 1-based).
        chain_path: Path to hg19ToHg38.over.chain.gz.

    Returns:
        Dataframe with pos_hg38 column added; rows failing liftover removed.
    """
    try:
        from pyliftover import LiftOver
    except ImportError:
        logger.error("pyliftover not installed. Run: pip install pyliftover")
        sys.exit(1)

    logger.info(f"Loading liftover chain: {chain_path}")
    lo = LiftOver(chain_path)

    pos_hg38 = []
    n_fail = 0

    for _, row in df.iterrows():
        chrom = f"chr{row['CHR_ID']}"
        pos_hg37 = int(row["CHR_POS"])
        # pyliftover is 0-based; CHR_POS is 1-based
        result = lo.convert_coordinate(chrom, pos_hg37 - 1)
        if result and len(result) > 0:
            pos_hg38.append(result[0][1] + 1)  # back to 1-based
        else:
            pos_hg38.append(None)
            n_fail += 1

    df = df.copy()
    df["pos_hg38"] = pos_hg38

    n_before = len(df)
    df = df[df["pos_hg38"].notna()].copy()
    df["pos_hg38"] = df["pos_hg38"].astype(int)

    logger.info(
        f"Liftover: {n_before:,} input, {len(df):,} success, {n_fail:,} failed (dropped)"
    )
    return df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Parse GWAS Catalog for EUR-only MS loci and liftover to GRCh38. "
            "All paths resolved from config/config.json."
        )
    )
    parser.add_argument(
        "--config",
        default="config/config.json",
        help="Path to config.json (default: config/config.json)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing output file",
    )
    return parser.parse_args()


def main() -> None:
    """Filter GWAS Catalog to EUR-only MS loci, liftover to GRCh38, write TSV."""
    args = parse_args()
    run_start = time.time()

    config = load_config(args.config)
    inputs = config.get("inputs", {})
    outputs = config.get("outputs", {})
    params = config.get("params", {})

    # Resolve paths from config
    catalog_path = inputs.get("gwas_catalog_raw", "")
    chain_path = inputs.get("liftover_chain", "")
    locus_def_dir = outputs.get("locus_def_dir", "results/2-locus_definition")
    output_path = os.path.join(locus_def_dir, "gwas_catalog_ms_eur_hg38.tsv")
    p_threshold = float(params.get("gwas_catalog_p_threshold", 5e-8))
    gcs_dest = params.get(
        "gwas_catalog_gcs_path",
        "gs://fc-secure-b43840eb-548f-464d-bece-31ac7a969abd/reference/gwas_catalog_ms_eur_hg38.tsv",
    )

    logger.info("=" * 60)
    logger.info("02d: Parse GWAS Catalog -- EUR MS Loci")
    logger.info("=" * 60)
    logger.info(f"  Timestamp : {datetime.now().isoformat()}")
    logger.info(f"  Config    : {args.config}")
    logger.info(f"  Catalog   : {catalog_path}")
    logger.info(f"  Chain     : {chain_path}")
    logger.info(f"  Output    : {output_path}")
    logger.info(f"  P thresh  : {p_threshold}")
    logger.info("")

    # Pre-flight
    if not catalog_path:
        logger.error("inputs.gwas_catalog_raw not set in config.json")
        sys.exit(1)
    if not os.path.exists(catalog_path):
        logger.error(f"Catalog file not found: {catalog_path}")
        sys.exit(1)
    if not chain_path:
        logger.error("inputs.liftover_chain not set in config.json")
        sys.exit(1)

    if os.path.exists(output_path) and not args.force:
        logger.info(f"Output already exists (use --force to overwrite): {output_path}")
        sys.exit(0)

    os.makedirs(locus_def_dir, exist_ok=True)

    # Ensure chain file (downloads if missing)
    ensure_chain_file(chain_path)

    # Read and filter catalog in chunks to manage memory (~500MB file)
    logger.info(f"Reading catalog in chunks of {CHUNK_SIZE:,} rows...")
    t0 = time.time()

    chunks_filtered = []
    n_total = 0

    for i, chunk in enumerate(
        pd.read_csv(
            catalog_path,
            sep="\t",
            dtype=str,
            chunksize=CHUNK_SIZE,
            on_bad_lines="warn",
            low_memory=False,
        )
    ):
        n_total += len(chunk)
        filtered = filter_catalog(chunk, p_threshold)
        if len(filtered) > 0:
            chunks_filtered.append(filtered)
        if (i + 1) % 10 == 0:
            logger.info(
                f"  Processed {n_total:,} rows... "
                f"({len(chunks_filtered)} chunks with hits)"
            )

    logger.info(
        f"Catalog read complete: {n_total:,} total rows in {_fmt_elapsed(time.time() - t0)}"
    )

    if not chunks_filtered:
        logger.error(
            "No rows passed filtering. "
            "Check catalog path, column names, and that the trait phrase 'multiple sclerosis' is present."
        )
        sys.exit(1)

    df = pd.concat(chunks_filtered, ignore_index=True)
    logger.info(f"Rows after all filters (pre-liftover): {len(df):,}")

    # Liftover GRCh37 -> GRCh38
    t1 = time.time()
    df = liftover_positions(df, chain_path)
    logger.info(f"Liftover complete in {_fmt_elapsed(time.time() - t1)}")

    # Select and rename output columns
    df["chrom"] = "chr" + df["CHR_ID"].astype(str)
    df = df.rename(
        columns={
            "CHR_POS": "pos_hg37",
            "SNPS": "rsid",
            "DISEASE/TRAIT": "trait",
            "P-VALUE": "p_value",
            "PUBMEDID": "pubmed_id",
            "INITIAL SAMPLE SIZE": "sample_size",
        }
    )

    out_cols = [
        "chrom", "pos_hg37", "pos_hg38", "rsid",
        "trait", "p_value", "pubmed_id", "sample_size",
    ]
    out_cols = [c for c in out_cols if c in df.columns]
    df = df[out_cols].sort_values(["chrom", "pos_hg38"]).reset_index(drop=True)

    # Deduplicate on chrom + pos_hg38 + rsid (keep lowest p per locus)
    n_before_dedup = len(df)
    df = (
        df.sort_values("p_value")
        .drop_duplicates(subset=["chrom", "pos_hg38", "rsid"], keep="first")
        .sort_values(["chrom", "pos_hg38"])
        .reset_index(drop=True)
    )
    logger.info(f"Deduplicated: {n_before_dedup:,} -> {len(df):,} unique loci")

    # Write output
    df.to_csv(output_path, sep="\t", index=False)
    logger.info(f"Written: {output_path}  ({len(df):,} EUR MS loci, GRCh38)")

    # Summary
    logger.info("")
    logger.info("Summary:")
    logger.info(f"  Total catalog rows processed : {n_total:,}")
    logger.info(f"  EUR-only MS loci (GRCh38)    : {len(df):,}")
    logger.info(f"  Chromosomes represented       : {sorted(df['chrom'].unique())}")
    logger.info(f"  Unique studies (PubMed IDs)   : {df['pubmed_id'].nunique()}")
    logger.info("")
    logger.info(f"Total runtime: {_fmt_elapsed(time.time() - run_start)}")
    logger.info("Done.")
    logger.info("")
    logger.info("Next step: upload to GCS workspace bucket with:")
    logger.info(f"  gsutil cp {output_path} {gcs_dest}")


if __name__ == "__main__":
    main()
