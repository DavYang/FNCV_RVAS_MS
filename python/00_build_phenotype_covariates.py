#!/usr/bin/env python3
"""Build REGENIE-formatted phenotype and covariate files for Phase 0.

Generates two space-delimited output files for REGENIE Step 1:
  - MS_phenotype.txt  : FID IID MS  (1=case, 0=control, NA=missing)
  - MS_covariates.txt : FID IID Age Sex PC1..PC10  (no missing values)

Data sources:
  - EUR sample IDs + PCs : ancestry_preds.tsv (CDR bucket)
  - Sex at birth         : genomic_metrics.tsv (CDR bucket)
  - Year of birth        : BigQuery person table
  - MS case status       : BigQuery condition_occurrence table

Usage:
    python python/00_build_phenotype_covariates.py --config config/config.json
    python python/00_build_phenotype_covariates.py --config config/config.json --test-mode
    python python/00_build_phenotype_covariates.py --config config/config.json --force
"""

import argparse
import ast
import json
import logging
import os
import subprocess
import sys
import tempfile
import time
from typing import List, Set, Tuple

import numpy as np
import pandas as pd
from google.cloud import bigquery

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


logger = setup_logger("phenotype_covariates")


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
# GCS helpers
# ---------------------------------------------------------------------------

def gsutil_cp(src: str, dst: str, google_project: str) -> None:
    """Copy a single file from/to GCS using gsutil with requester-pays."""
    cmd = ["gsutil", "-u", google_project, "cp", src, dst]
    logger.info(f"  gsutil cp {src} -> {dst}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"gsutil cp failed (rc={result.returncode}):\n{result.stderr}"
        )


def gcs_exists(path: str, google_project: str) -> bool:
    """Return True if a GCS object exists."""
    cmd = ["gsutil", "-u", google_project, "ls", path]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.returncode == 0


def gsutil_upload(local_path: str, gcs_path: str, google_project: str) -> None:
    """Upload a local file to GCS."""
    cmd = ["gsutil", "-u", google_project, "cp", local_path, gcs_path]
    logger.info(f"  Uploading {local_path} -> {gcs_path}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"gsutil upload failed (rc={result.returncode}):\n{result.stderr}"
        )


# ---------------------------------------------------------------------------
# Step 1: EUR samples + PCs from ancestry_preds.tsv
# ---------------------------------------------------------------------------

def load_eur_pcs(
    ancestry_preds_path: str,
    google_project: str,
    test_mode: bool,
    test_n: int,
    random_seed: int,
) -> Tuple[pd.DataFrame, List[str]]:
    """Load EUR sample IDs and expand pca_features array into PC1..PC10.

    Returns:
        eur_df  : DataFrame with columns [research_id, PC1..PC10]
        pc_cols : List of PC column names ["PC1", ..., "PC10"]
    """
    logger.info("Step 1: Loading ancestry_preds.tsv...")
    with tempfile.NamedTemporaryFile(suffix=".tsv", delete=False) as tmp:
        tmp_path = tmp.name

    try:
        gsutil_cp(ancestry_preds_path, tmp_path, google_project)
        ancestry_df = pd.read_csv(tmp_path, sep="\t", low_memory=False)
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)

    logger.info(f"  Total samples in ancestry file: {len(ancestry_df):,}")
    logger.info(
        f"  Ancestry distribution:\n"
        + ancestry_df["ancestry_pred"].value_counts().to_string()
    )

    eur_df = ancestry_df[ancestry_df["ancestry_pred"] == "eur"].copy()
    eur_df["research_id"] = eur_df["research_id"].astype(str)
    logger.info(f"  EUR samples: {len(eur_df):,}")

    if "pca_features" not in eur_df.columns:
        raise KeyError(
            f"'pca_features' column not found. Available columns: {list(eur_df.columns)}"
        )

    # Parse pca_features: may be a Python list object or a string representation
    def _parse_pcs(val) -> list:
        if isinstance(val, (list, np.ndarray)):
            return list(val)
        if isinstance(val, str):
            return ast.literal_eval(val)
        return [np.nan] * 10

    logger.info("  Parsing pca_features array into PC1..PC10...")
    pc_matrix = np.array(eur_df["pca_features"].apply(_parse_pcs).tolist())

    if pc_matrix.shape[1] < 10:
        raise ValueError(
            f"Expected >= 10 PCs in pca_features, found {pc_matrix.shape[1]}"
        )

    pc_cols = [f"PC{i}" for i in range(1, 11)]
    for i, col in enumerate(pc_cols):
        eur_df[col] = pc_matrix[:, i]

    eur_df = eur_df[["research_id"] + pc_cols].copy()

    if test_mode:
        eur_df = eur_df.sample(
            n=min(test_n, len(eur_df)), random_state=random_seed
        ).copy()
        logger.info(f"  [TEST MODE] Subsampled to {len(eur_df):,} EUR samples")

    logger.info(f"  Final EUR sample set: {len(eur_df):,}")
    return eur_df, pc_cols


# ---------------------------------------------------------------------------
# Step 2: Sex from genomic_metrics.tsv
# ---------------------------------------------------------------------------

def load_sex(genomic_metrics_path: str, google_project: str) -> pd.DataFrame:
    """Load sex_at_birth from genomic_metrics.tsv and recode to 1=Male, 2=Female.

    Returns:
        DataFrame with columns [research_id, Sex]
    """
    logger.info("Step 2: Loading genomic_metrics.tsv for sex_at_birth...")
    with tempfile.NamedTemporaryFile(suffix=".tsv", delete=False) as tmp:
        tmp_path = tmp.name

    try:
        gsutil_cp(genomic_metrics_path, tmp_path, google_project)
        metrics_df = pd.read_csv(tmp_path, sep="\t", low_memory=False)
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)

    logger.info(f"  Genomic metrics rows: {len(metrics_df):,}")

    if "sex_at_birth" not in metrics_df.columns:
        raise KeyError(
            f"'sex_at_birth' column not found. Available: {list(metrics_df.columns)}"
        )

    sex_map = {"Male": 1, "Female": 2}
    metrics_df["Sex"] = metrics_df["sex_at_birth"].map(sex_map)
    metrics_df["research_id"] = metrics_df["research_id"].astype(str)

    n_unknown = metrics_df["Sex"].isna().sum()
    logger.info(
        f"  Sex distribution (1=M, 2=F): "
        + metrics_df["Sex"].value_counts(dropna=False).to_dict().__str__()
    )
    logger.info(f"  Samples with unknown/other sex (will be excluded): {n_unknown:,}")

    return metrics_df[["research_id", "Sex"]].copy()


# ---------------------------------------------------------------------------
# Step 3: Age from BigQuery person table
# ---------------------------------------------------------------------------

def load_age(
    workspace_cdr: str,
    google_project: str,
    age_reference_year: int,
) -> pd.DataFrame:
    """Query year_of_birth from BigQuery person table and compute Age.

    Returns:
        DataFrame with columns [research_id, Age]
    """
    logger.info("Step 3: Querying person table from BigQuery for year_of_birth...")
    client = bigquery.Client(project=google_project)

    query = f"""
    SELECT
        CAST(p.person_id AS STRING) AS research_id,
        p.year_of_birth
    FROM
        `{workspace_cdr}.person` AS p
    """

    age_df = client.query(query).to_dataframe()
    age_df["research_id"] = age_df["research_id"].astype(str)
    age_df["Age"] = age_reference_year - age_df["year_of_birth"]

    logger.info(f"  Person table rows: {len(age_df):,}")
    logger.info(
        f"  Age range: {age_df['Age'].min():.0f} - {age_df['Age'].max():.0f} "
        f"(mean={age_df['Age'].mean():.1f})"
    )

    return age_df[["research_id", "Age"]].copy()


# ---------------------------------------------------------------------------
# Step 4: MS case/control status from BigQuery condition_occurrence
# ---------------------------------------------------------------------------

def load_ms_cases(
    workspace_cdr: str,
    google_project: str,
    omop_concept_ids: List[int],
    icd10_codes: List[str],
    eur_ids: Set[str],
) -> pd.DataFrame:
    """Query condition_occurrence for MS cases and build phenotype series.

    Cases  = EUR samples with MS diagnosis (OMOP concept or ICD-10 source code).
    Controls = all remaining EUR samples.

    Returns:
        DataFrame with columns [research_id, MS]  (1=case, 0=control)
    """
    logger.info("Step 4: Querying condition_occurrence for MS cases...")
    client = bigquery.Client(project=google_project)

    concept_id_list = ", ".join(str(c) for c in omop_concept_ids)
    icd_code_list = ", ".join(f"'{c}'" for c in icd10_codes)

    query = f"""
    SELECT DISTINCT
        CAST(co.person_id AS STRING) AS research_id
    FROM
        `{workspace_cdr}.condition_occurrence` AS co
    WHERE
        co.condition_concept_id IN ({concept_id_list})
        OR co.condition_source_value IN ({icd_code_list})
    """

    cases_df = client.query(query).to_dataframe()
    cases_df["research_id"] = cases_df["research_id"].astype(str)

    logger.info(f"  Total MS cases in CDR (all ancestries): {len(cases_df):,}")

    eur_cases = set(cases_df["research_id"]) & eur_ids
    eur_controls = eur_ids - eur_cases

    logger.info(f"  EUR MS cases    : {len(eur_cases):,}")
    logger.info(f"  EUR controls    : {len(eur_controls):,}")
    logger.info(f"  EUR total       : {len(eur_ids):,}")

    records = [
        {"research_id": rid, "MS": 1 if rid in eur_cases else 0}
        for rid in eur_ids
    ]
    pheno_df = pd.DataFrame(records)
    pheno_df["research_id"] = pheno_df["research_id"].astype(str)

    return pheno_df


# ---------------------------------------------------------------------------
# Step 5: Merge + format output DataFrames
# ---------------------------------------------------------------------------

def build_output_dfs(
    eur_df: pd.DataFrame,
    sex_df: pd.DataFrame,
    age_df: pd.DataFrame,
    pheno_df: pd.DataFrame,
    pc_cols: List[str],
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Merge all sources and build phenotype + covariate DataFrames.

    Returns:
        pheno_out : DataFrame [FID, IID, MS]
        covar_out : DataFrame [FID, IID, Age, Sex, PC1..PC10]
    """
    logger.info("Step 5: Merging all data sources...")

    merged = eur_df.merge(sex_df, on="research_id", how="inner")
    merged = merged.merge(age_df, on="research_id", how="inner")
    merged = merged.merge(pheno_df, on="research_id", how="inner")
    logger.info(f"  After merge: {len(merged):,} samples")

    # Drop samples with missing Age, Sex, or any PC
    required_cols = ["Age", "Sex"] + pc_cols
    n_before = len(merged)
    merged = merged.dropna(subset=required_cols)
    n_dropped = n_before - len(merged)
    if n_dropped > 0:
        logger.info(f"  Dropped {n_dropped:,} samples with missing Age/Sex/PC values")
    logger.info(f"  Final sample count: {len(merged):,}")

    merged["Sex"] = merged["Sex"].astype(int)
    merged["Age"] = merged["Age"].astype(int)

    # Add FID = IID = research_id
    merged.insert(0, "FID", merged["research_id"])
    merged.insert(1, "IID", merged["research_id"])

    pheno_out = merged[["FID", "IID", "MS"]].copy()
    covar_out = merged[["FID", "IID", "Age", "Sex"] + pc_cols].copy()

    return pheno_out, covar_out


# ---------------------------------------------------------------------------
# Step 6: Validation
# ---------------------------------------------------------------------------

def validate_outputs(
    pheno_out: pd.DataFrame,
    covar_out: pd.DataFrame,
    pc_cols: List[str],
) -> None:
    """Run validation checks on output DataFrames. Raises on any failure."""
    logger.info("Step 6: Running validation checks...")
    errors = []

    # No duplicate IIDs
    n_dup_pheno = pheno_out["IID"].duplicated().sum()
    n_dup_covar = covar_out["IID"].duplicated().sum()
    _check("No duplicate IIDs in phenotype file", n_dup_pheno == 0, errors)
    _check("No duplicate IIDs in covariate file", n_dup_covar == 0, errors)

    # Phenotype values in {0, 1}
    actual_vals = set(pheno_out["MS"].dropna().unique())
    _check(
        f"Phenotype values are subset of {{0, 1}} (found {actual_vals})",
        actual_vals.issubset({0, 1}),
        errors,
    )

    # No missing covariates
    covar_data_cols = ["Age", "Sex"] + pc_cols
    n_missing_covar = covar_out[covar_data_cols].isna().sum().sum()
    _check(
        f"No missing values in covariate columns (found {n_missing_covar})",
        n_missing_covar == 0,
        errors,
    )

    # FID == IID
    _check(
        "FID == IID in phenotype file",
        (pheno_out["FID"] == pheno_out["IID"]).all(),
        errors,
    )

    # Sample summary
    n_cases = (pheno_out["MS"] == 1).sum()
    n_controls = (pheno_out["MS"] == 0).sum()
    n_missing = pheno_out["MS"].isna().sum()
    logger.info(f"  Cases    : {n_cases:,}")
    logger.info(f"  Controls : {n_controls:,}")
    logger.info(f"  Missing  : {n_missing:,}")
    logger.info(f"  Total    : {len(pheno_out):,}")

    # PC sanity check
    pc_max = covar_out[pc_cols].abs().max().max()
    logger.info(f"  Max absolute PC value: {pc_max:.6f}")

    if errors:
        for e in errors:
            logger.error(f"  FAIL: {e}")
        raise ValueError(f"Validation failed with {len(errors)} error(s)")

    logger.info("  All validation checks passed.")


def _check(label: str, condition: bool, errors: list) -> None:
    status = "PASS" if condition else "FAIL"
    logger.info(f"  [{status}] {label}")
    if not condition:
        errors.append(label)


# ---------------------------------------------------------------------------
# Step 7: Write files + upload to GCS
# ---------------------------------------------------------------------------

def write_and_upload(
    pheno_out: pd.DataFrame,
    covar_out: pd.DataFrame,
    local_dir: str,
    gcs_dir: str,
    google_project: str,
    test_mode: bool,
) -> Tuple[str, str]:
    """Write space-delimited output files locally and upload to GCS.

    Returns:
        (pheno_gcs_path, covar_gcs_path)
    """
    logger.info("Step 7: Writing output files...")
    os.makedirs(local_dir, exist_ok=True)

    suffix = "_TEST" if test_mode else ""
    pheno_fname = f"MS_phenotype{suffix}.txt"
    covar_fname = f"MS_covariates{suffix}.txt"

    pheno_local = os.path.join(local_dir, pheno_fname)
    covar_local = os.path.join(local_dir, covar_fname)

    pheno_out.to_csv(pheno_local, sep=" ", index=False, na_rep="NA")
    covar_out.to_csv(covar_local, sep=" ", index=False, na_rep="NA")

    logger.info(f"  Written: {pheno_local}  ({os.path.getsize(pheno_local):,} bytes)")
    logger.info(f"  Written: {covar_local}  ({os.path.getsize(covar_local):,} bytes)")

    # Upload to GCS (skip for test mode)
    if not test_mode:
        pheno_gcs = f"{gcs_dir}/{pheno_fname}"
        covar_gcs = f"{gcs_dir}/{covar_fname}"
        gsutil_upload(pheno_local, pheno_gcs, google_project)
        gsutil_upload(covar_local, covar_gcs, google_project)
        logger.info(f"  GCS: {pheno_gcs}")
        logger.info(f"  GCS: {covar_gcs}")
        return pheno_gcs, covar_gcs
    else:
        logger.info("  [TEST MODE] Skipping GCS upload.")
        return pheno_local, covar_local


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build REGENIE phenotype and covariate files for Phase 0."
    )
    parser.add_argument(
        "--config",
        default="config/config.json",
        help="Path to config JSON file (default: config/config.json)",
    )
    parser.add_argument(
        "--test-mode",
        action="store_true",
        help="Run on a small subsample (TEST_N samples) without uploading to GCS",
    )
    parser.add_argument(
        "--test-n",
        type=int,
        default=500,
        help="Number of EUR samples to use in test mode (default: 500)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing output files in GCS",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    t_start = time.time()

    logger.info("=" * 60)
    logger.info("Phase 0: Phenotype & Covariate File Generation")
    logger.info("=" * 60)
    logger.info(f"Config     : {args.config}")
    logger.info(f"Test mode  : {args.test_mode} (n={args.test_n if args.test_mode else 'ALL'})")
    logger.info(f"Force      : {args.force}")

    # Load config
    config = load_config(args.config)
    inputs = config["inputs"]
    outputs = config["outputs"]
    params = config["params"]

    # Required environment variables
    workspace_bucket = os.environ.get("WORKSPACE_BUCKET", "")
    google_project = os.environ.get("GOOGLE_PROJECT", "")
    workspace_cdr = os.environ.get("WORKSPACE_CDR", "")

    if not workspace_bucket:
        logger.error("WORKSPACE_BUCKET environment variable not set")
        sys.exit(1)
    if not google_project:
        logger.error("GOOGLE_PROJECT environment variable not set")
        sys.exit(1)
    if not workspace_cdr:
        logger.error("WORKSPACE_CDR environment variable not set")
        sys.exit(1)

    if not workspace_bucket.startswith("gs://"):
        workspace_bucket = f"gs://{workspace_bucket}"

    logger.info(f"Workspace  : {workspace_bucket}")
    logger.info(f"Project    : {google_project}")
    logger.info(f"CDR        : {workspace_cdr}")

    # Output paths
    local_out_dir = outputs["phenotype_dir"]
    gcs_out_dir = f"{workspace_bucket}/{outputs['phenotype_dir']}"

    # Resume check (skip if outputs already exist and --force not set)
    if not args.test_mode and not args.force:
        pheno_gcs_check = f"{gcs_out_dir}/MS_phenotype.txt"
        covar_gcs_check = f"{gcs_out_dir}/MS_covariates.txt"
        if gcs_exists(pheno_gcs_check, google_project) and gcs_exists(
            covar_gcs_check, google_project
        ):
            logger.info(
                "Output files already exist in GCS. Use --force to overwrite."
            )
            logger.info(f"  {pheno_gcs_check}")
            logger.info(f"  {covar_gcs_check}")
            sys.exit(0)

    # Step 1: EUR samples + PCs
    eur_df, pc_cols = load_eur_pcs(
        ancestry_preds_path=inputs["ancestry_pred"],
        google_project=google_project,
        test_mode=args.test_mode,
        test_n=args.test_n,
        random_seed=params["random_seed"],
    )
    eur_ids = set(eur_df["research_id"])

    # Step 2: Sex
    sex_df = load_sex(
        genomic_metrics_path=inputs["genomic_metrics"],
        google_project=google_project,
    )

    # Step 3: Age
    age_df = load_age(
        workspace_cdr=workspace_cdr,
        google_project=google_project,
        age_reference_year=params["age_reference_year"],
    )

    # Step 4: MS case/control status
    pheno_df = load_ms_cases(
        workspace_cdr=workspace_cdr,
        google_project=google_project,
        omop_concept_ids=params["ms_omop_concept_ids"],
        icd10_codes=params["ms_icd10_codes"],
        eur_ids=eur_ids,
    )

    # Step 5: Merge + format
    pheno_out, covar_out = build_output_dfs(
        eur_df=eur_df,
        sex_df=sex_df,
        age_df=age_df,
        pheno_df=pheno_df,
        pc_cols=pc_cols,
    )

    # Step 6: Validate
    validate_outputs(pheno_out, covar_out, pc_cols)

    # Step 7: Write + upload
    pheno_path, covar_path = write_and_upload(
        pheno_out=pheno_out,
        covar_out=covar_out,
        local_dir=local_out_dir,
        gcs_dir=gcs_out_dir,
        google_project=google_project,
        test_mode=args.test_mode,
    )

    elapsed = time.time() - t_start
    logger.info("=" * 60)
    logger.info(f"Phase 0 complete in {elapsed:.1f}s")
    logger.info(f"  Phenotype : {pheno_path}")
    logger.info(f"  Covariates: {covar_path}")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
