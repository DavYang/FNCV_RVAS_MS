#!/bin/bash
set -e

# ---------------------------------------------------------------------------
# QC per-chromosome PLINK files and merge into a single genome-wide fileset
# for Regenie Step 1 null model fitting.
# ---------------------------------------------------------------------------
# Usage:
#   nohup bash bash/test_02_qc_merge_background_snps.sh > /dev/null 2>&1 &
#
# Monitor:
#   tail -f logs/test_02_qc_merge_background_snps.log
#
# Input:  results/1-bg_snp/plink_no-qc/test/chrN_background.{bed,bim,fam}
# Output: results/1-bg_snp/plink_qc_test/chrN_background_qc.{bed,bim,fam}
#         results/1-bg_snp/plink_qc_test/all_background_qc.{bed,bim,fam}
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_DIR="/home/jupyter/FNCV_RVAS_MS"

INPUT_DIR="${PROJECT_DIR}/results/1-bg_snp/plink_no-qc/test"
OUTPUT_DIR="${PROJECT_DIR}/results/1-bg_snp/plink_qc/test"
LOG_DIR="${PROJECT_DIR}/logs"
LOG_FILE="${LOG_DIR}/02_qc_merge_background_snps.log"
INPUT_PREFIX="${INPUT_DIR}/chr21_background"
OUTPUT_PREFIX="${OUTPUT_DIR}/chr21_background"
# ---------------------------------------------------------------------------
# QC thresholds
# ---------------------------------------------------------------------------
GENO=0.05
HWE=1e-6
MAF=0.01
MIND=0.01

plink2 \
    --bfile "${INPUT_PREFIX}" \
    --maf "${MAF}" \
    --geno "${GENO}" \
    --mind "${MIND}" \
    --hwe "${HWE}" \
    --max-alleles 2 \
    --make-bed \
    --out "${OUTPUT_PREFIX}" \
    --threads 50 \
    --memory 16000



# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}"
