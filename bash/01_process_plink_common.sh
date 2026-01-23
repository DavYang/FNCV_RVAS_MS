#!/bin/bash
set -e

# --- Usage Help ---
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <path_to_input_plink_prefix> <path_to_ancestry_tsv> <output_dir>"
    echo "Example: $0 gs://bucket/path/to/plink_data /path/to/ancestry.tsv results/common_variants"
    exit 1
fi

# --- Inputs ---
INPUT_PLINK="$1"
ANCESTRY_FILE="$2"
OUT_DIR="$3"

# --- Configuration ---
# QC Thresholds
MAF=0.01
GENO=0.05
MIND=0.05
HWE=1e-6
TARGET_VARIANTS=500000

# --- Setup ---
mkdir -p "$OUT_DIR"
mkdir -p "$OUT_DIR/tmp"

LOG_FILE="${OUT_DIR}/process_plink_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "=========================================================="
echo "   Processing Common Variants (PLINK Strategy)"
echo "=========================================================="
echo "Input PLINK:  ${INPUT_PLINK}"
echo "Ancestry TSV: ${ANCESTRY_FILE}"
echo "Output Dir:   ${OUT_DIR}"
echo "QC Params:    MAF=${MAF}, GENO=${GENO}, MIND=${MIND}, HWE=${HWE}"
echo "=========================================================="

# 1. Extract European Samples
echo "[1/4] Extracting EUR sample IDs..."
EUR_SAMPLES_FILE="${OUT_DIR}/tmp/eur_samples.txt"

# Expecting TSV with header, columns: research_id, ancestry_pred, ...
# Adjust column numbers ($1, $2) if your TSV format differs
awk -F'\t' '$2 == "eur" {print $1, $1}' "${ANCESTRY_FILE}" > "${EUR_SAMPLES_FILE}"

EUR_COUNT=$(wc -l < "${EUR_SAMPLES_FILE}")
echo "      Found ${EUR_COUNT} European samples."

# 2 & 3. Filter Samples AND Perform QC
echo "[2/4] Filtering samples and running QC..."
QC_OUT="${OUT_DIR}/eur_common_qc"

# Note: Using plink2 for speed. If input is strictly binary .bed/.bim/.fam, plink 1.9 syntax might be needed 
# but plink2 usually handles it. If input is .pgen, plink2 is required.
# --make-bed creates standard binary PLINK files.

plink2 \
    --bfile "${INPUT_PLINK}" \
    --keep "${EUR_SAMPLES_FILE}" \
    --maf "${MAF}" \
    --geno "${GENO}" \
    --mind "${MIND}" \
    --hwe "${HWE}" \
    --make-bed \
    --out "${QC_OUT}" \
    --silent

# Check results
if [ -f "${QC_OUT}.bim" ]; then
    VAR_COUNT_QC=$(wc -l < "${QC_OUT}.bim")
    echo "      QC Complete. Variants remaining: ${VAR_COUNT_QC}"
else
    echo "ERROR: QC output not found. Check PLINK logs."
    exit 1
fi

# 4. Downsample to 500k
echo "[3/4] Downsampling to ${TARGET_VARIANTS} variants..."
FINAL_OUT="${OUT_DIR}/eur_common_500k_final"

plink2 \
    --bfile "${QC_OUT}" \
    --max-alleles 2 \
    --thin-count "${TARGET_VARIANTS}" \
    --make-bed \
    --out "${FINAL_OUT}" \
    --silent

# Final Stats
if [ -f "${FINAL_OUT}.bim" ]; then
    FINAL_VAR_COUNT=$(wc -l < "${FINAL_OUT}.bim")
    echo "=========================================================="
    echo "   Pipeline Completed Successfully"
    echo "=========================================================="
    echo "QC'ed Dataset:      ${QC_OUT} (${VAR_COUNT_QC} variants)"
    echo "Final 500k Dataset: ${FINAL_OUT} (${FINAL_VAR_COUNT} variants)"
    echo "Log File:           ${LOG_FILE}"
else
    echo "ERROR: Final downsampling failed."
    exit 1
fi
