#!/bin/bash
set -eo pipefail

# ---------------------------------------------------------------------------
# 01_build_background_snps.sh
#
# Single entry point for the background SNP PLINK pipeline (Regenie Step 1).
#
# Phase 1: Hail samples loci from ACAF splitMT rows Table (no genotypes),
#           exports plink2-compatible range/keep files, then exits.
# Phase 2: plink2 extracts sampled loci + EUR samples from AoU pre-built
#           per-chromosome PLINK files. No Spark/JVM involved.
#
# Usage (run with nohup to survive disconnects):
#   nohup bash/01_build_background_snps.sh > /dev/null 2>&1 &
#   tail -f logs/01_background_snps_*.log
# ---------------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
PYTHON_SCRIPT="${PROJECT_DIR}/python/01a_sample_loci.py"
CONFIG_FILE="${PROJECT_DIR}/config/config.json"
LOG_DIR="${PROJECT_DIR}/logs"

mkdir -p "$LOG_DIR"

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="${LOG_DIR}/01_background_snps_${TIMESTAMP}.log"

# Redirect all output to log file and stdout
exec > >(tee -a "${LOG_FILE}") 2>&1

# ---------------------------------------------------------------------------
# Pre-flight checks
# ---------------------------------------------------------------------------
echo "============================================================"
echo "Background SNP PLINK Pipeline (plink2-based)"
echo "============================================================"
echo "Started at: $(date)"
echo ""

if [ -z "$WORKSPACE_BUCKET" ]; then
    echo "ERROR: WORKSPACE_BUCKET environment variable not set"
    exit 1
fi
echo "Workspace bucket: $WORKSPACE_BUCKET"

if [ -z "$GOOGLE_PROJECT" ]; then
    echo "WARNING: GOOGLE_PROJECT not set. Requester Pays buckets may fail."
fi

if [ ! -f "$PYTHON_SCRIPT" ]; then
    echo "ERROR: Python script not found: $PYTHON_SCRIPT"
    exit 1
fi

if [ ! -f "$CONFIG_FILE" ]; then
    echo "ERROR: Config file not found: $CONFIG_FILE"
    exit 1
fi

if ! command -v plink2 &> /dev/null; then
    echo "ERROR: plink2 not found in PATH"
    exit 1
fi
echo "plink2: $(which plink2)"

# ---------------------------------------------------------------------------
# Read config
# ---------------------------------------------------------------------------
TARGET_SNPS=$(python3 -c "import json; c=json.load(open('${CONFIG_FILE}')); print(c['sampling']['target_total_snps'])")
TEST_MODE=$(python3 -c "import json; c=json.load(open('${CONFIG_FILE}')); print(c['params'].get('test_mode', False))")
TEST_CHR=$(python3 -c "import json; c=json.load(open('${CONFIG_FILE}')); print(c['params'].get('test_chromosome', 'None'))")
PLINK_BED_BASE=$(python3 -c "import json; c=json.load(open('${CONFIG_FILE}')); print(c['inputs']['plink_bed_base'])")

echo "Target SNPs: ${TARGET_SNPS}"
echo "PLINK source: ${PLINK_BED_BASE}"

if [ "$TEST_MODE" = "True" ]; then
    echo "Mode: TEST (${TEST_CHR} only)"
else
    echo "Mode: PRODUCTION (chr1-22)"
fi

# Format SNP label
if [ "$TARGET_SNPS" -ge 1000000 ]; then
    SNP_LABEL="$((TARGET_SNPS / 1000000))M"
elif [ "$TARGET_SNPS" -ge 1000 ]; then
    SNP_LABEL="$((TARGET_SNPS / 1000))K"
else
    SNP_LABEL="$TARGET_SNPS"
fi

# Output directory on GCS
OUTPUT_DIR="${WORKSPACE_BUCKET}/results/FNCV_RVAS_MS/${SNP_LABEL}_background_snps_${TIMESTAMP}"

echo "Output directory: ${OUTPUT_DIR}"
echo "Log file: ${LOG_FILE}"
echo "------------------------------------------------------------"
echo ""

# ---------------------------------------------------------------------------
# Phase 1: Hail loci sampling + plink2 input file export
# ---------------------------------------------------------------------------
echo "============================================================"
echo "PHASE 1: Hail loci sampling"
echo "============================================================"
echo ""

cd "$PROJECT_DIR"
echo "Starting Phase 1 (Hail)..."
if ! python3 "$PYTHON_SCRIPT" --output-dir "$OUTPUT_DIR"; then
    echo "ERROR: Phase 1 (Hail) failed with exit code $?"
    exit 1
fi

echo ""
echo "Phase 1 complete. Hail has exited. JVM memory freed."
echo ""

# ---------------------------------------------------------------------------
# Phase 2: plink2 extraction per chromosome
# ---------------------------------------------------------------------------
echo "============================================================"
echo "PHASE 2: plink2 extraction from pre-built PLINK files"
echo "============================================================"
echo ""

LOCAL_TMP="${PROJECT_DIR}/tmp_plink2_work"
mkdir -p "$LOCAL_TMP"

PLINK2_INPUTS="${OUTPUT_DIR}/plink2_inputs"
PLINK_OUTPUT_DIR="${OUTPUT_DIR}/null_model_data"
EUR_KEEP_FILE="${LOCAL_TMP}/eur_samples.keep"

# Download EUR keep-file once
echo "Downloading EUR keep-file ..."
gsutil -u "${GOOGLE_PROJECT}" -q cp "${PLINK2_INPUTS}/eur_samples.keep" "${EUR_KEEP_FILE}"
EUR_COUNT=$(wc -l < "${EUR_KEEP_FILE}")
echo "EUR samples: ${EUR_COUNT}"
echo ""

# Determine which chromosomes to process
if [ "$TEST_MODE" = "True" ]; then
    # Extract chromosome number from test_chromosome (e.g., "chr21" -> "21")
    CHR_NUM="${TEST_CHR#chr}"
    CHROMOSOMES="${CHR_NUM}"
else
    CHROMOSOMES=$(seq 1 22)
fi

TOTAL_VARIANTS=0
FAILED_CHROMS=""

for chr in ${CHROMOSOMES}; do
    echo "------------------------------------------------------------"
    echo "Processing chr${chr} ..."

    # Check if output already exists (resume support)
    OUTPUT_EXISTS=$(gsutil -u "${GOOGLE_PROJECT}" -q stat "${PLINK_OUTPUT_DIR}/background_snps_chr${chr}.bed" 2>&1 || true)
    if echo "$OUTPUT_EXISTS" | grep -q "No URLs matched"; then
        : # File does not exist, proceed
    else
        if ! echo "$OUTPUT_EXISTS" | grep -q "CommandException"; then
            echo "  chr${chr}: Output already exists, skipping"
            continue
        fi
    fi

    # Download pre-built PLINK files first (need .bim to detect chr format)
    echo "  Downloading pre-built PLINK files ..."
    GCS_PLINK="${PLINK_BED_BASE}/chr${chr}"
    gsutil -u "${GOOGLE_PROJECT}" -q cp "${GCS_PLINK}.bed" "${LOCAL_TMP}/chr${chr}.bed"
    gsutil -u "${GOOGLE_PROJECT}" -q cp "${GCS_PLINK}.bim" "${LOCAL_TMP}/chr${chr}.bim"
    gsutil -u "${GOOGLE_PROJECT}" -q cp "${GCS_PLINK}.fam" "${LOCAL_TMP}/chr${chr}.fam"

    # Auto-detect chromosome naming in .bim (chr-prefixed vs bare number)
    BIM_CHR=$(head -1 "${LOCAL_TMP}/chr${chr}.bim" | cut -f1)
    if echo "$BIM_CHR" | grep -q "^chr"; then
        RANGE_SUFFIX="_variants.range"
        echo "  .bim uses chr-prefixed chromosomes (${BIM_CHR})"
    else
        RANGE_SUFFIX="_variants_bare.range"
        echo "  .bim uses bare chromosome numbers (${BIM_CHR})"
    fi

    # Download the matching range file
    RANGE_FILE_GCS="${PLINK2_INPUTS}/chr${chr}${RANGE_SUFFIX}"
    RANGE_FILE_LOCAL="${LOCAL_TMP}/chr${chr}_variants.range"
    echo "  Downloading variant range file ..."
    if ! gsutil -u "${GOOGLE_PROJECT}" -q cp "${RANGE_FILE_GCS}" "${RANGE_FILE_LOCAL}" 2>/dev/null; then
        echo "  WARNING: No range file for chr${chr}, skipping"
        rm -f "${LOCAL_TMP}/chr${chr}"*
        continue
    fi

    VARIANT_COUNT=$(wc -l < "${RANGE_FILE_LOCAL}")
    echo "  Variants to extract: ${VARIANT_COUNT}"

    if [ "$VARIANT_COUNT" -eq 0 ]; then
        echo "  WARNING: Empty range file for chr${chr}, skipping"
        rm -f "${RANGE_FILE_LOCAL}" "${LOCAL_TMP}/chr${chr}"*
        continue
    fi

    # Run plink2 extraction
    echo "  Running plink2 --extract range --keep --make-bed ..."
    LOCAL_OUT="${LOCAL_TMP}/background_snps_chr${chr}"

    if ! plink2 \
        --bfile "${LOCAL_TMP}/chr${chr}" \
        --extract range "${RANGE_FILE_LOCAL}" \
        --keep "${EUR_KEEP_FILE}" \
        --make-bed \
        --memory 4000 \
        --out "${LOCAL_OUT}"; then
        echo "  ERROR: plink2 failed for chr${chr}"
        FAILED_CHROMS="${FAILED_CHROMS} chr${chr}"
        rm -f "${LOCAL_TMP}/chr${chr}"* "${LOCAL_OUT}"*
        continue
    fi

    # Check output
    if [ ! -f "${LOCAL_OUT}.bed" ]; then
        echo "  ERROR: plink2 output not found for chr${chr}"
        FAILED_CHROMS="${FAILED_CHROMS} chr${chr}"
        rm -f "${LOCAL_TMP}/chr${chr}"* "${LOCAL_OUT}"*
        continue
    fi

    EXTRACTED=$(wc -l < "${LOCAL_OUT}.bim")
    TOTAL_VARIANTS=$((TOTAL_VARIANTS + EXTRACTED))
    echo "  Extracted ${EXTRACTED} variants for ${EUR_COUNT} EUR samples"

    # Upload results
    echo "  Uploading to GCS ..."
    gsutil -u "${GOOGLE_PROJECT}" -q cp "${LOCAL_OUT}.bed" "${PLINK_OUTPUT_DIR}/background_snps_chr${chr}.bed"
    gsutil -u "${GOOGLE_PROJECT}" -q cp "${LOCAL_OUT}.bim" "${PLINK_OUTPUT_DIR}/background_snps_chr${chr}.bim"
    gsutil -u "${GOOGLE_PROJECT}" -q cp "${LOCAL_OUT}.fam" "${PLINK_OUTPUT_DIR}/background_snps_chr${chr}.fam"

    # Cleanup local files for this chromosome
    rm -f "${LOCAL_TMP}/chr${chr}"* "${LOCAL_OUT}"* "${RANGE_FILE_LOCAL}"

    echo "  chr${chr}: Done"
done

# Cleanup
rm -rf "${LOCAL_TMP}"

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo ""
echo "============================================================"
echo "PIPELINE COMPLETE"
echo "============================================================"
echo "Total variants extracted: ${TOTAL_VARIANTS}"
echo ""

if [ -n "$FAILED_CHROMS" ]; then
    echo "WARNING: Failed chromosomes:${FAILED_CHROMS}"
    echo ""
fi

echo "Output locations:"
echo "  PLINK files: ${PLINK_OUTPUT_DIR}/"
echo "  Sampled loci: ${OUTPUT_DIR}/sampled_loci.ht"
echo "  plink2 inputs: ${PLINK2_INPUTS}/"
echo "  Log file: ${LOG_FILE}"
echo ""

if [ "$TEST_MODE" = "True" ]; then
    echo "=== TEST COMPLETE: ${TEST_CHR} ==="
else
    echo "=== PRODUCTION COMPLETE: chr1-22 ==="
fi

echo "Finished at: $(date)"
echo "============================================================"
