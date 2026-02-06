#!/bin/bash
set -e

# --- Configuration ---
SCRIPT_NAME="00_sample_common_SNPs"
BASE_DIR=$(pwd)
PYTHON_DIR="${BASE_DIR}/python"
LOG_DIR="${BASE_DIR}/logs"
CONFIG_DIR="${BASE_DIR}/config"
LOCAL_OUTPUT_DIR="${BASE_DIR}/results"

# --- Pre-flight checks ---
# Ensure WORKSPACE_BUCKET is set
if [ -z "${WORKSPACE_BUCKET}" ]; then
    echo "ERROR: WORKSPACE_BUCKET environment variable is not set!"
    echo "This script must be run in an AoU workspace with bucket access."
    exit 1
fi

# Check if the python script exists
if [ ! -f "${PYTHON_DIR}/${SCRIPT_NAME}.py" ]; then
    echo "Error: Python script not found at ${PYTHON_DIR}/${SCRIPT_NAME}.py"
    exit 1
fi

# Check if config exists
if [ ! -f "${CONFIG_DIR}/config.json" ]; then
    echo "Error: Config file not found at ${CONFIG_DIR}/config.json"
    exit 1
fi

# --- Setup ---
mkdir -p "$LOG_DIR"
mkdir -p "$LOCAL_OUTPUT_DIR"

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="${LOG_DIR}/snp_sampling_${TIMESTAMP}.log"

# Read target SNP count from config
TARGET_SNPS=$(grep "target_total_snps" "${CONFIG_DIR}/config.json" | head -1 | grep -o '[0-9]*')

if [ "$TARGET_SNPS" -ge 1000000 ]; then
    SNP_LABEL=$((TARGET_SNPS / 1000000))"M"
elif [ "$TARGET_SNPS" -ge 1000 ]; then
    SNP_LABEL=$((TARGET_SNPS / 1000))"K"
else
    SNP_LABEL="$TARGET_SNPS"
fi

BASE_FILENAME="${SNP_LABEL}_background_snps"

echo "==========================================="
echo "Date: ${TIMESTAMP}"
echo "Script: ${PYTHON_DIR}/${SCRIPT_NAME}.py"
echo "Log File: ${LOG_FILE}"
echo "Config: ${CONFIG_DIR}/config.json"
echo "Workspace Bucket: ${WORKSPACE_BUCKET}"
echo "==========================================="

echo "Target: ${TARGET_SNPS} SNPs (${BASE_FILENAME})"

# Check for test mode
if grep -q '"test_mode": true' "${CONFIG_DIR}/config.json"; then
    TEST_CHR=$(grep '"test_chromosome"' "${CONFIG_DIR}/config.json" | cut -d'"' -f4 | cut -d'"' -f1)
    echo ""
    echo "TEST MODE: Processing only ${TEST_CHR}"
    echo "Target for ${TEST_CHR}: proportional to chromosome size from total ${TARGET_SNPS}"
    echo ""
    echo "Steps:"
    echo "  1. Find VCF shards for ${TEST_CHR}"
    echo "  2. Import VCF, filter to EUR samples"
    echo "  3. Sample SNPs (proportional to chromosome size)"
    echo "  4. Export ${TEST_CHR} VCF"
    echo "  5. Generate test summary"
    echo ""
    echo "Expected output:"
    echo "  ${BASE_FILENAME}_<timestamp>/"
    echo "  +-- ${TEST_CHR}/"
    echo "  |   +-- ${TEST_CHR}_background_snps.vcf.bgz"
    echo "  |   +-- ${TEST_CHR}_sampling_summary.json"
    echo "  +-- ${BASE_FILENAME}_test_summary.json"
else
    echo ""
    echo "PRODUCTION MODE: Processing all chromosomes (chr1-22, chrX, chrY)"
    echo "Each chromosome target: proportional to chromosome size from total ${TARGET_SNPS}"
    echo ""
    echo "Steps:"
    echo "  1. Calculate proportional SNP targets per chromosome"
    echo "  2. For each chromosome:"
    echo "     a. Find matching VCF shards"
    echo "     b. Import VCF, filter to EUR samples"
    echo "     c. Sample target SNPs"
    echo "     d. Export chromosome VCF + summary"
    echo "  3. Create overall summary"
    echo ""
    echo "Expected output:"
    echo "  ${BASE_FILENAME}_<timestamp>/"
    echo "  +-- chr1/  ... chr22/ chrX/ chrY/"
    echo "  |   +-- <chr>_background_snps.vcf.bgz"
    echo "  |   +-- <chr>_sampling_summary.json"
    echo "  +-- ${BASE_FILENAME}_overall_summary.json"
fi
echo ""

# --- Execution ---
nohup python3 "${PYTHON_DIR}/${SCRIPT_NAME}.py" > "${LOG_FILE}" 2>&1 &
PID=$!
echo "Job submitted! PID: ${PID}"
echo ""
echo "Monitor progress:"
echo "  tail -f ${LOG_FILE}"
echo ""
echo "Copy results after completion:"
echo "  gsutil ls gs://\${WORKSPACE_BUCKET}/results/FNCV_RVAS_MS/${BASE_FILENAME}_*/"
echo "  gsutil -m cp -r gs://\${WORKSPACE_BUCKET}/results/FNCV_RVAS_MS/${BASE_FILENAME}_<timestamp>/ ${LOCAL_OUTPUT_DIR}/"
echo "==========================================="