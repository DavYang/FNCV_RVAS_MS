#!/bin/bash
set -e

# --- Configuration ---
SCRIPT_NAME="00_sample_common_SNPs_mt"
BASE_DIR=$(pwd)
PYTHON_DIR="${BASE_DIR}/python"
LOG_DIR="${BASE_DIR}/logs"
CONFIG_DIR="${BASE_DIR}/config"

# --- Pre-flight checks ---
if [ -z "${WORKSPACE_BUCKET}" ]; then
    echo "ERROR: WORKSPACE_BUCKET environment variable is not set!"
    echo "This script must be run in an AoU workspace with bucket access."
    exit 1
fi

if [ ! -f "${PYTHON_DIR}/${SCRIPT_NAME}.py" ]; then
    echo "Error: Python script not found at ${PYTHON_DIR}/${SCRIPT_NAME}.py"
    exit 1
fi

if [ ! -f "${CONFIG_DIR}/config.json" ]; then
    echo "Error: Config file not found at ${CONFIG_DIR}/config.json"
    exit 1
fi

# --- Setup ---
mkdir -p "$LOG_DIR"

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="${LOG_DIR}/snp_sampling_mt_${TIMESTAMP}.log"

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
echo "Approach: MT-based (ACAF splitMT)"
echo "Target: ${TARGET_SNPS} SNPs (${BASE_FILENAME})"
echo "Chromosomes: Autosomes only (chr1-22)"

# Check for test mode
if grep -q '"test_mode": true' "${CONFIG_DIR}/config.json"; then
    TEST_CHR=$(grep '"test_chromosome"' "${CONFIG_DIR}/config.json" | cut -d'"' -f4 | cut -d'"' -f1)
    echo ""
    echo "TEST MODE: Processing only ${TEST_CHR}"
    echo "Target for ${TEST_CHR}: proportional to chromosome size from total ${TARGET_SNPS}"
    echo ""
    echo "Steps:"
    echo "  1. Read ACAF splitMT (partition-pruned to ${TEST_CHR})"
    echo "  2. Filter to EUR samples, select GT only"
    echo "  3. Count variant pool"
    echo "  4. Sample proportional SNPs"
    echo "  5. Write per-chromosome MT to GCS"
    echo ""
    echo "Expected output:"
    echo "  gs://\${WORKSPACE_BUCKET}/results/FNCV_RVAS_MS/${BASE_FILENAME}_<timestamp>/"
    echo "  +-- ${TEST_CHR}/"
    echo "  |   +-- ${TEST_CHR}_background_snps.mt/"
    echo "  |   +-- ${TEST_CHR}_sampling_summary.json"
    echo "  +-- ${BASE_FILENAME}_test_summary.json"
else
    echo ""
    echo "PRODUCTION MODE: Processing all autosomes (chr1-22)"
    echo ""
    echo "Steps:"
    echo "  1. Calculate proportional targets per chromosome"
    echo "  2. For each autosome: partition-pruned read -> EUR filter -> sample -> write MT"
    echo "  3. Write overall summary"
    echo ""
    echo "Expected output:"
    echo "  gs://\${WORKSPACE_BUCKET}/results/FNCV_RVAS_MS/${BASE_FILENAME}_<timestamp>/"
    echo "  +-- chr1/ ... chr22/"
    echo "  |   +-- <chr>_background_snps.mt/"
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
echo "==========================================="
