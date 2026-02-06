#!/bin/bash
set -e

# --- Configuration ---
SCRIPT_NAME="00_sample_common_SNPs"
BASE_DIR=$(pwd)
PYTHON_DIR="${BASE_DIR}/python"
LOG_DIR="${BASE_DIR}/logs"
CONFIG_DIR="${BASE_DIR}/config"
LOCAL_OUTPUT_DIR="${BASE_DIR}/results"

# --- Setup ---
# Ensure directories exist
mkdir -p "$LOG_DIR"
mkdir -p "$PYTHON_DIR"
mkdir -p "$CONFIG_DIR"
mkdir -p "$LOCAL_OUTPUT_DIR"

# Generate a timestamp for the log file and export
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="${LOG_DIR}/snp_sampling_${TIMESTAMP}.log"

# Read target SNP count from config to determine naming
TARGET_SNPS=$(grep "target_total_snps" "${CONFIG_DIR}/config.json" | head -1 | grep -o '[0-9]*')

if [ "$TARGET_SNPS" -ge 1000000 ]; then
    SNP_LABEL=$((TARGET_SNPS / 1000000))"M"
elif [ "$TARGET_SNPS" -ge 1000 ]; then
    SNP_LABEL=$((TARGET_SNPS / 1000))"K"
else
    SNP_LABEL="$TARGET_SNPS"
fi

BASE_FILENAME="${SNP_LABEL}_background_snps"

echo "Date: ${TIMESTAMP}"
echo "Script: ${PYTHON_DIR}/${SCRIPT_NAME}.py"
echo "Log File: ${LOG_FILE}"
echo "Config: ${CONFIG_DIR}/config.json"
echo "==========================================="

# --- Execution ---
# Check if the python script exists
if [ ! -f "${PYTHON_DIR}/${SCRIPT_NAME}.py" ]; then
    echo "Error: Python script not found at ${PYTHON_DIR}/${SCRIPT_NAME}.py"
    exit 1
fi

echo "Starting chromosome-based SNP sampling from sharded VCFs..."
echo "Target: ${TARGET_SNPS} SNPs (${BASE_FILENAME})"

# Check for test mode
if grep -q '"test_mode": true' "${CONFIG_DIR}/config.json"; then
    TEST_CHR=$(grep '"test_chromosome"' "${CONFIG_DIR}/config.json" | cut -d'"' -f4 | cut -d'"' -f1)
    echo "TEST MODE DETECTED: Processing only ${TEST_CHR}"
    echo "Target for ${TEST_CHR}: Based on chromosome size proportion of total ${TARGET_SNPS}"
else
    echo "PRODUCTION MODE: Processing all chromosomes (chr1-22, chrX, chrY)"
    echo "Each chromosome target: Based on proportional size distribution of total ${TARGET_SNPS}"
fi

echo "This will:"
if grep -q '"test_mode": true' "${CONFIG_DIR}/config.json"; then
    echo "  2. Sample random intervals for ${TEST_CHR}"
    echo "  3. Filter to EUR samples (cached)"
    echo "  4. Sample ${TEST_CHR} SNPs for ${TEST_CHR} (proportional to chromosome size)"
    echo "  5. Export ${TEST_CHR} VCF immediately"
    echo "  6. Generate test summary"
else
    echo "  1. Calculate proportional SNP targets per chromosome"
    echo "  2. Process each chromosome sequentially:"
    echo "     - Sample random intervals per chromosome"
    echo "     - Filter to EUR samples (cached after first chromosome)"
    echo "     - Sample target SNPs per chromosome (proportional to chromosome size)"
    echo "     - Export chromosome-specific VCF immediately"
    echo "  3. Generate per-chromosome summaries"
    echo "  4. Create overall summary"
fi
echo ""

# Run in background with nohup
nohup python3 "${PYTHON_DIR}/${SCRIPT_NAME}.py" > "${LOG_FILE}" 2>&1 &
PID=$!
echo "Job submitted successfully!"
echo "PID: ${PID}"
echo ""
echo "To monitor progress, run:"
echo "  tail -f ${LOG_FILE}"
echo ""
echo "To wait for completion and auto-copy results:"
echo "  wait_for_job_and_copy() { wait $PID; echo 'Job completed, copying results...'; nohup gsutil -m cp -r gs://\${WORKSPACE_BUCKET}/results/FNCV_RVAS_MS/${BASE_FILENAME}_${TIMESTAMP}/ ${LOCAL_OUTPUT_DIR}/ > ${LOG_DIR}/gsutil_copy_${TIMESTAMP}.log 2>&1 & echo 'Copy started in background, check ${LOG_DIR}/gsutil_copy_${TIMESTAMP}.log'; }"
echo "  wait_for_job_and_copy"
echo ""
echo "Expected output structure:"
if grep -q '"test_mode": true' "${CONFIG_DIR}/config.json"; then
    echo "  ${BASE_FILENAME}_${TIMESTAMP}/"
    echo "  ├── ${TEST_CHR}/"
    echo "  │   ├── ${TEST_CHR}_background_snps.vcf.bgz"
    echo "  │   ├── ${TEST_CHR}_background_snps.vcf.bgz.tbi"
    echo "  │   └── ${TEST_CHR}_sampling_summary.json"
    echo "  └── ${BASE_FILENAME}_test_summary.json"
else
    echo "  ${BASE_FILENAME}_${TIMESTAMP}/"
    echo "  ├── chr1/"
    echo "  │   ├── chr1_background_snps.vcf.bgz"
    echo "  │   ├── chr1_background_snps.vcf.bgz.tbi"
    echo "  │   └── chr1_sampling_summary.json"
    echo "  ├── chr2/"
    echo "  │   ├── chr2_background_snps.vcf.bgz"
    echo "  │   ├── chr2_background_snps.vcf.bgz.tbi"
    echo "  │   └── chr2_sampling_summary.json"
    echo "  └── ... (chr1-22, chrX, chrY)"
    echo "  └── ${BASE_FILENAME}_overall_summary.json"
fi
echo ""
echo "To wait for completion and auto-copy results:"
echo "  wait_for_job_and_copy() { wait $PID; echo 'Job completed, copying results...'; nohup gsutil -m cp -r gs://\${WORKSPACE_BUCKET}/results/FNCV_RVAS_MS/${BASE_FILENAME}_${TIMESTAMP}/ ${LOCAL_OUTPUT_DIR}/ > ${LOG_DIR}/gsutil_copy_${TIMESTAMP}.log 2>&1 & echo 'Copy started in background, check ${LOG_DIR}/gsutil_copy_${TIMESTAMP}.log'; }"
echo "  wait_for_job_and_copy"
echo ""
echo "Or copy manually after completion:"
echo "  nohup gsutil -m cp -r gs://\${WORKSPACE_BUCKET}/results/FNCV_RVAS_MS/${BASE_FILENAME}_${TIMESTAMP}/ ${LOCAL_OUTPUT_DIR}/ > ${LOG_DIR}/gsutil_copy_${TIMESTAMP}.log 2>&1 &"
echo "  # Monitor copy progress: tail -f ${LOG_DIR}/gsutil_copy_${TIMESTAMP}.log"
echo ""
echo "==========================================="