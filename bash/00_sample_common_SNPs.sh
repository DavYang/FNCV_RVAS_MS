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
LOG_FILE="${LOG_DIR}/100k_snp_sampling_${TIMESTAMP}.log"

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

echo "Starting SNP sampling from sharded VCFs..."
echo "This will:"
echo "  1. Load EUR sample IDs from ancestry data"
echo "  2. Select random genome-wide intervals"
echo "  3. Find corresponding VCF shards"
echo "  4. Import and filter VCFs to EUR samples"
echo "  5. Sample exactly 100,000 SNPs (configurable)"
echo "  6. Export as compressed VCF with tabix index"
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
echo "  wait_for_job_and_copy() { wait $PID; echo 'Job completed, copying results...'; nohup gsutil -m cp -r gs://\${WORKSPACE_BUCKET}/results/FNCV_RVAS_MS/100k_background_snps_${TIMESTAMP}/ ${LOCAL_OUTPUT_DIR}/ > ${LOG_DIR}/gsutil_copy_${TIMESTAMP}.log 2>&1 & echo 'Copy started in background, check ${LOG_DIR}/gsutil_copy_${TIMESTAMP}.log'; }"
echo "  wait_for_job_and_copy"
echo ""
echo "Or copy manually after completion:"
echo "  nohup gsutil -m cp -r gs://\${WORKSPACE_BUCKET}/results/FNCV_RVAS_MS/100k_background_snps_${TIMESTAMP}/ ${LOCAL_OUTPUT_DIR}/ > ${LOG_DIR}/gsutil_copy_${TIMESTAMP}.log 2>&1 &"
echo "  # Monitor copy progress: tail -f ${LOG_DIR}/gsutil_copy_${TIMESTAMP}.log"
echo ""
echo "==========================================="