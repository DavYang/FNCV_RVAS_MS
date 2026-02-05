#!/bin/bash
set -e

# --- Configuration ---
SCRIPT_NAME="00_sample_common_SNPs"
BASE_DIR=$(pwd)
PYTHON_DIR="${BASE_DIR}/python"
LOG_DIR="${BASE_DIR}/logs"
CONFIG_DIR="${BASE_DIR}/config"

# --- Setup ---
# Ensure directories exist
mkdir -p "$LOG_DIR"
mkdir -p "$PYTHON_DIR"
mkdir -p "$CONFIG_DIR"

# Generate a timestamp for the log file and export
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="${LOG_DIR}/100k_snp_sampling_${TIMESTAMP}.log"

echo "==========================================="
echo "   100K SNP Sampling Pipeline (VCF-based)"
echo "==========================================="
echo "Date: ${TIMESTAMP}"
echo "Script: ${PYTHON_DIR}/${SCRIPT_NAME}.py"
echo "Log File: ${LOG_FILE}"
echo "Config: ${CONFIG_DIR}/config.json"
echo "Target: 100,000 random SNPs genome-wide"
echo "==========================================="

# --- Execution ---
# Check if the python script exists
if [ ! -f "${PYTHON_DIR}/${SCRIPT_NAME}.py" ]; then
    echo "Error: Python script not found at ${PYTHON_DIR}/${SCRIPT_NAME}.py"
    exit 1
fi

echo "Starting 100K SNP sampling from sharded VCFs..."
echo "This will:"
echo "  1. Load EUR sample IDs from ancestry data"
echo "  2. Select random genome-wide intervals"
echo "  3. Find corresponding VCF shards"
echo "  4. Import and filter VCFs to EUR samples"
echo "  5. Sample exactly 100,000 SNPs"
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
echo "Expected output:"
echo "  - 100k_background_snps.vcf.bgz (compressed VCF)"
echo "  - 100k_background_snps.vcf.bgz.tbi (tabix index)"
echo "  - 100k_sampling_summary.json (processing summary)"
echo "==========================================="