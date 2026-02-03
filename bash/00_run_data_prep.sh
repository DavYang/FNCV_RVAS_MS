#!/bin/bash
set -e

# --- Configuration ---
SCRIPT_NAME="00_background_SNPs_mt"
MERGE_SCRIPT_NAME="merge_chromosome_mts"
BASE_DIR=$(pwd)
PYTHON_DIR="${BASE_DIR}/python"
LOG_DIR="${BASE_DIR}/logs"
CONFIG_DIR="${BASE_DIR}/config"

# --- Setup ---
# Ensure directories exist
mkdir -p "$LOG_DIR"
mkdir -p "$PYTHON_DIR"
mkdir -p "$CONFIG_DIR"

# Generate a timestamp for the log file
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="${LOG_DIR}/${SCRIPT_NAME}_${TIMESTAMP}.log"
MERGE_LOG_FILE="${LOG_DIR}/${MERGE_SCRIPT_NAME}_${TIMESTAMP}.log"

echo "==========================================="
echo "   Starting Background SNP Sampling Pipeline"
echo "==========================================="
echo "Date: ${TIMESTAMP}"
echo "Step 1 Script: ${PYTHON_DIR}/${SCRIPT_NAME}.py"
echo "Step 1 Log: ${LOG_FILE}"
echo "Step 2 Script: ${PYTHON_DIR}/${MERGE_SCRIPT_NAME}.py"
echo "Step 2 Log: ${MERGE_LOG_FILE}"
echo "Config: ${CONFIG_DIR}/config.json"
echo "==========================================="

# --- Step 1: Process Individual Chromosomes ---
echo "🧬 Step 1: Processing individual chromosomes..."

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

# Check if merge script exists
if [ ! -f "${PYTHON_DIR}/${MERGE_SCRIPT_NAME}.py" ]; then
    echo "Error: Merge script not found at ${PYTHON_DIR}/${MERGE_SCRIPT_NAME}.py"
    exit 1
fi

# Run chromosome processing
echo "Starting chromosome processing..."
python3 "${PYTHON_DIR}/${SCRIPT_NAME}.py" > "${LOG_FILE}" 2>&1
CHROMOSOME_EXIT_CODE=$?

if [ $CHROMOSOME_EXIT_CODE -ne 0 ]; then
    echo "❌ Step 1 failed! Check log: ${LOG_FILE}"
    echo "Last 20 lines of error:"
    tail -20 "${LOG_FILE}"
    exit 1
fi

echo "✅ Step 1 completed successfully!"
echo "📊 Chromosome processing summary:"
grep "Successfully processed" "${LOG_FILE}" || echo "Check log for details"

# --- Step 2: Merge Chromosome Results ---
echo ""
echo "🔗 Step 2: Merging chromosome results..."

# Run merge script
python3 "${PYTHON_DIR}/${MERGE_SCRIPT_NAME}.py" > "${MERGE_LOG_FILE}" 2>&1
MERGE_EXIT_CODE=$?

if [ $MERGE_EXIT_CODE -ne 0 ]; then
    echo "❌ Step 2 failed! Check log: ${MERGE_LOG_FILE}"
    echo "Last 20 lines of error:"
    tail -20 "${MERGE_LOG_FILE}"
    exit 1
fi

echo "✅ Step 2 completed successfully!"
echo "📊 Merge summary:"
grep "Final dataset:" "${MERGE_LOG_FILE}" || echo "Check log for details"

# --- Completion ---
echo ""
echo "==========================================="
echo "🎉 Pipeline completed successfully!"
echo "==========================================="
echo "Step 1 Log: ${LOG_FILE}"
echo "Step 2 Log: ${MERGE_LOG_FILE}"
echo ""
echo "To review results:"
echo "  - Check chromosome processing: tail -50 ${LOG_FILE}"
echo "  - Check merge results: tail -50 ${MERGE_LOG_FILE}"
echo "  - Find output files in: results/FNCV_RVAS_MS/background_snps/"
echo "==========================================="
