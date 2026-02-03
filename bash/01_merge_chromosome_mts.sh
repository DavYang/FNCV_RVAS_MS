#!/bin/bash
set -e

# --- Configuration ---
MERGE_SCRIPT_NAME="merge_chromosome_mts"
BASE_DIR=$(pwd)
PYTHON_DIR="${BASE_DIR}/python"
LOG_DIR="${BASE_DIR}/logs"

# --- Setup ---
# Ensure directories exist
mkdir -p "$LOG_DIR"
mkdir -p "$PYTHON_DIR"

# Generate a timestamp for the log file
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="${LOG_DIR}/${MERGE_SCRIPT_NAME}_${TIMESTAMP}.log"

echo "==========================================="
echo "   Starting Chromosome MatrixTable Merge"
echo "==========================================="
echo "Date: ${TIMESTAMP}"
echo "Script: ${PYTHON_DIR}/${MERGE_SCRIPT_NAME}.py"
echo "Log File: ${LOG_FILE}"
echo "==========================================="

# --- Execution ---
# Check if the merge script exists
if [ ! -f "${PYTHON_DIR}/${MERGE_SCRIPT_NAME}.py" ]; then
    echo "Error: Merge script not found at ${PYTHON_DIR}/${MERGE_SCRIPT_NAME}.py"
    exit 1
fi

# Run merge script in background
echo "Starting chromosome MatrixTable merge..."
nohup python3 "${PYTHON_DIR}/${MERGE_SCRIPT_NAME}.py" > "${LOG_FILE}" 2>&1 &
MERGE_PID=$!

echo "Merge job submitted!"
echo "PID: ${MERGE_PID}"
echo ""
echo "To monitor progress, run:"
echo "  tail -f ${LOG_FILE}"
echo ""
echo "To check job status, run:"
echo "  ps -p ${MERGE_PID}"
echo "==========================================="
