#!/bin/bash
set -e

# --- Configuration ---
SCRIPT_NAME="test_mt_stats"
BASE_DIR=$(pwd)
PYTHON_DIR="${BASE_DIR}/python"
LOG_DIR="${BASE_DIR}/logs"

# Ensure log directory exists
mkdir -p "$LOG_DIR"

# Generate a timestamp
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="${LOG_DIR}/${SCRIPT_NAME}_${TIMESTAMP}.log"

echo "==========================================="
echo "   Running MatrixTable Verification & Export Test"
echo "==========================================="
echo "Date: ${TIMESTAMP}"
echo "Script: ${PYTHON_DIR}/${SCRIPT_NAME}.py"
echo "Log File: ${LOG_FILE}"
echo "Workspace Bucket: ${WORKSPACE_BUCKET}"
echo "==========================================="

if [ -z "${WORKSPACE_BUCKET}" ]; then
    echo "ERROR: WORKSPACE_BUCKET environment variable is not set."
    exit 1
fi

# Run the python script (in foreground, since it's a quick test)
python3 "${PYTHON_DIR}/${SCRIPT_NAME}.py" | tee "${LOG_FILE}"

echo "==========================================="
echo "Done."