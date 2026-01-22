#!/bin/bash
set -e

# Ensure we are running from the project root
if [[ ! -d "python" || ! -d "config" ]]; then
    echo "Error: Must be run from the project root directory (containing python/ and config/)."
    exit 1
fi

SCRIPT_NAME="01_define_loci"
BASE_DIR=$(pwd)
PYTHON_DIR="${BASE_DIR}/python"
LOG_DIR="${BASE_DIR}/logs"
TOOLS_DIR="${BASE_DIR}/tools"

mkdir -p "$LOG_DIR"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="${LOG_DIR}/${SCRIPT_NAME}_${TIMESTAMP}.log"

echo "=== Starting Locus Definition & Signal Dissection ==="
echo "Log: ${LOG_FILE}"

# 1. Install Dependencies (GCTA)
echo "Checking dependencies..."
python3 "${PYTHON_DIR}/install_gcta.py"

# Verify GCTA was installed
if [[ ! -f "${TOOLS_DIR}/gcta64" ]]; then
    echo "Error: GCTA binary not found at ${TOOLS_DIR}/gcta64"
    exit 1
fi

# 2. Run Main Analysis Script
echo "Launching Python analysis script..."
nohup python3 "${PYTHON_DIR}/${SCRIPT_NAME}.py" > "${LOG_FILE}" 2>&1 &

PID=$!
echo "Job PID: ${PID}"
echo "Monitor: tail -f ${LOG_FILE}"