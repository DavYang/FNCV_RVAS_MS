#!/bin/bash
set -e

SCRIPT_NAME="01_define_loci"
BASE_DIR=$(pwd)
PYTHON_DIR="${BASE_DIR}/python"
LOG_DIR="${BASE_DIR}/logs"

mkdir -p "$LOG_DIR"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="${LOG_DIR}/${SCRIPT_NAME}_${TIMESTAMP}.log"

echo "=== Starting Locus Definition ==="
echo "Log: ${LOG_FILE}"

# First check/install GCTA if needed
python3 "${PYTHON_DIR}/install_gcta.py"

# Run main script
nohup python3 "${PYTHON_DIR}/${SCRIPT_NAME}.py" > "${LOG_FILE}" 2>&1 &

PID=$!
echo "Job PID: ${PID}"
echo "Monitor: tail -f ${LOG_FILE}"
