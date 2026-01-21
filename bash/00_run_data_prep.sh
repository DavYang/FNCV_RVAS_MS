#!/bin/bash
set -e

# --- Configuration ---
SCRIPT_NAME="00_data_prep"
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

echo "==========================================="
echo "   Starting Data Prep Pipeline"
echo "==========================================="
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

# Run in background with nohup
# > redirects standard output to log
# 2>&1 redirects errors to the same log
# & puts it in the background
nohup python3 "${PYTHON_DIR}/${SCRIPT_NAME}.py" > "${LOG_FILE}" 2>&1 &

PID=$!
echo "Job submitted successfully!"
echo "PID: ${PID}"
echo ""
echo "To monitor progress, run:"
echo "  tail -f ${LOG_FILE}"
echo "==========================================="
