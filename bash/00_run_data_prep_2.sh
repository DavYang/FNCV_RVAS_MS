#!/bin/bash
set -e

# --- Configuration ---
SCRIPT_NAME="00_data_prep_2"
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
LOG_FILE="${LOG_DIR}/${SCRIPT_NAME}_${TIMESTAMP}.log"

# Export timestamp for Python script
export EXPORT_TIMESTAMP="${TIMESTAMP}"

echo "==========================================="
echo "   Starting Data Prep Pipeline (Optimized)"
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

echo "Checking for existing checkpoint will happen inside the script."

# Run in background with nohup
nohup python3 "${PYTHON_DIR}/${SCRIPT_NAME}.py" > "${LOG_FILE}" 2>&1 &

PID=$!
echo "Job submitted successfully!"
echo "PID: ${PID}"
echo ""
echo "To monitor progress, run:"
echo "  tail -f ${LOG_FILE}"
echo "==========================================="