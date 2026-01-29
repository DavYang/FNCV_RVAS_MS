#!/bin/bash


BASE_DIR=$(pwd)
PYTHON_DIR="${BASE_DIR}/python"
LOG_DIR="${BASE_DIR}/logs"

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="${LOG_DIR}/test_${TIMESTAMP}.log"

echo "Checking for existing checkpoint will happen inside the script."

# Run in background with nohup
nohup python3 "${PYTHON_DIR}/test.py" > "${LOG_FILE}" 2>&1 &

PID=$!
echo "Job submitted successfully!"
echo "PID: ${PID}"
echo ""
echo "To monitor progress, run:"
echo "  tail -f ${LOG_FILE}"