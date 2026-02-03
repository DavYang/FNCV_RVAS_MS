#!/bin/bash
set -e

# --- Configuration ---
SCRIPT_NAME="00_background_SNPs_mt"
BASE_DIR=$(pwd)
PYTHON_DIR="${BASE_DIR}/python"
LOG_DIR="${BASE_DIR}/logs"
CONFIG_DIR="${BASE_DIR}/config"
ERROR_LOG_DIR="${BASE_DIR}/error_logs"

# --- Setup ---
# Ensure directories exist
mkdir -p "$LOG_DIR"
mkdir -p "$PYTHON_DIR"
mkdir -p "$CONFIG_DIR"
mkdir -p "$ERROR_LOG_DIR"

# Generate a timestamp for the log file
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="${LOG_DIR}/${SCRIPT_NAME}_${TIMESTAMP}.log"
ERROR_LOG_FILE="${ERROR_LOG_DIR}/${SCRIPT_NAME}_${TIMESTAMP}_errors.log"

echo "==========================================="
echo "   Starting Background SNP Sampling Pipeline"
echo "==========================================="
echo "Date: ${TIMESTAMP}"
echo "Script: ${PYTHON_DIR}/${SCRIPT_NAME}.py"
echo "Log File: ${LOG_FILE}"
echo "Error Log: ${ERROR_LOG_FILE}"
echo "Config: ${CONFIG_DIR}/config.json"
echo "==========================================="

# --- Step 1: Process Individual Chromosomes ---
echo "Step 1: Processing individual chromosomes..."

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

# Function to monitor progress
monitor_progress() {
    local pid=$1
    local log_file=$2
    local error_log=$3
    
    echo "Monitoring progress for PID: $pid"
    echo "Log file: $log_file"
    echo "Error log: $error_log"
    echo ""
    
    # Monitor the process
    while kill -0 $pid 2>/dev/null; do
        # Show last 5 lines of main log
        echo "--- Recent Progress ($(date)) ---"
        tail -5 "$log_file" 2>/dev/null || echo "Log file not accessible yet"
        echo ""
        
        # Check for errors in the last minute
        if [ -f "$log_file" ]; then
            # Extract recent errors (last minute)
            recent_errors=$(tail -100 "$log_file" | grep -i "error\|failed\|exception" | tail -3)
            if [ ! -z "$recent_errors" ]; then
                echo "!!! RECENT ERRORS DETECTED !!!"
                echo "$recent_errors"
                echo ""
                # Also save to error log
                echo "$(date): $recent_errors" >> "$error_log"
            fi
        fi
        
        sleep 30  # Check every 30 seconds
    done
    
    # Process finished, check exit status
    wait $pid
    exit_code=$?
    
    if [ $exit_code -eq 0 ]; then
        echo "✅ Chromosome processing completed successfully!"
    else
        echo "❌ Chromosome processing failed with exit code: $exit_code"
        echo "Last 10 lines of log:"
        tail -10 "$log_file"
        echo ""
        echo "Error details saved to: $error_log"
    fi
    
    return $exit_code
}

# Run chromosome processing in background
echo "Starting chromosome processing..."
nohup python3 "${PYTHON_DIR}/${SCRIPT_NAME}.py" > "${LOG_FILE}" 2>&1 &
CHROMOSOME_PID=$!

echo "Chromosome processing job submitted!"
echo "PID: ${CHROMOSOME_PID}"
echo ""
echo "To monitor progress manually, run:"
echo "  tail -f ${LOG_FILE}"
echo ""
echo "To check job status, run:"
echo "  ps -p ${CHROMOSOME_PID}"
echo ""
echo "Automated monitoring started..."
echo ""

# Start monitoring in background
monitor_progress $CHROMOSOME_PID "$LOG_FILE" "$ERROR_LOG_FILE" &
MONITOR_PID=$!

# Wait for monitoring to complete
wait $MONITOR_PID
CHROMOSOME_EXIT_CODE=$?

echo ""
echo "==========================================="

if [ $CHROMOSOME_EXIT_CODE -eq 0 ]; then
    echo "✅ STEP 1 COMPLETED SUCCESSFULLY"
    echo ""
    echo "After chromosome processing completes, run:"
    echo "  bash bash/01_merge_chromosome_mts.sh"
    echo ""
    echo "Check results in:"
    echo "  - Summary: ${LOG_FILE}"
    echo "  - Errors: ${ERROR_LOG_FILE}"
else
    echo "❌ STEP 1 FAILED"
    echo ""
    echo "Troubleshooting:"
    echo "  - Check main log: ${LOG_FILE}"
    echo "  - Check error log: ${ERROR_LOG_FILE}"
    echo "  - Look for 'SessionTimeoutError' or 'Connection refused' messages"
fi

echo "==========================================="

exit $CHROMOSOME_EXIT_CODE
