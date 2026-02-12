#!/bin/bash
set -eo pipefail

# ---------------------------------------------------------------------------
# 01_build_background_snps.sh
#
# Wrapper script for the background SNP PLINK pipeline (Regenie Step 1).
# Runs a single Hail Python script that:
#   Pass 1: Samples loci from ACAF splitMT rows Table (no genotypes).
#   Pass 2: Extracts EUR genotypes at sampled loci and exports per-chromosome
#           PLINK files via checkpoint + hl.export_plink.
#
# The script is designed to run in the background and survive terminal
# disconnects. All output goes to a timestamped log file.
#
# Usage:
#   bash bash/01_build_background_snps.sh &
#
# Monitor:
#   tail -f logs/01_background_snps_*.log
#   cat logs/01_background_snps.pid   # check PID
#   ps -p $(cat logs/01_background_snps.pid) -o pid,etime,cmd  # check status
#   kill $(cat logs/01_background_snps.pid)  # stop if needed
# ---------------------------------------------------------------------------

# Ignore hangup signals so the script survives terminal disconnects
trap '' HUP

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
PYTHON_SCRIPT="${PROJECT_DIR}/python/01_build_background_snps.py"
CONFIG_FILE="${PROJECT_DIR}/config/config.json"
LOG_DIR="${PROJECT_DIR}/logs"

mkdir -p "$LOG_DIR"

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="${LOG_DIR}/01_background_snps_${TIMESTAMP}.log"
PID_FILE="${LOG_DIR}/01_background_snps.pid"

# Redirect all output to log file AND stdout, detach stdin
exec < /dev/null
exec > >(tee -a "${LOG_FILE}") 2>&1

# Write PID file for monitoring
echo $$ > "${PID_FILE}"

# Record start time for elapsed calculation
START_SECONDS=$(date +%s)

# ---------------------------------------------------------------------------
# Pre-flight checks
# ---------------------------------------------------------------------------
echo "============================================================"
echo "Background SNP PLINK Pipeline"
echo "============================================================"
echo "Started at : $(date)"
echo "PID        : $$"
echo "PID file   : ${PID_FILE}"
echo "Log file   : ${LOG_FILE}"
echo "Hail log   : /tmp/hail_background_snps.log"
echo ""

if [ -z "$WORKSPACE_BUCKET" ]; then
    echo "ERROR: WORKSPACE_BUCKET environment variable not set"
    exit 1
fi
echo "Workspace  : $WORKSPACE_BUCKET"

if [ -z "$GOOGLE_PROJECT" ]; then
    echo "WARNING: GOOGLE_PROJECT not set. Requester Pays buckets may fail."
else
    echo "Project    : $GOOGLE_PROJECT"
fi

if [ ! -f "$PYTHON_SCRIPT" ]; then
    echo "ERROR: Python script not found: $PYTHON_SCRIPT"
    exit 1
fi

if [ ! -f "$CONFIG_FILE" ]; then
    echo "ERROR: Config file not found: $CONFIG_FILE"
    exit 1
fi

# ---------------------------------------------------------------------------
# Read config for display
# ---------------------------------------------------------------------------
TARGET_SNPS=$(python3 -c "import json; c=json.load(open('${CONFIG_FILE}')); print(c['sampling']['target_total_snps'])")
TEST_MODE=$(python3 -c "import json; c=json.load(open('${CONFIG_FILE}')); print(c['params'].get('test_mode', False))")
TEST_CHR=$(python3 -c "import json; c=json.load(open('${CONFIG_FILE}')); print(c['params'].get('test_chromosome', 'None'))")

echo ""
echo "Target SNPs: ${TARGET_SNPS}"
if [ "$TEST_MODE" = "True" ]; then
    echo "Mode       : TEST (${TEST_CHR} only)"
else
    echo "Mode       : PRODUCTION (chr1-22)"
fi
echo "------------------------------------------------------------"
echo ""

# ---------------------------------------------------------------------------
# Run pipeline
# ---------------------------------------------------------------------------
cd "$PROJECT_DIR"

EXIT_CODE=0
python3 "$PYTHON_SCRIPT" || EXIT_CODE=$?

END_SECONDS=$(date +%s)
ELAPSED=$((END_SECONDS - START_SECONDS))
ELAPSED_MIN=$((ELAPSED / 60))
ELAPSED_SEC=$((ELAPSED % 60))

echo ""
echo "============================================================"
if [ $EXIT_CODE -eq 0 ]; then
    echo "Pipeline SUCCEEDED"
else
    echo "Pipeline FAILED (exit code: ${EXIT_CODE})"
    echo ""
    echo "--- Last 50 lines of Hail log ---"
    tail -50 /tmp/hail_background_snps.log 2>/dev/null || echo "(Hail log not found)"
    echo "--- End Hail log ---"
fi
echo "Elapsed    : ${ELAPSED_MIN}m ${ELAPSED_SEC}s"
echo "Finished at: $(date)"
echo "============================================================"

# Clean up PID file
rm -f "${PID_FILE}"

exit $EXIT_CODE
