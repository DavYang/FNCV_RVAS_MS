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
# Usage (run with nohup to survive disconnects):
#   nohup bash bash/01_build_background_snps.sh > /dev/null 2>&1 &
#   tail -f logs/01_background_snps_*.log
# ---------------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
PYTHON_SCRIPT="${PROJECT_DIR}/python/01_build_background_snps.py"
CONFIG_FILE="${PROJECT_DIR}/config/config.json"
LOG_DIR="${PROJECT_DIR}/logs"

mkdir -p "$LOG_DIR"

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="${LOG_DIR}/01_background_snps_${TIMESTAMP}.log"

exec > >(tee -a "${LOG_FILE}") 2>&1

# ---------------------------------------------------------------------------
# Pre-flight checks
# ---------------------------------------------------------------------------
echo "============================================================"
echo "Background SNP PLINK Pipeline"
echo "============================================================"
echo "Started at: $(date)"
echo ""

if [ -z "$WORKSPACE_BUCKET" ]; then
    echo "ERROR: WORKSPACE_BUCKET environment variable not set"
    exit 1
fi
echo "Workspace bucket: $WORKSPACE_BUCKET"

if [ -z "$GOOGLE_PROJECT" ]; then
    echo "WARNING: GOOGLE_PROJECT not set. Requester Pays buckets may fail."
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

echo "Target SNPs: ${TARGET_SNPS}"
if [ "$TEST_MODE" = "True" ]; then
    echo "Mode: TEST (${TEST_CHR} only)"
else
    echo "Mode: PRODUCTION (chr1-22)"
fi
echo "Log file: ${LOG_FILE}"
echo "------------------------------------------------------------"
echo ""

# ---------------------------------------------------------------------------
# Run pipeline
# ---------------------------------------------------------------------------
cd "$PROJECT_DIR"

if ! python3 "$PYTHON_SCRIPT"; then
    echo ""
    echo "ERROR: Pipeline failed with exit code $?"
    echo "Check Hail log: /tmp/hail_background_snps.log"
    echo "Check script log: ${LOG_FILE}"
    exit 1
fi

echo ""
echo "Finished at: $(date)"
echo "============================================================"
