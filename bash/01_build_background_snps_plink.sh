#!/bin/bash
set -e

# ---------------------------------------------------------------------------
# 01_build_background_snps_plink.sh
#
# Two-pass pipeline to build background SNP PLINK files for Regenie Step 1.
# Pass 1: Sample loci from ACAF splitMT rows Table (no entry data).
# Pass 2: Extract EUR genotypes at sampled loci, export per-chr PLINK.
# ---------------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
PYTHON_SCRIPT="${PROJECT_DIR}/python/01_build_background_snps_plink.py"
CONFIG_FILE="${PROJECT_DIR}/config/config.json"
LOG_DIR="${PROJECT_DIR}/logs"

# Create log directory
mkdir -p "$LOG_DIR"

# Timestamp for log file
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="${LOG_DIR}/01_background_snps_plink_${TIMESTAMP}.log"

# ---------------------------------------------------------------------------
# Pre-flight checks
# ---------------------------------------------------------------------------
echo "============================================================"
echo "Background SNP PLINK Pipeline"
echo "============================================================"

# Check workspace bucket
if [ -z "$WORKSPACE_BUCKET" ]; then
    echo "ERROR: WORKSPACE_BUCKET environment variable not set"
    exit 1
fi
echo "Workspace bucket: $WORKSPACE_BUCKET"

# Check Python script exists
if [ ! -f "$PYTHON_SCRIPT" ]; then
    echo "ERROR: Python script not found: $PYTHON_SCRIPT"
    exit 1
fi

# Check config exists
if [ ! -f "$CONFIG_FILE" ]; then
    echo "ERROR: Config file not found: $CONFIG_FILE"
    exit 1
fi

# Read target SNPs from config
TARGET_SNPS=$(python3 -c "import json; c=json.load(open('${CONFIG_FILE}')); print(c['sampling']['target_total_snps'])")
echo "Target SNPs: ${TARGET_SNPS}"

# Check test mode
TEST_MODE=$(python3 -c "import json; c=json.load(open('${CONFIG_FILE}')); print(c['params'].get('test_mode', False))")
TEST_CHR=$(python3 -c "import json; c=json.load(open('${CONFIG_FILE}')); print(c['params'].get('test_chromosome', 'None'))")

if [ "$TEST_MODE" = "True" ]; then
    echo "Mode: TEST (${TEST_CHR} only)"
else
    echo "Mode: PRODUCTION (chr1-22)"
fi

# Format SNP label for output paths
if [ "$TARGET_SNPS" -ge 1000000 ]; then
    SNP_LABEL="$((TARGET_SNPS / 1000000))M"
elif [ "$TARGET_SNPS" -ge 1000 ]; then
    SNP_LABEL="$((TARGET_SNPS / 1000))K"
else
    SNP_LABEL="$TARGET_SNPS"
fi

echo "------------------------------------------------------------"
echo "Output structure:"
echo "  PLINK: \${WORKSPACE_BUCKET}/results/FNCV_RVAS_MS/${SNP_LABEL}_background_snps_*/null_model_data/"
echo "  Loci:  \${WORKSPACE_BUCKET}/results/FNCV_RVAS_MS/${SNP_LABEL}_background_snps_*/sampled_loci.ht"
echo "  Log:   ${LOG_FILE}"
echo "------------------------------------------------------------"
echo ""

# ---------------------------------------------------------------------------
# Run the pipeline
# ---------------------------------------------------------------------------
echo "Starting pipeline at $(date)"
echo "Log file: ${LOG_FILE}"
echo ""

cd "$PROJECT_DIR"
nohup python3 "$PYTHON_SCRIPT" > "$LOG_FILE" 2>&1 &
PID=$!

echo "Pipeline running in background (PID: $PID)"
echo ""
echo "Monitor progress:"
echo "  tail -f ${LOG_FILE}"
echo ""
echo "Check if running:"
echo "  ps -p $PID"
echo ""
echo "============================================================"
