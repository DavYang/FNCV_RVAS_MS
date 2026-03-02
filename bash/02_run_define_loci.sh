#!/bin/bash
set -e

# ---------------------------------------------------------------------------
# Phase 2a: Lead SNP extraction + GCTA-COJO signal dissection
# ---------------------------------------------------------------------------
# Usage:
#   nohup bash bash/01_run_define_loci.sh > /dev/null 2>&1 &
#
# Monitor:
#   tail -f logs/01_define_loci_<timestamp>.log
#
# Input:  GWAS HT (GCS), LD reference PLINK (local Phase 1 output)
# Output: results/2-locus_definition/target_loci.bed
#         results/2-locus_definition/all_independent_signals.tsv
#         results/2-locus_definition/cojo/{locus_id}.jma.cojo
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Ensure we are running from the project root
# ---------------------------------------------------------------------------
if [[ ! -d "python" || ! -d "config" ]]; then
    echo "Error: Must be run from the project root directory (containing python/ and config/)."
    exit 1
fi

BASE_DIR=$(pwd)
PYTHON_DIR="${BASE_DIR}/python"
LOG_DIR="${BASE_DIR}/logs"
TOOLS_DIR="${BASE_DIR}/tools"
CONFIG_FILE="${BASE_DIR}/config/config.json"

mkdir -p "${LOG_DIR}"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="${LOG_DIR}/01_define_loci_${TIMESTAMP}.log"

echo "============================================================"
echo "  Phase 2a: Locus Definition (GCTA-COJO)"
echo "============================================================"
echo "  Config  : ${CONFIG_FILE}"
echo "  Log     : ${LOG_FILE}"
echo "  Started : $(date)"
echo "============================================================"

# ---------------------------------------------------------------------------
# Step 1: Install / verify GCTA
# ---------------------------------------------------------------------------
echo "Checking GCTA dependency ..."
python3 "${PYTHON_DIR}/install_gcta.py"

if [[ ! -f "${TOOLS_DIR}/gcta64" ]]; then
    echo "Error: GCTA binary not found at ${TOOLS_DIR}/gcta64"
    exit 1
fi

GCTA_VER=$("${TOOLS_DIR}/gcta64" --version 2>&1 | head -1 || true)
echo "  GCTA: ${GCTA_VER}"

# ---------------------------------------------------------------------------
# Step 2: Launch Python pipeline
# ---------------------------------------------------------------------------
echo ""
echo "Launching 01_define_loci.py ..."
nohup python3 "${PYTHON_DIR}/01_define_loci.py" \
    --config "${CONFIG_FILE}" \
    > "${LOG_FILE}" 2>&1 &

PID=$!
echo "  PID     : ${PID}"
echo "  Monitor : tail -f ${LOG_FILE}"
echo ""
echo "Job running in background. To check status:"
echo "  tail -f ${LOG_FILE}"