#!/bin/bash
set -eo pipefail

# ---------------------------------------------------------------------------
# 02a_export_gwas_ma.sh
#
# Phase 2 Step 1: Export per-chromosome GCTA-COJO .ma summary statistics
# from the AllxAll GWAS Hail Table. Runs Hail on the VM (/opt/conda env).
#
# GCTA-COJO requires full-chromosome .ma files (all SNPs) to correctly
# estimate phenotypic variance. This script wraps python/02a_export_gwas_ma.py.
#
# Inputs  (GCS, read via Hail):
#   AllxAll GWAS HT (config: inputs.phenotype_gwas)
#
# Outputs (local, then uploaded to GCS):
#   results/2-locus_definition/ma/chrN.ma  (22 files)
#
# Usage:
#   nohup bash bash/02a_export_gwas_ma.sh > /dev/null 2>&1 &
#
#   # Force re-export (overwrite existing .ma files)
#   bash bash/02a_export_gwas_ma.sh --force
#
#   # Single chromosome (use = syntax)
#   bash bash/02a_export_gwas_ma.sh --chrom=chr21
#
#   # Test mode: process only test chromosome (chr21) for iterative testing
#   bash bash/02a_export_gwas_ma.sh --test
#
# Monitor:
#   tail -f logs/02a_export_gwas_ma_*.log
#   cat logs/02a_export_gwas_ma.pid
# ---------------------------------------------------------------------------

trap '' HUP

PROJECT_DIR="/home/jupyter/FNCV_RVAS_MS"

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
CONFIG_FILE="${PROJECT_DIR}/config/config.json"
PYTHON_SCRIPT="${PROJECT_DIR}/python/02a_export_gwas_ma.py"
OUTPUT_DIR="${PROJECT_DIR}/results/2-locus_definition/ma"
LOG_DIR="${PROJECT_DIR}/logs"

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="${LOG_DIR}/02a_export_gwas_ma_${TIMESTAMP}.log"
PID_FILE="${LOG_DIR}/02a_export_gwas_ma.pid"

# ---------------------------------------------------------------------------
# Parse arguments (pass through to Python)
# ---------------------------------------------------------------------------
EXTRA_ARGS=()
FORCE=0
TEST_MODE=0
TEST_CHR=$(python3 -c "import json; print(json.load(open('${CONFIG_FILE}'))['params']['test_chromosome'])")
for arg in "$@"; do
    case "$arg" in
        --force) FORCE=1; EXTRA_ARGS+=("--force") ;;
        --test)  TEST_MODE=1; EXTRA_ARGS+=("--chrom" "${TEST_CHR}") ;;
        --chrom=*) EXTRA_ARGS+=("--chrom" "${arg#--chrom=}") ;;
        *) echo "WARNING: Unknown argument '$arg' -- ignoring" ;;
    esac
done

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}"
echo $$ > "${PID_FILE}"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "${LOG_FILE}"
}

START_SECONDS=$(date +%s)

log "============================================================"
log "  02a_export_gwas_ma.sh -- Export .ma files for GCTA-COJO"
log "============================================================"
log "  Timestamp : ${TIMESTAMP}"
log "  PID       : $$"
log "  Config    : ${CONFIG_FILE}"
log "  Output    : ${OUTPUT_DIR}"
log "  Force     : ${FORCE}"
log "  Test mode : ${TEST_MODE} (${TEST_CHR})"
log ""

# ---------------------------------------------------------------------------
# Pre-flight checks
# ---------------------------------------------------------------------------
if [ ! -f "${CONFIG_FILE}" ]; then
    log "ERROR: Config file not found: ${CONFIG_FILE}"
    exit 1
fi

if [ ! -f "${PYTHON_SCRIPT}" ]; then
    log "ERROR: Python script not found: ${PYTHON_SCRIPT}"
    exit 1
fi

# Check Hail is available
if ! python3 -c "import hail" 2>/dev/null; then
    log "ERROR: Hail not available in current Python environment."
    log "  Expected at /opt/conda. Ensure correct environment is active."
    exit 1
fi

# Check GOOGLE_PROJECT is set (needed for requester-pays GCS access)
if [ -z "${GOOGLE_PROJECT}" ]; then
    log "WARNING: GOOGLE_PROJECT not set. Hail may fail on requester-pays buckets."
fi

# ---------------------------------------------------------------------------
# Run Python script
# ---------------------------------------------------------------------------
log "Running: python3 ${PYTHON_SCRIPT} --config ${CONFIG_FILE} --output-dir ${OUTPUT_DIR} ${EXTRA_ARGS[*]}"
log ""

python3 "${PYTHON_SCRIPT}" \
    --config "${CONFIG_FILE}" \
    --output-dir "${OUTPUT_DIR}" \
    "${EXTRA_ARGS[@]}" 2>&1 | tee -a "${LOG_FILE}"

PYTHON_EXIT=$?

if [ "${PYTHON_EXIT}" -ne 0 ]; then
    log ""
    log "ERROR: Python script exited with code ${PYTHON_EXIT}"
    rm -f "${PID_FILE}"
    exit "${PYTHON_EXIT}"
fi

# ---------------------------------------------------------------------------
# Upload to GCS
# ---------------------------------------------------------------------------
log ""
log "============================================================"
log "  Upload to GCS"
log "============================================================"
if [ -n "${WORKSPACE_BUCKET}" ]; then
    if [[ "${WORKSPACE_BUCKET}" != gs://* ]]; then
        WORKSPACE_BUCKET="gs://${WORKSPACE_BUCKET}"
    fi
    GCS_DEST="${WORKSPACE_BUCKET}/results/2-locus_definition/ma/"
    log "Uploading to: ${GCS_DEST}"
    gsutil -u "${GOOGLE_PROJECT}" -m cp "${OUTPUT_DIR}"/*.ma "${GCS_DEST}"
    log "Upload complete"
else
    log "WARNING: WORKSPACE_BUCKET not set -- skipping GCS upload."
    log "  To upload manually:"
    log "  gsutil -u \$GOOGLE_PROJECT -m cp ${OUTPUT_DIR}/*.ma \$WORKSPACE_BUCKET/results/2-locus_definition/ma/"
fi

# ---------------------------------------------------------------------------
# Final summary
# ---------------------------------------------------------------------------
END_SECONDS=$(date +%s)
ELAPSED=$(( END_SECONDS - START_SECONDS ))
ELAPSED_FMT=$(printf '%02dh %02dm %02ds' $((ELAPSED/3600)) $(( (ELAPSED%3600)/60 )) $((ELAPSED%60)))

log ""
log "============================================================"
log "  02a COMPLETED in ${ELAPSED_FMT}"
log "============================================================"
log ".ma files: ${OUTPUT_DIR}/chrN.ma"
log ""
for f in "${OUTPUT_DIR}"/*.ma; do
    if [ -f "$f" ]; then
        n=$(( $(wc -l < "$f") - 1 ))
        log "  $(basename "$f"): ${n} variants"
    fi
done
log ""
log "Next step: bash bash/02b_run_cojo.sh"
log "============================================================"

rm -f "${PID_FILE}"
