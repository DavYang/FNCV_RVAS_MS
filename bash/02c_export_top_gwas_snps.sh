#!/bin/bash
set -eo pipefail

# ---------------------------------------------------------------------------
# 02c_export_top_gwas_snps.sh
#
# Phase 2 Step 2: Filter per-chromosome .ma files (from 02a) to SNPs reaching
# the significance threshold (default p < 5e-6) and write a consolidated TSV
# (top_gwas_snps.tsv). Cross-referenced against published EUR-only MS GWAS
# Catalog loci in the windowed locus definition step (02e).
#
# Outputs:
#   results/2-locus_definition/top_gwas_snps.tsv  (local)
#   $WORKSPACE_BUCKET/results/2-locus_definition/top_gwas_snps.tsv  (GCS)
#
# No Hail/Spark required -- pure pandas, runs on the VM in seconds.
#
# Usage:
#   bash bash/02c_export_top_gwas_snps.sh
#   bash bash/02c_export_top_gwas_snps.sh --p-threshold=5e-6
#   bash bash/02c_export_top_gwas_snps.sh --force
#
# Monitor:
#   tail -f logs/02c_top_gwas_snps_*.log
# ---------------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
PYTHON_SCRIPT="${PROJECT_DIR}/python/02c_export_top_gwas_snps.py"
CONFIG_FILE="${PROJECT_DIR}/config/config.json"
LOG_DIR="${PROJECT_DIR}/logs"

mkdir -p "$LOG_DIR"

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="${LOG_DIR}/02c_top_gwas_snps_${TIMESTAMP}.log"

# ---------------------------------------------------------------------------
# Parse arguments (pass through to Python script)
# ---------------------------------------------------------------------------
EXTRA_ARGS=()
FORCE=""

for arg in "$@"; do
    case "$arg" in
        --force)
            FORCE="--force"
            EXTRA_ARGS+=("$arg")
            ;;
        --p-threshold=*)
            EXTRA_ARGS+=("--p-threshold" "${arg#*=}")
            ;;
        --output=*)
            EXTRA_ARGS+=("--output" "${arg#*=}")
            ;;
        --ma-dir=*)
            EXTRA_ARGS+=("--ma-dir" "${arg#*=}")
            ;;
        *)
            EXTRA_ARGS+=("$arg")
            ;;
    esac
done

# ---------------------------------------------------------------------------
# Logging helper
# ---------------------------------------------------------------------------
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "${LOG_FILE}"
}

log "============================================================"
log "  02c: Export Top GWAS SNPs"
log "============================================================"
log "  Config : ${CONFIG_FILE}"
log "  Log    : ${LOG_FILE}"
log "  Args   : ${EXTRA_ARGS[*]}"
log ""

# ---------------------------------------------------------------------------
# Pre-flight checks
# ---------------------------------------------------------------------------
if [ ! -f "${PYTHON_SCRIPT}" ]; then
    log "ERROR: Python script not found: ${PYTHON_SCRIPT}"
    exit 1
fi

if [ ! -f "${CONFIG_FILE}" ]; then
    log "ERROR: Config not found: ${CONFIG_FILE}"
    exit 1
fi

# Check that .ma files exist (spot-check chr1)
MA_DIR=$(python3 -c "
import json
c = json.load(open('${CONFIG_FILE}'))
print(c['outputs'].get('gwas_ma_dir', 'results/2-locus_definition/ma'))
")

if [ ! -f "${MA_DIR}/chr1.ma" ]; then
    log "ERROR: .ma files not found in ${MA_DIR}"
    log "  Run: bash bash/02a_export_gwas_ma.sh"
    exit 1
fi

log "Pre-flight checks PASSED"
log ""

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
log "Running: python3 ${PYTHON_SCRIPT} --config ${CONFIG_FILE} ${EXTRA_ARGS[*]}"
log ""

python3 "${PYTHON_SCRIPT}" \
    --config "${CONFIG_FILE}" \
    "${EXTRA_ARGS[@]}" 2>&1 | tee -a "${LOG_FILE}"

EXIT_CODE=$?

if [ "${EXIT_CODE}" -ne 0 ]; then
    log "ERROR: Python script exited with code ${EXIT_CODE}"
    exit "${EXIT_CODE}"
fi

# ---------------------------------------------------------------------------
# Upload top_gwas_snps.tsv to GCS
# ---------------------------------------------------------------------------
TOP_SNPS_FILE=$(python3 -c "
import json, os
c = json.load(open('${CONFIG_FILE}'))
locus_def_dir = c['outputs'].get('locus_def_dir', 'results/2-locus_definition')
print(os.path.join(locus_def_dir, 'top_gwas_snps.tsv'))
")

# Resolve GCS destination from config, expanding \$WORKSPACE_BUCKET
GCS_DEST_TEMPLATE=$(python3 -c "
import json
c = json.load(open('${CONFIG_FILE}'))
print(c.get('params', {}).get('top_gwas_snps_gcs_path', ''))
")
GCS_DEST="${GCS_DEST_TEMPLATE//\$WORKSPACE_BUCKET/${WORKSPACE_BUCKET}}"

log ""
log "============================================================"
log "  Upload to GCS"
log "============================================================"

if [ -z "${GCS_DEST}" ] || [ -z "${WORKSPACE_BUCKET}" ]; then
    log "WARNING: WORKSPACE_BUCKET not set or top_gwas_snps_gcs_path not configured."
    log "  Set WORKSPACE_BUCKET and params.top_gwas_snps_gcs_path in config/config.json."
    log "  To upload manually:"
    log "  gsutil cp ${TOP_SNPS_FILE} \$WORKSPACE_BUCKET/results/2-locus_definition/top_gwas_snps.tsv"
else
    log "Uploading: ${TOP_SNPS_FILE}"
    log "       to: ${GCS_DEST}"
    gsutil cp "${TOP_SNPS_FILE}" "${GCS_DEST}" 2>&1 | tee -a "${LOG_FILE}"
    log "Upload complete."
fi

log ""
log "02c complete"
log "  Output : ${TOP_SNPS_FILE}"
log "  Next   : bash bash/02e_define_loci_catalog.sh"
