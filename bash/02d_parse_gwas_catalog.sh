#!/bin/bash
set -eo pipefail

# ---------------------------------------------------------------------------
# 02d_parse_gwas_catalog.sh
#
# Parse the full GWAS Catalog associations TSV to extract EUR-only MS loci
# and liftover positions from GRCh37 to GRCh38.
#
# Runs on the HPC (not the AoU workbench). All paths are read from
# config/config.json:
#   inputs.gwas_catalog_raw  -- full catalog TSV
#   inputs.liftover_chain    -- hg19ToHg38.over.chain.gz (downloaded if absent)
#   outputs.locus_def_dir    -- where gwas_catalog_ms_eur_hg38.tsv is written
#   params.gwas_catalog_gcs_path -- GCS destination shown in the upload reminder
#
# Prerequisites:
#   - pip install pyliftover pandas
#   - inputs.gwas_catalog_raw must exist (already downloaded)
#
# Usage (run from project root):
#   bash bash/02d_parse_gwas_catalog.sh
#   bash bash/02d_parse_gwas_catalog.sh --force
#
# Background:
#   nohup bash bash/02d_parse_gwas_catalog.sh > /dev/null 2>&1 &
#   tail -f logs/02d_parse_gwas_catalog_*.log
#
# After completion, upload the output to GCS with:
#   gsutil cp results/2-locus_definition/gwas_catalog_ms_eur_hg38.tsv \
#     $(python3 -c "import json; c=json.load(open('config/config.json')); \
#       print(c['params']['gwas_catalog_gcs_path'])")
# ---------------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
PYTHON_SCRIPT="${PROJECT_DIR}/python/02d_parse_gwas_catalog.py"
CONFIG_FILE="${PROJECT_DIR}/config/config.json"
LOG_DIR="${PROJECT_DIR}/logs"

mkdir -p "$LOG_DIR"

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="${LOG_DIR}/02d_parse_gwas_catalog_${TIMESTAMP}.log"

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
EXTRA_ARGS=()

for arg in "$@"; do
    case "$arg" in
        --force)
            EXTRA_ARGS+=("--force")
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
log "  02d: Parse GWAS Catalog -- EUR MS Loci + Liftover"
log "============================================================"
log "  Config  : ${CONFIG_FILE}"
log "  Log     : ${LOG_FILE}"
log "  Args    : ${EXTRA_ARGS[*]:-none}"
log ""

# ---------------------------------------------------------------------------
# Pre-flight checks
# ---------------------------------------------------------------------------
if [ ! -f "${PYTHON_SCRIPT}" ]; then
    log "ERROR: Python script not found: ${PYTHON_SCRIPT}"
    log "  Expected: ${PYTHON_SCRIPT}"
    exit 1
fi

if [ ! -f "${CONFIG_FILE}" ]; then
    log "ERROR: Config not found: ${CONFIG_FILE}"
    log "  Run from project root: bash bash/02d_parse_gwas_catalog.sh"
    exit 1
fi

# Verify catalog input exists
CATALOG_PATH=$(python3 -c "
import json
c = json.load(open('${CONFIG_FILE}'))
print(c.get('inputs', {}).get('gwas_catalog_raw', ''))
")

if [ -z "${CATALOG_PATH}" ]; then
    log "ERROR: inputs.gwas_catalog_raw not set in config.json"
    exit 1
fi

if [ ! -f "${CATALOG_PATH}" ]; then
    log "ERROR: GWAS Catalog file not found: ${CATALOG_PATH}"
    log "  Set inputs.gwas_catalog_raw in config.json to the correct path."
    exit 1
fi

log "Pre-flight checks PASSED"
log "  Catalog : ${CATALOG_PATH}"
log ""

# ---------------------------------------------------------------------------
# Run Python script
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
# Print upload reminder
# ---------------------------------------------------------------------------
LOCUS_DEF_DIR=$(python3 -c "
import json
c = json.load(open('${CONFIG_FILE}'))
print(c.get('outputs', {}).get('locus_def_dir', 'results/2-locus_definition'))
")
GCS_DEST=$(python3 -c "
import json
c = json.load(open('${CONFIG_FILE}'))
print(c.get('params', {}).get('gwas_catalog_gcs_path', ''))
")

OUTPUT_FILE="${LOCUS_DEF_DIR}/gwas_catalog_ms_eur_hg38.tsv"

log ""
log "02d complete."
log "  Output : ${OUTPUT_FILE}"
log ""
if [ -n "${GCS_DEST}" ]; then
    log "Upload to AoU workspace bucket with:"
    log "  gsutil cp ${OUTPUT_FILE} ${GCS_DEST}"
else
    log "WARNING: params.gwas_catalog_gcs_path not set in config.json -- set it before uploading."
fi
