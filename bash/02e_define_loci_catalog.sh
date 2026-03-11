#!/bin/bash
set -eo pipefail

# ---------------------------------------------------------------------------
# 02e_define_loci_catalog.sh
#
# Cross-reference AoU top GWAS SNPs against EUR-only MS GWAS Catalog loci
# (GRCh38, produced by 02d_parse_gwas_catalog.sh on HPC), generate +-250kb
# windows around validated SNPs, merge overlapping windows, exclude MHC,
# and export target_loci.bed.
#
# All paths are read from config/config.json:
#   outputs.locus_def_dir          -- base output directory
#   params.gwas_catalog_hg38_path  -- local path to catalog file
#   params.gwas_catalog_gcs_path   -- GCS URI to download catalog from if absent
#
# Prerequisites:
#   - results/2-locus_definition/top_gwas_snps.tsv  (from 02c)
#   - gwas_catalog_ms_eur_hg38.tsv uploaded to $WORKSPACE_BUCKET/reference/
#     by running 02d_parse_gwas_catalog.sh on HPC
#
# Usage (run from project root):
#   bash bash/02e_define_loci_catalog.sh
#   bash bash/02e_define_loci_catalog.sh --force
#
# Background (recommended):
#   nohup bash bash/02e_define_loci_catalog.sh > /dev/null 2>&1 &
#   tail -f logs/02e_define_loci_catalog_*.log
# ---------------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
PYTHON_SCRIPT="${PROJECT_DIR}/python/02e_define_loci_catalog.py"
CONFIG_FILE="${PROJECT_DIR}/config/config.json"
LOG_DIR="${PROJECT_DIR}/logs"

mkdir -p "$LOG_DIR"

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="${LOG_DIR}/02e_define_loci_catalog_${TIMESTAMP}.log"

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
FORCE=""
EXTRA_ARGS=()

for arg in "$@"; do
    case "$arg" in
        --force)
            FORCE="--force"
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
log "  02e: Define Loci via GWAS Catalog Cross-Reference"
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
    exit 1
fi

if [ ! -f "${CONFIG_FILE}" ]; then
    log "ERROR: Config not found: ${CONFIG_FILE}"
    log "  Run from project root: bash bash/02e_define_loci_catalog.sh"
    exit 1
fi

# Resolve locus_def_dir and catalog path from config
LOCUS_DEF_DIR=$(python3 -c "
import json
c = json.load(open('${CONFIG_FILE}'))
print(c.get('outputs', {}).get('locus_def_dir', 'results/2-locus_definition'))
")

CATALOG_LOCAL=$(python3 -c "
import json, os
c = json.load(open('${CONFIG_FILE}'))
path = c.get('params', {}).get('gwas_catalog_hg38_path', '')
print(path)
")

# Fall back to locus_def_dir if config path is empty
if [ -z "${CATALOG_LOCAL}" ]; then
    CATALOG_LOCAL="${LOCUS_DEF_DIR}/gwas_catalog_ms_eur_hg38.tsv"
fi

TOP_SNPS_FILE="${LOCUS_DEF_DIR}/top_gwas_snps.tsv"

log "Pre-flight checks PASSED"
log "  Top SNPs : ${TOP_SNPS_FILE}"
log "  Catalog  : ${CATALOG_LOCAL}"
log ""

# ---------------------------------------------------------------------------
# Download top_gwas_snps.tsv from GCS if not already local
# ---------------------------------------------------------------------------
if [ ! -f "${TOP_SNPS_FILE}" ] || [ -n "${FORCE}" ]; then
    TOP_SNPS_GCS_TEMPLATE=$(python3 -c "
import json
c = json.load(open('${CONFIG_FILE}'))
print(c.get('params', {}).get('top_gwas_snps_gcs_path', ''))
" 2>/dev/null)
    TOP_SNPS_GCS="${TOP_SNPS_GCS_TEMPLATE//\$WORKSPACE_BUCKET/${WORKSPACE_BUCKET}}"

    if [ -z "${TOP_SNPS_GCS}" ] || [ -z "${WORKSPACE_BUCKET}" ]; then
        log "ERROR: top_gwas_snps.tsv not found locally and params.top_gwas_snps_gcs_path not set."
        log "  Run: bash bash/02c_export_top_gwas_snps.sh"
        exit 1
    fi

    log "Downloading top_gwas_snps.tsv from GCS: ${TOP_SNPS_GCS}"
    mkdir -p "${LOCUS_DEF_DIR}"
    gsutil cp "${TOP_SNPS_GCS}" "${TOP_SNPS_FILE}" 2>&1 | tee -a "${LOG_FILE}"

    if [ ! -f "${TOP_SNPS_FILE}" ]; then
        log "ERROR: Download failed -- file not found after gsutil cp: ${TOP_SNPS_FILE}"
        exit 1
    fi
    log "Downloaded: ${TOP_SNPS_FILE}"
else
    log "top_gwas_snps.tsv already local: ${TOP_SNPS_FILE} (use --force to re-download)"
fi

log ""

# ---------------------------------------------------------------------------
# Download GWAS Catalog from GCS if not already local
# ---------------------------------------------------------------------------
if [ ! -f "${CATALOG_LOCAL}" ] || [ -n "${FORCE}" ]; then
    CATALOG_GCS_TEMPLATE=$(python3 -c "
import json
c = json.load(open('${CONFIG_FILE}'))
print(c.get('params', {}).get('gwas_catalog_gcs_path', ''))
" 2>/dev/null)
    CATALOG_GCS="${CATALOG_GCS_TEMPLATE//\$WORKSPACE_BUCKET/${WORKSPACE_BUCKET}}"

    if [ -z "${CATALOG_GCS}" ]; then
        log "ERROR: Catalog not found locally and params.gwas_catalog_gcs_path is not set."
        log "  Run bash/02d_parse_gwas_catalog.sh on HPC, upload the output, then set:"
        log "    params.gwas_catalog_gcs_path in config/config.json"
        exit 1
    fi

    log "Downloading GWAS Catalog from GCS: ${CATALOG_GCS}"
    mkdir -p "$(dirname "${CATALOG_LOCAL}")"
    gsutil cp "${CATALOG_GCS}" "${CATALOG_LOCAL}" 2>&1 | tee -a "${LOG_FILE}"

    if [ ! -f "${CATALOG_LOCAL}" ]; then
        log "ERROR: Download failed -- file not found after gsutil cp: ${CATALOG_LOCAL}"
        exit 1
    fi
    log "Downloaded: ${CATALOG_LOCAL}"
else
    log "Catalog already local: ${CATALOG_LOCAL} (use --force to re-download)"
fi

log ""

# ---------------------------------------------------------------------------
# Run Python script
# ---------------------------------------------------------------------------
log "Running: python3 ${PYTHON_SCRIPT} --config ${CONFIG_FILE} --catalog ${CATALOG_LOCAL} ${EXTRA_ARGS[*]}"
log ""

python3 "${PYTHON_SCRIPT}" \
    --config "${CONFIG_FILE}" \
    --catalog "${CATALOG_LOCAL}" \
    "${EXTRA_ARGS[@]}" 2>&1 | tee -a "${LOG_FILE}"

EXIT_CODE=$?

if [ "${EXIT_CODE}" -ne 0 ]; then
    log "ERROR: Python script exited with code ${EXIT_CODE}"
    exit "${EXIT_CODE}"
fi

# ---------------------------------------------------------------------------
# Upload outputs to GCS
# ---------------------------------------------------------------------------
WORKSPACE_BUCKET=$(python3 -c "
import os
print(os.environ.get('WORKSPACE_BUCKET', ''))
" 2>/dev/null)

if [ -z "${WORKSPACE_BUCKET}" ]; then
    log "WARNING: WORKSPACE_BUCKET not set -- skipping GCS upload"
else
    GCS_OUT="${WORKSPACE_BUCKET}/results/2-locus_definition"
    log ""
    log "Uploading outputs to GCS: ${GCS_OUT}"

    for f in \
        "${LOCUS_DEF_DIR}/gwas_catalog_validated_snps.tsv" \
        "${LOCUS_DEF_DIR}/novel_snps.tsv" \
        "${LOCUS_DEF_DIR}/target_loci.bed" \
        "${LOCUS_DEF_DIR}/mhc_excluded_windows.tsv"
    do
        if [ -f "${f}" ]; then
            gsutil cp "${f}" "${GCS_OUT}/$(basename "${f}")" 2>&1 | tee -a "${LOG_FILE}"
            log "  Uploaded: $(basename "${f}")"
        fi
    done

    log "GCS upload complete."
fi

log ""
log "02e complete."
log "  Outputs:"
log "    ${LOCUS_DEF_DIR}/gwas_catalog_validated_snps.tsv"
log "    ${LOCUS_DEF_DIR}/novel_snps.tsv"
log "    ${LOCUS_DEF_DIR}/target_loci.bed"
