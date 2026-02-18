#!/bin/bash
set -e

# ---------------------------------------------------------------------------
# Copy per-chromosome PLINK files from GCS into a single local directory
# ---------------------------------------------------------------------------
# Usage:
#   nohup bash bash/copy_plink_from_gcs.sh > /dev/null 2>&1 &
#
# Monitor:
#   tail -f logs/copy_plink_from_gcs.log
#
# Outputs all .bed/.bim/.fam files into:
#   results/1-bg_snp/plink_no-qc/
# ---------------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "${SCRIPT_DIR}")"

GCS_BASE="gs://fc-secure-b43840eb-548f-464d-bece-31ac7a969abd/results/FNCV_RVAS_MS/background_snps_20260218"
LOCAL_DIR="/home/jupyter/FNCV_RVAS_MS/results/1-bg_snp/plink_no-qc"
LOG_DIR="${PROJECT_DIR}/logs"
LOG_FILE="${LOG_DIR}/copy_plink_from_gcs.log"

mkdir -p "${LOCAL_DIR}" "${LOG_DIR}"

log() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $*" | tee -a "${LOG_FILE}"
}

log "============================================================"
log "Copying PLINK files from GCS to ${LOCAL_DIR}"
log "============================================================"
log ""

FAILED=()
for chr_num in $(seq 1 22); do
    chr_name="chr${chr_num}"
    log "[${chr_name}] Copying .bed/.bim/.fam ..."
    if gsutil -m cp \
        "${GCS_BASE}/${chr_name}/${chr_name}_background.bed" \
        "${GCS_BASE}/${chr_name}/${chr_name}_background.bim" \
        "${GCS_BASE}/${chr_name}/${chr_name}_background.fam" \
        "${LOCAL_DIR}/" >> "${LOG_FILE}" 2>&1; then
        log "[${chr_name}] Done"
    else
        log "[${chr_name}] FAILED"
        FAILED+=("${chr_name}")
    fi
done

log ""
log "============================================================"
if [ ${#FAILED[@]} -eq 0 ]; then
    log "All PLINK files copied successfully to ${LOCAL_DIR}"
else
    log "FAILED chromosomes: ${FAILED[*]}"
fi
log "============================================================"
log "File count: $(ls -1 "${LOCAL_DIR}"/*.{bed,bim,fam} 2>/dev/null | wc -l) files"
ls -lh "${LOCAL_DIR}/" >> "${LOG_FILE}" 2>&1
