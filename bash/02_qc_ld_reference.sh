#!/bin/bash
set -eo pipefail

# ---------------------------------------------------------------------------
# 02_qc_ld_reference.sh
#
# Phase 2 Step 0: Build per-chromosome LD reference panels for GCTA-COJO
# from the raw (no-QC) PLINK files exported by Hail in Phase 1.
#
# QC Pipeline per chromosome (PLINK2):
#   1. Variant QC: --maf 0.01, --hwe 1e-6, --geno 0.05, --mind 0.05,
#      --snps-only just-acgt, --max-alleles 2
#   2. Remove duplicate variant IDs (--rm-dup exclude-all)
#
# Samples are already EUR-only (234K) from Hail export in 01a.
# No ancestry PCA or relatedness filtering needed.
#
# Inputs  (local on VM):
#   results/1-bg_snp/plink_no-qc/chrN_background.{bed,bim,fam}
#
# Outputs (local, then uploaded to GCS):
#   results/2-locus_definition/ld_ref/chrN_ld_ref.{bed,bim,fam}
#
# Usage:
#   nohup bash bash/02_qc_ld_reference.sh > /dev/null 2>&1 &
#
#   # Force re-run (overwrites existing LD ref files)
#   bash bash/02_qc_ld_reference.sh --force
#
#   # Test mode: process only test chromosome (chr21) for iterative testing
#   bash bash/02_qc_ld_reference.sh --test
#
# Monitor:
#   tail -f logs/02_qc_ld_ref_*.log
#   cat logs/02_qc_ld_ref.pid
# ---------------------------------------------------------------------------

trap '' HUP

PROJECT_DIR="/home/jupyter/FNCV_RVAS_MS"
CONFIG_FILE="${PROJECT_DIR}/config/config.json"

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
INPUT_DIR="${PROJECT_DIR}/results/1-bg_snp/plink_no-qc"
OUTPUT_DIR="${PROJECT_DIR}/results/2-locus_definition/ld_ref"
TMP_DIR="${PROJECT_DIR}/tmp/ld_ref_qc"
LOG_DIR="${PROJECT_DIR}/logs"

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="${LOG_DIR}/02_qc_ld_ref_${TIMESTAMP}.log"
PID_FILE="${LOG_DIR}/02_qc_ld_ref.pid"

# ---------------------------------------------------------------------------
# QC thresholds (read from config.json plink_qc for cross-step consistency)
# ---------------------------------------------------------------------------
MAF=$(python3 -c "import json; print(json.load(open('${CONFIG_FILE}'))['plink_qc']['maf'])")
HWE=$(python3 -c "import json; print(json.load(open('${CONFIG_FILE}'))['plink_qc']['hwe'])")
GENO=$(python3 -c "import json; print(json.load(open('${CONFIG_FILE}'))['plink_qc']['geno'])")
MIND=$(python3 -c "import json; print(json.load(open('${CONFIG_FILE}'))['plink_qc']['mind'])")
MAX_ALLELES=$(python3 -c "import json; print(json.load(open('${CONFIG_FILE}'))['plink_qc']['max_alleles'])")
SNPS_ONLY=$(python3 -c "import json; print(json.load(open('${CONFIG_FILE}'))['plink_qc']['snps_only'])")

# ---------------------------------------------------------------------------
# Resources
# ---------------------------------------------------------------------------
THREADS=8
MEMORY=16000

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
FORCE=0
TEST_MODE=0
TEST_CHR=$(python3 -c "import json; print(json.load(open('${CONFIG_FILE}'))['params']['test_chromosome'])")
for arg in "$@"; do
    case "$arg" in
        --force) FORCE=1 ;;
        --test)  TEST_MODE=1 ;;
        *) echo "Unknown argument: $arg"; exit 1 ;;
    esac
done

# ---------------------------------------------------------------------------
# GCS location for raw PLINK (crash recovery download)
# ---------------------------------------------------------------------------
GCS_PLINK_BASE="gs://fc-secure-b43840eb-548f-464d-bece-31ac7a969abd/results/FNCV_RVAS_MS/background_snps_20260226"

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
mkdir -p "${OUTPUT_DIR}" "${TMP_DIR}" "${LOG_DIR}"
echo $$ > "${PID_FILE}"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "${LOG_FILE}"
}

START_SECONDS=$(date +%s)

log "============================================================"
log "  02_qc_ld_reference.sh -- LD Reference QC for GCTA-COJO"
log "============================================================"
log "  Timestamp : ${TIMESTAMP}"
log "  PID       : $$"
log "  Force     : ${FORCE}"
log "  Input     : ${INPUT_DIR}"
log "  Output    : ${OUTPUT_DIR}"
log "  QC params : MAF=${MAF} HWE=${HWE} GENO=${GENO} MIND=${MIND}"
log "  SNPs-only : ${SNPS_ONLY}  Max-alleles: ${MAX_ALLELES}"
log "  Test mode : ${TEST_MODE} (${TEST_CHR})"
log "  Threads   : ${THREADS}  Memory: ${MEMORY} MB"
log ""

# ---------------------------------------------------------------------------
# Per-chromosome QC loop
# ---------------------------------------------------------------------------
FAILED=()
TOTAL_INPUT=0
TOTAL_QC1=0
TOTAL_FINAL=0

if [ "${TEST_MODE}" -eq 1 ]; then
    CHR_LIST=$(echo "${TEST_CHR}" | sed 's/chr//')
    log "TEST MODE: processing only ${TEST_CHR}"
else
    CHR_LIST=$(seq 1 22)
fi

for chr_num in ${CHR_LIST}; do
    chr_name="chr${chr_num}"
    INPUT_PREFIX="${INPUT_DIR}/${chr_name}_background"
    QC1_PREFIX="${TMP_DIR}/${chr_name}_qc1"
    FINAL_PREFIX="${OUTPUT_DIR}/${chr_name}_ld_ref"

    # Resume: skip if final output already exists
    if [ -f "${FINAL_PREFIX}.bed" ] && [ "${FORCE}" -eq 0 ]; then
        FINAL_COUNT=$(wc -l < "${FINAL_PREFIX}.bim")
        TOTAL_FINAL=$((TOTAL_FINAL + FINAL_COUNT))
        log "[${chr_name}] SKIP (ld_ref exists) -- ${FINAL_COUNT} variants"
        continue
    fi

    # Download from GCS if local input is missing
    if [ ! -f "${INPUT_PREFIX}.bed" ]; then
        log "[${chr_name}] Local input missing, downloading from GCS..."
        mkdir -p "${INPUT_DIR}"
        if gsutil -m cp \
            "${GCS_PLINK_BASE}/${chr_name}/${chr_name}_background.bed" \
            "${GCS_PLINK_BASE}/${chr_name}/${chr_name}_background.bim" \
            "${GCS_PLINK_BASE}/${chr_name}/${chr_name}_background.fam" \
            "${INPUT_DIR}/" 2>&1; then
            log "[${chr_name}] Downloaded from GCS"
        else
            log "[${chr_name}] SKIPPED -- GCS download failed"
            FAILED+=("${chr_name}")
            continue
        fi
    fi

    # Verify input exists
    if [ ! -f "${INPUT_PREFIX}.bed" ]; then
        log "[${chr_name}] SKIPPED -- input .bed not found: ${INPUT_PREFIX}.bed"
        FAILED+=("${chr_name}")
        continue
    fi

    INPUT_COUNT=$(wc -l < "${INPUT_PREFIX}.bim")
    TOTAL_INPUT=$((TOTAL_INPUT + INPUT_COUNT))

    # --- Step 1: Variant QC ---
    log "[${chr_name}] Step 1: Variant QC (${INPUT_COUNT} input variants)..."

    if ! plink2 \
        --bfile "${INPUT_PREFIX}" \
        --maf "${MAF}" \
        --hwe "${HWE}" \
        --geno "${GENO}" \
        --mind "${MIND}" \
        --snps-only "${SNPS_ONLY}" \
        --max-alleles "${MAX_ALLELES}" \
        --make-bed \
        --out "${QC1_PREFIX}" \
        --threads "${THREADS}" \
        --memory "${MEMORY}" 2>&1 | tee -a "${LOG_FILE}"; then
        log "[${chr_name}] Step 1 FAILED"
        FAILED+=("${chr_name}")
        continue
    fi

    QC1_COUNT=$(wc -l < "${QC1_PREFIX}.bim")
    REMOVED1=$((INPUT_COUNT - QC1_COUNT))
    TOTAL_QC1=$((TOTAL_QC1 + QC1_COUNT))
    log "[${chr_name}] Step 1 done: ${QC1_COUNT} variants (${REMOVED1} removed)"

    # --- Step 2: Remove duplicate variant IDs ---
    log "[${chr_name}] Step 2: Removing duplicate variant IDs..."

    if ! plink2 \
        --bfile "${QC1_PREFIX}" \
        --rm-dup exclude-all \
        --make-bed \
        --out "${FINAL_PREFIX}" \
        --threads "${THREADS}" \
        --memory "${MEMORY}" 2>&1 | tee -a "${LOG_FILE}"; then
        log "[${chr_name}] Step 2 FAILED"
        FAILED+=("${chr_name}")
        continue
    fi

    FINAL_COUNT=$(wc -l < "${FINAL_PREFIX}.bim")
    DUPS_REMOVED=$((QC1_COUNT - FINAL_COUNT))
    TOTAL_FINAL=$((TOTAL_FINAL + FINAL_COUNT))
    log "[${chr_name}] Step 2 done: ${FINAL_COUNT} variants (${DUPS_REMOVED} duplicates removed)"

    # Cleanup intermediate QC1 files
    rm -f "${QC1_PREFIX}".{bed,bim,fam,log,nosex}

    log "[${chr_name}] COMPLETE: ${INPUT_COUNT} -> ${FINAL_COUNT} variants"
    log ""
done

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
log ""
log "============================================================"
log "  Per-chromosome summary"
log "============================================================"
log "Total input variants   : ${TOTAL_INPUT}"
log "Total after QC step 1  : ${TOTAL_QC1}"
log "Total final (LD ref)   : ${TOTAL_FINAL}"

if [ ${#FAILED[@]} -gt 0 ]; then
    log "FAILED chromosomes: ${FAILED[*]}"
    log "Fix failures and re-run."
    exit 1
fi

# ---------------------------------------------------------------------------
# Cleanup temp directory
# ---------------------------------------------------------------------------
log ""
log "Cleaning up temp files..."
rm -rf "${TMP_DIR}"
log "Done."

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
    GCS_DEST="${WORKSPACE_BUCKET}/results/2-locus_definition/ld_ref/"
    log "Uploading to: ${GCS_DEST}"
    gsutil -u "${GOOGLE_PROJECT}" -m cp "${OUTPUT_DIR}"/*.{bed,bim,fam} "${GCS_DEST}"
    log "Upload complete"
else
    log "WARNING: WORKSPACE_BUCKET not set -- skipping GCS upload."
    log "  To upload manually:"
    log "  gsutil -u \$GOOGLE_PROJECT -m cp ${OUTPUT_DIR}/*.{bed,bim,fam} \$WORKSPACE_BUCKET/results/2-locus_definition/ld_ref/"
fi

# ---------------------------------------------------------------------------
# Final summary
# ---------------------------------------------------------------------------
END_SECONDS=$(date +%s)
ELAPSED=$(( END_SECONDS - START_SECONDS ))
ELAPSED_FMT=$(printf '%02dh %02dm %02ds' $((ELAPSED/3600)) $(( (ELAPSED%3600)/60 )) $((ELAPSED%60)))

log ""
log "============================================================"
log "  Pipeline COMPLETED in ${ELAPSED_FMT}"
log "============================================================"
log "LD ref files : ${OUTPUT_DIR}/chrN_ld_ref.{bed,bim,fam}"
log "Total final  : ${TOTAL_FINAL} variants"
log ""
for chr_num in ${CHR_LIST}; do
    f="${OUTPUT_DIR}/chr${chr_num}_ld_ref.bim"
    if [ -f "$f" ]; then
        n=$(wc -l < "$f")
        log "  chr${chr_num}: ${n} variants"
    fi
done
log ""
log "Next steps:"
log "  bash bash/02a_export_gwas_ma.sh   (export .ma files)"
log "  bash bash/02b_run_cojo.sh         (GCTA-COJO)"
log "============================================================"

rm -f "${PID_FILE}"
