#!/bin/bash
set -e

# ---------------------------------------------------------------------------
# QC and merge per-chromosome PLINK files for downstream use.
# ---------------------------------------------------------------------------
# Pipeline per chromosome:
#   1. Variant QC  (--maf, --geno, --hwe, --mind, --snps-only, --max-alleles 2)
#   2. Remove duplicate variant IDs (--rm-dup exclude-all)
#
# After all 22 chromosomes are QCed:
#   3. Merge QC'd per-chrom files into a genome-wide fileset
#   4. Thin down to TARGET_SNPS if total exceeds it (for REGENIE Step 1)
#
# The per-chromosome QC'd files also serve as the LD reference panels
# for GCTA-COJO in Phase 2, making a separate LD-ref QC step unnecessary.
# ---------------------------------------------------------------------------
# Usage:
#   nohup bash bash/01b_qc_merge_background_snps.sh > /dev/null 2>&1 &
#
#   # Force re-run (deletes QC'd + merged outputs)
#   bash bash/01b_qc_merge_background_snps.sh --force
#
# Monitor:
#   tail -f logs/01b_qc_merge_*.log
#
# Input:  results/1-bg_snp/plink_no-qc/chrN_background.{bed,bim,fam}
# Output:
#   Per-chrom QC'd: results/1-bg_snp/plink_qc/chrN_background_qc.{bed,bim,fam}
#   Step 1 merged:  results/1-bg_snp/plink_step1/step1_500k.{bed,bim,fam}
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_DIR="/home/jupyter/FNCV_RVAS_MS"
CONFIG_FILE="${PROJECT_DIR}/config/config.json"

INPUT_DIR="${PROJECT_DIR}/results/1-bg_snp/plink_no-qc"
QC_DIR="${PROJECT_DIR}/results/1-bg_snp/plink_qc"
STEP1_DIR="${PROJECT_DIR}/results/1-bg_snp/plink_step1"
LOG_DIR="${PROJECT_DIR}/logs"

# GCS location of raw (no-QC) per-chromosome PLINK files from Hail export
GCS_PLINK_BASE="gs://fc-secure-b43840eb-548f-464d-bece-31ac7a969abd/results/FNCV_RVAS_MS/background_snps_20260226"

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="${LOG_DIR}/01b_qc_merge_${TIMESTAMP}.log"
MERGE_LIST="${QC_DIR}/merge_list.txt"

# ---------------------------------------------------------------------------
# QC thresholds (read from config.json plink_qc section for cross-step consistency)
# ---------------------------------------------------------------------------
MAF=$(python3 -c "import json; print(json.load(open('${CONFIG_FILE}'))['plink_qc']['maf'])")
HWE=$(python3 -c "import json; print(json.load(open('${CONFIG_FILE}'))['plink_qc']['hwe'])")
GENO=$(python3 -c "import json; print(json.load(open('${CONFIG_FILE}'))['plink_qc']['geno'])")
MIND=$(python3 -c "import json; print(json.load(open('${CONFIG_FILE}'))['plink_qc']['mind'])")
MAX_ALLELES=$(python3 -c "import json; print(json.load(open('${CONFIG_FILE}'))['plink_qc']['max_alleles'])")
SNPS_ONLY=$(python3 -c "import json; print(json.load(open('${CONFIG_FILE}'))['plink_qc']['snps_only'])")
HWE_SAMPLE_TERM=0        # explicit sample-size term to suppress PLINK2 warning

# ---------------------------------------------------------------------------
# Final target (for REGENIE Step 1 thinning)
# ---------------------------------------------------------------------------
TARGET_SNPS=500000

RANDOM_SEED=$(python3 -c "import json; c=json.load(open('${CONFIG_FILE}')); print(c.get('sampling',{}).get('random_seed', c.get('params',{}).get('random_seed', 42)))")


THREADS=50
MEMORY=16000

# ---------------------------------------------------------------------------
# Parse flags
# ---------------------------------------------------------------------------
FORCE=0
for arg in "$@"; do
    case "$arg" in
        --force) FORCE=1 ;;
        *) echo "WARNING: Unknown argument '$arg' -- ignoring" ;;
    esac
done

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
mkdir -p "${QC_DIR}" "${STEP1_DIR}" "${LOG_DIR}"

exec < /dev/null
exec > >(tee -a "${LOG_FILE}") 2>&1

START_SECONDS=$(date +%s)

log() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $*"
}

log "============================================================"
log "  QC + Merge Background SNPs"
log "============================================================"
log "Input dir    : ${INPUT_DIR}"
log "QC dir       : ${QC_DIR}"
log "Step1 dir    : ${STEP1_DIR}"
log "QC params    : MAF=${MAF}, GENO=${GENO}, HWE=${HWE} (st=${HWE_SAMPLE_TERM}), MIND=${MIND}, snps-only=${SNPS_ONLY}, max-alleles=${MAX_ALLELES}"
log "Target SNPs  : ${TARGET_SNPS}"
log "Force rerun  : ${FORCE}"
log "============================================================"
log ""

# ---------------------------------------------------------------------------
# If --force, clean QC'd and merged outputs for full re-run
# ---------------------------------------------------------------------------
if [ "${FORCE}" -eq 1 ]; then
    log "[--force] Cleaning QC'd and merged outputs..."
    rm -f "${QC_DIR}"/*.{bed,bim,fam,log,nosex} 2>/dev/null || true
    rm -f "${STEP1_DIR}"/step1_*.{bed,bim,fam,log} 2>/dev/null || true
    rm -f "${MERGE_LIST}" 2>/dev/null || true
fi

# ---------------------------------------------------------------------------
# Step 1+2: Per-chromosome variant QC + duplicate removal
# ---------------------------------------------------------------------------
FAILED=()
TOTAL_PRE=0
TOTAL_QC=0

rm -f "${MERGE_LIST}"

for chr_num in $(seq 1 22); do
    chr_name="chr${chr_num}"
    INPUT_PREFIX="${INPUT_DIR}/${chr_name}_background"
    QC_PREFIX="${QC_DIR}/${chr_name}_background_qc"

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
            log "[${chr_name}] SKIPPED - GCS download failed"
            FAILED+=("${chr_name}")
            continue
        fi
    fi

    # Verify input exists
    if [ ! -f "${INPUT_PREFIX}.bed" ]; then
        log "[${chr_name}] SKIPPED - input .bed not found: ${INPUT_PREFIX}.bed"
        FAILED+=("${chr_name}")
        continue
    fi

    PRE_COUNT=$(wc -l < "${INPUT_PREFIX}.bim")
    TOTAL_PRE=$((TOTAL_PRE + PRE_COUNT))

    # Resume: skip if QC output already exists
    if [ -f "${QC_PREFIX}.bed" ] && [ "${FORCE}" -eq 0 ]; then
        QC_COUNT=$(wc -l < "${QC_PREFIX}.bim")
        TOTAL_QC=$((TOTAL_QC + QC_COUNT))
        log "[${chr_name}] SKIP (QC exists) - ${QC_COUNT} variants"
        echo "${QC_PREFIX}" >> "${MERGE_LIST}"
        continue
    fi

    # --- Step 1: Variant QC ---
    log "[${chr_name}] QC starting - ${PRE_COUNT} variants pre-QC"
    QC_TMP="${QC_DIR}/${chr_name}_qc_tmp"

    if ! plink2 \
        --bfile "${INPUT_PREFIX}" \
        --maf "${MAF}" \
        --geno "${GENO}" \
        --hwe "${HWE}" "${HWE_SAMPLE_TERM}" \
        --mind "${MIND}" \
        --snps-only "${SNPS_ONLY}" \
        --max-alleles "${MAX_ALLELES}" \
        --make-bed \
        --out "${QC_TMP}" \
        --threads "${THREADS}" \
        --memory "${MEMORY}" 2>&1; then
        log "[${chr_name}] QC FAILED"
        FAILED+=("${chr_name}")
        continue
    fi

    QC_TMP_COUNT=$(wc -l < "${QC_TMP}.bim")
    REMOVED=$((PRE_COUNT - QC_TMP_COUNT))
    log "[${chr_name}] QC done - ${QC_TMP_COUNT} variants (${REMOVED} removed)"

    # --- Step 2: Remove duplicate variant IDs ---
    log "[${chr_name}] Removing duplicate variant IDs..."

    if ! plink2 \
        --bfile "${QC_TMP}" \
        --rm-dup exclude-all \
        --make-bed \
        --out "${QC_PREFIX}" \
        --threads "${THREADS}" \
        --memory "${MEMORY}" 2>&1; then
        log "[${chr_name}] DEDUP FAILED"
        FAILED+=("${chr_name}")
        continue
    fi

    QC_COUNT=$(wc -l < "${QC_PREFIX}.bim")
    DUPS_REMOVED=$((QC_TMP_COUNT - QC_COUNT))
    TOTAL_QC=$((TOTAL_QC + QC_COUNT))
    log "[${chr_name}] Final: ${QC_COUNT} variants (${DUPS_REMOVED} duplicates removed)"
    echo "${QC_PREFIX}" >> "${MERGE_LIST}"

    # Cleanup intermediate tmp files
    rm -f "${QC_TMP}".{bed,bim,fam,log,nosex}
done

log ""
log "============================================================"
log "  Per-chromosome summary"
log "============================================================"
log "Total pre-QC variants  : ${TOTAL_PRE}"
log "Total post-QC variants : ${TOTAL_QC}"

if [ ${#FAILED[@]} -gt 0 ]; then
    log "FAILED chromosomes: ${FAILED[*]}"
    log "Cannot proceed to merge. Fix failures and re-run."
    exit 1
fi

log "All 22 chromosomes passed QC"
log ""

# ---------------------------------------------------------------------------
# Step 3: Merge QC'd per-chrom files
# ---------------------------------------------------------------------------
MERGED_PREFIX="${QC_DIR}/all_background_qc"

if [ -f "${MERGED_PREFIX}.bed" ] && [ "${FORCE}" -eq 0 ]; then
    log "Merged file already exists: ${MERGED_PREFIX}.bed"
    log "Use --force to re-merge."
else
    N_FILES=$(wc -l < "${MERGE_LIST}")
    log "============================================================"
    log "  Merging ${N_FILES} QC'd chromosome files"
    log "============================================================"

    if plink2 \
        --pmerge-list "${MERGE_LIST}" bfile \
        --make-bed \
        --out "${MERGED_PREFIX}" \
        --threads "${THREADS}" \
        --memory 32000 2>&1; then

        MERGED_VARIANTS=$(wc -l < "${MERGED_PREFIX}.bim")
        MERGED_SAMPLES=$(wc -l < "${MERGED_PREFIX}.fam")
        log "Merge complete: ${MERGED_VARIANTS} variants x ${MERGED_SAMPLES} samples"
    else
        log "MERGE FAILED - check plink2 log: ${MERGED_PREFIX}.log"
        exit 1
    fi
fi

# ---------------------------------------------------------------------------
# Step 4: Thin to TARGET_SNPS if needed, write final Step 1 fileset
# ---------------------------------------------------------------------------
FINAL_PREFIX="${STEP1_DIR}/step1_500k"

if [ -f "${FINAL_PREFIX}.bed" ] && [ "${FORCE}" -eq 0 ]; then
    log "Final Step 1 file already exists: ${FINAL_PREFIX}.bed"
    log "Use --force to re-generate."
else
    MERGED_N=$(wc -l < "${MERGED_PREFIX}.bim")

    if [ "${MERGED_N}" -le "${TARGET_SNPS}" ]; then
        log "Merged set (${MERGED_N}) <= target (${TARGET_SNPS}), no thinning needed"
        # Just copy/link the merged files
        cp "${MERGED_PREFIX}.bed" "${FINAL_PREFIX}.bed"
        cp "${MERGED_PREFIX}.bim" "${FINAL_PREFIX}.bim"
        cp "${MERGED_PREFIX}.fam" "${FINAL_PREFIX}.fam"
    else
        log "============================================================"
        log "  Thinning ${MERGED_N} -> ${TARGET_SNPS} variants"
        log "============================================================"

        if plink2 \
            --bfile "${MERGED_PREFIX}" \
            --thin-count "${TARGET_SNPS}" \
            --seed "${RANDOM_SEED}" \
            --make-bed \
            --out "${FINAL_PREFIX}" \
            --threads "${THREADS}" \
            --memory 32000 2>&1; then
            log "Thinning complete"
        else
            log "THINNING FAILED - check plink2 log: ${FINAL_PREFIX}.log"
            exit 1
        fi
    fi

    FINAL_VARIANTS=$(wc -l < "${FINAL_PREFIX}.bim")
    FINAL_SAMPLES=$(wc -l < "${FINAL_PREFIX}.fam")
    log "Final Step 1 fileset: ${FINAL_VARIANTS} variants x ${FINAL_SAMPLES} samples"
fi

# ---------------------------------------------------------------------------
# Step 5: Upload to GCS
# ---------------------------------------------------------------------------
log ""
log "============================================================"
log "  Upload to GCS"
log "============================================================"
if [ -n "${WORKSPACE_BUCKET}" ]; then
    if [[ "${WORKSPACE_BUCKET}" != gs://* ]]; then
        WORKSPACE_BUCKET="gs://${WORKSPACE_BUCKET}"
    fi

    # Upload per-chrom QC files (also used as LD reference for GCTA-COJO)
    GCS_QC_DEST="${WORKSPACE_BUCKET}/results/1-bg_snp/plink_qc/"
    log "Uploading per-chrom QC files to: ${GCS_QC_DEST}"
    gsutil -u "${GOOGLE_PROJECT}" -m cp \
        "${QC_DIR}"/chr*_background_qc.{bed,bim,fam} "${GCS_QC_DEST}"
    log "Per-chrom QC upload complete"

    # Upload Step 1 thinned fileset
    GCS_DEST="${WORKSPACE_BUCKET}/results/1-bg_snp/plink_step1/"
    log "Uploading Step 1 fileset to: ${GCS_DEST}"
    gsutil -u "${GOOGLE_PROJECT}" -m cp "${FINAL_PREFIX}".* "${GCS_DEST}"
    log "Step 1 upload complete"
else
    log "WARNING: WORKSPACE_BUCKET not set -- skipping GCS upload."
    log "  To upload manually:"
    log "  gsutil -u \$GOOGLE_PROJECT -m cp ${QC_DIR}/chr*_background_qc.{bed,bim,fam} \$WORKSPACE_BUCKET/results/1-bg_snp/plink_qc/"
    log "  gsutil -u \$GOOGLE_PROJECT -m cp ${FINAL_PREFIX}.* \$WORKSPACE_BUCKET/results/1-bg_snp/plink_step1/"
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
log "Per-chrom QC : ${QC_DIR}/chrN_background_qc.{bed,bim,fam}"
log "Merged QC    : ${MERGED_PREFIX}.{bed,bim,fam}"
log "Step 1 final : ${FINAL_PREFIX}.{bed,bim,fam}"
log ""
ls -lh "${FINAL_PREFIX}".{bed,bim,fam} 2>/dev/null | while read -r line; do
    log "  ${line}"
done
log ""
log "Next step: bash bash/01c_run_regenie_step1.sh"
log "============================================================"
