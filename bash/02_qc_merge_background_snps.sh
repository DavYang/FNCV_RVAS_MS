#!/bin/bash
set -e

# ---------------------------------------------------------------------------
# QC per-chromosome PLINK files and merge into a single genome-wide fileset
# for Regenie Step 1 null model fitting.
# ---------------------------------------------------------------------------
# Usage:
#   nohup bash bash/02_qc_merge_background_snps.sh > /dev/null 2>&1 &
#
# Monitor:
#   tail -f logs/02_qc_merge_background_snps.log
#
# Input:  results/1-bg_snp/plink_no-qc/chrN_background.{bed,bim,fam}
# Output: results/1-bg_snp/plink_qc/chrN_background_qc.{bed,bim,fam}
#         results/1-bg_snp/plink_qc/all_background_qc.{bed,bim,fam}
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_DIR="/home/jupyter/FNCV_RVAS_MS"

INPUT_DIR="${PROJECT_DIR}/results/1-bg_snp/plink_no-qc"
OUTPUT_DIR="${PROJECT_DIR}/results/1-bg_snp/plink_qc"
LOG_DIR="${PROJECT_DIR}/logs"
LOG_FILE="${LOG_DIR}/02_qc_merge_background_snps.log"
MERGE_LIST="${OUTPUT_DIR}/merge_list.txt"

# ---------------------------------------------------------------------------
# QC thresholds
# ---------------------------------------------------------------------------
GENO=0.05
MAF=0.01
HWE=1e-6
MIND=0.1

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}"

log() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $*" | tee -a "${LOG_FILE}"
}

log "============================================================"
log "  QC + Merge Background SNPs for Regenie Step 1"
log "============================================================"
log "Input dir  : ${INPUT_DIR}"
log "Output dir : ${OUTPUT_DIR}"
log "QC params  : MAF=${MAF}, GENO=${GENO}, HWE=${HWE}, MIND=${MIND}"
log "============================================================"
log ""

# ---------------------------------------------------------------------------
# Step 1: Per-chromosome QC
# ---------------------------------------------------------------------------
FAILED=()
TOTAL_VARIANTS_PRE=0
TOTAL_VARIANTS_POST=0

rm -f "${MERGE_LIST}"

for chr_num in $(seq 1 22); do
    chr_name="chr${chr_num}"
    INPUT_PREFIX="${INPUT_DIR}/${chr_name}_background"
    OUTPUT_PREFIX="${OUTPUT_DIR}/${chr_name}_background_qc"

    # Verify input exists
    if [ ! -f "${INPUT_PREFIX}.bed" ]; then
        log "[${chr_name}] SKIPPED - input .bed not found: ${INPUT_PREFIX}.bed"
        FAILED+=("${chr_name}")
        continue
    fi

    # Count pre-QC variants
    PRE_COUNT=$(wc -l < "${INPUT_PREFIX}.bim")
    TOTAL_VARIANTS_PRE=$((TOTAL_VARIANTS_PRE + PRE_COUNT))

    # Resume: skip if QC output already exists
    if [ -f "${OUTPUT_PREFIX}.bed" ]; then
        POST_COUNT=$(wc -l < "${OUTPUT_PREFIX}.bim")
        TOTAL_VARIANTS_POST=$((TOTAL_VARIANTS_POST + POST_COUNT))
        log "[${chr_name}] SKIP (exists) - pre: ${PRE_COUNT}, post: ${POST_COUNT} variants"
        echo "${OUTPUT_PREFIX}" >> "${MERGE_LIST}"
        continue
    fi

    log "[${chr_name}] QC starting - ${PRE_COUNT} variants pre-QC"

    if plink2 \
        --bfile "${INPUT_PREFIX}" \
        --maf "${MAF}" \
        --geno "${GENO}" \
        --hwe "${HWE}" \
        --mind "${MIND}" \
        --max-alleles 2 \
        --make-bed \
        --out "${OUTPUT_PREFIX}" \
        --threads 50 \
        --memory 16000 >> "${LOG_FILE}" 2>&1; then

        POST_COUNT=$(wc -l < "${OUTPUT_PREFIX}.bim")
        TOTAL_VARIANTS_POST=$((TOTAL_VARIANTS_POST + POST_COUNT))
        REMOVED=$((PRE_COUNT - POST_COUNT))
        log "[${chr_name}] QC done - ${POST_COUNT} variants retained (${REMOVED} removed)"
        echo "${OUTPUT_PREFIX}" >> "${MERGE_LIST}"
    else
        log "[${chr_name}] QC FAILED"
        FAILED+=("${chr_name}")
    fi
done

log ""
log "============================================================"
log "  Per-chromosome QC summary"
log "============================================================"
log "Total variants pre-QC  : ${TOTAL_VARIANTS_PRE}"
log "Total variants post-QC : ${TOTAL_VARIANTS_POST}"
log "Total removed          : $((TOTAL_VARIANTS_PRE - TOTAL_VARIANTS_POST))"

if [ ${#FAILED[@]} -gt 0 ]; then
    log "FAILED chromosomes: ${FAILED[*]}"
    log "Cannot proceed to merge. Fix failures and re-run."
    exit 1
fi

log "All 22 chromosomes passed QC"
log ""

# ---------------------------------------------------------------------------
# Step 2: Merge all QC'd chromosomes into a single PLINK fileset
# ---------------------------------------------------------------------------
MERGED_PREFIX="${OUTPUT_DIR}/all_background_qc"

if [ -f "${MERGED_PREFIX}.bed" ]; then
    log "Merged file already exists: ${MERGED_PREFIX}.bed"
    log "Delete it to re-merge."
else
    N_FILES=$(wc -l < "${MERGE_LIST}")
    log "============================================================"
    log "  Merging ${N_FILES} chromosome files"
    log "============================================================"

    if plink2 \
        --pmerge-list "${MERGE_LIST}" bfile \
        --make-bed \
        --out "${MERGED_PREFIX}" \
        --threads 50 \
        --memory 32000 >> "${LOG_FILE}" 2>&1; then

        MERGED_VARIANTS=$(wc -l < "${MERGED_PREFIX}.bim")
        MERGED_SAMPLES=$(wc -l < "${MERGED_PREFIX}.fam")
        log "Merge complete: ${MERGED_VARIANTS} variants x ${MERGED_SAMPLES} samples"
    else
        log "MERGE FAILED - check plink2 log: ${MERGED_PREFIX}.log"
        exit 1
    fi
fi

# ---------------------------------------------------------------------------
# Step 3: Post-merge QC on the genome-wide merged file
# ---------------------------------------------------------------------------
FINAL_PREFIX="${OUTPUT_DIR}/all_background_final_qc"

if [ -f "${FINAL_PREFIX}.bed" ]; then
    log "Post-merge QC file already exists: ${FINAL_PREFIX}.bed"
    log "Delete it to re-run post-merge QC."
else
    MERGED_PRE=$(wc -l < "${MERGED_PREFIX}.bim")
    log "============================================================"
    log "  Post-merge QC (genome-wide)"
    log "============================================================"
    log "Input  : ${MERGED_PREFIX} (${MERGED_PRE} variants)"
    log "Output : ${FINAL_PREFIX}"

    if plink2 \
        --bfile "${MERGED_PREFIX}" \
        --maf "${MAF}" \
        --geno "${GENO}" \
        --hwe "${HWE}" \
        --mind "${MIND}" \
        --max-alleles 2 \
        --make-bed \
        --out "${FINAL_PREFIX}" \
        --threads 50 \
        --memory 32000 >> "${LOG_FILE}" 2>&1; then

        FINAL_VARIANTS=$(wc -l < "${FINAL_PREFIX}.bim")
        FINAL_SAMPLES=$(wc -l < "${FINAL_PREFIX}.fam")
        FINAL_REMOVED=$((MERGED_PRE - FINAL_VARIANTS))
        log "Post-merge QC done: ${FINAL_VARIANTS} variants x ${FINAL_SAMPLES} samples (${FINAL_REMOVED} removed)"
    else
        log "POST-MERGE QC FAILED - check plink2 log: ${FINAL_PREFIX}.log"
        exit 1
    fi
fi

# ---------------------------------------------------------------------------
# Final summary
# ---------------------------------------------------------------------------
log ""
log "============================================================"
log "  Pipeline complete"
log "============================================================"
log "Per-chrom QC files : ${OUTPUT_DIR}/chrN_background_qc.{bed,bim,fam}"
log "Merged file        : ${MERGED_PREFIX}.{bed,bim,fam}"
log "Final QC file      : ${FINAL_PREFIX}.{bed,bim,fam}"
log ""
ls -lh "${FINAL_PREFIX}".{bed,bim,fam} 2>/dev/null | while read -r line; do
    log "  ${line}"
done
log "============================================================"
