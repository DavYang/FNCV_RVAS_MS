#!/bin/bash
set -eo pipefail

# ---------------------------------------------------------------------------
# 02b_downsample_background_snps.sh
#
# Downsample the QC'd genome-wide background PLINK fileset to 500K SNPs for
# use as REGENIE Step 1 input. REGENIE recommends <=1M variants; using 500K
# provides a safe margin with minimal loss of null model accuracy.
#
# Strategy:
#   1. Apply LD pruning (--indep-pairwise 1000 100 0.9) to remove redundant
#      variants, then thin the pruned set to exactly TARGET_SNPS with
#      --thin-count. This preserves spread across the genome better than
#      random thinning alone.
#   2. Exclude the MHC region (chr6:25-35 Mb) prior to pruning to avoid
#      inflating LD structure in that region.
#
# Input  : results/1-bg_snp/plink_qc/all_background_final_qc.{bed,bim,fam}
# Output : results/1-bg_snp/plink_step1/step1_500k.{bed,bim,fam}
#
# Usage:
#   # Production run (nohup)
#   nohup bash bash/02b_downsample_background_snps.sh > /dev/null 2>&1 &
#
#   # Force re-run even if output exists
#   bash bash/02b_downsample_background_snps.sh --force
#
# Monitor:
#   tail -f logs/02b_downsample_background_snps_*.log
# ---------------------------------------------------------------------------

PROJECT_DIR="/home/jupyter/FNCV_RVAS_MS"

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
INPUT_PREFIX="${PROJECT_DIR}/results/1-bg_snp/plink_qc/all_background_final_qc"
OUT_DIR="${PROJECT_DIR}/results/1-bg_snp/plink_step1"
OUT_PREFIX="${OUT_DIR}/step1_500k"
TMP_DIR="${PROJECT_DIR}/tmp/downsample"
LOG_DIR="${PROJECT_DIR}/logs"

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
TARGET_SNPS=500000
THREADS=50

# Pre-LD QC thresholds
QC_MIND=0.05
QC_GENO=0.01
QC_MAF=0.01
QC_HWE="1e-6 midp"
# Relatedness threshold for up to 2nd degree relatives
QC_KING=0.0884

# LD pruning window/step/r2 — loose r2 threshold keeps more genome coverage
LD_WINDOW=1000
LD_STEP=100
LD_R2=0.9

# Parse MHC exclusion range from config.json
MHC_INTERVAL=$(python3 -c "import json; print(json.load(open('${PROJECT_DIR}/config/config.json'))['params']['mhc_interval'])")
MHC_CHROM=$(echo "$MHC_INTERVAL" | awk -F'[:-]' '{print $1}' | sed 's/chr//')
MHC_START=$(echo "$MHC_INTERVAL" | awk -F'[:-]' '{print $2}')
MHC_END=$(echo "$MHC_INTERVAL" | awk -F'[:-]' '{print $3}')

RANDOM_SEED=$(python3 -c "import json; print(json.load(open('${PROJECT_DIR}/config/config.json'))['params']['random_seed'])")

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
mkdir -p "${OUT_DIR}" "${TMP_DIR}" "${LOG_DIR}"

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="${LOG_DIR}/02b_downsample_background_snps_${TIMESTAMP}.log"

exec < /dev/null
exec > >(tee -a "${LOG_FILE}") 2>&1

START_SECONDS=$(date +%s)

# ---------------------------------------------------------------------------
# Parse flags
# ---------------------------------------------------------------------------
FORCE=0
for arg in "$@"; do
    case "$arg" in
        --force) FORCE=1 ;;
        *) echo "WARNING: Unknown argument '$arg' — ignoring" ;;
    esac
done

# ---------------------------------------------------------------------------
# Pre-flight
# ---------------------------------------------------------------------------
echo "============================================================"
echo "  Downsample Background SNPs for REGENIE Step 1"
echo "============================================================"
echo "Started at    : $(date)"
echo "Log file      : ${LOG_FILE}"
echo "Input         : ${INPUT_PREFIX}"
echo "Output        : ${OUT_PREFIX}"
echo "Pre-LD QC     : --mind ${QC_MIND} --geno ${QC_GENO} --maf ${QC_MAF} --hwe ${QC_HWE}"
echo "Relatedness   : --king-cutoff ${QC_KING} (up to 2nd degree)"
echo "Target SNPs   : ${TARGET_SNPS}"
echo "LD pruning    : window=${LD_WINDOW} step=${LD_STEP} r2=${LD_R2}"
echo "MHC exclusion : chr${MHC_CHROM}:${MHC_START}-${MHC_END}"
echo "Threads       : ${THREADS}"
echo ""

for ext in bed bim fam; do
    if [ ! -f "${INPUT_PREFIX}.${ext}" ]; then
        echo "ERROR: Missing input file: ${INPUT_PREFIX}.${ext}"
        exit 1
    fi
done

N_INPUT=$(wc -l < "${INPUT_PREFIX}.bim")
echo "Input variants: ${N_INPUT}"

# Resume check
if [ -f "${OUT_PREFIX}.bed" ] && [ "${FORCE}" -eq 0 ]; then
    N_OUT=$(wc -l < "${OUT_PREFIX}.bim")
    echo "Output already exists: ${OUT_PREFIX}.bed  (${N_OUT} variants)"
    echo "Use --force to re-run."
    exit 0
fi

if [ -f "${OUT_PREFIX}.bed" ] && [ "${FORCE}" -eq 1 ]; then
    echo "[--force] Removing existing output..."
    rm -f "${OUT_PREFIX}".{bed,bim,fam,log}
fi

# ---------------------------------------------------------------------------
# Step 1: Pre-LD QC and Relatedness Filtering
# ---------------------------------------------------------------------------
echo "------------------------------------------------------------"
echo "  Step 1: Pre-LD QC (mind, geno, maf, hwe, relatedness)"
echo "------------------------------------------------------------"

QC_PREFIX="${TMP_DIR}/qc_filtered"

plink2 \
    --bfile "${INPUT_PREFIX}" \
    --mind "${QC_MIND}" \
    --geno "${QC_GENO}" \
    --maf "${QC_MAF}" \
    --hwe ${QC_HWE} \
    --king-cutoff "${QC_KING}" \
    --make-bed \
    --threads "${THREADS}" \
    --memory 32000 \
    --out "${QC_PREFIX}"

N_QC_VARIANTS=$(wc -l < "${QC_PREFIX}.bim")
N_QC_SAMPLES=$(wc -l < "${QC_PREFIX}.fam")
echo "  Variants passing QC : ${N_QC_VARIANTS}"
echo "  Samples passing QC  : ${N_QC_SAMPLES}"

# ---------------------------------------------------------------------------
# Step 2: Exclude MHC region and write a SNP exclusion list
# ---------------------------------------------------------------------------
echo ""
echo "------------------------------------------------------------"
echo "  Step 2: Exclude MHC region (chr${MHC_CHROM}:${MHC_START}-${MHC_END})"
echo "------------------------------------------------------------"

MHC_EXCLUDE="${TMP_DIR}/mhc_snps.txt"

awk -v chrom="${MHC_CHROM}" \
    -v start="${MHC_START}" \
    -v end="${MHC_END}" \
    '$1 == chrom && $4 >= start && $4 <= end { print $2 }' \
    "${QC_PREFIX}.bim" > "${MHC_EXCLUDE}"

N_MHC=$(wc -l < "${MHC_EXCLUDE}")
echo "  MHC SNPs excluded: ${N_MHC}"

# ---------------------------------------------------------------------------
# Step 3: LD pruning on non-MHC variants
# ---------------------------------------------------------------------------
echo ""
echo "------------------------------------------------------------"
echo "  Step 3: LD pruning (r2 < ${LD_R2})"
echo "------------------------------------------------------------"

PRUNE_PREFIX="${TMP_DIR}/pruned"

plink2 \
    --bfile "${QC_PREFIX}" \
    --exclude "${MHC_EXCLUDE}" \
    --indep-pairwise "${LD_WINDOW}" "${LD_STEP}" "${LD_R2}" \
    --threads "${THREADS}" \
    --memory 32000 \
    --seed "${RANDOM_SEED}" \
    --out "${PRUNE_PREFIX}"

N_PRUNED_IN=$(wc -l < "${PRUNE_PREFIX}.prune.in")
N_PRUNED_OUT=$(wc -l < "${PRUNE_PREFIX}.prune.out")
echo "  Variants kept after LD pruning : ${N_PRUNED_IN}"
echo "  Variants removed by LD pruning : ${N_PRUNED_OUT}"

if [ "${N_PRUNED_IN}" -lt "${TARGET_SNPS}" ]; then
    echo ""
    echo "WARNING: Pruned set (${N_PRUNED_IN}) < target (${TARGET_SNPS})."
    echo "  Using all pruned variants without thinning."
    TARGET_SNPS=${N_PRUNED_IN}
fi

# ---------------------------------------------------------------------------
# Step 4: Thin pruned set to TARGET_SNPS and write final PLINK fileset
# ---------------------------------------------------------------------------
echo ""
echo "------------------------------------------------------------"
echo "  Step 4: Thin to ${TARGET_SNPS} SNPs and write PLINK fileset"
echo "------------------------------------------------------------"

plink2 \
    --bfile "${QC_PREFIX}" \
    --extract "${PRUNE_PREFIX}.prune.in" \
    --thin-count "${TARGET_SNPS}" \
    --seed "${RANDOM_SEED}" \
    --make-bed \
    --threads "${THREADS}" \
    --memory 32000 \
    --out "${OUT_PREFIX}"

N_FINAL=$(wc -l < "${OUT_PREFIX}.bim")
N_SAMPLES=$(wc -l < "${OUT_PREFIX}.fam")
echo ""
echo "  Final variants : ${N_FINAL}"
echo "  Samples        : ${N_SAMPLES}"
echo "  Output         : ${OUT_PREFIX}.{bed,bim,fam}"

# ---------------------------------------------------------------------------
# Step 5: Upload to GCS
# ---------------------------------------------------------------------------
echo ""
echo "------------------------------------------------------------"
echo "  Step 5: Upload to GCS"
echo "------------------------------------------------------------"
if [ -n "${WORKSPACE_BUCKET}" ]; then
    if [[ "${WORKSPACE_BUCKET}" != gs://* ]]; then
        WORKSPACE_BUCKET="gs://${WORKSPACE_BUCKET}"
    fi
    GCS_DEST="${WORKSPACE_BUCKET}/results/1-bg_snp/plink_step1/"
    echo "Uploading to GCS: ${GCS_DEST}"
    gsutil -u "${GOOGLE_PROJECT}" -m cp "${OUT_PREFIX}".* "${GCS_DEST}"
    echo "  Upload complete."
else
    echo "WARNING: WORKSPACE_BUCKET not set — skipping GCS upload."
    echo "  To upload manually:"
    echo "  gsutil -u \$GOOGLE_PROJECT -m cp ${OUT_PREFIX}.* \$WORKSPACE_BUCKET/results/1-bg_snp/plink_step1/"
fi

# ---------------------------------------------------------------------------
# Cleanup temp files
# ---------------------------------------------------------------------------
echo ""
echo "Cleaning up temp files..."
rm -f "${QC_PREFIX}".{bed,bim,fam,log} \
      "${PRUNE_PREFIX}".prune.in \
      "${PRUNE_PREFIX}".prune.out \
      "${PRUNE_PREFIX}".log \
      "${MHC_EXCLUDE}"
echo "  Done."

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
END_SECONDS=$(date +%s)
ELAPSED=$(( END_SECONDS - START_SECONDS ))
ELAPSED_FMT=$(printf '%02dh %02dm %02ds' $((ELAPSED/3600)) $(( (ELAPSED%3600)/60 )) $((ELAPSED%60)))

echo ""
echo "============================================================"
echo "Downsample COMPLETED in ${ELAPSED_FMT}"
echo ""
echo "  Input  : ${INPUT_PREFIX}  (${N_INPUT} variants)"
echo "  Output : ${OUT_PREFIX}  (${N_FINAL} variants, ${N_SAMPLES} samples)"
echo ""
echo "Pass to REGENIE Step 1 via:"
echo "  bash bash/03_run_regenie_step1.sh"
echo "============================================================"

exit 0
