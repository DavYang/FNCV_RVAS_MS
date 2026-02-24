#!/bin/bash
set -eo pipefail

# ---------------------------------------------------------------------------
# 03_run_regenie_step1.sh
#
# Phase 1: Run REGENIE Step 1 (Leave-One-Chromosome-Out whole-genome regression)
# to fit the null model for binary trait (MS) correcting for age, sex, and
# top 10 ancestry PCs.
#
# Inputs  (all local on VM):
#   results/1-bg_snp/plink_step1/step1_500k.{bed,bim,fam}  (500K downsampled; see 02b_downsample_background_snps.sh)
#   results/0-phenotype/MS_phenotype.txt
#   results/0-phenotype/MS_covariates.txt
#
# Outputs (local, then uploaded to GCS):
#   results/1-bg_snp/regenie/MS_null_model_pred.list
#   results/1-bg_snp/regenie/MS_null_model_1.loco  (one per phenotype)
#   results/1-bg_snp/regenie/MS_null_model.log
#
# Usage:
#   # Production run (nohup)
#   nohup bash bash/03_run_regenie_step1.sh > /dev/null 2>&1 &
#
#   # Force re-run even if outputs exist
#   bash bash/03_run_regenie_step1.sh --force
#
# Monitor:
#   tail -f logs/03_regenie_step1_*.log
#   cat logs/03_regenie_step1.pid
#   ps -p $(cat logs/03_regenie_step1.pid) -o pid,etime,cmd
#   kill $(cat logs/03_regenie_step1.pid)
# ---------------------------------------------------------------------------

trap '' HUP

PROJECT_DIR="/home/jupyter/FNCV_RVAS_MS"

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PLINK_PREFIX="${PROJECT_DIR}/results/1-bg_snp/plink_step1/step1_500k"
PHENO_FILE="${PROJECT_DIR}/results/0-phenotype/MS_phenotype.txt"
COVAR_FILE="${PROJECT_DIR}/results/0-phenotype/MS_covariates.txt"
OUT_DIR="${PROJECT_DIR}/results/1-bg_snp/regenie"
OUT_PREFIX="${OUT_DIR}/MS_null_model"
TMP_DIR="${PROJECT_DIR}/tmp/regenie_step1"
LOG_DIR="${PROJECT_DIR}/logs"

# ---------------------------------------------------------------------------
# REGENIE parameters
# ---------------------------------------------------------------------------
BSIZE=1000
THREADS=32
MEM=110000   # MB — leave headroom on 120GB node

# Covariate columns to use (must match header in MS_covariates.txt exactly)
COVAR_COLS="Age,Sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"

# ---------------------------------------------------------------------------
# GCS upload helper
# ---------------------------------------------------------------------------
_upload_results() {
    if [ -z "${WORKSPACE_BUCKET}" ]; then
        echo "WARNING: WORKSPACE_BUCKET not set — skipping GCS upload"
        return
    fi
    if [[ "${WORKSPACE_BUCKET}" != gs://* ]]; then
        WORKSPACE_BUCKET="gs://${WORKSPACE_BUCKET}"
    fi
    GCS_DEST="${WORKSPACE_BUCKET}/results/1-bg_snp/regenie/"
    echo "Uploading results to GCS: ${GCS_DEST}"
    gsutil -u "${GOOGLE_PROJECT}" -m cp \
        "${OUT_DIR}"/MS_null_model* \
        "${GCS_DEST}"
    echo "  Upload complete."
    echo ""
    echo "GCS output path: ${GCS_DEST}"
}

# ---------------------------------------------------------------------------
# Setup logging
# ---------------------------------------------------------------------------
mkdir -p "${OUT_DIR}" "${TMP_DIR}" "${LOG_DIR}"

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="${LOG_DIR}/03_regenie_step1_${TIMESTAMP}.log"
PID_FILE="${LOG_DIR}/03_regenie_step1.pid"

exec < /dev/null
exec > >(tee -a "${LOG_FILE}") 2>&1

echo $$ > "${PID_FILE}"

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
# Pre-flight checks
# ---------------------------------------------------------------------------
echo "============================================================"
echo "Phase 1: REGENIE Step 1 — Null Model Fitting"
echo "============================================================"
echo "Started at  : $(date)"
echo "PID         : $$"
echo "PID file    : ${PID_FILE}"
echo "Log file    : ${LOG_FILE}"
echo ""

# Check REGENIE is available
if ! command -v regenie &> /dev/null; then
    echo "ERROR: 'regenie' not found in PATH"
    exit 1
fi
REGENIE_VERSION=$(regenie --version 2>&1 | head -1 || true)
echo "REGENIE     : ${REGENIE_VERSION}"

# Validate PLINK fileset
for ext in bed bim fam; do
    if [ ! -f "${PLINK_PREFIX}.${ext}" ]; then
        echo "ERROR: Missing PLINK file: ${PLINK_PREFIX}.${ext}"
        exit 1
    fi
done
N_VARIANTS=$(wc -l < "${PLINK_PREFIX}.bim")
N_SAMPLES=$(wc -l < "${PLINK_PREFIX}.fam")
echo "PLINK input : ${PLINK_PREFIX}"
echo "  Variants  : ${N_VARIANTS}"
echo "  Samples   : ${N_SAMPLES}"

# Validate phenotype and covariate files
if [ ! -f "${PHENO_FILE}" ]; then
    echo "ERROR: Phenotype file not found: ${PHENO_FILE}"
    exit 1
fi
if [ ! -f "${COVAR_FILE}" ]; then
    echo "ERROR: Covariate file not found: ${COVAR_FILE}"
    exit 1
fi
N_PHENO=$(( $(wc -l < "${PHENO_FILE}") - 1 ))
N_COVAR=$(( $(wc -l < "${COVAR_FILE}") - 1 ))
echo "Phenotype   : ${PHENO_FILE}  (${N_PHENO} samples)"
echo "Covariates  : ${COVAR_FILE}  (${N_COVAR} samples)"
echo "Covar cols  : ${COVAR_COLS}"
echo "Params      : --bsize ${BSIZE}  --threads ${THREADS}  --lowmem"
echo "Output dir  : ${OUT_DIR}"
echo ""

# Resume check
PRED_LIST="${OUT_PREFIX}_pred.list"
if [ -f "${PRED_LIST}" ] && [ "${FORCE}" -eq 0 ]; then
    echo "Output already exists: ${PRED_LIST}"
    echo "Use --force to re-run."

    # Still attempt GCS upload in case it was skipped previously
    if [ -n "${WORKSPACE_BUCKET}" ]; then
        echo ""
        echo "Uploading existing results to GCS..."
        _upload_results
    fi

    rm -f "${PID_FILE}"
    exit 0
fi

if [ -f "${PRED_LIST}" ] && [ "${FORCE}" -eq 1 ]; then
    echo "[--force] Removing existing output files..."
    rm -f "${OUT_PREFIX}"*.loco "${OUT_PREFIX}"*.list "${OUT_PREFIX}".log
fi

# ---------------------------------------------------------------------------
# Run REGENIE Step 1
# ---------------------------------------------------------------------------
echo "============================================================"
echo "  Running REGENIE Step 1"
echo "============================================================"
echo "Command:"
echo "  regenie \\"
echo "    --step 1 \\"
echo "    --bed ${PLINK_PREFIX} \\"
echo "    --phenoFile ${PHENO_FILE} \\"
echo "    --covarFile ${COVAR_FILE} \\"
echo "    --covarColList ${COVAR_COLS} \\"
echo "    --bt \\"
echo "    --bsize ${BSIZE} \\"
echo "    --lowmem \\"
echo "    --lowmem-prefix ${TMP_DIR}/regenie_tmp \\"
echo "    --threads ${THREADS} \\"
echo "    --out ${OUT_PREFIX}"
echo ""

regenie \
    --step 1 \
    --bed "${PLINK_PREFIX}" \
    --phenoFile "${PHENO_FILE}" \
    --covarFile "${COVAR_FILE}" \
    --covarColList "${COVAR_COLS}" \
    --bt \
    --bsize "${BSIZE}" \
    --lowmem \
    --lowmem-prefix "${TMP_DIR}/regenie_tmp" \
    --threads "${THREADS}" \
    --out "${OUT_PREFIX}"

REGENIE_EXIT=$?

if [ ${REGENIE_EXIT} -ne 0 ]; then
    echo ""
    echo "ERROR: REGENIE Step 1 failed (exit code: ${REGENIE_EXIT})"
    echo "Check REGENIE log: ${OUT_PREFIX}.log"
    rm -f "${PID_FILE}"
    exit ${REGENIE_EXIT}
fi

# ---------------------------------------------------------------------------
# Validate output
# ---------------------------------------------------------------------------
echo ""
echo "============================================================"
echo "  Validating output"
echo "============================================================"

if [ ! -f "${PRED_LIST}" ]; then
    echo "ERROR: Expected output not found: ${PRED_LIST}"
    echo "Check REGENIE log: ${OUT_PREFIX}.log"
    rm -f "${PID_FILE}"
    exit 1
fi

N_LOCO=$(wc -l < "${PRED_LIST}")
echo "Prediction list : ${PRED_LIST}  (${N_LOCO} entries)"
echo ""
echo "Output files:"
ls -lh "${OUT_DIR}"/MS_null_model* 2>/dev/null | while read -r line; do
    echo "  ${line}"
done

# Clean up lowmem temp files
echo ""
echo "Cleaning up lowmem temp files in ${TMP_DIR}..."
rm -f "${TMP_DIR}"/regenie_tmp*
echo "  Done."

# ---------------------------------------------------------------------------
# Upload to GCS
# ---------------------------------------------------------------------------
echo ""
echo "============================================================"
echo "  Uploading to GCS"
echo "============================================================"

if [ -z "${WORKSPACE_BUCKET}" ]; then
    echo "WARNING: WORKSPACE_BUCKET not set — skipping GCS upload"
    echo "  To upload manually:"
    echo "  gsutil -u \$GOOGLE_PROJECT -m cp ${OUT_DIR}/MS_null_model* \$WORKSPACE_BUCKET/results/1-bg_snp/regenie/"
else
    _upload_results
fi

# ---------------------------------------------------------------------------
# Final summary
# ---------------------------------------------------------------------------
END_SECONDS=$(date +%s)
ELAPSED=$(( END_SECONDS - START_SECONDS ))
ELAPSED_FMT=$(printf '%02dh %02dm %02ds' $((ELAPSED/3600)) $(( (ELAPSED%3600)/60 )) $((ELAPSED%60)))

echo ""
echo "============================================================"
echo "Phase 1 REGENIE Step 1 COMPLETED in ${ELAPSED_FMT}"
echo ""
echo "Null model outputs:"
echo "  Local : ${OUT_PREFIX}_pred.list"
echo "  Local : ${OUT_DIR}/MS_null_model_1.loco  (+ per-chromosome .loco files)"
if [ -n "${WORKSPACE_BUCKET}" ]; then
    echo "  GCS   : ${WORKSPACE_BUCKET}/results/1-bg_snp/regenie/"
fi
echo ""
echo "REGENIE Step 2 inputs:"
echo "  --pred  ${OUT_PREFIX}_pred.list"
echo "============================================================"

rm -f "${PID_FILE}"
exit 0
