#!/bin/bash
set -eo pipefail

# ---------------------------------------------------------------------------
# 00_run_phenotype_covariates.sh
#
# Bash wrapper for Phase 0: Phenotype & Covariate File Generation.
# Calls python/00_build_phenotype_covariates.py to generate REGENIE-formatted
# output files for the full EUR dataset:
#   - MS_phenotype.txt  : FID IID MS  (1=case, 0=control)
#   - MS_covariates.txt : FID IID Age Sex PC1..PC10
#
# Data sources:
#   - EUR IDs + PCs : ancestry_preds.tsv (CDR bucket, requester-pays)
#   - Sex           : genomic_metrics.tsv (CDR bucket, requester-pays)
#   - Age           : BigQuery person table (via WORKSPACE_CDR)
#   - MS cases      : BigQuery condition_occurrence (OMOP 374919 / ICD-10 G35)
#
# Usage:
#   # Production run
#   nohup bash bash/00_run_phenotype_covariates.sh > /dev/null 2>&1 &
#
#   # Test run (500 EUR samples, no GCS upload)
#   bash bash/00_run_phenotype_covariates.sh --test-mode
#
#   # Force overwrite of existing GCS outputs
#   bash bash/00_run_phenotype_covariates.sh --force
#
#   # Restore files from GCS to local disk (crash recovery, no regeneration)
#   bash bash/00_run_phenotype_covariates.sh --restore
#
# Monitor:
#   tail -f logs/00_phenotype_covariates_*.log
#   cat logs/00_phenotype_covariates.pid
#   ps -p $(cat logs/00_phenotype_covariates.pid) -o pid,etime,cmd
#   kill $(cat logs/00_phenotype_covariates.pid)
# ---------------------------------------------------------------------------

trap '' HUP

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
PYTHON_SCRIPT="${PROJECT_DIR}/python/00_build_phenotype_covariates.py"
CONFIG_FILE="${PROJECT_DIR}/config/config.json"
LOG_DIR="${PROJECT_DIR}/logs"

mkdir -p "$LOG_DIR"

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="${LOG_DIR}/00_phenotype_covariates_${TIMESTAMP}.log"
PID_FILE="${LOG_DIR}/00_phenotype_covariates.pid"

exec < /dev/null
exec > >(tee -a "${LOG_FILE}") 2>&1

echo $$ > "${PID_FILE}"

START_SECONDS=$(date +%s)

# ---------------------------------------------------------------------------
# Pre-flight checks
# ---------------------------------------------------------------------------
echo "============================================================"
echo "Phase 0: Phenotype & Covariate File Generation"
echo "============================================================"
echo "Started at  : $(date)"
echo "PID         : $$"
echo "PID file    : ${PID_FILE}"
echo "Log file    : ${LOG_FILE}"
echo ""

# ---------------------------------------------------------------------------
# Parse flags early so --restore can skip BigQuery/CDR checks
# ---------------------------------------------------------------------------
EXTRA_ARGS=""
RESTORE_ONLY=0
for arg in "$@"; do
    case "$arg" in
        --test-mode) EXTRA_ARGS="$EXTRA_ARGS --test-mode" ;;
        --force)     EXTRA_ARGS="$EXTRA_ARGS --force" ;;
        --test-n=*)  EXTRA_ARGS="$EXTRA_ARGS $arg" ;;
        --restore)   RESTORE_ONLY=1 ;;
        *)
            echo "WARNING: Unknown argument '$arg' — passing through to Python script"
            EXTRA_ARGS="$EXTRA_ARGS $arg"
            ;;
    esac
done

# ---------------------------------------------------------------------------
# Environment checks (WORKSPACE_BUCKET + GOOGLE_PROJECT always required)
# ---------------------------------------------------------------------------
if [ -z "$WORKSPACE_BUCKET" ]; then
    echo "ERROR: WORKSPACE_BUCKET environment variable not set"
    exit 1
fi
if [[ "$WORKSPACE_BUCKET" != gs://* ]]; then
    WORKSPACE_BUCKET="gs://${WORKSPACE_BUCKET}"
fi
echo "Workspace   : $WORKSPACE_BUCKET"

if [ -z "$GOOGLE_PROJECT" ]; then
    echo "ERROR: GOOGLE_PROJECT environment variable not set"
    exit 1
fi
echo "Project     : $GOOGLE_PROJECT"

if [ ! -f "$CONFIG_FILE" ]; then
    echo "ERROR: Config file not found: ${CONFIG_FILE}"
    exit 1
fi
echo "Config      : $CONFIG_FILE"

# ---------------------------------------------------------------------------
# Restore mode: download existing outputs from GCS without regenerating
# (skips CDR, Python script, and BigQuery checks)
# ---------------------------------------------------------------------------
LOCAL_PHENO_DIR=$(python3 -c "import json; print(json.load(open('${CONFIG_FILE}'))['outputs']['phenotype_dir'])")
GCS_PHENO_DIR="${WORKSPACE_BUCKET}/${LOCAL_PHENO_DIR}"

if [ "${RESTORE_ONLY}" -eq 1 ]; then
    echo "[--restore] Downloading existing phenotype/covariate files from GCS..."
    mkdir -p "${LOCAL_PHENO_DIR}"
    gsutil -u "${GOOGLE_PROJECT}" cp "${GCS_PHENO_DIR}/MS_phenotype.txt" "${LOCAL_PHENO_DIR}/" \
        && echo "  Downloaded MS_phenotype.txt" \
        || { echo "ERROR: Failed to download MS_phenotype.txt from GCS"; exit 1; }
    gsutil -u "${GOOGLE_PROJECT}" cp "${GCS_PHENO_DIR}/MS_covariates.txt" "${LOCAL_PHENO_DIR}/" \
        && echo "  Downloaded MS_covariates.txt" \
        || { echo "ERROR: Failed to download MS_covariates.txt from GCS"; exit 1; }
    echo ""
    echo "Restored to: ${LOCAL_PHENO_DIR}/"
    ls -lh "${LOCAL_PHENO_DIR}"/MS_*.txt
    rm -f "${PID_FILE}"
    exit 0
fi

# ---------------------------------------------------------------------------
# Additional checks only needed for full regeneration (not --restore)
# ---------------------------------------------------------------------------
if [ -z "$WORKSPACE_CDR" ]; then
    echo "ERROR: WORKSPACE_CDR environment variable not set"
    exit 1
fi
echo "CDR dataset : $WORKSPACE_CDR"

if [ ! -f "$PYTHON_SCRIPT" ]; then
    echo "ERROR: Python script not found: ${PYTHON_SCRIPT}"
    exit 1
fi
echo "Script      : $PYTHON_SCRIPT"
echo ""

# ---------------------------------------------------------------------------
# Run Python script
# ---------------------------------------------------------------------------
echo "Running: python ${PYTHON_SCRIPT} --config ${CONFIG_FILE}${EXTRA_ARGS}"
echo "------------------------------------------------------------"

python "${PYTHON_SCRIPT}" \
    --config "${CONFIG_FILE}" \
    ${EXTRA_ARGS}

EXIT_CODE=$?

END_SECONDS=$(date +%s)
ELAPSED=$(( END_SECONDS - START_SECONDS ))
ELAPSED_FMT=$(printf '%02dh %02dm %02ds' $((ELAPSED/3600)) $(( (ELAPSED%3600)/60 )) $((ELAPSED%60)))

echo ""
echo "============================================================"
if [ $EXIT_CODE -eq 0 ]; then
    echo "Phase 0 COMPLETED successfully in ${ELAPSED_FMT}"
    echo ""
    echo "Output files:"
    echo "  GCS : ${WORKSPACE_BUCKET}/results/0-phenotype/MS_phenotype.txt"
    echo "  GCS : ${WORKSPACE_BUCKET}/results/0-phenotype/MS_covariates.txt"
    echo ""
    echo "Next step:"
    echo "  bash bash/01b_qc_merge_background_snps.sh"
    echo "  bash bash/01c_run_regenie_step1.sh"
else
    echo "Phase 0 FAILED (exit code: ${EXIT_CODE}) after ${ELAPSED_FMT}"
    echo "Check log: ${LOG_FILE}"
fi
echo "============================================================"

rm -f "${PID_FILE}"
exit $EXIT_CODE
