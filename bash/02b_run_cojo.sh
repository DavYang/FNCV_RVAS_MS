#!/bin/bash
set -eo pipefail

# ---------------------------------------------------------------------------
# 02b_run_cojo.sh
#
# Phase 2 Step 2: Run GCTA-COJO per chromosome for locus definition.
# Uses per-chromosome LD reference panels (from 01b_qc_merge_background_snps.sh) and
# per-chromosome .ma summary statistics (from 02a_export_gwas_ma.sh).
#
# Inputs  (local on VM):
#   results/1-bg_snp/plink_qc/chrN_background_qc.{bed,bim,fam}
#   results/2-locus_definition/ma/chrN.ma
#   tools/gcta64
#
# Outputs (local, then uploaded to GCS):
#   results/2-locus_definition/cojo/chrN/{locus}.jma.cojo
#   results/2-locus_definition/all_independent_signals.tsv
#   results/2-locus_definition/target_loci.bed
#
# Usage:
#   nohup bash bash/02b_run_cojo.sh > /dev/null 2>&1 &
#
#   # Force re-run (overwrite existing COJO results)
#   bash bash/02b_run_cojo.sh --force
#
#   # Single chromosome (use = syntax)
#   bash bash/02b_run_cojo.sh --chrom=chr6
#
#   # Test mode: process only test chromosome (chr21) for iterative testing
#   bash bash/02b_run_cojo.sh --test
#
# Monitor:
#   tail -f logs/02b_run_cojo_*.log
#   cat logs/02b_run_cojo.pid
# ---------------------------------------------------------------------------

trap '' HUP

PROJECT_DIR="/home/jupyter/FNCV_RVAS_MS"

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
CONFIG_FILE="${PROJECT_DIR}/config/config.json"
PYTHON_SCRIPT="${PROJECT_DIR}/python/02b_define_loci.py"
OUTPUT_DIR="${PROJECT_DIR}/results/2-locus_definition"
LD_REF_DIR="${PROJECT_DIR}/results/1-bg_snp/plink_qc"
MA_DIR="${PROJECT_DIR}/results/2-locus_definition/ma"
GCTA_BIN="${PROJECT_DIR}/tools/gcta64"
LOG_DIR="${PROJECT_DIR}/logs"

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="${LOG_DIR}/02b_run_cojo_${TIMESTAMP}.log"
PID_FILE="${LOG_DIR}/02b_run_cojo.pid"

# ---------------------------------------------------------------------------
# Parse arguments (pass through to Python)
# ---------------------------------------------------------------------------
EXTRA_ARGS=()
FORCE=0
TEST_MODE=0
TEST_CHR=$(python3 -c "import json; print(json.load(open('${CONFIG_FILE}'))['params']['test_chromosome'])")
for arg in "$@"; do
    case "$arg" in
        --force) FORCE=1; EXTRA_ARGS+=("--force") ;;
        --test)  TEST_MODE=1; EXTRA_ARGS+=("--chrom" "${TEST_CHR}") ;;
        --chrom=*) EXTRA_ARGS+=("--chrom" "${arg#--chrom=}") ;;
        *) echo "WARNING: Unknown argument '$arg' -- ignoring" ;;
    esac
done

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}"
echo $$ > "${PID_FILE}"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "${LOG_FILE}"
}

START_SECONDS=$(date +%s)

log "============================================================"
log "  02b_run_cojo.sh -- GCTA-COJO Locus Definition"
log "============================================================"
log "  Timestamp : ${TIMESTAMP}"
log "  PID       : $$"
log "  Config    : ${CONFIG_FILE}"
log "  LD ref    : ${LD_REF_DIR}"
log "  .ma dir   : ${MA_DIR}"
log "  GCTA      : ${GCTA_BIN}"
log "  Output    : ${OUTPUT_DIR}"
log "  Force     : ${FORCE}"
log "  Test mode : ${TEST_MODE} (${TEST_CHR})"
log ""

# ---------------------------------------------------------------------------
# Pre-flight checks
# ---------------------------------------------------------------------------
PREFLIGHT_OK=1

if [ ! -f "${CONFIG_FILE}" ]; then
    log "ERROR: Config file not found: ${CONFIG_FILE}"
    PREFLIGHT_OK=0
fi

if [ ! -f "${PYTHON_SCRIPT}" ]; then
    log "ERROR: Python script not found: ${PYTHON_SCRIPT}"
    PREFLIGHT_OK=0
fi

# Check GCTA binary
if [ ! -x "${GCTA_BIN}" ]; then
    log "GCTA binary not found or not executable: ${GCTA_BIN}"
    log "Attempting to install GCTA ..."
    if python3 "${PROJECT_DIR}/python/install_gcta.py"; then
        log "GCTA installed successfully"
    else
        log "ERROR: GCTA installation failed"
        PREFLIGHT_OK=0
    fi
fi

# Check LD reference files exist (use test chrom in test mode, chr1 otherwise)
CHECK_CHR="chr1"
if [ "${TEST_MODE}" -eq 1 ]; then
    CHECK_CHR="${TEST_CHR}"
fi
if [ ! -f "${LD_REF_DIR}/${CHECK_CHR}_background_qc.bed" ]; then
    log "ERROR: LD reference not found: ${LD_REF_DIR}/${CHECK_CHR}_background_qc.bed"
    log "  Run: bash bash/01b_qc_merge_background_snps.sh"
    PREFLIGHT_OK=0
fi

# Check .ma files exist
if [ ! -f "${MA_DIR}/${CHECK_CHR}.ma" ]; then
    log "ERROR: .ma file not found: ${MA_DIR}/${CHECK_CHR}.ma"
    log "  Run: bash bash/02a_export_gwas_ma.sh $([ ${TEST_MODE} -eq 1 ] && echo '--test')"
    PREFLIGHT_OK=0
fi

# Check plink2 is available
if ! command -v plink2 &>/dev/null; then
    log "ERROR: plink2 not found in PATH"
    PREFLIGHT_OK=0
fi

# Check Hail is available (needed for loading GWAS HT in Python script)
if ! python3 -c "import hail" 2>/dev/null; then
    log "ERROR: Hail not available in current Python environment."
    PREFLIGHT_OK=0
fi

if [ "${PREFLIGHT_OK}" -eq 0 ]; then
    log ""
    log "Pre-flight checks FAILED. Fix the above errors and re-run."
    rm -f "${PID_FILE}"
    exit 1
fi

log "Pre-flight checks PASSED"
log ""

# ---------------------------------------------------------------------------
# GCS download of inputs if missing locally
# ---------------------------------------------------------------------------
if [ -n "${WORKSPACE_BUCKET}" ]; then
    if [[ "${WORKSPACE_BUCKET}" != gs://* ]]; then
        WORKSPACE_BUCKET="gs://${WORKSPACE_BUCKET}"
    fi

    # Download LD ref files from GCS if missing
    MISSING_LD=0
    for chr_num in $(seq 1 22); do
        if [ ! -f "${LD_REF_DIR}/chr${chr_num}_background_qc.bed" ]; then
            MISSING_LD=1
            break
        fi
    done

    if [ "${MISSING_LD}" -eq 1 ]; then
        log "Some LD ref files missing locally. Downloading from GCS..."
        GCS_LD="${WORKSPACE_BUCKET}/results/1-bg_snp/plink_qc/"
        mkdir -p "${LD_REF_DIR}"
        gsutil -u "${GOOGLE_PROJECT}" -m cp "${GCS_LD}*.bed" "${GCS_LD}*.bim" "${GCS_LD}*.fam" "${LD_REF_DIR}/" 2>&1 | tee -a "${LOG_FILE}" || true
        log "GCS download of LD ref complete"
    fi

    # Download .ma files from GCS if missing
    MISSING_MA=0
    for chr_num in $(seq 1 22); do
        if [ ! -f "${MA_DIR}/chr${chr_num}.ma" ]; then
            MISSING_MA=1
            break
        fi
    done

    if [ "${MISSING_MA}" -eq 1 ]; then
        log "Some .ma files missing locally. Downloading from GCS..."
        GCS_MA="${WORKSPACE_BUCKET}/results/2-locus_definition/ma/"
        mkdir -p "${MA_DIR}"
        gsutil -u "${GOOGLE_PROJECT}" -m cp "${GCS_MA}*.ma" "${MA_DIR}/" 2>&1 | tee -a "${LOG_FILE}" || true
        log "GCS download of .ma files complete"
    fi
fi

# ---------------------------------------------------------------------------
# Run Python script
# ---------------------------------------------------------------------------
log "Running: python3 ${PYTHON_SCRIPT} --config ${CONFIG_FILE} --output-dir ${OUTPUT_DIR} ${EXTRA_ARGS[*]}"
log ""

python3 "${PYTHON_SCRIPT}" \
    --config "${CONFIG_FILE}" \
    --output-dir "${OUTPUT_DIR}" \
    "${EXTRA_ARGS[@]}" 2>&1 | tee -a "${LOG_FILE}"

PYTHON_EXIT=$?

if [ "${PYTHON_EXIT}" -ne 0 ]; then
    log ""
    log "ERROR: Python script exited with code ${PYTHON_EXIT}"
    rm -f "${PID_FILE}"
    exit "${PYTHON_EXIT}"
fi

# ---------------------------------------------------------------------------
# Upload to GCS
# ---------------------------------------------------------------------------
log ""
log "============================================================"
log "  Upload to GCS"
log "============================================================"
if [ -n "${WORKSPACE_BUCKET}" ]; then
    GCS_DEST="${WORKSPACE_BUCKET}/results/2-locus_definition/"
    log "Uploading to: ${GCS_DEST}"

    # Upload signals and BED
    for f in "${OUTPUT_DIR}/all_independent_signals.tsv" "${OUTPUT_DIR}/target_loci.bed"; do
        if [ -f "$f" ]; then
            gsutil -u "${GOOGLE_PROJECT}" cp "$f" "${GCS_DEST}"
            log "  Uploaded: $(basename "$f")"
        fi
    done

    # Upload per-chrom cojo directories
    if [ -d "${OUTPUT_DIR}/cojo" ]; then
        gsutil -u "${GOOGLE_PROJECT}" -m cp -r "${OUTPUT_DIR}/cojo" "${GCS_DEST}"
        log "  Uploaded: cojo/"
    fi

    log "Upload complete"
else
    log "WARNING: WORKSPACE_BUCKET not set -- skipping GCS upload."
fi

# ---------------------------------------------------------------------------
# Final summary
# ---------------------------------------------------------------------------
END_SECONDS=$(date +%s)
ELAPSED=$(( END_SECONDS - START_SECONDS ))
ELAPSED_FMT=$(printf '%02dh %02dm %02ds' $((ELAPSED/3600)) $(( (ELAPSED%3600)/60 )) $((ELAPSED%60)))

log ""
log "============================================================"
log "  02b COMPLETED in ${ELAPSED_FMT}"
log "============================================================"

if [ -f "${OUTPUT_DIR}/all_independent_signals.tsv" ]; then
    N_SIGNALS=$(( $(wc -l < "${OUTPUT_DIR}/all_independent_signals.tsv") - 1 ))
    log "  Independent signals: ${N_SIGNALS}"
fi

if [ -f "${OUTPUT_DIR}/target_loci.bed" ]; then
    N_LOCI=$(( $(wc -l < "${OUTPUT_DIR}/target_loci.bed") - 1 ))
    log "  Target loci: ${N_LOCI}"
fi

log ""
log "Output files:"
log "  ${OUTPUT_DIR}/all_independent_signals.tsv"
log "  ${OUTPUT_DIR}/target_loci.bed"
log "  ${OUTPUT_DIR}/cojo/"
log "============================================================"

rm -f "${PID_FILE}"
