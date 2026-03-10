#!/bin/bash
set -eo pipefail

# ---------------------------------------------------------------------------
# 02c_export_top_gwas_snps.sh
#
# Export top GWAS SNPs from per-chromosome .ma files (produced by 02a).
# Filters to SNPs below a configurable p-value threshold and writes a
# single consolidated TSV for cross-referencing with published MS loci.
#
# No Hail/Spark required -- pure pandas, runs on the VM in seconds.
#
# Usage:
#   bash bash/02c_export_top_gwas_snps.sh
#   bash bash/02c_export_top_gwas_snps.sh --p-threshold=5e-6
#   bash bash/02c_export_top_gwas_snps.sh --force
#
# Monitor:
#   tail -f logs/02c_top_gwas_snps_*.log
# ---------------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
PYTHON_SCRIPT="${PROJECT_DIR}/python/02c_export_top_gwas_snps.py"
CONFIG_FILE="${PROJECT_DIR}/config/config.json"
LOG_DIR="${PROJECT_DIR}/logs"

mkdir -p "$LOG_DIR"

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="${LOG_DIR}/02c_top_gwas_snps_${TIMESTAMP}.log"

# ---------------------------------------------------------------------------
# Parse arguments (pass through to Python script)
# ---------------------------------------------------------------------------
EXTRA_ARGS=()
FORCE=""

for arg in "$@"; do
    case "$arg" in
        --force)
            FORCE="--force"
            EXTRA_ARGS+=("$arg")
            ;;
        --p-threshold=*)
            EXTRA_ARGS+=("--p-threshold" "${arg#*=}")
            ;;
        --output=*)
            EXTRA_ARGS+=("--output" "${arg#*=}")
            ;;
        --ma-dir=*)
            EXTRA_ARGS+=("--ma-dir" "${arg#*=}")
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
log "  02c: Export Top GWAS SNPs"
log "============================================================"
log "  Config : ${CONFIG_FILE}"
log "  Log    : ${LOG_FILE}"
log "  Args   : ${EXTRA_ARGS[*]}"
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
    exit 1
fi

# Check that .ma files exist (spot-check chr1)
MA_DIR=$(python3 -c "
import json
c = json.load(open('${CONFIG_FILE}'))
print(c['outputs'].get('gwas_ma_dir', 'results/2-locus_definition/ma'))
")

if [ ! -f "${MA_DIR}/chr1.ma" ]; then
    log "ERROR: .ma files not found in ${MA_DIR}"
    log "  Run: bash bash/02a_export_gwas_ma.sh"
    exit 1
fi

log "Pre-flight checks PASSED"
log ""

# ---------------------------------------------------------------------------
# Run
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

log ""
log "02c complete"
