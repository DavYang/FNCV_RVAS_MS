#!/bin/bash
set -eo pipefail

# ---------------------------------------------------------------------------
# 01_build_background_snps.sh
#
# Wrapper script for the background SNP PLINK pipeline (Regenie Step 1).
# Loops over autosomes (chr1-22), calling the Python script ONCE per
# chromosome so each gets a fresh Hail/JVM session. This avoids JVM
# state accumulation and OOM crashes seen in multi-chromosome runs.
#
# Features:
#   - Computes per-chromosome SNP targets proportionally (same formula as Python).
#   - Fault tolerant: failure on one chromosome does not block others.
#   - Resume support: skips chromosomes with existing summary.json.
#   - Writes pipeline_summary.json at the end with overall status.
#
# Usage:
#   nohup bash bash/01_build_background_snps.sh > /dev/null 2>&1 &
#
# Resume a previous run (point at existing output dir):
#   RESUME_OUTPUT_DIR=gs://bucket/results/FNCV_RVAS_MS/500K_background_snps \
#     nohup bash bash/01_build_background_snps.sh > /dev/null 2>&1 &
#
# Monitor:
#   tail -f logs/01_background_snps_*.log           # main pipeline log
#   tail -f logs/01_background_snps_chr21_*.log      # per-chromosome log
#   cat logs/01_background_snps.pid
#   ps -p $(cat logs/01_background_snps.pid) -o pid,etime,cmd
#   kill $(cat logs/01_background_snps.pid)
# ---------------------------------------------------------------------------

trap '' HUP

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
PYTHON_SCRIPT="${PROJECT_DIR}/python/01_build_background_snps.py"
CONFIG_FILE="${PROJECT_DIR}/config/config.json"
LOG_DIR="${PROJECT_DIR}/logs"

mkdir -p "$LOG_DIR"

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="${LOG_DIR}/01_background_snps_${TIMESTAMP}.log"
PID_FILE="${LOG_DIR}/01_background_snps.pid"

exec < /dev/null
exec > >(tee -a "${LOG_FILE}") 2>&1

echo $$ > "${PID_FILE}"

START_SECONDS=$(date +%s)

# ---------------------------------------------------------------------------
# Pre-flight checks
# ---------------------------------------------------------------------------
echo "============================================================"
echo "Background SNP PLINK Pipeline (per-chromosome mode)"
echo "============================================================"
echo "Started at : $(date)"
echo "PID        : $$"
echo "PID file   : ${PID_FILE}"
echo "Log file   : ${LOG_FILE}"
echo ""

if [ -z "$WORKSPACE_BUCKET" ]; then
    echo "ERROR: WORKSPACE_BUCKET environment variable not set"
    exit 1
fi
# Ensure gs:// prefix
if [[ "$WORKSPACE_BUCKET" != gs://* ]]; then
    WORKSPACE_BUCKET="gs://${WORKSPACE_BUCKET}"
fi
echo "Workspace  : $WORKSPACE_BUCKET"

if [ -z "$GOOGLE_PROJECT" ]; then
    echo "WARNING: GOOGLE_PROJECT not set. Requester Pays buckets may fail."
else
    echo "Project    : $GOOGLE_PROJECT"
fi

if [ ! -f "$PYTHON_SCRIPT" ]; then
    echo "ERROR: Python script not found: $PYTHON_SCRIPT"
    exit 1
fi

if [ ! -f "$CONFIG_FILE" ]; then
    echo "ERROR: Config file not found: $CONFIG_FILE"
    exit 1
fi

# ---------------------------------------------------------------------------
# Read config
# ---------------------------------------------------------------------------
TARGET_SNPS=$(python3 -c "import json; c=json.load(open('${CONFIG_FILE}')); print(c['sampling']['target_total_snps'])")
TEST_MODE=$(python3 -c "import json; c=json.load(open('${CONFIG_FILE}')); print(c['params'].get('test_mode', False))")
TEST_CHR=$(python3 -c "import json; c=json.load(open('${CONFIG_FILE}')); print(c['params'].get('test_chromosome', 'None'))")

echo ""
echo "Target SNPs: ${TARGET_SNPS}"
if [ "$TEST_MODE" = "True" ]; then
    echo "Mode       : TEST (${TEST_CHR} only)"
    CHROMOSOMES=("$TEST_CHR")
else
    echo "Mode       : PRODUCTION (chr1-22)"
    CHROMOSOMES=(chr{1..22})
fi

# ---------------------------------------------------------------------------
# Build output directory name (stable across re-runs for resume support)
# ---------------------------------------------------------------------------
if [ "$TARGET_SNPS" -ge 1000000 ]; then
    SNP_LABEL="$((TARGET_SNPS / 1000000))M"
elif [ "$TARGET_SNPS" -ge 1000 ]; then
    SNP_LABEL="$((TARGET_SNPS / 1000))K"
else
    SNP_LABEL="${TARGET_SNPS}"
fi

# Allow explicit output dir override via env var or use a stable (no date) default
if [ -n "$RESUME_OUTPUT_DIR" ]; then
    OUTPUT_DIR="$RESUME_OUTPUT_DIR"
else
    OUTPUT_DIR="${WORKSPACE_BUCKET}/results/FNCV_RVAS_MS/${SNP_LABEL}_background_snps"
fi

echo "Output dir : ${OUTPUT_DIR}"
echo "------------------------------------------------------------"
echo ""

# ---------------------------------------------------------------------------
# Compute per-chromosome targets (proportional to chromosome size)
# GRCh38 sizes in bp
# ---------------------------------------------------------------------------
declare -A CHR_SIZES
CHR_SIZES=(
    [chr1]=248956422  [chr2]=242193529  [chr3]=198295559
    [chr4]=190214555  [chr5]=181538259  [chr6]=170805979
    [chr7]=159345973  [chr8]=145138636  [chr9]=138394717
    [chr10]=133797422 [chr11]=135086622 [chr12]=133275309
    [chr13]=114364328 [chr14]=107043718 [chr15]=101991189
    [chr16]=90338345  [chr17]=83257441  [chr18]=80373285
    [chr19]=58617616  [chr20]=64444167  [chr21]=46709983
    [chr22]=50818468
)

TOTAL_GENOME_SIZE=0
for chr_name in "${CHROMOSOMES[@]}"; do
    TOTAL_GENOME_SIZE=$((TOTAL_GENOME_SIZE + CHR_SIZES[$chr_name]))
done

declare -A CHR_TARGETS
ALLOCATED=0
for chr_name in "${CHROMOSOMES[@]}"; do
    size=${CHR_SIZES[$chr_name]}
    # Integer proportional allocation: target * size / total_size
    chr_target=$(python3 -c "print(max(1, int(${TARGET_SNPS} * ${size} / ${TOTAL_GENOME_SIZE})))")
    CHR_TARGETS[$chr_name]=$chr_target
    ALLOCATED=$((ALLOCATED + chr_target))
done

# Distribute remainder to largest chromosomes
REMAINDER=$((TARGET_SNPS - ALLOCATED))
if [ "$REMAINDER" -gt 0 ]; then
    for chr_name in chr1 chr2 chr3 chr4 chr5; do
        if [ "$REMAINDER" -le 0 ]; then break; fi
        CHR_TARGETS[$chr_name]=$((CHR_TARGETS[$chr_name] + 1))
        REMAINDER=$((REMAINDER - 1))
    done
fi

echo "Per-chromosome targets:"
for chr_name in "${CHROMOSOMES[@]}"; do
    printf "  %-8s %10d\n" "$chr_name" "${CHR_TARGETS[$chr_name]}"
done
echo ""

# ---------------------------------------------------------------------------
# Loop over chromosomes, one fresh Python/Hail process each
# ---------------------------------------------------------------------------
SUCCEEDED=()
FAILED=()
SKIPPED=()

TOTAL_CHROMS=${#CHROMOSOMES[@]}
CURRENT=0

for chr_name in "${CHROMOSOMES[@]}"; do
    CURRENT=$((CURRENT + 1))
    chr_target=${CHR_TARGETS[$chr_name]}
    chr_start=$(date +%s)

    echo "============================================================"
    echo "[${CURRENT}/${TOTAL_CHROMS}] Processing ${chr_name} (target: ${chr_target})"
    echo "============================================================"

    # Check if already completed (resume support)
    SUMMARY_FILE="${OUTPUT_DIR}/${chr_name}/summary.json"
    # Use gsutil to check existence since we can't use hfs in bash
    if gsutil -q stat "${SUMMARY_FILE}" 2>/dev/null; then
        echo "${chr_name}: summary.json exists, skipping (already completed)"
        SKIPPED+=("$chr_name")
        echo ""
        continue
    fi

    # Per-chromosome log captures ALL Python/Spark/JVM output
    CHR_LOG="${LOG_DIR}/01_background_snps_${chr_name}_${TIMESTAMP}.log"
    HAIL_LOG="/tmp/hail_${chr_name}.log"

    echo "  Python log : ${CHR_LOG}"
    echo "  Hail log   : ${HAIL_LOG}"
    echo "  Started at : $(date)"

    CHR_EXIT=0
    python3 "$PYTHON_SCRIPT" \
        --chrom "$chr_name" \
        --target "$chr_target" \
        --output-dir "$OUTPUT_DIR" \
        --config "$CONFIG_FILE" \
        </dev/null >"${CHR_LOG}" 2>&1 \
    || CHR_EXIT=$?

    chr_end=$(date +%s)
    chr_elapsed=$((chr_end - chr_start))
    chr_min=$((chr_elapsed / 60))
    chr_sec=$((chr_elapsed % 60))

    if [ $CHR_EXIT -eq 0 ]; then
        echo "${chr_name}: SUCCEEDED in ${chr_min}m ${chr_sec}s"
        SUCCEEDED+=("$chr_name")
    else
        echo "${chr_name}: FAILED (exit code ${CHR_EXIT}) after ${chr_min}m ${chr_sec}s"

        # Diagnose the failure mode
        if [ $CHR_EXIT -gt 128 ]; then
            SIGNAL=$((CHR_EXIT - 128))
            echo "  >>> Process killed by signal ${SIGNAL}"
            case $SIGNAL in
                9)  echo "  >>> SIGKILL: likely OOM-killer or external kill" ;;
                11) echo "  >>> SIGSEGV: JVM segmentation fault (memory corruption)" ;;
                6)  echo "  >>> SIGABRT: JVM abort (assertion failure or OOM)" ;;
                *)  echo "  >>> Signal ${SIGNAL}: unexpected termination" ;;
            esac
            # Check dmesg for OOM kills (may need root, best effort)
            dmesg -T 2>/dev/null | grep -i "oom\|killed process" | tail -5 \
                && echo "  >>> (dmesg OOM entries above)" \
                || true
        elif [ $CHR_EXIT -eq 1 ]; then
            echo "  >>> Python exception (check log for traceback)"
        else
            echo "  >>> Unexpected exit code ${CHR_EXIT}"
        fi

        # Show tail of per-chromosome Python log
        echo "--- Last 40 lines of Python log (${CHR_LOG}) ---"
        tail -40 "${CHR_LOG}" 2>/dev/null || echo "(Python log not found)"
        echo "--- End Python log ---"

        # Show tail of Hail JVM log
        echo "--- Last 20 lines of Hail log (${HAIL_LOG}) ---"
        tail -20 "${HAIL_LOG}" 2>/dev/null || echo "(Hail log not found)"
        echo "--- End Hail log ---"

        FAILED+=("$chr_name")
    fi
    echo ""

    # Memory mitigation: flush filesystem caches and pause between chromosomes
    # to allow OS to fully reclaim memory from the terminated JVM process
    sync
    echo 3 | sudo tee /proc/sys/vm/drop_caches > /dev/null 2>&1 || true
    echo "  [memory cleanup] Sleeping 60s between chromosomes ..."
    sleep 60
done

# ---------------------------------------------------------------------------
# Final summary
# ---------------------------------------------------------------------------
END_SECONDS=$(date +%s)
ELAPSED=$((END_SECONDS - START_SECONDS))
ELAPSED_MIN=$((ELAPSED / 60))
ELAPSED_SEC=$((ELAPSED % 60))

echo ""
echo "============================================================"
echo "PIPELINE SUMMARY"
echo "============================================================"
echo "Total time : ${ELAPSED_MIN}m ${ELAPSED_SEC}s"
echo "Succeeded  : ${#SUCCEEDED[@]} - ${SUCCEEDED[*]:-none}"
echo "Failed     : ${#FAILED[@]} - ${FAILED[*]:-none}"
echo "Skipped    : ${#SKIPPED[@]} - ${SKIPPED[*]:-none}"
echo "Output dir : ${OUTPUT_DIR}"
echo "Finished at: $(date)"
echo "============================================================"

# Write pipeline summary JSON to GCS
_json_array() {
    local arr=("$@")
    if [ ${#arr[@]} -eq 0 ]; then
        echo '[]'
    else
        printf '['; printf '"%s",' "${arr[@]}" | sed 's/,$//'; printf ']'
    fi
}

PIPELINE_SUMMARY=$(cat <<EOF
{
  "total_chromosomes": ${TOTAL_CHROMS},
  "succeeded": ${#SUCCEEDED[@]},
  "failed": ${#FAILED[@]},
  "skipped": ${#SKIPPED[@]},
  "succeeded_list": $(_json_array "${SUCCEEDED[@]}"),
  "failed_list": $(_json_array "${FAILED[@]}"),
  "skipped_list": $(_json_array "${SKIPPED[@]}"),
  "output_dir": "${OUTPUT_DIR}",
  "total_time_seconds": ${ELAPSED},
  "timestamp": "$(date -Iseconds)"
}
EOF
)

echo "$PIPELINE_SUMMARY" | gsutil -q cp - "${OUTPUT_DIR}/pipeline_summary.json" 2>/dev/null \
    || echo "WARNING: Failed to write pipeline_summary.json to GCS"

rm -f "${PID_FILE}"

if [ ${#FAILED[@]} -gt 0 ]; then
    exit 1
else
    exit 0
fi
