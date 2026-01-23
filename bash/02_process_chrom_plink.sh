#!/bin/bash
set -e

# --- Configuration ---
# Hardcoded Paths
GCS_INPUT_BASE="gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/acaf_threshold/plink_bed"
GCS_ANCESTRY_FILE="gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv"

# Output to Workspace Bucket (must be set in environment)
if [ -z "${WORKSPACE_BUCKET}" ]; then
    echo "ERROR: WORKSPACE_BUCKET environment variable is not set."
    exit 1
fi
GCS_OUTPUT_DIR="${WORKSPACE_BUCKET}/data/common_variants_plink"

# QC Thresholds
GENO=0.05
MIND=0.05
HWE=1e-6

# --- Local Setup ---
LOCAL_WORK_DIR="tmp_work_dir"
mkdir -p "$LOCAL_WORK_DIR"
LOCAL_ANCESTRY_FILE="${LOCAL_WORK_DIR}/ancestry_preds.tsv"
EUR_SAMPLES_FILE="${LOCAL_WORK_DIR}/eur_samples.txt"

# Logging setup
LOG_FILE="process_chrom_plink_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "=========================================================="
echo "   Processing Common Variants (Cloud-Optimized)"
echo "=========================================================="
echo "Input Base:   ${GCS_INPUT_BASE}"
echo "Ancestry:     ${GCS_ANCESTRY_FILE}"
echo "Output (GCS): ${GCS_OUTPUT_DIR}"
echo "Work Dir:     ${LOCAL_WORK_DIR}"
echo "Project:      ${GOOGLE_PROJECT}"
echo "----------------------------------------------------------"
echo "QC Parameters:"
echo "  --geno ${GENO}  (Exclude variants missing in >${GENO//./}% of samples)"
echo "  --mind ${MIND}  (Exclude samples missing >${MIND//./}% of variants)"
echo "  --hwe  ${HWE}   (Exclude variants failing HWE)"
echo "=========================================================="

if [ -z "${GOOGLE_PROJECT}" ]; then
    echo "WARNING: GOOGLE_PROJECT is not set. Requester Pays buckets may fail."
fi

# 0. Download Ancestry File
echo "[0/5] Downloading ancestry file..."
gsutil -u "${GOOGLE_PROJECT}" cp "${GCS_ANCESTRY_FILE}" "${LOCAL_ANCESTRY_FILE}"

# 1. Extract European Samples
echo "[1/5] Extracting EUR sample IDs..."
awk -F'\t' '$2 == "eur" {print $1, $1}' "${LOCAL_ANCESTRY_FILE}" > "${EUR_SAMPLES_FILE}"
EUR_COUNT=$(wc -l < "${EUR_SAMPLES_FILE}")
echo "      Found ${EUR_COUNT} European samples."

# 2. Process Each Chromosome (Download -> QC -> Upload -> Delete)
echo "[2/5] Processing Chromosomes 1-22..."

for chr in {1..22}; do
    echo "------------------------------------------------"
    echo "Processing Chromosome ${chr}..."
    
    # Paths
    GCS_IN="${GCS_INPUT_BASE}/chr${chr}"
    LOCAL_IN="${LOCAL_WORK_DIR}/chr${chr}"
    LOCAL_QC="${LOCAL_WORK_DIR}/chr${chr}_qc"
    GCS_QC="${GCS_OUTPUT_DIR}/chroms/chr${chr}_qc"
    
    # A. Download Raw
    echo "  > Downloading raw data..."
    gsutil -u "${GOOGLE_PROJECT}" -q cp "${GCS_IN}.bed" "${LOCAL_IN}.bed"
    gsutil -u "${GOOGLE_PROJECT}" -q cp "${GCS_IN}.bim" "${LOCAL_IN}.bim"
    gsutil -u "${GOOGLE_PROJECT}" -q cp "${GCS_IN}.fam" "${LOCAL_IN}.fam"
    
    # B. Run QC
    echo "  > Running QC (EUR subset + GENO/MIND/HWE)..."
    plink2 \
        --bfile "${LOCAL_IN}" \
        --keep "${EUR_SAMPLES_FILE}" \
        --geno "${GENO}" \
        --mind "${MIND}" \
        --hwe "${HWE}" \
        --make-bed \
        --out "${LOCAL_QC}" \
        --silent
        
    # C. Upload QC'd Result
    echo "  > Uploading QC results to GCS..."
    gsutil -u "${GOOGLE_PROJECT}" -q cp "${LOCAL_QC}.bed" "${GCS_QC}.bed"
    gsutil -u "${GOOGLE_PROJECT}" -q cp "${LOCAL_QC}.bim" "${GCS_QC}.bim"
    gsutil -u "${GOOGLE_PROJECT}" -q cp "${LOCAL_QC}.fam" "${GCS_QC}.fam"
    
    # D. Cleanup Everything
    echo "  > Cleaning up local files..."
    rm "${LOCAL_IN}"* "${LOCAL_QC}"*
    
    echo "  > Done."
done

# 3. Download All QC'd Files for Merge
echo "=========================================================="
echo "[3/5] Downloading all QC'd chromosomes for merging..."
mkdir -p "${LOCAL_WORK_DIR}/merge_stage"
gsutil -u "${GOOGLE_PROJECT}" -m cp -r "${GCS_OUTPUT_DIR}/chroms/*" "${LOCAL_WORK_DIR}/merge_stage/"

# 4. Create Merge List
MERGE_LIST="${LOCAL_WORK_DIR}/merge_list.txt"
rm -f "$MERGE_LIST"
for chr in {1..22}; do
    echo "${LOCAL_WORK_DIR}/merge_stage/chr${chr}_qc" >> "$MERGE_LIST"
done

# 5. Merge
echo "[4/5] Merging datasets..."
MERGED_PREFIX="${LOCAL_WORK_DIR}/eur_common_merged"
plink2 \
    --pmerge-list "${MERGE_LIST}" bfile \
    --make-bed \
    --out "${MERGED_PREFIX}"

# 6. Final Upload
echo "[5/5] Uploading final merged dataset to GCS..."
gsutil -u "${GOOGLE_PROJECT}" cp "${MERGED_PREFIX}.bed" "${GCS_OUTPUT_DIR}/final/"
gsutil -u "${GOOGLE_PROJECT}" cp "${MERGED_PREFIX}.bim" "${GCS_OUTPUT_DIR}/final/"
gsutil -u "${GOOGLE_PROJECT}" cp "${MERGED_PREFIX}.fam" "${GCS_OUTPUT_DIR}/final/"
gsutil -u "${GOOGLE_PROJECT}" cp "${LOG_FILE}" "${GCS_OUTPUT_DIR}/"

# Cleanup
echo "Cleaning up workspace..."
rm -rf "${LOCAL_WORK_DIR}"

echo "DONE! Final data in: ${GCS_OUTPUT_DIR}/final/"
