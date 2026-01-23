#!/bin/bash
set -e

# --- Usage Help ---
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <gcs_input_base> <path_to_ancestry_tsv> <gcs_output_dir>"
    echo "Example: $0 gs://fc-aou.../plink_bed /path/to/ancestry.tsv gs://my-bucket/data/common_variants"
    exit 1
fi

# --- Inputs ---
GCS_INPUT_BASE="$1"
ANCESTRY_FILE="$2"
GCS_OUTPUT_DIR="$3"

# --- Configuration ---
# QC Thresholds
GENO=0.05   # Filter VARIANTS with >5% missingness (proxy for low quality)
MIND=0.05   # Filter SAMPLES with >5% missingness (low quality samples)
HWE=1e-6    # Filter variants failing HWE test (p < 1e-6)

# --- Local Setup ---
LOCAL_WORK_DIR="tmp_work_dir"
mkdir -p "$LOCAL_WORK_DIR"
EUR_SAMPLES_FILE="${LOCAL_WORK_DIR}/eur_samples.txt"

# Logging setup
LOG_FILE="process_chrom_plink_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "=========================================================="
echo "   Processing Common Variants (Cloud-Optimized)"
echo "=========================================================="
echo "Input:        ${GCS_INPUT_BASE}"
echo "Output (GCS): ${GCS_OUTPUT_DIR}"
echo "Work Dir:     ${LOCAL_WORK_DIR}"
echo "----------------------------------------------------------"
echo "QC Parameters:"
echo "  --geno ${GENO}  (Exclude variants missing in >${GENO//./}% of samples)"
echo "  --mind ${MIND}  (Exclude samples missing >${MIND//./}% of variants)"
echo "  --hwe  ${HWE}   (Exclude variants failing HWE)"
echo "=========================================================="

# 1. Extract European Samples
echo "[1/5] Extracting EUR sample IDs..."
awk -F'\t' '$2 == "eur" {print $1, $1}' "${ANCESTRY_FILE}" > "${EUR_SAMPLES_FILE}"
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
    gsutil -q cp "${GCS_IN}.bed" "${LOCAL_IN}.bed"
    gsutil -q cp "${GCS_IN}.bim" "${LOCAL_IN}.bim"
    gsutil -q cp "${GCS_IN}.fam" "${LOCAL_IN}.fam"
    
    # B. Run QC
    # Note: --geno filters low-quality variants (high missingness)
    #       --mind filters low-quality samples (high missingness)
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
    gsutil -q cp "${LOCAL_QC}.bed" "${GCS_QC}.bed"
    gsutil -q cp "${LOCAL_QC}.bim" "${GCS_QC}.bim"
    gsutil -q cp "${LOCAL_QC}.fam" "${GCS_QC}.fam"
    
    # D. Cleanup Everything
    echo "  > Cleaning up local files..."
    rm "${LOCAL_IN}"* "${LOCAL_QC}"*
    
    echo "  > Done."
done

# 3. Download All QC'd Files for Merge
echo "=========================================================="
echo "[3/5] Downloading all QC'd chromosomes for merging..."
mkdir -p "${LOCAL_WORK_DIR}/merge_stage"
gsutil -m cp -r "${GCS_OUTPUT_DIR}/chroms/*" "${LOCAL_WORK_DIR}/merge_stage/"

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
gsutil cp "${MERGED_PREFIX}.bed" "${GCS_OUTPUT_DIR}/final/"
gsutil cp "${MERGED_PREFIX}.bim" "${GCS_OUTPUT_DIR}/final/"
gsutil cp "${MERGED_PREFIX}.fam" "${GCS_OUTPUT_DIR}/final/"
gsutil cp "${LOG_FILE}" "${GCS_OUTPUT_DIR}/"

# Cleanup
echo "Cleaning up workspace..."
rm -rf "${LOCAL_WORK_DIR}"

echo "DONE! Final data in: ${GCS_OUTPUT_DIR}/final/"
