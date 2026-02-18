#!/bin/bash
set -e

# ---------------------------------------------------------------------------
# Copy per-chromosome PLINK files from GCS into a single local directory
# ---------------------------------------------------------------------------
# Usage:
#   bash bash/copy_plink_from_gcs.sh
#
# Outputs all .bed/.bim/.fam files into:
#   results/1-bg_snp/plink/
# ---------------------------------------------------------------------------

GCS_BASE="gs://fc-secure-b43840eb-548f-464d-bece-31ac7a969abd/results/FNCV_RVAS_MS/1M_background_snps"
LOCAL_DIR="/home/jupyter/FNCV_RVAS_MS/results/1-bg_snp/plink_no-qc"

mkdir -p "${LOCAL_DIR}"

echo "============================================================"
echo "Copying PLINK files from GCS to ${LOCAL_DIR}"
echo "============================================================"
echo ""

for chr_num in $(seq 1 22); do
    chr_name="chr${chr_num}"
    echo "[${chr_name}] Copying .bed/.bim/.fam ..."
    gsutil -m cp \
        "${GCS_BASE}/${chr_name}/${chr_name}_background.bed" \
        "${GCS_BASE}/${chr_name}/${chr_name}_background.bim" \
        "${GCS_BASE}/${chr_name}/${chr_name}_background.fam" \
        "${LOCAL_DIR}/"
    echo "[${chr_name}] Done"
done

echo ""
echo "============================================================"
echo "All PLINK files copied to ${LOCAL_DIR}"
echo "============================================================"
echo ""
echo "File count: $(ls -1 "${LOCAL_DIR}"/*.{bed,bim,fam} 2>/dev/null | wc -l) files"
ls -lh "${LOCAL_DIR}/"
