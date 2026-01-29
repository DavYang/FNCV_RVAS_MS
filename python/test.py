#!/usr/bin/env python3
import os
import sys
import hail as hl
import pandas as pd
from datetime import datetime
from utils import load_config, init_hail, setup_logger

logger = setup_logger("smoke_test")

def main():
    # 1. Load config
    config = load_config("config/config.json")
    
    # 2. Init Hail with your 32g settings
    init_hail("smoke_test", driver_mem="32g", reference="GRCh38")

    # 3. Path Setup
    workspace_bucket = os.environ.get('WORKSPACE_BUCKET')
    if not workspace_bucket:
        logger.error("WORKSPACE_BUCKET not set. Please run in an AoU terminal/notebook.")
        sys.exit(1)
    
    base_data_dir = workspace_bucket if workspace_bucket.startswith('gs://') else f"gs://{workspace_bucket}"
    test_output_dir = f"{base_data_dir}/smoke_test_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    
    logger.info(f"TESTING PERMISSIONS: Writing to {test_output_dir}")

    # 4. Use the AoU splitMT (Standard for WGS analysis)
    # Note: We use the common splitMT path directly if not in config
    SPLIT_MT_PATH = config['inputs'].get('split_mt', 'gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/acaf_threshold/splitMT/hail.mt')
    ANCESTRY_PATH = config['inputs']['ancestry_pred']
    
    # Load 10 EUR samples for the test
    logger.info("Accessing ancestry file...")
    with hl.hadoop_open(ANCESTRY_PATH, 'r') as f:
        ancestry_df = pd.read_csv(f, sep="\t")
    eur_sample_subset = list(ancestry_df[ancestry_df['ancestry_pred'] == 'eur']['research_id'].astype(str))[:10]

    # 5. Load and Slice Data
    logger.info(f"Reading MatrixTable: {SPLIT_MT_PATH}")
    mt = hl.read_matrix_table(SPLIT_MT_PATH)

    # Filter to a stable, variant-dense region on Chr22 (smaller than Chr1 for a faster test)
    logger.info("Extracting test interval on Chr22...")
    test_interval = [hl.parse_locus_interval('chr22:20000000-20050000', reference_genome='GRCh38')]
    
    mt_test = hl.filter_intervals(mt, test_interval)
    mt_test = mt_test.filter_cols(hl.literal(set(eur_sample_subset)).contains(mt_test.s))

    # Force 1 variant for the smoke test
    if mt_test.count_rows() > 0:
        mt_test = mt_test.head(1)
    else:
        logger.error("No variants found in test interval. Check genomic coordinates.")
        sys.exit(1)

    # 6. Export to PLINK
    # This checks if your workspace bucket is writable by the Spark workers
    plink_test_base = f"{test_output_dir}/smoke_test_output"
    logger.info(f"Exporting test PLINK files to: {plink_test_base}")
    
    hl.export_plink(
        mt_test,
        plink_test_base,
        ind_id=mt_test.s,
        fam_id=mt_test.s
    )

    logger.info("SUCCESS: Smoke test completed. PLINK files generated.")
    logger.info(f"Verify files: gsutil ls {test_output_dir}")

if __name__ == "__main__":
    main()