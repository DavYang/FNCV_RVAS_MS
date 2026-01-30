#!/usr/bin/env python3
import os
import sys
import hail as hl
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
    SPLIT_MT_PATH = config['inputs'].get('split_mt', 'gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/acaf_threshold/splitMT/hail.mt')
    
    # 5. Load and Slice Data
    logger.info(f"Reading MatrixTable: {SPLIT_MT_PATH}")
    mt = hl.read_matrix_table(SPLIT_MT_PATH)

    # Filter to a predetermined region on Chr1 for the smoke test
    logger.info("Extracting test interval on Chr1...")
    test_interval = [hl.parse_locus_interval('chr1:10000000-10100000', reference_genome='GRCh38')]
    mt_test = hl.filter_intervals(mt, test_interval)
    
    # Sample 1 variant immediately
    mt_test = mt_test.head(1)
    
    # Use first 10 samples from the dataset (fastest for smoke test)
    logger.info("Selecting first 10 samples...")
    sample_ids = mt_test.cols().s.take(10)
    
    # Filter to those 10 samples
    mt_test = mt_test.filter_cols(hl.literal(set(sample_ids)).contains(mt_test.s))

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