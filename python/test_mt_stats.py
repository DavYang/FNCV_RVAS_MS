#!/usr/bin/env python3
import sys
import os
import hail as hl
from utils import load_config, init_hail, setup_logger

logger = setup_logger("test_mt")

def main():
    config = load_config("config/config.json")
    
    # Initialize Hail
    init_hail("test_mt_stats", driver_mem="8g", reference="GRCh38")
    
    # Path to check
    checkpoint_path = f"{config['outputs']['intermediate_dir']}/wgs_eur_intervals.mt"
    
    logger.info(f"Checking for MatrixTable at: {checkpoint_path}")
    
    if not hl.hadoop_exists(checkpoint_path):
        logger.error("❌ MatrixTable NOT found at path!")
        sys.exit(1)
        
    logger.info("✅ MatrixTable directory exists.")
    
    try:
        logger.info("Attempting to read MatrixTable...")
        mt = hl.read_matrix_table(checkpoint_path)
        
        # Gather basic stats
        n_rows = mt.count_rows()
        n_cols = mt.count_cols()
        n_partitions = mt.n_partitions()
        
        logger.info("=======================================")
        logger.info("       MatrixTable Statistics")
        logger.info("=======================================")
        logger.info(f"Path:       {checkpoint_path}")
        logger.info(f"Variants:   {n_rows:,}")
        logger.info(f"Samples:    {n_cols:,}")
        logger.info(f"Partitions: {n_partitions:,}")
        
        # Check entry fields
        entry_fields = list(mt.entry)
        logger.info(f"Entry fields: {entry_fields}")
        
        # Check row fields (info, locus, alleles)
        row_fields = list(mt.row)
        logger.info(f"Row fields:   {row_fields[:5]}...")  # First 5 fields
        
        logger.info("=======================================")
        
        # --- TEST EXPORT TO BUCKET ---
        workspace_bucket = os.environ.get('WORKSPACE_BUCKET')
        if not workspace_bucket:
            logger.error("❌ WORKSPACE_BUCKET environment variable not found!")
            sys.exit(1)
            
        logger.info(f"Detected Workspace Bucket: {workspace_bucket}")
        test_export_path = f"{workspace_bucket}/tmp/test_plink_export/test_sample"
        
        logger.info(f"Testing PLINK export to: {test_export_path}")
        
        # REPARTITION TO FIX PERFORMANCE BLOCK
        logger.info("Repartitioning to 100 partitions before sampling/export...")
        mt = mt.repartition(100)
        
        # Downsample drastically for test (0.001% of rows)
        mt_small = mt.sample_rows(0.00001, seed=42)
        
        # Force a count to ensure repartition is processed
        n_small = mt_small.count_rows()
        logger.info(f"Exporting small sample ({n_small} variants)...")
        
        if n_small == 0:
             logger.warning("Sample too small, taking first 100 rows instead.")
             mt_small = mt.head(100)
        
        hl.export_plink(
            mt_small,
            test_export_path,
            ind_id=mt_small.s,
            fam_id=mt_small.s
        )
        
        # Verify files exist
        extensions = ['.bed', '.bim', '.fam']
        missing = []
        for ext in extensions:
            path = test_export_path + ext
            if hl.hadoop_exists(path):
                logger.info(f"✅ Created: {path}")
            else:
                logger.error(f"❌ FAILED to create: {path}")
                missing.append(path)
                
        if missing:
            logger.error("Export test FAILED.")
            sys.exit(1)
            
        logger.info("✅ PLINK export test successful.")
        
        logger.info("=======================================")
        logger.info("Test completed successfully.")
        
    except Exception as e:
        logger.error(f"Test failed with error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()