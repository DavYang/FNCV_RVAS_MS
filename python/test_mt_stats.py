#!/usr/bin/env python3
import sys
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
        logger.info("Test completed successfully.")
        
    except Exception as e:
        logger.error(f"Failed to read MatrixTable: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()