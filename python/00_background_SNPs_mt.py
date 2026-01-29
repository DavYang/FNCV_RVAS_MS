#!/usr/bin/env python3
import os
import sys
import random
import hail as hl
import pandas as pd
from datetime import datetime
from utils import load_config, init_hail, setup_logger

logger = setup_logger("data_prep")

def perform_plink_export(mt, output_dir, config, timestamp):
    """Export MatrixTable to PLINK with dated directory and optimizations."""
    
    # Export base path
    plink_base = f"{output_dir}/eur_common_snps_sampled"
    
    logger.info(f"Exporting to PLINK: {plink_base}")
    
    # OPTIMIZATION: Repartition before export to prevent memory issues
    # More partitions = smaller tasks = less memory per executor
    logger.info("Repartitioning to 500 partitions for export...")
    mt = mt.repartition(500)
    
    hl.export_plink(
        mt,
        plink_base,
        ind_id=mt.s,
        fam_id=mt.s
    )
    
    # Write export metadata
    metadata_file = f"{output_dir}/export_metadata.txt"
    try:
        with hl.hadoop_open(metadata_file, 'w') as f:
            f.write(f"Export timestamp: {timestamp}\n")
            f.write(f"Variants: {mt.count_rows()}\n")
            f.write(f"Samples: {mt.count_cols()}\n")
            f.write(f"Source: AoU splitMT (EUR samples)\n")
            f.write(f"Config: {config}\n")
    except Exception as e:
        logger.warning(f"Failed to write metadata file to {metadata_file}: {e}")
    
    logger.info(f"PLINK export complete: {plink_base}")
    logger.info(f"Output directory: {output_dir}")

def perform_filtering(config, params):
    """Load pre-processed MT and filter to EUR samples."""
    ANCESTRY_PATH = config['inputs']['ancestry_pred']
    SPLIT_MT_PATH = config['inputs'].get('split_mt', 
                                          'gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/splitMT/splitMT.mt')
    
    # 1. Load EUR sample list
    logger.info("Loading EUR sample list...")
    with hl.hadoop_open(ANCESTRY_PATH, 'r') as f:
        ancestry_df = pd.read_csv(f, sep="\t")
    eur_ids = set(ancestry_df[ancestry_df['ancestry_pred'] == 'eur']['research_id'].astype(str))
    logger.info(f"Found {len(eur_ids)} EUR samples")

    # 2. Read pre-processed MatrixTable (MUCH faster than importing VCFs)
    logger.info(f"Reading pre-processed MatrixTable from: {SPLIT_MT_PATH}")
    mt = hl.read_matrix_table(SPLIT_MT_PATH)
    logger.info(f"Loaded MT - Rows: {mt.count_rows()}, Cols: {mt.count_cols()}")

    # 3. Filter to EUR samples
    logger.info(f"Filtering to {len(eur_ids)} EUR samples...")
    mt = mt.filter_cols(hl.literal(eur_ids).contains(mt.s))
    logger.info(f"After sample filter - Cols: {mt.count_cols()}")
    
    # 4. Optional: Filter to common variants (uncomment if desired)
    # This can dramatically reduce export time
    # logger.info("Filtering to common variants (MAF > 1%, call rate > 95%)...")
    # mt = hl.variant_qc(mt)
    # mt = mt.filter_rows(
    #     (mt.variant_qc.AF[1] > 0.01) &           # MAF > 1%
    #     (mt.variant_qc.call_rate > 0.95)         # Call rate > 95%
    # )
    # logger.info(f"After variant filter - Rows: {mt.count_rows()}")
    
    # 5. Optional: Sample variants randomly (uncomment if you want a subset)
    variant_sample_fraction = params.get('variant_sample_fraction', None)
    if variant_sample_fraction:
        logger.info(f"Randomly sampling {variant_sample_fraction*100}% of variants...")
        mt = mt.sample_rows(variant_sample_fraction, seed=42)
        logger.info(f"After sampling - Rows: {mt.count_rows()}")
    
    # 6. Repartitioning for efficient processing
    logger.info("Coalescing to 100 partitions...")
    mt = mt.naive_coalesce(100)
    
    logger.info(f"Filtering complete. Rows: {mt.count_rows()}, Cols: {mt.count_cols()}")
    return mt

def main():
    config = load_config("config/config.json")
    
    # OPTIMIZATION: Increase memory allocation for large datasets
    init_hail("data_prep", driver_mem="32g", executor_mem="16g", reference="GRCh38")

    params = config['params']
    
    # Get timestamp from environment or generate new one
    timestamp = os.environ.get('EXPORT_TIMESTAMP', datetime.now().strftime("%Y%m%d_%H%M%S"))
    
    # Determine Output Directory - Always use workspace bucket
    workspace_bucket = os.environ.get('WORKSPACE_BUCKET')
    
    if not workspace_bucket:
        logger.error("WORKSPACE_BUCKET environment variable not set!")
        logger.error("This script requires running in an AoU workspace with bucket access.")
        sys.exit(1)
    
    # Handle case where workspace_bucket already includes gs:// prefix
    if workspace_bucket.startswith('gs://'):
        base_data_dir = workspace_bucket
    else:
        base_data_dir = f"gs://{workspace_bucket}"
    
    logger.info(f"Detected AoU Workspace Bucket: {base_data_dir}")

    if params.get('dated_exports', True):
        dated_output_dir = f"{base_data_dir}/{timestamp}"
    else:
        dated_output_dir = base_data_dir

    logger.info(f"Target Output Directory: {dated_output_dir}")
    
    # Define checkpoint path for filtered MatrixTable
    checkpoint_path = f"{dated_output_dir}/filtered_mt.mt"
    
    # Check if checkpoint exists from previous run
    if hl.hadoop_exists(checkpoint_path):
        logger.info(f"Found existing checkpoint: {checkpoint_path}")
        logger.info("Loading filtered MatrixTable from checkpoint...")
        mt = hl.read_matrix_table(checkpoint_path)
        logger.info(f"Loaded MT - Rows: {mt.count_rows()}, Cols: {mt.count_cols()}")
    else:
        # Run filtering pipeline
        logger.info("No checkpoint found. Running filtering pipeline...")
        mt = perform_filtering(config, params)
        
        # Save checkpoint before attempting export
        logger.info(f"Saving checkpoint to: {checkpoint_path}")
        logger.info("This allows resuming from this point if export fails.")
        mt = mt.checkpoint(checkpoint_path, overwrite=True)
        logger.info("Checkpoint saved successfully.")

    # Export with optimizations
    logger.info("Proceeding to PLINK export...")
    perform_plink_export(mt, dated_output_dir, config, timestamp)
    
    logger.info("Job completed successfully.")

if __name__ == "__main__":
    main()