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
    """Export MatrixTable to PLINK with dated directory."""
    
    # Note: os.makedirs is not needed for gs:// paths and won't work
    if not output_dir.startswith("gs://"):
        os.makedirs(output_dir, exist_ok=True)
    
    # Symlinks don't work on GS, so we skip that logic if using GS
    if not output_dir.startswith("gs://"):
        # Update latest symlink (Local only)
        latest_link = f"{os.path.dirname(output_dir)}/latest"
        
        if os.path.lexists(latest_link):
            try:
                if os.path.islink(latest_link):
                    os.remove(latest_link)
                else:
                    logger.warning(f"{latest_link} exists but is not a symlink. Skipping update.")
            except Exception as e:
                logger.warning(f"Could not remove existing symlink {latest_link}: {e}")
                
        try:
            os.symlink(os.path.basename(output_dir), latest_link)
        except Exception as e:
            logger.warning(f"Could not create symlink {latest_link}: {e}")
    else:
        logger.info(f"Skipping symlink creation for GS path: {output_dir}")
    
    # Export base path
    plink_base = f"{output_dir}/eur_common_snps_500k"
    
    logger.info(f"Exporting to PLINK: {plink_base}")
    
    hl.export_plink(
        mt,
        plink_base,
        ind_id=mt.s,
        fam_id=mt.s
    )
    
    # Write export metadata (Supports both local and GS)
    metadata_file = f"{output_dir}/export_metadata.txt"
    try:
        with hl.hadoop_open(metadata_file, 'w') as f:
            f.write(f"Export timestamp: {timestamp}\n")
            f.write(f"Variants: {mt.count_rows()}\n")
            f.write(f"Samples: {mt.count_cols()}\n")
            f.write(f"Source: ACAF common variants MatrixTable\n")
            f.write(f"Config: {config}\n")
    except Exception as e:
        logger.warning(f"Failed to write metadata file to {metadata_file}: {e}")
    
    logger.info(f"PLINK export complete: {plink_base}")
    logger.info(f"Output directory: {output_dir}")

def perform_filtering(config, params):
    """Run the initial filtering pipeline."""
    ANCESTRY_PATH = config['inputs']['ancestry_pred']
    WGS_MT_PATH = config['inputs']['wgs_mt']
    
    # Samples
    logger.info("Loading EUR sample list...")
    with hl.hadoop_open(ANCESTRY_PATH, 'r') as f:
        ancestry_df = pd.read_csv(f, sep="\t")
    eur_ids = set(ancestry_df[ancestry_df['ancestry_pred'] == 'eur']['research_id'].astype(str))

    # Intervals - Using contiguous blocks for faster processing
    logger.info("Generating contiguous genomic blocks (Optimized for speed)...")
    rg = hl.get_reference('GRCh38')
    lengths = rg.lengths
    chroms = [f"chr{i}" for i in range(1, 23)]
    hl_intervals = []
    
    # --- CONTIGUOUS BLOCK STRATEGY (FASTER) ---
    # Divide genome into systematic contiguous blocks for better partitioning
    # Total coverage: ~200Mb in contiguous segments
    
    n_blocks = 40
    large_window_size = 5_000_000  # 5 Mb
    
    logger.info(f"Strategy: {n_blocks} contiguous blocks of {large_window_size/1e6} Mb each")
    
    # Calculate total autosomal genome size
    total_size = sum(lengths[chrom] for chrom in chroms)
    block_size = total_size // n_blocks
    
    # Create contiguous blocks across chromosomes
    for chrom in chroms:
        if len(hl_intervals) >= n_blocks:
            break
            
        chrom_len = lengths[chrom]
        current_pos = 1
        
        while current_pos < chrom_len and len(hl_intervals) < n_blocks:
            start = current_pos
            end = min(current_pos + large_window_size, chrom_len)
            
            # Create interval
            interval = hl.locus_interval(chrom, start, end, reference_genome='GRCh38')
            hl_intervals.append(interval)
            
            current_pos = end + 1

    # 1. Read WGS MatrixTable
    logger.info("Reading WGS MatrixTable...")
    wgs_mt = hl.read_matrix_table(WGS_MT_PATH)
    
    # 2. Filter to Intervals (Optimized)
    logger.info(f"Filtering MatrixTable to {len(hl_intervals)} large intervals using filter_intervals...")
    # Use filter_intervals which is much faster than filter_rows with lookup
    mt = hl.filter_intervals(wgs_mt, hl_intervals)
    
    # 3. Filter Samples
    logger.info(f"Filtering to {len(eur_ids)} EUR samples...")
    mt = mt.filter_cols(hl.literal(eur_ids).contains(mt.s))
    
    # 4. Keep all entry fields (GT, AD, DP, GQ, etc.)
    # No select_entries() call - keeping all fields for complete data
    
    # 5. Optimize partitioning
    # Fix the "9000 partitions" issue here. Since we filtered down to ~200Mb of genome,
    # we don't need the thousands of partitions from the original WGS MT.
    logger.info("Coalescing filtered data to 100 partitions...")
    mt = mt.naive_coalesce(100)
    
    logger.info(f"Filtering complete. Rows: {mt.count_rows()}, Cols: {mt.count_cols()}")
    return mt

def process_mt_for_export(mt, params):
    """Process MatrixTable for export."""
    
    # HARDCODED: SKIP VARIANT QC
    logger.info("Skipping Hail Variant QC (QC will be performed in PLINK)")

    # NO DOWNSAMPLING - Export all variants
    n_total = mt.count_rows()
    logger.info(f"Total variants to export: {n_total}")
    logger.info("Keeping all variants (no downsampling)")
        
    return mt

def main():
    config = load_config("config/config.json")
    
    # Initialize with high memory
    init_hail("data_prep", driver_mem="16g", reference="GRCh38")

    params = config['params']
    
    # Get timestamp from environment or generate new one
    timestamp = os.environ.get('EXPORT_TIMESTAMP', datetime.now().strftime("%Y%m%d_%H%M%S"))
    
    # Determine Output Directory (Workspace Bucket vs Local)
    workspace_bucket = os.environ.get('WORKSPACE_BUCKET')
    data_dir_suffix = config['outputs'].get('data_dir_suffix', 'data')
    
    if workspace_bucket:
        # Handle case where workspace_bucket already includes gs:// prefix
        if workspace_bucket.startswith('gs://'):
            base_data_dir = f"{workspace_bucket}/{data_dir_suffix}"
        else:
            base_data_dir = f"gs://{workspace_bucket}/{data_dir_suffix}"
        logger.info(f"Detected AoU Workspace Bucket: {workspace_bucket}")
    else:
        logger.warning("WORKSPACE_BUCKET not found. Falling back to local 'data' directory (WARNING: Disk space low!)")
        base_data_dir = "data"

    if params.get('dated_exports', True):
        dated_output_dir = f"{base_data_dir}/{timestamp}"
    else:
        dated_output_dir = base_data_dir

    logger.info(f"Target Output Directory: {dated_output_dir}")

    # Always run full filtering pipeline (no checkpoint)
    logger.info("Running full filtering pipeline...")
    mt = perform_filtering(config, params)

    # Process MT (SKIP QC, NO DOWNSAMPLING)
    mt = process_mt_for_export(mt, params)
    
    # Export
    perform_plink_export(mt, dated_output_dir, config, timestamp)
    
    logger.info("Job completed successfully.")

if __name__ == "__main__":
    main()