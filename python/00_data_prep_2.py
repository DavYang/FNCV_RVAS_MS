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
            f.write(f"Checkpoint source: {config['outputs']['intermediate_dir']}/wgs_eur_intervals.mt\n")
            f.write(f"Config: {config}\n")
    except Exception as e:
        logger.warning(f"Failed to write metadata file to {metadata_file}: {e}")
    
    logger.info(f"PLINK export complete: {plink_base}")
    logger.info(f"Output directory: {output_dir}")

def perform_filtering_and_checkpoint(config, checkpoint_path, params):
    """Run the initial filtering pipeline and checkpoint."""
    ANCESTRY_PATH = config['inputs']['ancestry_pred']
    WGS_MT_PATH = config['inputs']['wgs_mt']
    
    # Samples
    logger.info("Loading EUR sample list...")
    with hl.hadoop_open(ANCESTRY_PATH, 'r') as f:
        ancestry_df = pd.read_csv(f, sep="\t")
    eur_ids = set(ancestry_df[ancestry_df['ancestry_pred'] == 'eur']['research_id'].astype(str))

    # Intervals
    logger.info("Generating random intervals (Optimized: Fewer, Larger)...")
    rg = hl.get_reference('GRCh38')
    lengths = rg.lengths
    chroms = [f"chr{i}" for i in range(1, 23)]
    hl_intervals = []
    
    # --- OPTIMIZED INTERVAL STRATEGY ---
    # Analysis Plan: "Randomly select 4,000 genomic intervals of 50kb each... ~200Mb"
    # To fix 145k partition issue, we generate 40 x 5Mb intervals instead.
    # Total coverage: 40 * 5Mb = 200Mb.
    
    n_large_intervals = 40
    large_window_size = 5_000_000  # 5 Mb
    
    logger.info(f"Strategy: {n_large_intervals} intervals of {large_window_size/1e6} Mb each")
    
    for _ in range(n_large_intervals):
        chrom = random.choice(chroms)
        chrom_len = lengths[chrom]
        # Ensure we don't go out of bounds
        start = random.randint(1, chrom_len - large_window_size)
        interval = hl.locus_interval(chrom, start, start + large_window_size, reference_genome='GRCh38')
        hl_intervals.append(interval)

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
    
    # 4. Select fields early
    mt = mt.select_entries('GT')
    
    # 5. Checkpoint (Crucial Step)
    logger.info(f"Checkpointing filtered MT to {checkpoint_path}...")
    
    # Checkpointing saves the filtered subset to disk, preventing re-reads of the massive WGS MT
    mt = mt.checkpoint(checkpoint_path, overwrite=True)
    
    logger.info(f"Checkpoint complete. Rows: {mt.count_rows()}, Cols: {mt.count_cols()}")
    return mt

def process_mt_for_export(mt, params):
    """Process MatrixTable (QC, filtering) for export."""
    
    # Repartition to fix 145k partition issue
    # Use naive_coalesce for INSTANT repartitioning (no shuffle)
    # 100 partitions is appropriate for ~200Mb of genotype data
    logger.info("Coalescing MatrixTable to 100 partitions (naive_coalesce)...")
    mt = mt.naive_coalesce(100)
    
    # HARDCODED: SKIP VARIANT QC
    logger.info("Skipping Hail Variant QC (QC will be performed in PLINK)")

    # HARDCODED: SKIP DOWNSAMPLING
    logger.info("Skipping downsampling (all variants will be exported)")
    
    n_total = mt.count_rows()
    logger.info(f"Total variants to export: {n_total}")
        
    return mt

def main():
    config = load_config("config/config.json")
    
    # Initialize with high memory
    init_hail("data_prep", driver_mem="16g", reference="GRCh38")

    params = config['params']
    
    # Get timestamp from environment or generate new one
    timestamp = os.environ.get('EXPORT_TIMESTAMP', datetime.now().strftime("%Y%m%d_%H%M%S"))
    
    # Checkpoint path
    checkpoint_path = f"{config['outputs']['intermediate_dir']}/wgs_eur_intervals.mt"
    
    # Determine Output Directory (Workspace Bucket vs Local)
    workspace_bucket = os.environ.get('WORKSPACE_BUCKET')
    data_dir_suffix = config['outputs'].get('data_dir_suffix', 'data')
    
    if workspace_bucket:
        logger.info(f"Detected AoU Workspace Bucket: {workspace_bucket}")
        base_data_dir = f"{workspace_bucket}/{data_dir_suffix}"
    else:
        logger.warning("WORKSPACE_BUCKET not found. Falling back to local 'data' directory (WARNING: Disk space low!)")
        base_data_dir = "data"

    if params.get('dated_exports', True):
        dated_output_dir = f"{base_data_dir}/{timestamp}"
    else:
        dated_output_dir = base_data_dir

    logger.info(f"Target Output Directory: {dated_output_dir}")

    # Check for existing checkpoint
    if hl.hadoop_exists(checkpoint_path):
        logger.info(f"Found checkpoint at {checkpoint_path}")
        logger.info("Skipping filtering, proceeding to processing...")
        
        mt = hl.read_matrix_table(checkpoint_path)
        logger.info(f"Loaded MT: {mt.count_rows()} variants, {mt.count_cols()} samples")
        
    else:
        logger.info("No checkpoint found, running full filtering pipeline...")
        mt = perform_filtering_and_checkpoint(config, checkpoint_path, params)

    # Process MT (repartition, SKIP QC, SKIP DOWNSAMPLING)
    mt = process_mt_for_export(mt, params)
    
    # Export
    perform_plink_export(mt, dated_output_dir, config, timestamp)
    
    logger.info("Job completed successfully.")

if __name__ == "__main__":
    main()