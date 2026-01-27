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
    plink_base = f"{output_dir}/eur_common_snps_sampled"
    
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
            f.write(f"Source: ACAF Threshold VCFs (Sampled)\n")
            f.write(f"Config: {config}\n")
    except Exception as e:
        logger.warning(f"Failed to write metadata file to {metadata_file}: {e}")
    
    logger.info(f"PLINK export complete: {plink_base}")
    logger.info(f"Output directory: {output_dir}")

def perform_filtering(config, params):
    """Run the VCF sampling and filtering pipeline."""
    ANCESTRY_PATH = config['inputs']['ancestry_pred']
    VCF_BASE = config['inputs']['wgs_vcf_base']
    
    # 1. Samples
    logger.info("Loading EUR sample list...")
    with hl.hadoop_open(ANCESTRY_PATH, 'r') as f:
        ancestry_df = pd.read_csv(f, sep="\t")
    eur_ids = set(ancestry_df[ancestry_df['ancestry_pred'] == 'eur']['research_id'].astype(str))

    # 2. Select Random VCF Shards
    total_shards = params.get('total_shards', 24744)
    n_sample = params.get('n_shards_to_sample', 500)
    seed = params.get('random_seed', 42)
    
    logger.info(f"Sampling {n_sample} random VCF shards from pool of {total_shards} using seed {seed}...")
    
    # Generate random shard indices
    # We use sample to get unique indices
    random.seed(seed)
    random_indices = random.sample(range(total_shards), n_sample)
    
    # Construct paths: 0000000000.vcf.bgz format
    vcf_paths = [f"{VCF_BASE}/{i:010d}.vcf.bgz" for i in random_indices]
    
    # Log first few for verification
    if len(vcf_paths) > 0:
        logger.info(f"Example VCFs to load: {vcf_paths[:3]}")

    # 3. Import VCFs
    logger.info("Importing VCFs...")
    # force_bgz is required for AoU files; reference_genome ensures proper contig handling
    # array_elements_required=False handles fields that might vary in length
    mt = hl.import_vcf(
        vcf_paths, 
        reference_genome='GRCh38', 
        force_bgz=True,
        array_elements_required=False
    )
    
    # 4. Split Multi-allelic Variants
    # The 'splitMT' you used before had this done already. For raw VCFs, we must do it manually
    # to ensure compatibility with PLINK and downstream tools.
    logger.info("Splitting multi-allelic variants (split_multi_hts)...")
    mt = hl.split_multi_hts(mt)

    # 5. Filter Samples
    logger.info(f"Filtering to {len(eur_ids)} EUR samples...")
    mt = mt.filter_cols(hl.literal(eur_ids).contains(mt.s))
    
    # 6. Repartitioning
    # Importing many small files might result in awkward partitioning. 
    # Coalesce to a reasonable number (e.g., 100) for downstream processing.
    logger.info("Coalescing to 100 partitions...")
    mt = mt.naive_coalesce(100)
    
    logger.info(f"Filtering complete. Rows: {mt.count_rows()}, Cols: {mt.count_cols()}")
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
    
    # We want the 'data_dir_suffix' (e.g., FNCV_RVAS_MS)
    # utils.load_config may have already injected the bucket prefix into config['outputs'] if the env var was set.
    # To be safe and avoid duplication, we check if the path already starts with gs://
    
    raw_suffix = config['outputs'].get('data_dir_suffix', 'data')
    
    if workspace_bucket:
        logger.info(f"Detected AoU Workspace Bucket: {workspace_bucket}")
        
        if raw_suffix.startswith("gs://"):
            # It's already a full path (likely modified by utils.load_config)
            base_data_dir = raw_suffix
        else:
            # It's just a suffix, need to prepend bucket
            # Ensure we don't double up if workspace_bucket ends with / or suffix starts with /
            # But simpler here: assume standard format
            if workspace_bucket.startswith("gs://"):
                 base_data_dir = f"{workspace_bucket}/{raw_suffix}"
            else:
                 base_data_dir = f"gs://{workspace_bucket}/{raw_suffix}"
    else:
        logger.warning("WORKSPACE_BUCKET not found. Falling back to local 'data' directory (WARNING: Disk space low!)")
        # If no workspace bucket, we assume we are local. 
        # If load_config didn't inject anything, raw_suffix is just "FNCV_RVAS_MS" or "data"
        base_data_dir = raw_suffix

    if params.get('dated_exports', True):
        dated_output_dir = f"{base_data_dir}/background_EUR_common_snps_sampled_{timestamp}"
    else:
        dated_output_dir = base_data_dir

    logger.info(f"Target Output Directory: {dated_output_dir}")

    # Run filtering pipeline
    logger.info("Running VCF sampling pipeline...")
    mt = perform_filtering(config, params)

    # Export
    perform_plink_export(mt, dated_output_dir, config, timestamp)
    
    logger.info("Job completed successfully.")

if __name__ == "__main__":
    main()
