#!/usr/bin/env python3
import os
import sys
import random
import hail as hl
import pandas as pd
import json
import time
from datetime import datetime
from utils import load_config, init_hail, setup_logger

logger = setup_logger("100k_snp_sampling")

def load_ancestry_data(config, logger):
    """Load ancestry predictions and filter to European ancestry samples."""
    logger.info("Loading ancestry data...")
    
    # Load ancestry TSV as Hail Table
    ancestry_ht = hl.import_table(
        config['inputs']['ancestry_pred'],
        impute=True,
        types={'research_id': hl.tstr, 'ancestry_pred': hl.tstr}
    )
    
    # Filter to European ancestry samples
    eur_samples = ancestry_ht.filter(ancestry_ht.ancestry_pred == 'eur')
    eur_sample_ids = eur_samples.aggregate(hl.agg.collect(eur_samples.research_id))
    
    logger.info(f"Found {len(eur_sample_ids)} European ancestry samples")
    return eur_sample_ids

def parse_interval_list(interval_list_path, config, logger):
    """Parse interval list file and extract genomic intervals."""
    logger.info(f"Parsing interval list: {interval_list_path}")
    
    intervals = []
    hla_region = config['sampling']['hla_region']
    
    try:
        with hl.hadoop_open(interval_list_path, 'r') as f:
            lines = f.readlines()
        
        for line in lines:
            line = line.strip()
            if line.startswith('@') or not line:
                continue
            
            # Parse interval format: chr:start-end
            parts = line.split('\t')
            if len(parts) >= 3:
                chrom = parts[0]  # Keep as-is (should include 'chr' prefix)
                start = int(parts[1])
                end = int(parts[2])
                
                # Skip HLA region
                if chrom == f"chr{hla_region['chrom']}" and start >= hla_region['start'] and end <= hla_region['end']:
                    continue
                
                intervals.append({
                    'chrom': chrom,
                    'start': start,
                    'end': end
                })
        
        logger.info(f"Parsed {len(intervals)} intervals from {interval_list_path}")
        return intervals
        
    except Exception as e:
        logger.error(f"Failed to parse interval list {interval_list_path}: {e}")
        return []

def sample_intervals_iteratively(config, logger):
    """Sample intervals iteratively until target variant count is reached."""
    target_variants = config['sampling']['target_total_snps']
    intervals_per_iteration = 1000  # Start with 1000 intervals
    max_iterations = 10
    
    # Set random seed for reproducibility
    random.seed(config['sampling']['random_seed'])
    
    interval_list_base = config['vcf_processing']['interval_list_base']
    total_interval_files = config['vcf_processing']['total_shards']
    
    # Track all collected intervals and VCF paths
    all_intervals = []
    all_vcf_paths = []
    processed_files = set()  # Avoid reprocessing same files
    
    iteration = 0
    current_variants = 0
    
    while current_variants < target_variants and iteration < max_iterations:
        iteration += 1
        logger.info(f"\n=== Interval Sampling Iteration {iteration}/{max_iterations} ===")
        logger.info(f"Current variants: {current_variants:,}, Target: {target_variants:,}")
        
        # Sample new intervals
        new_intervals = _sample_new_intervals(
            interval_list_base, total_interval_files, intervals_per_iteration, 
            processed_files, config, logger
        )
        
        if not new_intervals:
            logger.warning("No more intervals available to sample")
            break
            
        all_intervals.extend(new_intervals)
        logger.info(f"Added {len(new_intervals)} intervals (total: {len(all_intervals)})")
        
        # Find VCF shards for new intervals
        new_vcf_paths = find_shards_for_intervals(new_intervals, config, logger)
        all_vcf_paths.extend(new_vcf_paths)
        logger.info(f"Added {len(new_vcf_paths)} VCF shards (total: {len(all_vcf_paths)})")
        
        # Import all VCFs and count variants
        logger.info("Importing and filtering VCFs...")
        mt = _import_and_filter_vcfs(all_vcf_paths, config, logger)
        current_variants = mt.count_rows()
        logger.info(f"Total variants after filtering: {current_variants:,}")
        
        # Adjust next iteration size based on variant density
        if iteration == 1:
            variants_per_interval = current_variants / len(all_intervals)
            estimated_needed = (target_variants - current_variants) / variants_per_interval
            intervals_per_iteration = min(5000, max(500, int(estimated_needed * 1.2)))  # 20% buffer
            logger.info(f"Adjusted intervals per iteration to {intervals_per_iteration} based on variant density")
        
        if current_variants >= target_variants:
            logger.info(f"Target reached! {current_variants:,} variants collected")
            break
    
    if current_variants < target_variants:
        logger.warning(f"Could not reach target after {max_iterations} iterations. Final: {current_variants:,} variants")
    
    return all_vcf_paths, mt

def _sample_new_intervals(interval_list_base, total_files, target_count, processed_files, config, logger):
    """Sample new intervals from unprocessed files."""
    # Get unprocessed files
    available_files = [f for f in range(total_files) if f not in processed_files]
    random.shuffle(available_files)
    
    new_intervals = []
    
    for file_idx in available_files:
        if len(new_intervals) >= target_count:
            break
            
        interval_list_path = f"{interval_list_base}/{file_idx:010d}.interval_list"
        intervals = parse_interval_list(interval_list_path, config, logger)
        
        # Take intervals from this file
        remaining_needed = target_count - len(new_intervals)
        if len(intervals) <= remaining_needed:
            new_intervals.extend(intervals)
        else:
            sampled_from_file = random.sample(intervals, remaining_needed)
            new_intervals.extend(sampled_from_file)
        
        processed_files.add(file_idx)
    
    return new_intervals

def _import_and_filter_vcfs(vcf_paths, config, logger):
    """Import VCFs and apply standard filters."""
    # Import VCFs
    mt = hl.import_vcf(
        vcf_paths,
        reference_genome='GRCh38',
        force_bgz=True,
        array_elements_required=False
    )
    
    # Load EUR samples (cached if possible)
    if not hasattr(_import_and_filter_vcfs, 'eur_samples'):
        _import_and_filter_vcfs.eur_samples = load_ancestry_data(config, logger)
    
    # Filter to EUR samples
    mt = mt.filter_cols(hl.literal(_import_and_filter_vcfs.eur_samples).contains(mt.s))
    
    # Split multi-allelic variants
    mt = hl.split_multi_hts(mt)
    
    return mt

def find_shards_for_intervals(intervals, config, logger):
    """Find VCF shards that contain the selected intervals."""
    logger.info("Finding VCF shards for selected intervals...")
    
    sharded_vcf_base = config['vcf_processing']['sharded_vcf_base']
    total_shards = config['vcf_processing']['total_shards']
    
    # Since we're using the actual VCF files (not sharded by chromosome),
    # we'll sample shards directly from the total pool
    # Each shard corresponds to a VCF file like 0000000000.vcf.bgz
    
    # Calculate how many shards we need based on interval distribution
    # Aim for ~50-100 shards to get good coverage while keeping data manageable
    n_shards_to_sample = min(100, max(50, len(intervals) // 10))
    
    # Set random seed for reproducibility
    random.seed(config['sampling']['random_seed'])
    
    # Sample shard indices from the total pool
    sampled_shard_indices = random.sample(range(total_shards), n_shards_to_sample)
    
    # Generate VCF paths
    vcf_paths = []
    for shard_idx in sampled_shard_indices:
        vcf_path = f"{sharded_vcf_base}/{shard_idx:010d}.vcf.bgz"
        vcf_paths.append(vcf_path)
    
    logger.info(f"Selected {len(vcf_paths)} VCF shards from {total_shards} total shards")
    logger.info(f"Shard range: {min(sampled_shard_indices):010d} to {max(sampled_shard_indices):010d}")
    
    return vcf_paths


def export_100k_snps(mt, output_dir, config, logger):
    """Sample exactly 100K SNPs with iterative approach and export as VCF."""
    target_snps = config['sampling']['target_total_snps']
    random_seed = config['sampling']['random_seed']
    max_iterations = 10  # Prevent infinite loops
    
    logger.info(f"Target: {target_snps:,} SNPs")
    
    # Count total variants first
    logger.info("Counting total variants in filtered data...")
    total_variants = mt.count_rows()
    logger.info(f"Total variants available: {total_variants:,}")
    
    if total_variants <= target_snps:
        logger.warning(f"Only {total_variants:,} variants available, less than target {target_snps:,}")
        mt_sampled = mt
        estimated_sampled = total_variants
        iteration = 1
    else:
        # Iterative sampling approach
        mt_sampled = None
        estimated_sampled = 0
        iteration = 0
        
        while estimated_sampled < target_snps and iteration < max_iterations:
            iteration += 1
            logger.info(f"Sampling iteration {iteration}/{max_iterations}")
            
            # Calculate sampling fraction for this iteration
            remaining_needed = target_snps - estimated_sampled
            remaining_variants = total_variants - estimated_sampled
            
            if remaining_variants <= 0:
                break
                
            # Sample a bit more than needed to account for sampling variance
            sampling_fraction = min(1.0, (remaining_needed * 1.1) / remaining_variants)
            logger.info(f"Sampling fraction: {sampling_fraction:.6f}")
            
            # Sample variants
            logger.info("Performing variant sampling...")
            mt_current = mt.sample_rows(sampling_fraction, seed=random_seed + iteration)
            
            # If first iteration, use as base
            if mt_sampled is None:
                mt_sampled = mt_current
            else:
                # Union with previously sampled variants
                mt_sampled = mt_sampled.union_rows(mt_current)
            
            # Update counts
            estimated_sampled = mt_sampled.count_rows()
            logger.info(f"Total sampled variants so far: {estimated_sampled:,}")
            
            if estimated_sampled >= target_snps:
                logger.info(f"Target reached! {estimated_sampled:,} variants sampled")
                break
        
        if estimated_sampled < target_snps:
            logger.warning(f"Could not reach target after {max_iterations} iterations. Final count: {estimated_sampled:,}")
    
    # Export as compressed VCF
    target_snps_k = target_snps // 1000  # Convert to thousands for filename
    output_vcf = f"{output_dir}/{target_snps_k}k_background_snps.vcf.bgz"
    logger.info(f"Exporting to compressed VCF: {output_vcf}")
    
    export_start = time.time()
    hl.export_vcf(
        mt_sampled,
        output_vcf,
        parallel='header_per_shard',
        tabix=True
    )
    export_time = time.time() - export_start
    
    logger.info(f"VCF export completed in {export_time:.1f} seconds")
    
    # Create summary
    summary = {
        'target_snps': target_snps,
        'total_variants_before_sampling': total_variants,
        'final_sampled_variants': estimated_sampled,
        'iterations_used': iteration,
        'export_path': output_vcf,
        'export_time_seconds': export_time,
        'timestamp': datetime.now().isoformat()
    }
    
    summary_path = f"{output_dir}/{target_snps_k}k_sampling_summary.json"
    with hl.hadoop_open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)
    
    logger.info(f"Sampling summary saved to: {summary_path}")
    return summary


def main():
    config = load_config("config/config.json")
    
    # Initialize Hail with optimized settings
    init_hail("100k_snp_sampling", driver_mem="32g", reference="GRCh38")
    
    # Get timestamp for output directory
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Determine Output Directory
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
    
    # Create output directory
    target_snps_k = config['sampling']['target_total_snps'] // 1000
    output_dir = f"{base_data_dir}/results/FNCV_RVAS_MS/{target_snps_k}k_background_snps_{timestamp}"
    logger.info(f"Output Directory: {output_dir}")
    
    try:
        # Step 1: Use variant-driven interval sampling
        vcf_paths, mt = sample_intervals_iteratively(config, logger)
        
        # Step 2: Repartition for efficient processing
        logger.info("Repartitioning to 200 partitions...")
        mt = mt.naive_coalesce(200)
        
        logger.info(f"Filtered data - Rows: {mt.count_rows():,}, Cols: {mt.count_cols():,}")
        
        # Step 8: Sample exactly 100K SNPs and export
        summary = export_100k_snps(mt, output_dir, config, logger)
        
        logger.info("100K SNP sampling completed successfully!")
        logger.info(f"Results saved to: {output_dir}")
        logger.info(f"VCF file: {summary['export_path']}")
        
    except Exception as e:
        logger.error(f"Pipeline failed: {e}")
        raise

if __name__ == "__main__":
    main()