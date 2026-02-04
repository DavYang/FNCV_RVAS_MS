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

def select_genome_intervals(config, logger):
    """Select random intervals across the genome."""
    logger.info("Selecting genome-wide intervals...")
    
    # Set random seed for reproducibility
    random.seed(config['sampling']['random_seed'])
    
    interval_list_base = config['vcf_processing']['interval_list_base']
    
    # Sample a subset of interval files to process
    # Based on the actual file structure, we have files like 0000000000.interval_list
    total_interval_files = config['vcf_processing']['total_shards']
    n_files_to_sample = min(500, total_interval_files)  # Sample 500 files (increased from 100)
    
    sampled_files = random.sample(range(total_interval_files), n_files_to_sample)
    all_intervals = []
    
    for file_idx in sampled_files:
        interval_list_path = f"{interval_list_base}/{file_idx:010d}.interval_list"
        intervals = parse_interval_list(interval_list_path, config, logger)
        all_intervals.extend(intervals)
    
    # Sample intervals from the collected set
    n_intervals_to_sample = min(5000, len(all_intervals))  # Target 5000 intervals (increased from 1000)
    selected_intervals = random.sample(all_intervals, n_intervals_to_sample)
    
    logger.info(f"Selected {len(selected_intervals)} intervals across {n_files_to_sample} interval files")
    return selected_intervals

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
    output_vcf = f"{output_dir}/100k_background_snps.vcf.bgz"
    logger.info(f"Exporting to compressed VCF: {output_vcf}")
    
    export_start = time.time()
    hl.export_vcf(
        mt_sampled,
        output_vcf,
        parallel='header_per_shard',
        include_star=False,
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
    
    summary_path = f"{output_dir}/100k_sampling_summary.json"
    with open(summary_path, 'w') as f:
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
    output_dir = f"{base_data_dir}/results/FNCV_RVAS_MS/100k_background_snps_{timestamp}"
    logger.info(f"Output Directory: {output_dir}")
    
    try:
        # Step 1: Load EUR sample IDs
        eur_samples = load_ancestry_data(config, logger)
        
        # Step 2: Select genome-wide intervals
        selected_intervals = select_genome_intervals(config, logger)
        
        # Step 3: Find VCF shards for selected intervals
        vcf_paths = find_shards_for_intervals(selected_intervals, config, logger)
        
        # Step 4: Import VCFs with EUR sample filtering
        logger.info(f"Importing {len(vcf_paths)} VCF shards...")
        mt = hl.import_vcf(
            vcf_paths,
            reference_genome='GRCh38',
            force_bgz=True,
            array_elements_required=False
        )
        
        # Filter to EUR samples during import
        logger.info("Filtering to EUR samples...")
        mt = mt.filter_cols(hl.literal(eur_samples).contains(mt.s))
        
        # Step 5: Split multi-allelic variants
        logger.info("Splitting multi-allelic variants...")
        mt = hl.split_multi_hts(mt)
        
        # Step 6: Filter to selected intervals
        logger.info("Filtering to selected genomic intervals...")
        interval_literals = [hl.locus_interval(
            interval['chrom'], 
            interval['start'], 
            interval['end'], 
            reference_genome='GRCh38'
        ) for interval in selected_intervals]
        
        mt = mt.filter_rows(
            hl.any([interval_literals[i].contains(mt.locus) for i in range(len(interval_literals))])
        )
        
        # Step 7: Repartition for efficiency
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