import os
import sys
import random
import hail as hl
import pandas as pd
import json
import time
from datetime import datetime
from utils import load_config, init_hail, setup_logger

logger = setup_logger("SNP_sampling")

# Global cache for EUR samples
eur_samples_cache = None

def get_eur_samples(config, logger):
    """Load EUR samples with caching to avoid repeated loading."""
    global eur_samples_cache
    if eur_samples_cache is None:
        logger.info("Loading EUR ancestry samples (first time only)...")
        eur_samples_cache = load_ancestry_data(config, logger)
        logger.info(f"Cached {len(eur_samples_cache)} EUR samples for reuse")
    return eur_samples_cache

def calculate_chromosome_targets(target_snps, logger):
    """Calculate proportional SNP targets per chromosome based on GRCh38 sizes."""
    # GRCh38 approximate chromosome sizes (in base pairs)
    chr_sizes = {
        'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559,
        'chr4': 190214555, 'chr5': 181538259, 'chr6': 170805979,
        'chr7': 159345973, 'chr8': 145138636, 'chr9': 138394717,
        'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
        'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189,
        'chr16': 90338345, 'chr17': 83257441, 'chr18': 80373285,
        'chr19': 58617616, 'chr20': 64444167, 'chr21': 46709983,
        'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415
    }
    
    total_genome_size = sum(chr_sizes.values())
    chr_targets = {}
    
    for chrom, size in chr_sizes.items():
        proportion = size / total_genome_size
        target = int(target_snps * proportion)
        chr_targets[chrom] = max(1, target)  # Ensure at least 1 SNP per chromosome
    
    # Adjust for rounding to match exact target
    current_total = sum(chr_targets.values())
    if current_total != target_snps:
        # Add remaining SNPs to largest chromosomes
        remaining = target_snps - current_total
        sorted_chroms = sorted(chr_targets.items(), key=lambda item: item[1], reverse=True)
        for i, (chrom, target) in enumerate(sorted_chroms):
            if i < remaining:
                chr_targets[chrom] += 1
            else:
                break
    
    logger.info(f"Chromosome SNP targets: {chr_targets}")
    logger.info(f"Total target SNPs: {sum(chr_targets.values())}")
    
    return chr_targets

def sample_chromosome_intervals(chrom, target_snps, config, logger):
    """Sample intervals for a specific chromosome until target SNP count is reached."""
    random.seed(config['sampling']['random_seed'])
    
    interval_list_base = config['vcf_processing']['interval_list_base']
    total_interval_files = config['vcf_processing']['total_shards']
    
    logger.info(f"{chrom}: Starting interval sampling (target: {target_snps} variants)")
    logger.info(f"{chrom}: Interval list base: {interval_list_base}")
    logger.info(f"{chrom}: Total interval files: {total_interval_files}")
    
    # Get intervals for this chromosome
    all_intervals = []
    processed_files = set()
    
    # Sample intervals until we have enough variants
    current_variants = 0
    intervals_sampled = 0
    max_intervals = 100  # Safety limit
    
    while current_variants < target_snps and intervals_sampled < max_intervals:
        logger.info(f"{chrom}: Iteration {intervals_sampled + 1}/{max_intervals} - Current variants: {current_variants}/{target_snps}")
        
        # Sample one interval for this chromosome
        available_files = [f for f in range(total_interval_files) if f not in processed_files]
        random.shuffle(available_files)
        
        if not available_files:
            logger.warning(f"No more interval files available for {chrom}")
            break
        
        # Sample one file
        file_idx = available_files[0]
        interval_list_path = f"{interval_list_base}/{file_idx:010d}.interval_list"
        logger.info(f"{chrom}: Checking file {file_idx}: {interval_list_path}")
        
        # Parse intervals and stop early when we find a matching chromosome interval
        sampled_interval = parse_interval_list_and_sample(interval_list_path, chrom, config, logger)
        
        if sampled_interval:
            all_intervals.append(sampled_interval)
            intervals_sampled += 1
            logger.info(f"{chrom:} Found interval {intervals_sampled}: {sampled_interval['chrom']}:{sampled_interval['start']}-{sampled_interval['end']}")
            
            # Import VCF for this interval and count variants
            vcf_paths = find_shards_for_intervals([sampled_interval], config, logger)
            if vcf_paths:
                logger.info(f"{chrom}: Importing VCF for interval {intervals_sampled}...")
                mt = _import_and_filter_vcfs(vcf_paths, config, logger)
                interval_variants = mt.count_rows()
                current_variants += interval_variants
                logger.info(f"{chrom}: Sampled interval {intervals_sampled}, added {interval_variants} variants (total: {current_variants}/{target_snps})")
            else:
                logger.warning(f"{chrom}: No VCF paths found for interval {intervals_sampled}")
        else:
            logger.info(f"{chrom}: No {chrom} intervals found in file {file_idx}")
        
        processed_files.add(file_idx)
    
    logger.info(f"{chrom}: Sampled {len(all_intervals)} intervals, collected {current_variants} variants")
    return all_intervals, current_variants

def process_chromosome(chrom, target_snps, output_dir, config, logger, base_dir):
    """Process a single chromosome: sample intervals, filter, sample SNPs, and export."""
    logger.info(f"\n=== Processing {chrom} (target: {target_snps} SNPs) ===")
    
    # Step 1: Sample intervals for this chromosome
    intervals, available_variants = sample_chromosome_intervals(chrom, target_snps, config, logger)
    
    if available_variants == 0:
        logger.warning(f"No variants found for {chrom}")
        return None
    
    # Step 2: Import all VCFs for this chromosome
    logger.info(f"{chrom}: Importing VCFs for {len(intervals)} intervals...")
    vcf_paths = find_shards_for_intervals(intervals, config, logger)
    mt = _import_and_filter_vcfs(vcf_paths, config, logger)
    
    # Step 3: Sample target SNPs for this chromosome
    logger.info(f"{chrom}: Sampling {target_snps} SNPs from {mt.count_rows()} available variants...")
    mt_sampled = sample_snps_from_matrix(mt, target_snps, config, logger)
    
    # Step 4: Export chromosome VCF
    chr_output_dir = f"{output_dir}/{chrom}"
    hl.hadoop_mkdir(chr_output_dir)
    
    chr_summary = export_chromosome_vcf(mt_sampled, chr_output_dir, chrom, target_snps, config, logger, base_dir)
    
    logger.info(f"{chrom}: Completed - exported {chr_summary['final_sampled_variants']} SNPs")
    
    return chr_summary

def sample_snps_from_matrix(mt, target_snps, config, logger):
    """Sample exact number of SNPs from a MatrixTable."""
    random_seed = config['sampling']['random_seed']
    max_iterations = 10
    
    total_variants = mt.count_rows()
    logger.info(f"Available variants: {total_variants}, Target: {target_snps}")
    
    if total_variants <= target_snps:
        logger.warning(f"Only {total_variants} variants available, using all")
        return mt
    
    # Iterative sampling
    mt_sampled = None
    estimated_sampled = 0
    iteration = 0
    
    while estimated_sampled < target_snps and iteration < max_iterations:
        iteration += 1
        remaining_needed = target_snps - estimated_sampled
        remaining_variants = total_variants - estimated_sampled
        
        if remaining_variants <= 0:
            break
        
        sampling_fraction = min(1.0, (remaining_needed * 1.1) / remaining_variants)
        logger.info(f"Sampling iteration {iteration}: fraction={sampling_fraction:.6f}")
        
        mt_current = mt.sample_rows(sampling_fraction, seed=random_seed + iteration)
        
        if mt_sampled is None:
            mt_sampled = mt_current
        else:
            mt_sampled = mt_sampled.union_rows(mt_current)
        
        estimated_sampled = mt_sampled.count_rows()
        logger.info(f"Sampled {estimated_sampled} variants so far")
        
        if estimated_sampled >= target_snps:
            break
    
    logger.info(f"Final sampled variants: {estimated_sampled}")
    return mt_sampled

def export_chromosome_vcf(mt, output_dir, chrom, target_snps, config, logger, base_dir):
    """Export a single chromosome VCF with smart partitioning."""
    # Smart partitioning based on variant count
    total_variants = mt.count_rows()
    optimal_partitions = max(5, min(50, total_variants // 1000))
    
    logger.info(f"{chrom}: Coalescing to {optimal_partitions} partitions...")
    mt = mt.naive_coalesce(optimal_partitions)
    
    output_vcf = f"{output_dir}/{chrom}_background_snps.vcf.bgz"
    logger.info(f"{chrom}: Exporting to {output_vcf}")
    
    export_start = time.time()
    hl.export_vcf(
        mt,
        output_vcf,
        parallel='separate_header',
        tabix=True
    )
    export_time = time.time() - export_start
    
    # Create chromosome summary
    summary = {
        'chromosome': chrom,
        'target_snps': target_snps,
        'total_variants_before_sampling': total_variants,
        'final_sampled_variants': mt.count_rows(),
        'export_path': output_vcf,
        'export_time_seconds': export_time,
        'timestamp': datetime.now().isoformat()
    }
    
    summary_path = f"{output_dir}/{chrom}_sampling_summary.json"
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)
    
    logger.info(f"{chrom}: Export completed in {export_time:.1f} seconds")
    return summary

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

def parse_interval_list_and_sample(interval_list_path, target_chrom, config, logger):
    """Parse interval list and return one random interval for target chromosome, or None if not found."""
    logger.info(f"Parsing interval list for {target_chrom}: {interval_list_path}")
    
    target_intervals = []
    hla_region = config['sampling']['hla_region']
    total_lines = 0
    target_lines = 0
    
    try:
        with hl.hadoop_open(interval_list_path, 'r') as f:
            lines = f.readlines()
        
        for line in lines:
            total_lines += 1
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
                
                # Only collect intervals for target chromosome
                if chrom == target_chrom:
                    target_lines += 1
                    target_intervals.append({
                        'chrom': chrom,
                        'start': start,
                        'end': end
                    })
        
        logger.info(f"Processed {total_lines} lines, found {target_lines} {target_chrom} intervals")
        
        if target_intervals:
            # Return one random interval from this file
            sampled_interval = random.choice(target_intervals)
            logger.info(f"Selected {target_chrom} interval: {sampled_interval['start']:,}-{sampled_interval['end']:,}")
            return sampled_interval
        else:
            logger.info(f"No {target_chrom} intervals found in {interval_list_path}")
            return None
        
    except Exception as e:
        logger.error(f"Failed to parse interval list {interval_list_path}: {e}")
        return None

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
    
    # Get EUR samples (cached)
    eur_samples = get_eur_samples(config, logger)
    
    # Filter to EUR samples
    mt = mt.filter_cols(hl.literal(eur_samples).contains(mt.s))
    
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


def export_snps(mt, output_dir, config, logger, base_dir):
    """Sample exactly target SNPs with iterative approach and export as VCF."""
    target_snps = config['sampling']['target_total_snps']
    random_seed = config['sampling']['random_seed']
    max_iterations = 10  # Prevent infinite loops
    
    # Generate filename based on target SNP count
    if target_snps >= 1000000:
        snp_label = f"{target_snps//1000000}M"
    elif target_snps >= 1000:
        snp_label = f"{target_snps//1000}K"
    else:
        snp_label = f"{target_snps}"
    
    base_filename = f"{snp_label}_background_snps"
    
    logger.info(f"Target: {target_snps:,} SNPs ({snp_label})")
    
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
    # Smart partitioning based on variant count
    optimal_partitions = max(10, min(100, estimated_sampled // 1000))  # 1000 variants per partition
    logger.info(f"Coalescing to {optimal_partitions} partitions for VCF export...")
    mt_sampled = mt_sampled.naive_coalesce(optimal_partitions)
    
    output_vcf = f"{output_dir}/{base_filename}.vcf.bgz"
    logger.info(f"Exporting to compressed VCF: {output_vcf}")
    
    export_start = time.time()
    hl.export_vcf(
        mt_sampled,
        output_vcf,
        parallel='separate_header',
        tabix=True
    )
    export_time = time.time() - export_start
    
    # Post-process to create single VCF if needed
    if optimal_partitions > 1:
        logger.info("Multiple partitions detected, copying shards locally for consolidation...")
        
        # Create local directory in results folder
        import os
        local_parts_dir = f"{output_dir.replace('gs://', '/tmp/gs_')}/vcf_parts"
        os.makedirs(local_parts_dir, exist_ok=True)
        logger.info(f"Created local parts directory: {local_parts_dir}")
        
        # Copy all part files locally
        import subprocess
        try:
            # List part files in GS
            list_cmd = ["gsutil", "ls", f"{output_vcf}/part-*.bgz"]
            result = subprocess.run(list_cmd, capture_output=True, text=True, check=True)
            part_files = result.stdout.strip().split('\n')
            
            logger.info(f"Found {len(part_files)} part files to copy")
            
            # Copy each part file locally
            for gs_file in part_files:
                if gs_file.strip():
                    filename = os.path.basename(gs_file)
                    local_file = os.path.join(local_parts_dir, filename)
                    copy_cmd = ["gsutil", "cp", gs_file, local_file]
                    subprocess.run(copy_cmd, check=True)
                    logger.debug(f"Copied {filename} locally")
            
            # Consolidate locally
            logger.info("Consolidating VCF parts locally...")
            consolidated_vcf = f"{output_dir}/{base_filename}_consolidated.vcf.bgz"
            
            concat_cmd = [
                "bcftools", "concat", 
                "-Oz", "-o", consolidated_vcf,
                f"{local_parts_dir}/part-*.bgz"
            ]
            subprocess.run(concat_cmd, check=True)
            
            # Create index
            subprocess.run(["tabix", "-p", "vcf", consolidated_vcf], check=True)
            
            logger.info(f"Consolidated VCF created: {consolidated_vcf}")
            summary['export_path'] = consolidated_vcf
            
            # Copy consolidated VCF to local results directory
            local_results_dir = f"{base_dir}/results"
            os.makedirs(local_results_dir, exist_ok=True)
            local_consolidated_vcf = f"{local_results_dir}/{base_filename}_consolidated.vcf.bgz"
            
            logger.info(f"Copying consolidated VCF to local results: {local_consolidated_vcf}")
            subprocess.run(["gsutil", "cp", consolidated_vcf, local_consolidated_vcf], check=True)
            subprocess.run(["gsutil", "cp", f"{consolidated_vcf}.tbi", f"{local_consolidated_vcf}.tbi"], check=True)
            
            # Clean up temp directory
            import shutil
            shutil.rmtree(local_parts_dir)
            logger.info("Cleaned up temporary files")
            
        except (subprocess.CalledProcessError, Exception) as e:
            logger.warning(f"Could not consolidate VCF: {e}")
            logger.info("Keeping multi-part VCF structure")
            # Clean up temp directory on error
            if 'local_parts_dir' in locals() and os.path.exists(local_parts_dir):
                shutil.rmtree(local_parts_dir)
    
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
    
    summary_path = f"{output_dir}/{base_filename}_sampling_summary.json"
    
    # Create output directory if it doesn't exist
    try:
        hl.hadoop_mkdir(output_dir)
        logger.info(f"Created output directory: {output_dir}")
    except Exception as e:
        logger.debug(f"Directory may already exist or creation failed: {e}")
    
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)
    
    logger.info(f"Sampling summary saved to: {summary_path}")
    return summary


def main():
    config = load_config("config/config.json")
    
    # Initialize Hail with optimized settings
    init_hail("snp_sampling", driver_mem="32g", reference="GRCh38")
    
    # Get base directory
    base_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Get timestamp for output directory
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Generate filename based on target SNP count
    target_snps = config['sampling']['target_total_snps']
    if target_snps >= 1000000:
        snp_label = f"{target_snps//1000000}M"
    elif target_snps >= 1000:
        snp_label = f"{target_snps//1000}K"
    else:
        snp_label = f"{target_snps}"
    
    base_filename = f"{snp_label}_background_snps"
    
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
    output_dir = f"{base_data_dir}/results/FNCV_RVAS_MS/{base_filename}_{timestamp}"
    logger.info(f"Output Directory: {output_dir}")
    
    try:
        # Check for test mode
        test_mode = config['params'].get('test_mode', False)
        test_chromosome = config['params'].get('test_chromosome', None)
        
        if test_mode and test_chromosome:
            logger.info(f"=== TEST MODE: Processing only {test_chromosome} ===")
            # Test mode: process single chromosome with proportional target
            chr_targets = calculate_chromosome_targets(target_snps, logger)
            test_chr_target = chr_targets.get(test_chromosome, 0)
            if test_chr_target == 0:
                logger.error(f"Test chromosome {test_chromosome} not found in targets!")
                sys.exit(1)
            chr_targets = {test_chromosome: test_chr_target}
            logger.info(f"Target for {test_chromosome}: {chr_targets[test_chromosome]} SNPs (based on chromosome size)")
            
            # Process single chromosome
            chr_summary = process_chromosome(test_chromosome, chr_targets[test_chromosome], output_dir, config, logger, base_dir)
            
            if chr_summary:
                # Create test summary
                test_summary = {
                    'test_mode': True,
                    'test_chromosome': test_chromosome,
                    'target_snps': chr_targets[test_chromosome],
                    'final_sampled_variants': chr_summary['final_sampled_variants'],
                    'output_directory': output_dir,
                    'timestamp': datetime.now().isoformat()
                }
                
                summary_path = f"{output_dir}/{base_filename}_test_summary.json"
                with open(summary_path, 'w') as f:
                    json.dump(test_summary, f, indent=2)
                
                logger.info(f"=== TEST MODE COMPLETED: {test_chromosome} ===")
                logger.info(f"Sampled {chr_summary['final_sampled_variants']} SNPs")
                logger.info(f"Results saved to: {output_dir}")
            else:
                logger.error(f"TEST MODE FAILED: No variants found for {test_chromosome}")
        
        else:
            # Full production mode: Process all chromosomes
            # Step 1: Calculate chromosome targets
            chr_targets = calculate_chromosome_targets(target_snps, logger)
            
            # Step 2: Process each chromosome sequentially
            all_summaries = []
            total_sampled = 0
            
            for chrom, target in chr_targets.items():
                chr_summary = process_chromosome(chrom, target, output_dir, config, logger, base_dir)
                if chr_summary:
                    all_summaries.append(chr_summary)
                    total_sampled += chr_summary['final_sampled_variants']
            
            # Step 3: Create overall summary
            overall_summary = {
                'target_snps': target_snps,
                'chromosome_targets': chr_targets,
                'total_sampled_variants': total_sampled,
                'chromosomes_processed': len(all_summaries),
                'output_directory': output_dir,
                'timestamp': datetime.now().isoformat()
            }
            
            summary_path = f"{output_dir}/{base_filename}_overall_summary.json"
            with open(summary_path, 'w') as f:
                json.dump(overall_summary, f, indent=2)
            
            logger.info("=== Chromosome-based SNP sampling completed successfully! ===")
            logger.info(f"Total chromosomes processed: {len(all_summaries)}")
            logger.info(f"Total SNPs sampled: {total_sampled:,}")
            logger.info(f"Results saved to: {output_dir}")
        
    except Exception as e:
        logger.error(f"Pipeline failed: {e}")
        raise

if __name__ == "__main__":
    main()