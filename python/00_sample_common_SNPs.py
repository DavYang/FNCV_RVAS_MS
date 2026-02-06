import os
import sys
import random
import hail as hl
import json
import time
from datetime import datetime
from utils import load_config, init_hail, setup_logger

logger = setup_logger("SNP_sampling")

# Global cache for EUR sample IDs and Hail Set
eur_samples_cache = None
eur_samples_set_cache = None

def get_eur_samples(config, logger):
    """Load EUR sample IDs with caching to avoid repeated loading."""
    global eur_samples_cache
    if eur_samples_cache is None:
        logger.info("Loading EUR ancestry samples (first time only)...")
        eur_samples_cache = load_ancestry_data(config, logger)
        logger.info(f"Cached {len(eur_samples_cache)} EUR samples for reuse")
    return eur_samples_cache

def get_eur_samples_set(config, logger):
    """Get EUR samples as a Hail Set for efficient filtering."""
    global eur_samples_set_cache
    if eur_samples_set_cache is None:
        eur_ids = get_eur_samples(config, logger)
        eur_samples_set_cache = hl.literal(set(eur_ids))
        logger.info("Created Hail Set for EUR sample filtering")
    return eur_samples_set_cache

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

def _discover_shards_lazy(chrom, config, logger):
    """Lazily discover shard files containing intervals for a given chromosome.
    
    Yields shard numbers one at a time. In test mode, the specified test shard
    is yielded first. Additional shards are discovered on-demand by scanning
    interval list files only as the caller requests more.
    """
    interval_list_base = config['vcf_processing']['interval_list_base'].rstrip('/')
    total_interval_files = config['vcf_processing']['total_shards']
    
    test_mode = config['params'].get('test_mode', False)
    test_interval_file = config['params'].get('test_interval_file', None)
    
    skip_set = set()
    
    # In test mode, yield the priority shard first (no scanning needed)
    if test_mode and test_interval_file:
        test_shard = int(test_interval_file.split('.')[0])
        logger.info(f"{chrom}: TEST MODE - yielding priority shard: {test_shard:010d}")
        skip_set.add(test_shard)
        yield test_shard
    
    # Lazily scan remaining shards in shuffled order
    rng = random.Random(config['sampling']['random_seed'])
    all_shard_indices = [i for i in range(total_interval_files) if i not in skip_set]
    rng.shuffle(all_shard_indices)
    
    scanned = 0
    for file_idx in all_shard_indices:
        interval_list_path = f"{interval_list_base}/{file_idx:010d}.interval_list"
        try:
            with hl.hadoop_open(interval_list_path, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('@') or not line:
                        continue
                    parts = line.split()
                    if len(parts) >= 1 and parts[0] == chrom:
                        logger.info(f"{chrom}: Discovered matching shard {file_idx:010d} "
                                    f"(scanned {scanned + 1} files)")
                        yield file_idx
                    break  # Only need first data line
        except Exception as e:
            logger.debug(f"{chrom}: Error reading shard {file_idx:010d}: {e}")
        scanned += 1
        if scanned % 500 == 0:
            logger.info(f"{chrom}: Shard scan progress: "
                        f"{scanned}/{len(all_shard_indices)} files checked")


def process_chromosome(chrom, target_snps, output_dir, config, logger, base_dir):
    """Process a single chromosome: import one shard at a time, filter to
    chromosome intervals, sample SNPs, and export.
    
    Strategy to avoid Spark DAG explosion (union_rows OOM):
      1. Import shard VCF, filter to chromosome.
      2. Count variants in the shard (single Spark job).
      3. If the shard has enough variants, sample and export directly.
      4. If not, checkpoint the shard to a temp table, then move to next shard.
         Checkpointing materializes the data and breaks the DAG.
      5. After accumulating enough variants across shards, sample and export.
    """
    logger.info(f"\n=== Processing {chrom} (target: {target_snps} SNPs) ===")
    
    vcf_base = config['vcf_processing']['sharded_vcf_base'].rstrip('/')
    max_shards = config.get('vcf_processing', {}).get('max_shards_per_chrom', 100)
    seed = config['sampling']['random_seed']
    rng = random.Random(seed)
    
    shard_generator = _discover_shards_lazy(chrom, config, logger)
    
    # Collect VCF paths and their interval counts to decide how many shards we need
    # We import one shard at a time, count, and stop when we have enough
    shard_vcf_paths = []
    total_variants = 0
    shards_used = 0
    
    for shard_num in shard_generator:
        if total_variants >= target_snps or shards_used >= max_shards:
            break
        
        vcf_path = f"{vcf_base}/{shard_num:010d}.vcf.bgz"
        logger.info(f"{chrom}: Importing shard {shard_num:010d}: {vcf_path}")
        
        # Import and count variants for this shard (single Spark job)
        mt_shard = hl.import_vcf(
            vcf_path,
            reference_genome='GRCh38',
            force_bgz=True,
            array_elements_required=False
        )
        # Quick count: filter to chrom rows only, no sample filtering needed for counting
        shard_count = mt_shard.filter_rows(
            mt_shard.locus.contig == chrom
        ).count_rows()
        
        shards_used += 1
        logger.info(f"{chrom}: Shard {shard_num:010d} has {shard_count:,} "
                    f"{chrom} variants")
        
        if shard_count == 0:
            logger.info(f"{chrom}: Shard empty, trying next...")
            continue
        
        shard_vcf_paths.append(vcf_path)
        total_variants += shard_count
        logger.info(f"{chrom}: Running total: {total_variants:,}/{target_snps:,} "
                    f"variants from {shards_used} shard(s)")
    
    if not shard_vcf_paths or total_variants == 0:
        logger.warning(f"{chrom}: No variants found after {shards_used} shard(s)")
        return None
    
    logger.info(f"{chrom}: Collected {len(shard_vcf_paths)} shard path(s) with "
                f"{total_variants:,} total variants. Now importing for sampling...")
    
    # Single import of all needed shards at once (Hail handles multi-VCF efficiently)
    mt = hl.import_vcf(
        shard_vcf_paths,
        reference_genome='GRCh38',
        force_bgz=True,
        array_elements_required=False
    )
    
    # Filter to EUR samples
    eur_set = get_eur_samples_set(config, logger)
    mt = mt.filter_cols(eur_set.contains(mt.s))
    
    # Split multi-allelic variants
    mt = hl.split_multi_hts(mt)
    
    # Filter to target chromosome
    mt = mt.filter_rows(mt.locus.contig == chrom)
    
    # Count after EUR filtering + split (may differ from raw count)
    final_pool = mt.count_rows()
    logger.info(f"{chrom}: Pool after EUR filtering + split: {final_pool:,} variants")
    
    # Sample target SNPs from pool
    if final_pool > target_snps:
        sampling_fraction = target_snps / final_pool
        logger.info(f"{chrom}: Sampling {target_snps:,} from {final_pool:,} "
                    f"(fraction: {sampling_fraction:.6f})")
        mt_sampled = mt.sample_rows(sampling_fraction, seed=seed)
        sampled_count = mt_sampled.count_rows()
        
        # Trim if overshot
        if sampled_count > target_snps:
            logger.info(f"{chrom}: Trimming {sampled_count:,} -> {target_snps:,}")
            mt_sampled = mt_sampled.head(target_snps)
            sampled_count = target_snps
        
        logger.info(f"{chrom}: Sampled {sampled_count:,} variants")
    else:
        logger.info(f"{chrom}: Using all {final_pool:,} variants "
                    f"(fewer than target {target_snps:,})")
        mt_sampled = mt
        sampled_count = final_pool
    
    # Export chromosome VCF
    chr_output_dir = f"{output_dir}/{chrom}"
    chr_summary = export_chromosome_vcf(mt_sampled, chr_output_dir, chrom, target_snps,
                                        config, logger, base_dir)
    chr_summary['shards_used'] = shards_used
    chr_summary['total_pool_variants'] = final_pool
    
    logger.info(f"{chrom}: Completed - exported {chr_summary['final_sampled_variants']:,} SNPs "
                f"from {shards_used} shard(s)")
    return chr_summary


def export_chromosome_vcf(mt, output_dir, chrom, target_snps, config, logger, base_dir):
    """Export a single chromosome VCF as one .vcf.bgz file."""
    # Count once and reuse
    total_variants = mt.count_rows()
    logger.info(f"{chrom}: Exporting {total_variants:,} variants")
    
    # Coalesce to single partition for a single-file VCF output
    mt = mt.naive_coalesce(1)
    
    output_vcf = f"{output_dir}/{chrom}_background_snps.vcf.bgz"
    logger.info(f"{chrom}: Exporting to {output_vcf}")
    
    export_start = time.time()
    hl.export_vcf(
        mt,
        output_vcf,
        tabix=True
    )
    export_time = time.time() - export_start
    
    # Create chromosome summary
    summary = {
        'chromosome': chrom,
        'target_snps': target_snps,
        'final_sampled_variants': total_variants,
        'export_path': output_vcf,
        'export_time_seconds': round(export_time, 1),
        'timestamp': datetime.now().isoformat()
    }
    
    # Write summary to cloud storage
    summary_path = f"{output_dir}/{chrom}_sampling_summary.json"
    with hl.hadoop_open(summary_path, 'w') as f:
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
            chr_targets = calculate_chromosome_targets(target_snps, logger)
            test_chr_target = chr_targets.get(test_chromosome, 0)
            if test_chr_target == 0:
                logger.error(f"Test chromosome {test_chromosome} not found in targets!")
                sys.exit(1)
            chr_targets = {test_chromosome: test_chr_target}
            logger.info(f"Target for {test_chromosome}: {chr_targets[test_chromosome]:,} SNPs (based on chromosome size)")
            
            chr_summary = process_chromosome(test_chromosome, chr_targets[test_chromosome],
                                            output_dir, config, logger, base_dir)
            
            if chr_summary:
                test_summary = {
                    'test_mode': True,
                    'test_chromosome': test_chromosome,
                    'target_snps': chr_targets[test_chromosome],
                    'final_sampled_variants': chr_summary['final_sampled_variants'],
                    'output_directory': output_dir,
                    'timestamp': datetime.now().isoformat()
                }
                
                summary_path = f"{output_dir}/{base_filename}_test_summary.json"
                with hl.hadoop_open(summary_path, 'w') as f:
                    json.dump(test_summary, f, indent=2)
                
                logger.info(f"=== TEST MODE COMPLETED: {test_chromosome} ===")
                logger.info(f"Sampled {chr_summary['final_sampled_variants']:,} SNPs")
                logger.info(f"Results saved to: {output_dir}")
            else:
                logger.error(f"TEST MODE FAILED: No variants found for {test_chromosome}")
        
        else:
            # Full production mode: Process all chromosomes
            chr_targets = calculate_chromosome_targets(target_snps, logger)
            
            all_summaries = []
            total_sampled = 0
            
            for chrom, target in chr_targets.items():
                chr_summary = process_chromosome(chrom, target, output_dir, config, logger, base_dir)
                if chr_summary:
                    all_summaries.append(chr_summary)
                    total_sampled += chr_summary['final_sampled_variants']
            
            overall_summary = {
                'target_snps': target_snps,
                'chromosome_targets': chr_targets,
                'total_sampled_variants': total_sampled,
                'chromosomes_processed': len(all_summaries),
                'output_directory': output_dir,
                'timestamp': datetime.now().isoformat()
            }
            
            summary_path = f"{output_dir}/{base_filename}_overall_summary.json"
            with hl.hadoop_open(summary_path, 'w') as f:
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