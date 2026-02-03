#!/usr/bin/env python3
import os
import sys
import random
import hail as hl
import time
import pandas as pd
from datetime import datetime
from utils import load_config, init_hail, setup_logger


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


def calculate_chromosome_intervals_and_targets(target_variants=1000000):
    """Calculate intervals and variant targets per chromosome based on GRCh38 lengths."""
    
    chr_lengths = {
        'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555,
        'chr5': 181538259, 'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636,
        'chr9': 138394717, 'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
        'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189, 'chr16': 90338345,
        'chr17': 83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167,
        'chr21': 46709983, 'chr22': 50818468
    }
    
    total_autosome_length = sum(chr_lengths.values())
    targets_per_chr = {}
    
    # Calculate targets per chromosome
    for chrom, length in chr_lengths.items():
        proportion = length / total_autosome_length
        targets_per_chr[chrom] = int(target_variants * proportion)
    
    # Generate intervals for each chromosome (10Mb intervals for manageable chunks)
    intervals_per_chr = {}
    for chrom, length in chr_lengths.items():
        intervals = []
        # Use 10Mb intervals to reduce memory load
        interval_size = 10000000  # 10Mb
        
        for start in range(1, length, interval_size):
            end = min(start + interval_size, length)
            interval = hl.Interval(
                start=hl.Locus(chrom, start, reference_genome='GRCh38'),
                end=hl.Locus(chrom, end, reference_genome='GRCh38'),
                includes_start=True,
                includes_end=False
            )
            intervals.append(interval)
        
        intervals_per_chr[chrom] = intervals
    
    return targets_per_chr, intervals_per_chr


def check_session_health(logger):
    """Check if Hail session is healthy."""
    try:
        # Try a simple operation to check session health
        hl.utils.range_table(1, 1).count()
        return True
    except Exception as e:
        logger.warning(f"Session health check failed: {e}")
        return False


def initialize_fresh_session(logger):
    """Initialize a fresh Hail session with conservative settings."""
    try:
        # Stop any existing session first
        try:
            hl.stop()
        except:
            pass
        
        # Initialize fresh session with conservative memory
        hl.init(
            log='/tmp/hail_background_snps.log',
            spark_conf={
                'spark.driver.memory': '6g',   # Conservative for per-chromosome sessions
                'spark.executor.memory': '6g', # Conservative for per-chromosome sessions
                'spark.network.timeout': '1200s',
                'spark.executor.heartbeatInterval': '60s',
                'spark.yarn.am.waitTime': '600s',
                'spark.yarn.applicationMaster.waitTime': '600s',
                'spark.yarn.maxAppAttempts': '3'
            }
        )
        hl.default_reference('GRCh38')
        logger.info("Fresh Hail session initialized")
        return True
    except Exception as e:
        logger.error(f"Failed to initialize Hail session: {e}")
        return False


def load_and_filter_eur_data(config, logger):
    """Load MatrixTable and filter to EUR samples for current session."""
    try:
        # Load ancestry data
        logger.info("Loading ancestry data...")
        ancestry_ht = hl.import_table(
            config['inputs']['ancestry_pred'],
            impute=True,
            types={'research_id': hl.tstr, 'ancestry_pred': hl.tstr}
        )
        
        # Filter to European ancestry samples
        eur_samples = ancestry_ht.filter(ancestry_ht.ancestry_pred == 'eur')
        eur_sample_ids = eur_samples.aggregate(hl.agg.collect(eur_samples.research_id))
        logger.info(f"Found {len(eur_sample_ids)} European ancestry samples")
        
        # Load MatrixTable
        logger.info("Loading ACAF MatrixTable...")
        mt_path = config['inputs']['wgs_matrix_table']
        mt = hl.read_matrix_table(mt_path)
        logger.info(f"Loaded MatrixTable with {mt.count_rows():,} variants and {mt.count_cols():,} samples")
        
        # Filter to EUR samples
        logger.info("Filtering to EUR samples...")
        mt_eur = mt.filter_cols(hl.literal(eur_sample_ids).contains(mt.s))
        logger.info(f"Filtered to {mt_eur.count_cols():,} EUR samples")
        
        return mt_eur, eur_sample_ids
        
    except Exception as e:
        logger.error(f"Failed to load and filter EUR data: {e}")
        raise


def process_chromosome_with_session(mt_eur, chrom, intervals, target_variants_chr, logger):
    """Process chromosome in current session with immediate export."""
    
    logger.info(f"Chr {chrom}: Processing {len(intervals)} intervals...")
    
    processed_intervals = []
    total_filtered_variants_chr = 0
    
    for interval_idx, interval in enumerate(intervals):
        logger.info(f"Chr {chrom}: Processing interval {interval_idx + 1}/{len(intervals)}")
        
        # Load only the specific interval
        mt_interval = mt_eur.filter_rows(
            (mt_eur.locus.contig == chrom) & 
            (mt_eur.locus.position >= interval.start.position) & 
            (mt_eur.locus.position < interval.end.position)
        )
        
        # Quick count
        interval_variants = mt_interval.count_rows()
        
        if interval_variants == 0:
            logger.info(f"Chr {chrom} interval {interval_idx + 1}: No variants")
            continue
        
        filtered_variants = mt_interval.count_rows()
        total_filtered_variants_chr += filtered_variants
        
        logger.info(f"Chr {chrom} interval {interval_idx + 1}: {filtered_variants:,} filtered variants")
        
        if filtered_variants > 0:
            processed_intervals.append(mt_interval)
        
        # Clear memory for this interval
        del mt_interval
        
        # Small delay between intervals
        time.sleep(1)
    
    logger.info(f"Chr {chrom}: Total {total_filtered_variants_chr:,} filtered variants from {len(processed_intervals)} intervals")
    
    # Combine all intervals for this chromosome
    if processed_intervals:
        logger.info(f"Chr {chrom}: Combining {len(processed_intervals)} intervals...")
        mt_chrom_combined = hl.MatrixTable.union_rows(*processed_intervals)
        
        # Calculate final sampling for this chromosome
        if total_filtered_variants_chr > 0:
            sampling_fraction = target_variants_chr / total_filtered_variants_chr
            logger.info(f"Chr {chrom}: Final sampling fraction = {sampling_fraction:.6f}")
            
            # Sample variants from combined chromosome data
            mt_sampled = mt_chrom_combined.sample_rows(sampling_fraction, seed=42)
            actual_variants = mt_sampled.count_rows()
            logger.info(f"Chr {chrom}: Sampled {actual_variants:,} variants")
            
            # Clear memory aggressively
            del mt_chrom_combined, processed_intervals
            
            return mt_sampled
        else:
            logger.warning(f"Chr {chrom}: No filtered variants for sampling")
            return None
    else:
        logger.warning(f"Chr {chrom}: No intervals with data")
        return None


def process_all_chromosomes_with_sessions(config, logger):
    """Process all chromosomes with fresh Hail session per chromosome."""
    
    # Calculate targets and intervals per chromosome
    targets_per_chr, intervals_per_chr = calculate_chromosome_intervals_and_targets(
        config['params']['target_variants']
    )
    
    logger.info(f"Chromosome targets: {targets_per_chr}")
    
    # Count total intervals for progress tracking
    total_intervals = sum(len(intervals) for intervals in intervals_per_chr.values())
    logger.info(f"Total intervals to process: {total_intervals}")
    
    # Create output directory using WORKSPACE_BUCKET
    workspace_bucket = os.environ.get('WORKSPACE_BUCKET', 'gs://default-bucket')
    output_dir = f"{workspace_bucket}/results/FNCV_RVAS_MS/background_snps"
    logger.info(f"Output directory: {output_dir}")
    
    # Sort chromosomes by size (smallest first) for stability
    chromosome_order = sorted(targets_per_chr.keys(), 
                            key=lambda x: targets_per_chr[x], 
                            reverse=False)
    
    logger.info(f"Processing order (smallest to largest): {chromosome_order}")
    logger.info("Per-chromosome session processing with EUR reloading!")
    
    # Process each chromosome with fresh session
    processed_chromosomes = []
    failed_chromosomes = {}
    
    for chrom_idx, chrom in enumerate(chromosome_order):
        intervals = intervals_per_chr[chrom]
        logger.info(f"Starting chromosome {chrom} ({len(intervals)} intervals) - {chrom_idx + 1}/{len(chromosome_order)}")
        
        try:
            # Initialize fresh session for this chromosome
            if not initialize_fresh_session(logger):
                raise Exception("Failed to initialize Hail session")
            
            # Load and filter EUR data for this session
            mt_eur, eur_sample_ids = load_and_filter_eur_data(config, logger)
            
            # Process chromosome in current session
            mt_chrom_sampled = process_chromosome_with_session(
                mt_eur, chrom, intervals, targets_per_chr[chrom], logger
            )
            
            if mt_chrom_sampled is not None:
                processed_chromosomes.append(chrom)
                
                # Export immediately to WORKSPACE_BUCKET
                chromosome_output_path = f"{output_dir}/chr{chrom}_background_snps.mt"
                mt_chrom_sampled.write(chromosome_output_path, overwrite=True)
                logger.info(f"Exported {chrom} to {chromosome_output_path}")
                
                # Create summary for this chromosome
                chr_summary = {
                    'chromosome': chrom,
                    'target_variants': targets_per_chr[chrom],
                    'actual_variants': mt_chrom_sampled.count_rows(),
                    'output_path': chromosome_output_path,
                    'timestamp': datetime.now().isoformat()
                }
                
                # Save chromosome summary
                chr_summary_path = f"{output_dir}/chr{chrom}_summary.json"
                with open(chr_summary_path, 'w') as f:
                    json.dump(chr_summary, f, indent=2)
                
                # Clear memory aggressively
                del mt_chrom_sampled, mt_eur
                logger.info(f"Chr {chrom}: Memory cleared and exported")
                
                # Force garbage collection
                import gc
                gc.collect()
                
            else:
                logger.warning(f"Chr {chrom}: No data processed")
                failed_chromosomes[chrom] = "No data processed"
        
        except Exception as e:
            logger.error(f"Chromosome {chrom} failed: {e}")
            failed_chromosomes[chrom] = str(e)
            # Continue to next chromosome instead of stopping
            continue
        
        finally:
            # Always stop session after each chromosome
            try:
                logger.info(f"Chr {chrom}: Stopping Hail session...")
                hl.stop()
            except:
                pass
        
        # Delay between chromosomes to let resources stabilize
        logger.info(f"Chr {chrom}: Waiting 5 seconds before next chromosome...")
        time.sleep(5)
    
    logger.info(f"Successfully processed {len(processed_chromosomes)} chromosomes: {processed_chromosomes}")
    if failed_chromosomes:
        logger.info(f"Failed chromosomes: {list(failed_chromosomes.keys())}")
        for chrom, error in failed_chromosomes.items():
            logger.error(f"  {chrom}: {error}")
    
    return processed_chromosomes, output_dir


def sample_background_snps(config, logger):
    """Main function to sample 1M random common SNPs from EUR ancestry samples."""
    logger.info("Starting background SNP sampling...")
    
    try:
        # Process chromosomes with per-chromosome sessions
        processed_chromosomes, output_dir = process_all_chromosomes_with_sessions(config, logger)
        
        # Generate summary report
        summary = {
            'timestamp': datetime.now().isoformat(),
            'target_variants': config['params']['target_variants'],
            'ancestry': 'EUR',
            'processing_method': 'per_chromosome_sessions',
            'chromosomes_processed': processed_chromosomes,
            'chromosome_count': len(processed_chromosomes),
            'output_directory': output_dir,
            'individual_chromosome_mts': [f"chr{chrom}_background_snps.mt" for chrom in processed_chromosomes],
            'next_step': 'Run merge_chromosome_mts.py to combine results'
        }
        
        # Save summary report to WORKSPACE_BUCKET
        summary_path = f"{output_dir}/chromosome_processing_summary.json"
        import json
        with open(summary_path, 'w') as f:
            json.dump(summary, f, indent=2)
        
        logger.info(f"Chromosome processing completed!")
        logger.info(f"Summary saved to {summary_path}")
        logger.info(f"Individual chromosome MTs saved to {output_dir}")
        logger.info(f"Next step: Run 'python python/merge_chromosome_mts.py' to combine results")
        
        return summary
        
    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        raise


def main():
    """Main execution function."""
    # Setup logging
    logger = setup_logger("background_snps")
    
    try:
        # Load configuration
        config = load_config()
        logger.info("Configuration loaded successfully")
        
        # Run background SNP sampling
        summary = sample_background_snps(config, logger)
        
        logger.info("Analysis completed successfully!")
        return summary
        
    except Exception as e:
        logger.error(f"Analysis failed: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()