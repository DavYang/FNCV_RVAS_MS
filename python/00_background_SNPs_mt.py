#!/usr/bin/env python3
import os
import sys
import random
import hail as hl
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


def process_chromosome_intervals(mt, chrom, intervals, eur_sample_ids, target_variants_chr, logger):
    """Process chromosome using interval-based loading for memory efficiency."""
    
    logger.info(f"Processing chromosome {chrom} with {len(intervals)} intervals...")
    
    processed_intervals = []
    total_variants_chr = 0
    total_filtered_variants_chr = 0
    
    try:
        for interval_idx, interval in enumerate(intervals):
            logger.info(f"Chr {chrom}: Processing interval {interval_idx + 1}/{len(intervals)}")
            
            # Load only the specific interval (much more memory efficient!)
            mt_interval = mt.filter_intervals([interval])
            
            # Quick count
            interval_variants = mt_interval.count_rows()
            total_variants_chr += interval_variants
            
            if interval_variants == 0:
                logger.info(f"Chr {chrom} interval {interval_idx + 1}: No variants")
                continue
            
            # Filter to EUR samples for this interval
            mt_interval_eur = mt_interval.filter_cols(
                hl.literal(eur_sample_ids).contains(mt_interval.s)
            )
            
            filtered_variants = mt_interval_eur.count_rows()
            total_filtered_variants_chr += filtered_variants
            
            logger.info(f"Chr {chrom} interval {interval_idx + 1}: {filtered_variants:,} filtered variants")
            
            if filtered_variants > 0:
                processed_intervals.append(mt_interval_eur)
            
            # Clear memory for this interval
            del mt_interval, mt_interval_eur
        
        logger.info(f"Chr {chrom}: Total {total_variants_chr:,} variants, {total_filtered_variants_chr:,} filtered")
        
        # Combine all intervals for this chromosome
        if processed_intervals:
            mt_chrom_combined = hl.MatrixTable.union_rows(*processed_intervals)
            
            # Calculate final sampling for this chromosome
            if total_filtered_variants_chr > 0:
                sampling_fraction = target_variants_chr / total_filtered_variants_chr
                logger.info(f"Chr {chrom}: Final sampling fraction = {sampling_fraction:.6f}")
                
                # Sample variants from combined chromosome data
                mt_sampled = mt_chrom_combined.sample_rows(sampling_fraction, seed=42)
                actual_variants = mt_sampled.count_rows()
                logger.info(f"Chr {chrom}: Sampled {actual_variants:,} variants")
                
                # Clear memory
                del mt_chrom_combined, processed_intervals
                
                return mt_sampled
            else:
                logger.warning(f"Chr {chrom}: No filtered variants for sampling")
                return None
        else:
            logger.warning(f"Chr {chrom}: No intervals with data")
            return None
            
    except Exception as e:
        logger.error(f"Error processing chromosome {chrom}: {str(e)}")
        return None


def process_all_chromosomes(mt, eur_sample_ids, config, logger):
    """Process all chromosomes using interval-based loading with single session."""
    
    # Calculate targets and intervals per chromosome
    targets_per_chr, intervals_per_chr = calculate_chromosome_intervals_and_targets(
        config['params']['target_variants']
    )
    
    logger.info(f"Chromosome targets: {targets_per_chr}")
    
    # Count total intervals for progress tracking
    total_intervals = sum(len(intervals) for intervals in intervals_per_chr.values())
    logger.info(f"Total intervals to process: {total_intervals}")
    
    # Create output directory for individual chromosome MTs
    output_dir = os.path.join(
        config['outputs']['base_dir'],
        config['outputs']['data_dir_suffix'],
        "background_snps"
    )
    os.makedirs(output_dir, exist_ok=True)
    
    # Sort chromosomes by size (smallest first) to test stability
    chromosome_order = sorted(targets_per_chr.keys(), 
                            key=lambda x: targets_per_chr[x], 
                            reverse=False)
    
    logger.info(f"Processing order (smallest to largest): {chromosome_order}")
    
    # Process each chromosome with single session
    processed_chromosomes = []
    
    for chrom in chromosome_order:
        intervals = intervals_per_chr[chrom]
        logger.info(f"Starting chromosome {chrom} ({len(intervals)} intervals)...")
        
        try:
            mt_chrom_sampled = process_chromosome_intervals(
                mt, chrom, intervals, eur_sample_ids, 
                targets_per_chr[chrom], logger
            )
            
            if mt_chrom_sampled is not None:
                processed_chromosomes.append(chrom)
                
                # Export individual chromosome MT
                chromosome_output_path = os.path.join(output_dir, f"chr_{chrom}_sampled.mt")
                mt_chrom_sampled.write(chromosome_output_path, overwrite=True)
                logger.info(f"Exported {chrom} to {chromosome_output_path}")
                
                # Clear memory
                del mt_chrom_sampled
                logger.info(f"Chr {chrom}: Memory cleared")
            else:
                logger.warning(f"Chr {chrom}: No data processed")
        
        except Exception as e:
            logger.error(f"Error processing chromosome {chrom}: {e}")
            # Continue to next chromosome instead of stopping
            continue
    
    logger.info(f"Successfully processed {len(processed_chromosomes)} chromosomes: {processed_chromosomes}")
    return processed_chromosomes


def combine_chromosome_results(processed_chromosomes, config, logger):
    """Combine all chromosome results into final dataset."""
    
    logger.info("Combining chromosome results...")
    
    if not processed_chromosomes:
        raise ValueError("No chromosomes processed successfully")
    
    # Combine all chromosome MatrixTables
    final_mt = hl.MatrixTable.union_rows(*processed_chromosomes)
    
    # Final counts
    total_variants = final_mt.count_rows()
    total_samples = final_mt.count_cols()
    
    logger.info(f"Final dataset: {total_variants:,} variants, {total_samples:,} samples")
    
    return final_mt


def sample_background_snps(config, logger):
    """Main function to sample 1M random common SNPs from EUR ancestry samples."""
    logger.info("Starting background SNP sampling...")
    
    # Initialize Hail once for the entire processing
    hl.init(
        log='/tmp/hail_background_snps.log',
        spark_conf={
            'spark.driver.memory': '12g',  # More memory for single session
            'spark.executor.memory': '12g',
            'spark.network.timeout': '800s',
            'spark.executor.heartbeatInterval': '60s'
        }
    )
    hl.default_reference('GRCh38')
    
    try:
        # Load ancestry data
        eur_sample_ids = load_ancestry_data(config, logger)
        
        # Load pre-filtered ACAF MatrixTable
        logger.info("Loading ACAF MatrixTable...")
        mt_path = config['inputs']['wgs_matrix_table']
        
        # Read the existing MatrixTable directly
        mt = hl.read_matrix_table(mt_path)
        
        logger.info(f"Loaded MatrixTable with {mt.count_rows():,} variants and {mt.count_cols():,} samples")
        
        # Process chromosomes with single session
        processed_chromosomes = process_all_chromosomes(mt, eur_sample_ids, config, logger)
        
        # Generate summary report
        output_dir = os.path.join(
            config['outputs']['base_dir'],
            config['outputs']['data_dir_suffix'],
            "background_snps"
        )
        
        summary = {
            'timestamp': datetime.now().isoformat(),
            'total_samples': len(eur_sample_ids),
            'target_variants': config['params']['target_variants'],
            'ancestry': 'EUR',
            'processing_method': 'single_session_individual_exports',
            'chromosomes_processed': processed_chromosomes,
            'chromosome_count': len(processed_chromosomes),
            'output_directory': output_dir,
            'individual_chromosome_mts': [f"chr_{chrom}_sampled.mt" for chrom in processed_chromosomes],
            'next_step': 'Run merge_chromosome_mts.py to combine results'
        }
        
        # Save summary report
        summary_path = os.path.join(output_dir, "chromosome_processing_summary.json")
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
    finally:
        # Always stop Hail session
        logger.info("Stopping Hail session...")
        hl.stop()


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