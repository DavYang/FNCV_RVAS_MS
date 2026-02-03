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


def process_chromosome_single_session(mt, chrom, intervals, eur_sample_ids, target_variants_chr, logger):
    """Process chromosome with single session for all intervals."""
    
    try:
        # Single session for entire chromosome
        hl.init(
            log=f'/tmp/hail_chromosome_{chrom}.log',
            spark_conf={
                'spark.driver.memory': '8g',    # Moderate memory for chromosome
                'spark.executor.memory': '8g',  # Moderate memory for chromosome
                'spark.network.timeout': '1200s',
                'spark.executor.heartbeatInterval': '60s',
                'spark.yarn.am.waitTime': '600s',
                'spark.yarn.applicationMaster.waitTime': '600s',
                'spark.yarn.maxAppAttempts': '3'
            }
        )
        hl.default_reference('GRCh38')
        
        logger.info(f"Chr {chrom}: Processing {len(intervals)} intervals with single session...")
        
        processed_intervals = []
        total_filtered_variants_chr = 0
        
        for interval_idx, interval in enumerate(intervals):
            logger.info(f"Chr {chrom}: Processing interval {interval_idx + 1}/{len(intervals)}")
            
            # Load only the specific interval (no new session)
            mt_interval = mt.filter_rows(
                (mt.locus.contig == chrom) & 
                (mt.locus.position >= interval.start.position) & 
                (mt.locus.position < interval.end.position)
            )
            
            # Quick count
            interval_variants = mt_interval.count_rows()
            
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
            
            # Small delay between intervals (no session restart)
            time.sleep(1)
        
        logger.info(f"Chr {chrom}: Total {total_filtered_variants_chr:,} filtered variants from {len(processed_intervals)} intervals")
        
        # Combine all intervals for this chromosome (same session)
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
        logger.error(f"Chr {chrom}: Failed to process - {str(e)}")
        return None
    finally:
        # Always cleanup session
        try:
            hl.stop()
            logger.info(f"Chr {chrom}: Session stopped")
        except:
            pass


def process_all_chromosomes_chromosome_sessions(mt, eur_sample_ids, config, logger, max_concurrent=2):
    """Process chromosomes with one session per chromosome (resource-efficient)."""
    
    import time
    import threading
    from concurrent.futures import ThreadPoolExecutor, as_completed
    
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
    logger.info(f"Max concurrent chromosomes: {max_concurrent}")
    
    def process_chromosome_worker(chrom):
        """Worker function for processing a single chromosome."""
        intervals = intervals_per_chr[chrom]
        logger.info(f"[Worker] Starting chromosome {chrom} ({len(intervals)} intervals)...")
        
        try:
            # Process chromosome with single session
            mt_chrom_sampled = process_chromosome_single_session(
                mt, chrom, intervals, eur_sample_ids, 
                targets_per_chr[chrom], logger
            )
            
            if mt_chrom_sampled is not None:
                # Export individual chromosome MT
                chromosome_output_path = os.path.join(output_dir, f"chr_{chrom}_sampled.mt")
                mt_chrom_sampled.write(chromosome_output_path, overwrite=True)
                logger.info(f"[Worker] Exported {chrom} to {chromosome_output_path}")
                
                # Clear memory
                del mt_chrom_sampled
                logger.info(f"[Worker] Chr {chrom}: Memory cleared")
                
                return chrom, True, None
            else:
                logger.warning(f"[Worker] Chr {chrom}: No data processed")
                return chrom, False, "No data processed"
        
        except Exception as e:
            logger.error(f"[Worker] Chromosome {chrom} failed: {e}")
            return chrom, False, str(e)
    
    # Process chromosomes with limited concurrency
    processed_chromosomes = []
    failed_chromosomes = {}
    
    with ThreadPoolExecutor(max_workers=max_concurrent) as executor:
        # Submit all chromosome jobs
        future_to_chrom = {
            executor.submit(process_chromosome_worker, chrom): chrom 
            for chrom in chromosome_order
        }
        
        # Process completed jobs
        for future in as_completed(future_to_chrom):
            chrom = future_to_chrom[future]
            try:
                result_chrom, success, error = future.result()
                if success:
                    processed_chromosomes.append(result_chrom)
                    logger.info(f"[Main] Chromosome {result_chrom} completed successfully")
                else:
                    failed_chromosomes[result_chrom] = error
                    logger.error(f"[Main] Chromosome {result_chrom} failed: {error}")
                
                # Small delay between completions to let resources stabilize
                time.sleep(5)
                
            except Exception as e:
                logger.error(f"[Main] Unexpected error for chromosome {chrom}: {e}")
                failed_chromosomes[chrom] = str(e)
    
    logger.info(f"Successfully processed {len(processed_chromosomes)} chromosomes: {processed_chromosomes}")
    if failed_chromosomes:
        logger.info(f"Failed chromosomes: {list(failed_chromosomes.keys())}")
        for chrom, error in failed_chromosomes.items():
            logger.error(f"  {chrom}: {error}")
    
    return processed_chromosomes


def sample_background_snps(config, logger):
    """Main function to sample 1M random common SNPs from EUR ancestry samples."""
    logger.info("Starting background SNP sampling...")
    
    # Initialize Hail once for the entire processing
    hl.init(
        log='/tmp/hail_background_snps.log',
        spark_conf={
            'spark.driver.memory': '8g',  # Reduced to fit cluster limits
            'spark.executor.memory': '8g',  # Reduced to fit cluster limits
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
        
        # Stop main session before per-interval processing
        logger.info("Stopping main session to prepare for per-interval processing...")
        hl.stop()
        
        # Process chromosomes with chromosome-level sessions
        processed_chromosomes = process_all_chromosomes_chromosome_sessions(mt, eur_sample_ids, config, logger)
        
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
            'processing_method': 'chromosome_level_sessions',
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
        # Always stop Hail session (if still running)
        try:
            logger.info("Stopping Hail session...")
            hl.stop()
        except:
            pass


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