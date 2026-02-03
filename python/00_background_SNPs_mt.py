#!/usr/bin/env python3
import os
import sys
import random
import hail as hl
import time
import pandas as pd
import json
import signal
import threading
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


class SessionTimeoutError(Exception):
    """Custom exception for session initialization timeout."""
    pass


def timeout_handler(signum, frame):
    """Handler for session timeout."""
    raise SessionTimeoutError("Session initialization timed out")


def check_session_health(logger):
    """Check if Hail session is healthy."""
    try:
        logger.debug("Performing session health check...")
        # Try a simple operation to check session health
        result = hl.utils.range_table(1, 1).count()
        logger.debug(f"Session health check passed: {result}")
        return True
    except Exception as e:
        logger.warning(f"Session health check failed: {e}")
        return False


def cleanup_resources(logger):
    """Clean up resources before initializing new session."""
    try:
        logger.info("Cleaning up resources before session initialization...")
        
        # Stop any existing Hail session
        try:
            logger.debug("Stopping existing Hail session...")
            hl.stop()
            logger.debug("Existing Hail session stopped")
        except Exception as e:
            logger.debug(f"No existing session to stop: {e}")
        
        # Force garbage collection
        import gc
        gc.collect()
        logger.debug("Garbage collection completed")
        
        # Small delay to let resources stabilize
        time.sleep(2)
        logger.info("Resource cleanup completed")
        
    except Exception as e:
        logger.warning(f"Resource cleanup failed: {e}")


def initialize_fresh_session(logger, timeout_seconds=600):
    """Initialize a fresh Hail session with conservative settings and timeout."""
    logger.info(f"Initializing fresh Hail session (timeout: {timeout_seconds}s)...")
    
    try:
        # Set up timeout
        signal.signal(signal.SIGALRM, timeout_handler)
        signal.alarm(timeout_seconds)
        
        logger.debug("Step 1: Cleaning up resources...")
        cleanup_resources(logger)
        
        logger.debug("Step 2: Configuring Spark settings...")
        spark_conf = {
            'spark.driver.memory': '6g',   # Increased from 4g for AM stability
            'spark.executor.memory': '6g', # Increased from 4g for AM stability
            'spark.driver.cores': '1',      # Reduced cores to reduce resource pressure
            'spark.executor.cores': '1',    # Reduced cores to reduce resource pressure
            'spark.network.timeout': '2400s', # Increased from 1200s
            'spark.executor.heartbeatInterval': '120s', # Increased from 60s
            'spark.yarn.am.waitTime': '1800s', # Increased from 600s (30 minutes)
            'spark.yarn.applicationMaster.waitTime': '1800s', # Increased from 600s
            'spark.yarn.maxAppAttempts': '1', # Reduced to prevent retry loops
            'spark.yarn.am.memory': '2g',    # Explicit AM memory
            'spark.yarn.am.cores': '1',     # Reduced AM cores
            'spark.driver.extraJavaOptions': '-XX:+UseG1GC -XX:MaxGCPauseMillis=200 -XX:+UnlockExperimentalVMOptions -XX:+UseCGroupMemoryLimitForHeap',
            'spark.executor.extraJavaOptions': '-XX:+UseG1GC -XX:MaxGCPauseMillis=200 -XX:+UnlockExperimentalVMOptions -XX:+UseCGroupMemoryLimitForHeap'
        }
        logger.debug(f"Spark configuration: {spark_conf}")
        
        logger.debug("Step 3: Starting Hail initialization...")
        start_time = time.time()
        
        hl.init(
            log='/tmp/hail_background_snps.log',
            spark_conf=spark_conf,
            quiet=False
        )
        
        init_time = time.time() - start_time
        logger.info(f"Hail initialization completed in {init_time:.1f} seconds")
        
        logger.debug("Step 4: Setting reference genome...")
        hl.default_reference('GRCh38')
        logger.debug("Reference genome set to GRCh38")
        
        logger.debug("Step 5: Verifying session health...")
        if check_session_health(logger):
            logger.info("Fresh Hail session initialized and verified")
            signal.alarm(0)  # Cancel timeout
            return True
        else:
            logger.error("Session initialization failed health check")
            signal.alarm(0)
            return False
            
    except SessionTimeoutError as e:
        logger.error(f"Session initialization timed out after {timeout_seconds} seconds")
        signal.alarm(0)
        return False
    except Exception as e:
        logger.error(f"Failed to initialize Hail session: {e}")
        logger.debug(f"Session initialization error details: {type(e).__name__}: {e}")
        signal.alarm(0)
        return False


def load_and_filter_eur_data(config, logger):
    """Load MatrixTable and filter to EUR samples for current session."""
    try:
        logger.info("Starting EUR data loading process...")
        
        # Step 1: Load ancestry data
        logger.debug("Step 1: Loading ancestry predictions...")
        start_time = time.time()
        
        ancestry_ht = hl.import_table(
            config['inputs']['ancestry_pred'],
            impute=True,
            types={'research_id': hl.tstr, 'ancestry_pred': hl.tstr}
        )
        
        ancestry_time = time.time() - start_time
        logger.info(f"Ancestry data loaded in {ancestry_time:.1f} seconds")
        
        # Step 2: Filter to European ancestry samples
        logger.debug("Step 2: Filtering to European ancestry samples...")
        start_time = time.time()
        
        eur_samples = ancestry_ht.filter(ancestry_ht.ancestry_pred == 'eur')
        eur_sample_ids = eur_samples.aggregate(hl.agg.collect(eur_samples.research_id))
        
        filter_time = time.time() - start_time
        logger.info(f"Found {len(eur_sample_ids):,} European ancestry samples in {filter_time:.1f} seconds")
        
        # Step 3: Load MatrixTable
        logger.debug("Step 3: Loading ACAF MatrixTable...")
        start_time = time.time()
        
        mt_path = config['inputs']['wgs_matrix_table']
        logger.debug(f"MatrixTable path: {mt_path}")
        
        mt = hl.read_matrix_table(mt_path)
        
        load_time = time.time() - start_time
        variant_count = mt.count_rows()
        sample_count = mt.count_cols()
        logger.info(f"MatrixTable loaded in {load_time:.1f} seconds: {variant_count:,} variants, {sample_count:,} samples")
        
        # Step 4: Filter to EUR samples
        logger.debug("Step 4: Filtering MatrixTable to EUR samples...")
        start_time = time.time()
        
        mt_eur = mt.filter_cols(hl.literal(eur_sample_ids).contains(mt.s))
        
        eur_filter_time = time.time() - start_time
        eur_sample_count = mt_eur.count_cols()
        logger.info(f"EUR filtering completed in {eur_filter_time:.1f} seconds: {eur_sample_count:,} EUR samples")
        
        # Cleanup intermediate objects
        del ancestry_ht, eur_samples, mt
        
        total_time = ancestry_time + filter_time + load_time + eur_filter_time
        logger.info(f"EUR data loading completed in {total_time:.1f} seconds total")
        
        return mt_eur, eur_sample_ids
        
    except Exception as e:
        logger.error(f"Failed to load and filter EUR data: {e}")
        logger.debug(f"EUR data loading error details: {type(e).__name__}: {e}")
        raise


def process_chromosome_with_session(mt_eur, chrom, intervals, target_variants_chr, logger):
    """Process chromosome in current session with immediate export and progress tracking."""
    
    logger.info(f"Chr {chrom}: Processing {len(intervals)} intervals...")
    chrom_start_time = time.time()
    
    processed_intervals = []
    total_filtered_variants_chr = 0
    
    for interval_idx, interval in enumerate(intervals):
        interval_start_time = time.time()
        logger.info(f"Chr {chrom}: Processing interval {interval_idx + 1}/{len(intervals)} (pos: {interval.start.position}-{interval.end.position})")
        
        try:
            # Load only the specific interval
            logger.debug(f"Chr {chrom} interval {interval_idx + 1}: Loading interval data...")
            mt_interval = mt_eur.filter_rows(
                (mt_eur.locus.contig == chrom) & 
                (mt_eur.locus.position >= interval.start.position) & 
                (mt_eur.locus.position < interval.end.position)
            )
            
            # Quick count
            logger.debug(f"Chr {chrom} interval {interval_idx + 1}: Counting variants...")
            interval_variants = mt_interval.count_rows()
            
            if interval_variants == 0:
                interval_time = time.time() - interval_start_time
                logger.info(f"Chr {chrom} interval {interval_idx + 1}: No variants ({interval_time:.1f}s)")
                continue
            
            # Count filtered variants (already EUR filtered)
            filtered_variants = mt_interval.count_rows()
            total_filtered_variants_chr += filtered_variants
            
            interval_time = time.time() - interval_start_time
            logger.info(f"Chr {chrom} interval {interval_idx + 1}: {filtered_variants:,} filtered variants ({interval_time:.1f}s)")
            
            if filtered_variants > 0:
                processed_intervals.append(mt_interval)
            
            # Clear memory for this interval
            del mt_interval
            
            # Small delay between intervals
            time.sleep(1)
            
        except Exception as e:
            logger.error(f"Chr {chrom} interval {interval_idx + 1} failed: {e}")
            continue
    
    chrom_time = time.time() - chrom_start_time
    logger.info(f"Chr {chrom}: Total {total_filtered_variants_chr:,} filtered variants from {len(processed_intervals)} intervals ({chrom_time:.1f}s)")
    
    # Combine all intervals for this chromosome
    if processed_intervals:
        logger.info(f"Chr {chrom}: Combining {len(processed_intervals)} intervals...")
        combine_start_time = time.time()
        
        try:
            mt_chrom_combined = hl.MatrixTable.union_rows(*processed_intervals)
            
            combine_time = time.time() - combine_start_time
            logger.info(f"Chr {chrom}: Intervals combined in {combine_time:.1f}s")
            
            # Calculate final sampling for this chromosome
            if total_filtered_variants_chr > 0:
                sampling_fraction = target_variants_chr / total_filtered_variants_chr
                logger.info(f"Chr {chrom}: Final sampling fraction = {sampling_fraction:.6f}")
                
                # Sample variants from combined chromosome data
                logger.debug(f"Chr {chrom}: Sampling variants...")
                sample_start_time = time.time()
                
                mt_sampled = mt_chrom_combined.sample_rows(sampling_fraction, seed=42)
                actual_variants = mt_sampled.count_rows()
                
                sample_time = time.time() - sample_start_time
                logger.info(f"Chr {chrom}: Sampled {actual_variants:,} variants in {sample_time:.1f}s")
                
                # Clear memory aggressively
                del mt_chrom_combined, processed_intervals
                
                return mt_sampled
            else:
                logger.warning(f"Chr {chrom}: No filtered variants for sampling")
                return None
                
        except Exception as e:
            logger.error(f"Chr {chrom}: Failed to combine intervals: {e}")
            return None
    else:
        logger.warning(f"Chr {chrom}: No intervals with data")
        return None


def process_all_chromosomes_with_sessions(config, logger):
    """Process all chromosomes with fresh Hail session per chromosome."""
    
    # Calculate targets and intervals per chromosome
    logger.info("Calculating chromosome targets and intervals...")
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
    
    # Check for test mode
    test_mode = config['params'].get('test_mode', False)
    if test_mode:
        test_chromosome = config['params'].get('test_chromosome', 'chr22')
        chromosome_order = [test_chromosome]
        logger.info(f"TEST MODE: Processing only {test_chromosome}")
    else:
        # Sort chromosomes by size (smallest first) for stability
        chromosome_order = sorted(targets_per_chr.keys(), 
                                key=lambda x: targets_per_chr[x], 
                                reverse=False)
    
    logger.info(f"Processing order: {chromosome_order}")
    logger.info("Per-chromosome session processing with EUR reloading!")
    
    # Process each chromosome with fresh session
    processed_chromosomes = []
    failed_chromosomes = {}
    pipeline_start_time = time.time()
    
    for chrom_idx, chrom in enumerate(chromosome_order):
        intervals = intervals_per_chr[chrom]
        chrom_start_time = time.time()
        logger.info(f"="*60)
        logger.info(f"Starting chromosome {chrom} ({len(intervals)} intervals) - {chrom_idx + 1}/{len(chromosome_order)}")
        logger.info(f"="*60)
        
        try:
            # Initialize fresh session for this chromosome
            logger.info(f"Chr {chrom}: Initializing session...")
            session_start_time = time.time()
            
            # Use longer timeout for test mode
            timeout_seconds = 900 if test_mode else 600
            logger.info(f"Chr {chrom}: Using {timeout_seconds}s timeout for session initialization")
            
            if not initialize_fresh_session(logger, timeout_seconds=timeout_seconds):
                raise Exception("Failed to initialize Hail session")
            
            session_time = time.time() - session_start_time
            logger.info(f"Chr {chrom}: Session initialized in {session_time:.1f}s")
            
            # Load and filter EUR data for this session
            logger.info(f"Chr {chrom}: Loading EUR data...")
            data_start_time = time.time()
            
            mt_eur, eur_sample_ids = load_and_filter_eur_data(config, logger)
            
            data_time = time.time() - data_start_time
            logger.info(f"Chr {chrom}: EUR data loaded in {data_time:.1f}s")
            
            # Process chromosome in current session
            logger.info(f"Chr {chrom}: Starting interval processing...")
            process_start_time = time.time()
            
            mt_chrom_sampled = process_chromosome_with_session(
                mt_eur, chrom, intervals, targets_per_chr[chrom], logger
            )
            
            process_time = time.time() - process_start_time
            logger.info(f"Chr {chrom}: Interval processing completed in {process_time:.1f}s")
            
            if mt_chrom_sampled is not None:
                processed_chromosomes.append(chrom)
                
                # Export immediately to WORKSPACE_BUCKET
                logger.info(f"Chr {chrom}: Exporting results...")
                export_start_time = time.time()
                
                chromosome_output_path = f"{output_dir}/chr{chrom}_background_snps.mt"
                mt_chrom_sampled.write(chromosome_output_path, overwrite=True)
                
                export_time = time.time() - export_start_time
                logger.info(f"Chr {chrom}: Exported to {chromosome_output_path} in {export_time:.1f}s")
                
                # Create summary for this chromosome
                chr_summary = {
                    'chromosome': chrom,
                    'target_variants': targets_per_chr[chrom],
                    'actual_variants': mt_chrom_sampled.count_rows(),
                    'output_path': chromosome_output_path,
                    'processing_time_seconds': round(process_time, 1),
                    'session_time_seconds': round(session_time, 1),
                    'data_load_time_seconds': round(data_time, 1),
                    'export_time_seconds': round(export_time, 1),
                    'timestamp': datetime.now().isoformat()
                }
                
                # Save chromosome summary
                chr_summary_path = f"{output_dir}/chr{chrom}_summary.json"
                with open(chr_summary_path, 'w') as f:
                    json.dump(chr_summary, f, indent=2)
                
                # Clear memory aggressively
                del mt_chrom_sampled, mt_eur
                
                chrom_total_time = time.time() - chrom_start_time
                logger.info(f"Chr {chrom}: COMPLETED in {chrom_total_time:.1f}s total")
                logger.info(f"Chr {chrom}: Memory cleared and exported")
                
                # Force garbage collection
                import gc
                gc.collect()
                
            else:
                logger.warning(f"Chr {chrom}: No data processed")
                failed_chromosomes[chrom] = "No data processed"
        
        except Exception as e:
            chrom_total_time = time.time() - chrom_start_time
            logger.error(f"Chr {chrom}: FAILED after {chrom_total_time:.1f}s - {e}")
            logger.debug(f"Chr {chrom} error details: {type(e).__name__}: {e}")
            failed_chromosomes[chrom] = str(e)
            # Continue to next chromosome instead of stopping
            continue
        
        finally:
            # Always stop session after each chromosome
            try:
                logger.info(f"Chr {chrom}: Stopping Hail session...")
                hl.stop()
                logger.debug(f"Chr {chrom}: Hail session stopped")
            except Exception as e:
                logger.debug(f"Chr {chrom}: Error stopping session: {e}")
        
        # Delay between chromosomes to let resources stabilize
        logger.info(f"Chr {chrom}: Waiting 3 seconds before next chromosome...")
        time.sleep(3)
    
    pipeline_time = time.time() - pipeline_start_time
    logger.info(f"="*60)
    logger.info(f"PIPELINE COMPLETED in {pipeline_time:.1f}s")
    logger.info(f"="*60)
    logger.info(f"Successfully processed {len(processed_chromosomes)}/{len(chromosome_order)} chromosomes: {processed_chromosomes}")
    
    if failed_chromosomes:
        logger.info(f"Failed chromosomes: {list(failed_chromosomes.keys())}")
        for chrom, error in failed_chromosomes.items():
            logger.error(f"  {chrom}: {error}")
    
    return processed_chromosomes, output_dir


def sample_background_snps(config, logger):
    """Main function to sample 1M random common SNPs from EUR ancestry samples."""
    logger.info("Starting background SNP sampling...")
    pipeline_start_time = time.time()
    
    try:
        # Process chromosomes with per-chromosome sessions
        processed_chromosomes, output_dir = process_all_chromosomes_with_sessions(config, logger)
        
        # Generate summary report
        pipeline_time = time.time() - pipeline_start_time
        summary = {
            'timestamp': datetime.now().isoformat(),
            'target_variants': config['params']['target_variants'],
            'ancestry': 'EUR',
            'processing_method': 'per_chromosome_sessions',
            'chromosomes_processed': processed_chromosomes,
            'chromosome_count': len(processed_chromosomes),
            'output_directory': output_dir,
            'individual_chromosome_mts': [f"chr{chrom}_background_snps.mt" for chrom in processed_chromosomes],
            'pipeline_time_seconds': round(pipeline_time, 1),
            'next_step': 'Run merge_chromosome_mts.py to combine results'
        }
        
        # Save summary report to WORKSPACE_BUCKET
        summary_path = f"{output_dir}/chromosome_processing_summary.json"
        with open(summary_path, 'w') as f:
            json.dump(summary, f, indent=2)
        
        logger.info(f"="*60)
        logger.info(f"PIPELINE SUMMARY")
        logger.info(f"="*60)
        logger.info(f"Chromosome processing completed in {pipeline_time:.1f}s")
        logger.info(f"Summary saved to {summary_path}")
        logger.info(f"Individual chromosome MTs saved to {output_dir}")
        logger.info(f"Next step: Run 'python python/merge_chromosome_mts.py' to combine results")
        
        return summary
        
    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        logger.debug(f"Analysis error details: {type(e).__name__}: {e}")
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