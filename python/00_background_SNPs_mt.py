#!/usr/bin/env python3
import os
import sys
import random
import hail as hl
import time
import pandas as pd
import json
import gc
from datetime import datetime
from utils import load_config, setup_logger


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


def initialize_hail_session(logger):
    """Initialize Hail session in local mode with optimized settings."""
    logger.info("Initializing Hail session in local mode...")
    
    try:
        spark_conf = {
            'spark.master': 'local[8]',           # Use 8 of 16 available cores
            'spark.driver.memory': '20g',         # Generous memory allocation for local mode
            'spark.driver.cores': '4',           # 4 cores for driver
            'spark.sql.shuffle.partitions': '8',  # Match core count for optimal performance
            'spark.network.timeout': '2400s',     # Keep generous timeout
            'spark.executor.heartbeatInterval': '120s', # Keep heartbeat settings
            'spark.driver.extraJavaOptions': '-XX:+UseG1GC -XX:MaxGCPauseMillis=200 -XX:+UnlockExperimentalVMOptions -XX:+UseCGroupMemoryLimitForHeap',
        }
        
        logger.info(f"Spark configuration: {spark_conf}")
        
        start_time = time.time()
        hl.init(
            log='/tmp/hail_background_snps.log',
            spark_conf=spark_conf,
            quiet=False
        )
        
        init_time = time.time() - start_time
        logger.info(f"Hail initialization completed in {init_time:.1f} seconds")
        
        # Set reference genome
        hl.default_reference('GRCh38')
        logger.info("Reference genome set to GRCh38")
        
        return True
        
    except Exception as e:
        logger.error(f"Failed to initialize Hail session: {e}")
        return False


def load_eur_data(config, logger):
    """Load ancestry data and filter MatrixTable to EUR samples."""
    logger.info("Loading EUR ancestry data...")
    
    try:
        # Load ancestry data
        start_time = time.time()
        ancestry_ht = hl.import_table(
            config['inputs']['ancestry_pred'],
            impute=True,
            types={'research_id': hl.tstr, 'ancestry_pred': hl.tstr}
        )
        
        # Filter to European ancestry samples
        eur_samples = ancestry_ht.filter(ancestry_ht.ancestry_pred == 'eur')
        eur_sample_ids = eur_samples.aggregate(hl.agg.collect(eur_samples.research_id))
        
        ancestry_time = time.time() - start_time
        logger.info(f"Found {len(eur_sample_ids):,} European ancestry samples in {ancestry_time:.1f} seconds")
        
        # Load MatrixTable
        logger.info("Loading MatrixTable...")
        start_time = time.time()
        
        mt_path = config['inputs']['wgs_matrix_table']
        mt = hl.read_matrix_table(mt_path)
        
        load_time = time.time() - start_time
        variant_count = mt.count_rows()
        sample_count = mt.count_cols()
        logger.info(f"MatrixTable loaded: {variant_count:,} variants, {sample_count:,} samples ({load_time:.1f}s)")
        
        # Filter to EUR samples
        logger.info("Filtering to EUR samples...")
        start_time = time.time()
        
        mt_eur = mt.filter_cols(hl.literal(eur_sample_ids).contains(mt.s))
        
        eur_filter_time = time.time() - start_time
        eur_sample_count = mt_eur.count_cols()
        logger.info(f"EUR filtering completed: {eur_sample_count:,} EUR samples ({eur_filter_time:.1f}s)")
        
        # Cleanup
        del ancestry_ht, eur_samples, mt
        
        total_time = ancestry_time + load_time + eur_filter_time
        logger.info(f"EUR data loading completed in {total_time:.1f} seconds total")
        
        return mt_eur, eur_sample_ids
        
    except Exception as e:
        logger.error(f"Failed to load EUR data: {e}")
        raise


def process_chromosome(mt_eur, chrom, intervals, target_variants_chr, output_dir, logger):
    """Process a single chromosome with its intervals."""
    logger.info(f"Processing {chrom} ({len(intervals)} intervals)...")
    
    chrom_start_time = time.time()
    processed_intervals = []
    total_filtered_variants_chr = 0
    
    for interval_idx, interval in enumerate(intervals):
        interval_start_time = time.time()
        logger.info(f"  {chrom} interval {interval_idx + 1}/{len(intervals)} (pos: {interval.start.position}-{interval.end.position})")
        
        try:
            # Filter to chromosome and interval
            mt_interval = mt_eur.filter_rows(
                (mt_eur.locus.contig == chrom) & 
                (mt_eur.locus.position >= interval.start.position) & 
                (mt_eur.locus.position < interval.end.position)
            )
            
            interval_variants = mt_interval.count_rows()
            
            if interval_variants == 0:
                logger.info(f"    No variants ({time.time() - interval_start_time:.1f}s)")
                continue
            
            total_filtered_variants_chr += interval_variants
            processed_intervals.append(mt_interval)
            
            logger.info(f"    {interval_variants:,} variants ({time.time() - interval_start_time:.1f}s)")
            
            # Cleanup interval
            del mt_interval
            
        except Exception as e:
            logger.error(f"    Interval {interval_idx + 1} failed: {e}")
            continue
    
    chrom_time = time.time() - chrom_start_time
    logger.info(f"  {chrom}: {total_filtered_variants_chr:,} total variants from {len(processed_intervals)} intervals ({chrom_time:.1f}s)")
    
    # Combine intervals and sample
    if processed_intervals and total_filtered_variants_chr > 0:
        logger.info(f"  {chrom}: Combining intervals and sampling...")
        
        try:
            # Combine intervals
            mt_chrom_combined = hl.MatrixTable.union_rows(*processed_intervals)
            
            # Calculate sampling fraction
            sampling_fraction = target_variants_chr / total_filtered_variants_chr
            logger.info(f"  {chrom}: Sampling fraction = {sampling_fraction:.6f}")
            
            # Sample variants
            mt_sampled = mt_chrom_combined.sample_rows(sampling_fraction, seed=42)
            actual_variants = mt_sampled.count_rows()
            
            logger.info(f"  {chrom}: Sampled {actual_variants:,} variants")
            
            # Export to GCS
            chromosome_output_path = f"{output_dir}/chr{chrom}_background_snps.mt"
            logger.info(f"  {chrom}: Exporting to {chromosome_output_path}...")
            
            export_start_time = time.time()
            mt_sampled.write(chromosome_output_path, overwrite=True)
            export_time = time.time() - export_start_time
            
            logger.info(f"  {chrom}: Export completed in {export_time:.1f}s")
            
            # Create chromosome summary
            chr_summary = {
                'chromosome': chrom,
                'target_variants': target_variants_chr,
                'actual_variants': actual_variants,
                'output_path': chromosome_output_path,
                'processing_time_seconds': round(chrom_time, 1),
                'export_time_seconds': round(export_time, 1),
                'timestamp': datetime.now().isoformat()
            }
            
            # Save chromosome summary
            chr_summary_path = f"{output_dir}/chr{chrom}_summary.json"
            with open(chr_summary_path, 'w') as f:
                json.dump(chr_summary, f, indent=2)
            
            # Cleanup
            del mt_chrom_combined, mt_sampled, processed_intervals
            gc.collect()
            
            return chr_summary
            
        except Exception as e:
            logger.error(f"  {chrom}: Failed to process: {e}")
            return None
    else:
        logger.warning(f"  {chrom}: No data to process")
        return None


def main():
    """Main execution function - single session local mode pipeline."""
    # Setup logging
    logger = setup_logger("background_snps")
    
    try:
        # Load configuration
        config = load_config()
        logger.info("Configuration loaded successfully")
        
        # Check for test mode
        test_mode = config['params'].get('test_mode', False)
        if test_mode:
            test_chromosome = config['params'].get('test_chromosome', 'chr22')
            logger.info(f"TEST MODE: Processing only {test_chromosome}")
        
        # Initialize Hail session
        if not initialize_hail_session(logger):
            raise Exception("Failed to initialize Hail session")
        
        # Load EUR data once
        mt_eur, eur_sample_ids = load_eur_data(config, logger)
        
        # Calculate chromosome targets and intervals
        logger.info("Calculating chromosome targets and intervals...")
        targets_per_chr, intervals_per_chr = calculate_chromosome_intervals_and_targets(
            config['params']['target_variants']
        )
        
        logger.info(f"Chromosome targets: {targets_per_chr}")
        
        # Create output directory
        workspace_bucket = os.environ.get('WORKSPACE_BUCKET', 'gs://default-bucket')
        output_dir = f"{workspace_bucket}/results/FNCV_RVAS_MS/background_snps"
        logger.info(f"Output directory: {output_dir}")
        
        # Determine chromosome order
        if test_mode:
            chromosome_order = [test_chromosome]
        else:
            # Sort by size (smallest first) for stability
            chromosome_order = sorted(targets_per_chr.keys(), 
                                    key=lambda x: targets_per_chr[x], 
                                    reverse=False)
        
        logger.info(f"Processing order: {chromosome_order}")
        
        # Process chromosomes
        processed_chromosomes = []
        failed_chromosomes = {}
        pipeline_start_time = time.time()
        
        for chrom_idx, chrom in enumerate(chromosome_order):
            logger.info(f"="*60)
            logger.info(f"Processing chromosome {chrom} - {chrom_idx + 1}/{len(chromosome_order)}")
            logger.info(f"="*60)
            
            chrom_start_time = time.time()
            
            try:
                intervals = intervals_per_chr[chrom]
                target_variants_chr = targets_per_chr[chrom]
                
                # Process chromosome
                chrom_summary = process_chromosome(
                    mt_eur, chrom, intervals, target_variants_chr, output_dir, logger
                )
                
                if chrom_summary:
                    processed_chromosomes.append(chrom)
                    logger.info(f"{chrom}: COMPLETED in {time.time() - chrom_start_time:.1f}s")
                else:
                    failed_chromosomes[chrom] = "Processing failed"
                    logger.error(f"{chrom}: FAILED")
                
                # Memory cleanup between chromosomes
                logger.info(f"{chrom}: Cleaning up memory...")
                gc.collect()
                time.sleep(3)  # Brief pause for memory to settle
                
            except Exception as e:
                failed_chromosomes[chrom] = str(e)
                logger.error(f"{chrom}: FAILED - {e}")
        
        # Stop Hail session
        logger.info("Stopping Hail session...")
        try:
            hl.stop()
        except:
            pass
        
        # Generate final summary
        pipeline_time = time.time() - pipeline_start_time
        summary = {
            'timestamp': datetime.now().isoformat(),
            'target_variants': config['params']['target_variants'],
            'ancestry': 'EUR',
            'processing_mode': 'local_single_session',
            'chromosomes_processed': processed_chromosomes,
            'chromosome_count': len(processed_chromosomes),
            'failed_chromosomes': failed_chromosomes,
            'output_directory': output_dir,
            'individual_chromosome_mts': [f"chr{chrom}_background_snps.mt" for chrom in processed_chromosomes],
            'pipeline_time_seconds': round(pipeline_time, 1),
            'next_step': 'Run merge_chromosome_mts.py to combine results'
        }
        
        # Save summary
        summary_path = f"{output_dir}/chromosome_processing_summary.json"
        with open(summary_path, 'w') as f:
            json.dump(summary, f, indent=2)
        
        # Final report
        logger.info(f"="*60)
        logger.info(f"PIPELINE COMPLETED in {pipeline_time:.1f}s")
        logger.info(f"="*60)
        logger.info(f"Successfully processed: {len(processed_chromosomes)}/{len(chromosome_order)} chromosomes")
        logger.info(f"Output directory: {output_dir}")
        logger.info(f"Summary saved: {summary_path}")
        
        if failed_chromosomes:
            logger.warning(f"Failed chromosomes: {list(failed_chromosomes.keys())}")
            for chrom, error in failed_chromosomes.items():
                logger.error(f"  {chrom}: {error}")
        
        logger.info(f"Next step: Run 'python python/merge_chromosome_mts.py' to combine results")
        
        return summary
        
    except Exception as e:
        logger.error(f"Pipeline failed: {e}")
        
        # Ensure Hail session is stopped on error
        try:
            hl.stop()
        except:
            pass
        
        sys.exit(1)


if __name__ == "__main__":
    main()