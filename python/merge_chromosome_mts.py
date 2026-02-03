#!/usr/bin/env python3
"""
Merge individual chromosome MatrixTables into final background SNP dataset.

This script loads all individual chromosome MTs created by 00_background_SNPs_mt.py
and merges them into a single MatrixTable with quality checks and final exports.
"""

import os
import json
import logging
from datetime import datetime
import hail as hl

# Import utility functions
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from utils import load_config, setup_logger


def load_chromosome_summary(output_dir, logger):
    """Load the chromosome processing summary to get list of processed chromosomes."""
    summary_path = os.path.join(output_dir, "chromosome_processing_summary.json")
    
    if not os.path.exists(summary_path):
        raise FileNotFoundError(f"Chromosome processing summary not found: {summary_path}")
    
    with open(summary_path, 'r') as f:
        summary = json.load(f)
    
    logger.info(f"Loaded chromosome processing summary from {summary_path}")
    logger.info(f"Processed chromosomes: {summary['chromosomes_processed']}")
    
    return summary


def merge_chromosome_mts(output_dir, processed_chromosomes, logger):
    """Merge all individual chromosome MatrixTables into final dataset."""
    
    logger.info("Starting chromosome MatrixTable merge...")
    
    # Load all chromosome MTs
    chromosome_mts = []
    missing_chromosomes = []
    
    for chrom in processed_chromosomes:
        mt_path = os.path.join(output_dir, f"chr_{chrom}_sampled.mt")
        
        if os.path.exists(mt_path):
            logger.info(f"Loading {chrom} from {mt_path}")
            mt_chrom = hl.read_matrix_table(mt_path)
            
            # Quick count
            variants = mt_chrom.count_rows()
            samples = mt_chrom.count_cols()
            logger.info(f"  {chrom}: {variants:,} variants, {samples:,} samples")
            
            chromosome_mts.append(mt_chrom)
        else:
            logger.warning(f"Chromosome MT not found: {mt_path}")
            missing_chromosomes.append(chrom)
    
    if not chromosome_mts:
        raise ValueError("No chromosome MatrixTables found to merge")
    
    if missing_chromosomes:
        logger.warning(f"Missing chromosomes: {missing_chromosomes}")
    
    logger.info(f"Loaded {len(chromosome_mts)} chromosome MatrixTables")
    
    # Merge all chromosome MatrixTables
    logger.info("Merging all chromosome MatrixTables...")
    final_mt = hl.MatrixTable.union_rows(*chromosome_mts)
    
    # Final counts
    total_variants = final_mt.count_rows()
    total_samples = final_mt.count_cols()
    
    logger.info(f"Merged dataset: {total_variants:,} variants, {total_samples:,} samples")
    
    return final_mt


def perform_quality_checks(final_mt, logger):
    """Perform quality checks on the final merged dataset."""
    
    logger.info("Performing quality checks on merged dataset...")
    
    # Check allele frequency distribution
    final_mt = final_mt.annotate_rows(
        af = hl.agg.mean(final_mt.GT.n_alt_alleles()) / (2 * final_mt.count_cols())
    )
    
    af_stats = final_mt.aggregate_rows(hl.struct(
        mean_af = hl.agg.mean(final_mt.af),
        median_af = hl.agg.approx_median(final_mt.af),
        min_af = hl.agg.min(final_mt.af),
        max_af = hl.agg.max(final_mt.af)
    ))
    
    logger.info(f"AF stats - Mean: {af_stats.mean_af:.4f}, Median: {af_stats.median_af:.4f}")
    logger.info(f"AF range - Min: {af_stats.min_af:.4f}, Max: {af_stats.max_af:.4f}")
    
    # Check chromosome distribution
    chr_counts = final_mt.aggregate_rows(hl.agg.counter(final_mt.locus.contig))
    logger.info(f"Chromosome distribution: {dict(chr_counts)}")
    
    return af_stats, chr_counts


def export_final_results(final_mt, output_dir, af_stats, chr_counts, logger):
    """Export final merged MatrixTable and summary."""
    
    logger.info("Exporting final results...")
    
    # Export final MatrixTable
    mt_output_path = os.path.join(output_dir, "background_1M_snps.mt")
    final_mt.write(mt_output_path, overwrite=True)
    logger.info(f"Final MatrixTable exported to {mt_output_path}")
    
    # Export as Hail Table for easier querying
    ht_output_path = os.path.join(output_dir, "background_1M_snps.ht")
    final_mt.rows().write(ht_output_path, overwrite=True)
    logger.info(f"Variant table exported to {ht_output_path}")
    
    # Generate final summary report
    final_summary = {
        'timestamp': datetime.now().isoformat(),
        'merge_completed': True,
        'total_samples': final_mt.count_cols(),
        'total_variants': final_mt.count_rows(),
        'target_variants': 1000000,
        'ancestry': 'EUR',
        'processing_method': 'single_session_individual_exports + merge',
        'af_stats': {
            'mean': af_stats.mean_af,
            'median': af_stats.median_af,
            'min': af_stats.min_af,
            'max': af_stats.max_af
        },
        'chromosome_distribution': dict(chr_counts),
        'outputs': {
            'final_matrix_table': mt_output_path,
            'variant_table': ht_output_path
        }
    }
    
    # Save final summary
    final_summary_path = os.path.join(output_dir, "final_merge_summary.json")
    with open(final_summary_path, 'w') as f:
        json.dump(final_summary, f, indent=2)
    
    logger.info(f"Final summary saved to {final_summary_path}")
    
    return final_summary


def main():
    """Main execution function."""
    # Setup logging
    logger = setup_logger("merge_chromosome_mts")
    
    try:
        # Load configuration
        config = load_config()
        
        # Define output directory
        output_dir = os.path.join(
            config['outputs']['base_dir'],
            config['outputs']['data_dir_suffix'],
            "background_snps"
        )
        
        logger.info(f"Starting chromosome MatrixTable merge...")
        logger.info(f"Output directory: {output_dir}")
        
        # Initialize Hail with more memory for merge operation
        logger.info("Initializing Hail for merge operation...")
        hl.init(
            log='/tmp/hail_merge.log',
            spark_conf={
                'spark.driver.memory': '20g',  # More memory for merge
                'spark.executor.memory': '20g',
                'spark.network.timeout': '1200s',
                'spark.executor.heartbeatInterval': '120s'
            }
        )
        hl.default_reference('GRCh38')
        
        # Load chromosome processing summary
        chromosome_summary = load_chromosome_summary(output_dir, logger)
        processed_chromosomes = chromosome_summary['chromosomes_processed']
        
        # Merge all chromosome MatrixTables
        final_mt = merge_chromosome_mts(output_dir, processed_chromosomes, logger)
        
        # Perform quality checks
        af_stats, chr_counts = perform_quality_checks(final_mt, logger)
        
        # Export final results
        final_summary = export_final_results(final_mt, output_dir, af_stats, chr_counts, logger)
        
        logger.info("Chromosome merge completed successfully!")
        logger.info(f"Final dataset: {final_summary['total_variants']:,} variants, {final_summary['total_samples']:,} samples")
        logger.info(f"Final outputs in: {output_dir}")
        
        return final_summary
        
    except Exception as e:
        logger.error(f"Merge failed: {e}")
        raise
    finally:
        # Always stop Hail session
        logger.info("Stopping Hail session...")
        hl.stop()


if __name__ == "__main__":
    main()
