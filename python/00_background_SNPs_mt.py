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


def sample_background_snps(config, logger):
    """Main function to sample 1M random common SNPs from EUR ancestry samples."""
    logger.info("Starting background SNP sampling...")
    
    # Initialize Hail
    init_hail(log_prefix="background_snps", driver_mem="8g", reference="GRCh38")
    
    try:
        # Load ancestry data
        eur_sample_ids = load_ancestry_data(config, logger)
        
        # Load pre-filtered ACAF MatrixTable
        logger.info("Loading ACAF MatrixTable...")
        mt_path = config['inputs']['wgs_matrix_table']
        
        # Read the existing MatrixTable directly
        mt = hl.read_matrix_table(mt_path)
        
        logger.info(f"Loaded MatrixTable with {mt.count_rows()} variants and {mt.count_cols()} samples")
        
        # Filter to EUR ancestry samples
        logger.info("Filtering to European ancestry samples...")
        mt = mt.filter_cols(hl.literal(eur_sample_ids).contains(mt.s))
        
        logger.info(f"After ancestry filtering: {mt.count_rows()} variants, {mt.count_cols()} samples")
        
        # Calculate sampling fraction for 1M variants
        total_variants = mt.count_rows()
        target_variants = config['params']['target_variants']
        sampling_fraction = target_variants / total_variants
        
        logger.info(f"Total variants: {total_variants}, Target: {target_variants}")
        logger.info(f"Sampling fraction: {sampling_fraction:.6f}")
        
        # Random sampling of variants
        logger.info("Sampling random variants...")
        mt_sampled = mt.sample_rows(sampling_fraction, seed=config['params']['random_seed'])
        
        actual_variants = mt_sampled.count_rows()
        logger.info(f"Sampled {actual_variants} variants")
        
        # Quality checks
        logger.info("Performing quality checks...")
        
        # Check allele frequency distribution
        mt_sampled = mt_sampled.annotate_rows(
            af = hl.agg.mean(mt_sampled.GT.n_alt_alleles()) / (2 * mt_sampled.count_cols())
        )
        
        af_stats = mt_sampled.aggregate_rows(hl.struct(
            mean_af = hl.agg.mean(mt_sampled.af),
            median_af = hl.agg.approx_median(mt_sampled.af),
            min_af = hl.agg.min(mt_sampled.af),
            max_af = hl.agg.max(mt_sampled.af)
        ))
        
        logger.info(f"AF stats - Mean: {af_stats.mean_af:.4f}, Median: {af_stats.median_af:.4f}")
        logger.info(f"AF range - Min: {af_stats.min_af:.4f}, Max: {af_stats.max_af:.4f}")
        
        # Export results
        logger.info("Exporting results...")
        
        # Create output directory
        output_dir = os.path.join(
            config['outputs']['base_dir'],
            config['outputs']['data_dir_suffix'],
            "background_snps"
        )
        os.makedirs(output_dir, exist_ok=True)
        
        # Export filtered MatrixTable
        mt_output_path = os.path.join(output_dir, "background_1M_snps.mt")
        mt_sampled.write(mt_output_path, overwrite=True)
        logger.info(f"Filtered MatrixTable exported to {mt_output_path}")
        
        # Also export as Hail Table for easier querying
        ht_output_path = os.path.join(output_dir, "background_1M_snps.ht")
        mt_sampled.rows().write(ht_output_path, overwrite=True)
        logger.info(f"Variant table exported to {ht_output_path}")
        
        # Generate summary report
        summary = {
            'timestamp': datetime.now().isoformat(),
            'total_samples': mt_sampled.count_cols(),
            'total_variants': actual_variants,
            'target_variants': target_variants,
            'sampling_fraction': sampling_fraction,
            'ancestry': 'EUR',
            'af_stats': {
                'mean': af_stats.mean_af,
                'median': af_stats.median_af,
                'min': af_stats.min_af,
                'max': af_stats.max_af
            },
            'outputs': {
                'matrix_table': mt_output_path,
                'variant_table': ht_output_path
            }
        }
        
        # Save summary report
        summary_path = os.path.join(output_dir, "sampling_summary.json")
        import json
        with open(summary_path, 'w') as f:
            json.dump(summary, f, indent=2)
        
        logger.info(f"Summary report saved to {summary_path}")
        logger.info("Background SNP sampling completed successfully!")
        
        return summary
        
    except Exception as e:
        logger.error(f"Error during background SNP sampling: {str(e)}")
        raise
    finally:
        # Stop Hail session
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