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


def calculate_chromosome_sample_targets(target_variants=1000000):
    """Calculate variant targets per chromosome based on GRCh38 lengths."""
    
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
    
    for chrom, length in chr_lengths.items():
        proportion = length / total_autosome_length
        targets_per_chr[chrom] = int(target_variants * proportion)
    
    return targets_per_chr


def process_chromosome(mt, chrom, eur_sample_ids, target_variants_chr, logger):
    """Process a single chromosome: filter → sample → save."""
    
    logger.info(f"Processing chromosome {chrom}...")
    
    # Step 1: Filter to chromosome
    mt_chrom = mt.filter_rows(mt.locus.contig == chrom)
    total_variants = mt_chrom.count_rows()
    logger.info(f"Chr {chrom}: {total_variants:,} total variants")
    
    # Step 2: Filter to EUR samples
    mt_chrom_eur = mt_chrom.filter_cols(
        hl.literal(eur_sample_ids).contains(mt_chrom.s)
    )
    filtered_variants = mt_chrom_eur.count_rows()
    filtered_samples = mt_chrom_eur.count_cols()
    logger.info(f"Chr {chrom}: {filtered_variants:,} variants, {filtered_samples:,} samples")
    
    # Step 3: Calculate sampling fraction
    if filtered_variants > 0:
        sampling_fraction = target_variants_chr / filtered_variants
        logger.info(f"Chr {chrom}: sampling fraction = {sampling_fraction:.6f}")
        
        # Step 4: Sample variants
        mt_sampled = mt_chrom_eur.sample_rows(sampling_fraction, seed=42)
        actual_variants = mt_sampled.count_rows()
        logger.info(f"Chr {chrom}: sampled {actual_variants:,} variants")
        
        return mt_sampled
    else:
        logger.warning(f"Chr {chrom}: No variants after filtering")
        return None


def process_all_chromosomes(mt, eur_sample_ids, config, logger):
    """Process all chromosomes sequentially."""
    
    # Calculate targets per chromosome
    targets_per_chr = calculate_chromosome_sample_targets(
        config['params']['target_variants']
    )
    
    logger.info(f"Chromosome targets: {targets_per_chr}")
    
    processed_chromosomes = []
    temp_dir = os.path.join(
        config['outputs']['base_dir'],
        config['outputs']['data_dir_suffix'],
        "temp"
    )
    os.makedirs(temp_dir, exist_ok=True)
    
    # Process each chromosome
    for chrom in range(1, 23):
        chrom_str = f"chr{chrom}"
        
        try:
            mt_chrom_sampled = process_chromosome(
                mt, chrom_str, eur_sample_ids, 
                targets_per_chr[chrom_str], logger
            )
            
            if mt_chrom_sampled is not None:
                processed_chromosomes.append(mt_chrom_sampled)
                
                # Checkpoint after each chromosome
                checkpoint_path = os.path.join(temp_dir, f"chr_{chrom_str}_sampled.mt")
                mt_chrom_sampled.write(checkpoint_path, overwrite=True)
                logger.info(f"Checkpoint saved: {checkpoint_path}")
        
        except Exception as e:
            logger.error(f"Error processing chromosome {chrom}: {e}")
            continue
    
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
    
    # Initialize Hail with increased memory
    hl.init(
        log='/tmp/hail_background_snps.log',
        spark_conf={
            'spark.driver.memory': '10g',
            'spark.executor.memory': '10g',
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
        
        # Process chromosomes sequentially
        processed_chromosomes = process_all_chromosomes(mt, eur_sample_ids, config, logger)
        
        # Combine all chromosome results
        final_mt = combine_chromosome_results(processed_chromosomes, config, logger)
        
        # Quality checks
        logger.info("Performing quality checks...")
        
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
        final_mt.write(mt_output_path, overwrite=True)
        logger.info(f"Filtered MatrixTable exported to {mt_output_path}")
        
        # Also export as Hail Table for easier querying
        ht_output_path = os.path.join(output_dir, "background_1M_snps.ht")
        final_mt.rows().write(ht_output_path, overwrite=True)
        logger.info(f"Variant table exported to {ht_output_path}")
        
        # Generate summary report
        summary = {
            'timestamp': datetime.now().isoformat(),
            'total_samples': final_mt.count_cols(),
            'total_variants': final_mt.count_rows(),
            'target_variants': config['params']['target_variants'],
            'ancestry': 'EUR',
            'processing_method': 'chromosome_sequential',
            'chromosomes_processed': len(processed_chromosomes),
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
        
        # Cleanup temporary files
        temp_dir = os.path.join(
            config['outputs']['base_dir'],
            config['outputs']['data_dir_suffix'],
            "temp"
        )
        if os.path.exists(temp_dir):
            import shutil
            shutil.rmtree(temp_dir)
            logger.info("Temporary files cleaned up")
        
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