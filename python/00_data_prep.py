#!/usr/bin/env python3
import os
import random
import hail as hl
import pandas as pd
from utils import load_config, init_hail, setup_logger

logger = setup_logger("data_prep")

def main():
    config = load_config("config/config.json")
    
    # Initialize with high memory
    init_hail("data_prep", driver_mem="16g", reference="GRCh38")

    # Paths
    ANCESTRY_PATH = config['inputs']['ancestry_pred']
    WGS_MT_PATH = config['inputs']['wgs_mt']
    OUT_DIR = config['outputs']['data_dir']
    OUTPUT_PATH = f"{OUT_DIR}/eur_common_snps_500k.mt"
    params = config['params']

    # Samples
    logger.info("Loading EUR sample list...")
    with hl.hadoop_open(ANCESTRY_PATH, 'r') as f:
        ancestry_df = pd.read_csv(f, sep="\t")
    eur_ids = set(ancestry_df[ancestry_df['ancestry_pred'] == 'eur']['research_id'].astype(str))

    # Intervals
    logger.info("Generating random intervals...")
    rg = hl.get_reference('GRCh38')
    lengths = rg.lengths
    chroms = [f"chr{i}" for i in range(1, 23)]
    hl_intervals = []
    
    # Generate random intervals
    for _ in range(params['n_intervals']):
        chrom = random.choice(chroms)
        chrom_len = lengths[chrom]
        start = random.randint(1, chrom_len - params['window_size'])
        # Ensure we don't go off the end
        if start + params['window_size'] > chrom_len:
            start = chrom_len - params['window_size']
        interval = hl.locus_interval(chrom, start, start + params['window_size'], reference_genome='GRCh38')
        hl_intervals.append(interval)

    # 1. Read and Push-Down Filter
    logger.info("Reading WGS MatrixTable and filtering intervals...")
    # NOTE: filter_intervals is optimized to push down into the reader if possible
    wgs_mt = hl.read_matrix_table(WGS_MT_PATH)
    
    # Filter to intervals FIRST (Massive speedup)
    mt = hl.filter_intervals(wgs_mt, hl_intervals)
    
    # 2. Filter Samples
    logger.info(f"Filtering to {len(eur_ids)} EUR samples...")
    mt = mt.filter_cols(hl.literal(eur_ids).contains(mt.s))
    
    # 3. Select fields early to minimize data movement
    mt = mt.select_entries('GT')
    
    # 4. Checkpoint (Crucial Step)
    # Write this intermediate result to a temporary location to avoid re-reading WGS
    # Use a hashed path or standard tmp location
    checkpoint_path = f"{config['outputs']['intermediate_dir']}/wgs_eur_intervals.mt"
    logger.info(f"Checkpointing filtered MT to {checkpoint_path}...")
    
    # Overwrite if exists
    mt = mt.checkpoint(checkpoint_path, overwrite=True)
    
    logger.info(f"Checkpoint complete. Rows: {mt.count_rows()}, Cols: {mt.count_cols()}")

    # 5. Variant QC & Final Filter
    logger.info("Running Variant QC...")
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows(
        (mt.variant_qc.AF[1] > params['maf_threshold']) & 
        (mt.variant_qc.AF[1] < (1 - params['maf_threshold'])) &
        (mt.variant_qc.call_rate > params['call_rate_threshold'])
    )
    
    n_total = mt.count_rows()
    logger.info(f"Variants passing QC: {n_total}")
    
    # 6. Downsample if needed
    target = params['target_variants']
    if n_total > target:
        fraction = target / n_total
        logger.info(f"Downsampling to {target} variants (fraction={fraction:.5f})")
        mt = mt.sample_rows(fraction, seed=42)
    
    # 7. Export
    plink_output_base = OUTPUT_PATH.replace(".mt", "")
    logger.info(f"Exporting directly to PLINK: {plink_output_base}")
    
    # Ensure output dir exists
    os.makedirs(os.path.dirname(plink_output_base), exist_ok=True)
    
    hl.export_plink(
        mt,
        plink_output_base,
        ind_id=mt.s,
        fam_id=mt.s
    )
    
    # Clean up checkpoint if desired? (Optional, maybe keep for debugging)
    # hl.hadoop_delete(checkpoint_path)
    
    logger.info("Job completed successfully.")

if __name__ == "__main__":
    main()
