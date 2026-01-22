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
        if start + params['window_size'] > chrom_len:
            start = chrom_len - params['window_size']
        interval = hl.locus_interval(chrom, start, start + params['window_size'], reference_genome='GRCh38')
        hl_intervals.append(interval)

    # 1. Create Interval Table (Fixes LineTooLong error)
    logger.info(f"Parallelizing {len(hl_intervals)} intervals to Hail Table...")
    interval_structs = [{'interval': i} for i in hl_intervals]
    interval_ht = hl.Table.parallelize(
        interval_structs, 
        schema=hl.tstruct(interval=hl.tinterval(hl.tlocus('GRCh38')))
    ).key_by('interval')

    # 2. Read and Filter
    logger.info("Reading WGS MatrixTable...")
    wgs_mt = hl.read_matrix_table(WGS_MT_PATH)
    
    logger.info("Filtering MatrixTable to intervals...")
    # Using filter_rows with a table join avoids passing a massive list literal
    mt = wgs_mt.filter_rows(hl.is_defined(interval_ht[wgs_mt.locus]))
    
    # 3. Filter Samples
    logger.info(f"Filtering to {len(eur_ids)} EUR samples...")
    mt = mt.filter_cols(hl.literal(eur_ids).contains(mt.s))
    
    # 4. Select fields early
    mt = mt.select_entries('GT')
    
    # 5. Checkpoint (Crucial Step)
    checkpoint_path = f"{config['outputs']['intermediate_dir']}/wgs_eur_intervals.mt"
    logger.info(f"Checkpointing filtered MT to {checkpoint_path}...")
    
    # Checkpointing saves the filtered subset to disk, preventing re-reads of the massive WGS MT
    mt = mt.checkpoint(checkpoint_path, overwrite=True)
    
    logger.info(f"Checkpoint complete. Rows: {mt.count_rows()}, Cols: {mt.count_cols()}")

    # 6. Variant QC & Final Filter
    logger.info("Running Variant QC...")
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows(
        (mt.variant_qc.AF[1] > params['maf_threshold']) & 
        (mt.variant_qc.AF[1] < (1 - params['maf_threshold'])) &
        (mt.variant_qc.call_rate > params['call_rate_threshold'])
    )
    
    n_total = mt.count_rows()
    logger.info(f"Variants passing QC: {n_total}")
    
    # 7. Downsample if needed
    target = params['target_variants']
    if n_total > target:
        fraction = target / n_total
        logger.info(f"Downsampling to {target} variants (fraction={fraction:.5f})")
        mt = mt.sample_rows(fraction, seed=42)
    
    # 8. Export
    plink_output_base = OUTPUT_PATH.replace(".mt", "")
    logger.info(f"Exporting directly to PLINK: {plink_output_base}")
    
    os.makedirs(os.path.dirname(plink_output_base), exist_ok=True)
    
    hl.export_plink(
        mt,
        plink_output_base,
        ind_id=mt.s,
        fam_id=mt.s
    )
    
    logger.info("Job completed successfully.")

if __name__ == "__main__":
    main()