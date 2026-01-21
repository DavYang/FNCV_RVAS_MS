#!/usr/bin/env python3
import os
import random
import hail as hl
import pandas as pd
from utils import load_config, init_hail, setup_logger

# Setup
logger = setup_logger("data_prep")

def main():
    # 1. Load Config
    config = load_config("config/config.json")
    
    # --- STABILITY FIX: Configure Spark Driver Memory ---
    hl.init(
        default_reference='GRCh38', 
        log='/tmp/hail_data_prep.log',
        spark_conf={
            'spark.driver.memory': '8g',
            'spark.executor.memory': '8g'
        }
    )

    # Extract paths from config
    WGS_MT_PATH = config['inputs']['wgs_mt']
    ANCESTRY_PATH = config['inputs']['ancestry_pred']
    
    # Construct output paths using config base dirs
    OUT_DIR = config['outputs']['data_dir']
    TMP_DIR = config['outputs']['intermediate_dir']
    
    OUTPUT_PATH = f"{OUT_DIR}/eur_common_snps_500k.mt"
    
    params = config['params']

    # 2. Get European Samples
    logger.info("Loading EUR sample list...")
    with hl.hadoop_open(ANCESTRY_PATH, 'r') as f:
        ancestry_df = pd.read_csv(f, sep="\t")
    eur_ids = set(ancestry_df[ancestry_df['ancestry_pred'] == 'eur']['research_id'].astype(str))
    
    # 3. Generate Random Intervals
    logger.info(f"Generating {params['n_intervals']} random genomic intervals...")
    rg = hl.get_reference('GRCh38')
    lengths = rg.lengths
    
    chroms = [f"chr{i}" for i in range(1, 23)]
    hl_intervals = []
    
    for _ in range(params['n_intervals']):
        chrom = random.choice(chroms)
        chrom_len = lengths[chrom]
        start = random.randint(1, chrom_len - params['window_size'])
        interval = hl.locus_interval(chrom, start, start + params['window_size'], reference_genome='GRCh38')
        hl_intervals.append(interval)

    # Convert to Table
    interval_structs = [{'interval': i} for i in hl_intervals]
    interval_ht = hl.Table.parallelize(
        interval_structs, 
        schema=hl.tstruct(interval=hl.tinterval(hl.tlocus('GRCh38')))
    ).key_by('interval')

    # 4. Filter WGS MatrixTable
    logger.info("Loading and filtering MatrixTable...")
    wgs_mt = hl.read_matrix_table(WGS_MT_PATH)
    mt = wgs_mt.filter_cols(hl.literal(eur_ids).contains(wgs_mt.s))
    
    # Aggressive Pruning: Keep only GT
    mt = mt.select_entries('GT')
    
    # Filter Rows
    mt = mt.filter_rows(hl.is_defined(interval_ht[mt.locus]))
    
    # --- STABILITY FIX: Reduce partitions further ---
    logger.info("Coalescing partitions to 500...")
    mt = mt.naive_coalesce(500)

    # 5. QC and Filtering (Using Config Thresholds)
    logger.info("Running QC...")
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows(
        (mt.variant_qc.AF[1] > params['maf_threshold']) & 
        (mt.variant_qc.AF[1] < (1 - params['maf_threshold'])) &
        (mt.variant_qc.call_rate > params['call_rate_threshold'])
    )

    # 6. Checkpoint REMOVED (Hardcoded removal as requested)
    # Proceeding directly to persistence and sampling

    # 7. Sampling
    logger.info("Persisting dataset in memory...")
    mt = mt.persist() 
    n_total = mt.count_rows()
    logger.info(f"Total variants passing filters: {n_total}")
    
    target = params['target_variants']
    if n_total > target:
        fraction = target / n_total
        logger.info(f"Downsampling to {target} (fraction={fraction:.4f})...")
        mt_final = mt.sample_rows(fraction, seed=42)
    else:
        mt_final = mt

    # 8. Final Save (MatrixTable)
    logger.info(f"Saving final dataset to {OUTPUT_PATH}...")
    mt_final.write(OUTPUT_PATH, overwrite=True)

    # 9. EXPORT TO PLINK
    plink_output_base = OUTPUT_PATH.replace(".mt", "")
    logger.info(f"Exporting to PLINK format: {plink_output_base}")
    hl.export_plink(
        mt_final,
        plink_output_base,
        ind_id=mt_final.s,
        fam_id=mt_final.s
    )

    logger.info("Job completed successfully.")

if __name__ == "__main__":
    main()
