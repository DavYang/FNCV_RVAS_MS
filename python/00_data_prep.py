#!/usr/bin/env python3
import os
import random
import hail as hl
import pandas as pd
from utils import load_config, init_hail, setup_logger

logger = setup_logger("data_prep")

def get_filtered_mt(config, interval_ht, eur_ids):
    """Re-creates the filtered MatrixTable query."""
    WGS_MT_PATH = config['inputs']['wgs_mt']
    params = config['params']
    
    wgs_mt = hl.read_matrix_table(WGS_MT_PATH)
    mt = wgs_mt.filter_cols(hl.literal(eur_ids).contains(wgs_mt.s))
    mt = mt.select_entries('GT')
    mt = mt.filter_rows(hl.is_defined(interval_ht[mt.locus]))
    mt = mt.naive_coalesce(500)
    
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows(
        (mt.variant_qc.AF[1] > params['maf_threshold']) & 
        (mt.variant_qc.AF[1] < (1 - params['maf_threshold'])) &
        (mt.variant_qc.call_rate > params['call_rate_threshold'])
    )
    return mt

def main():
    config = load_config("config/config.json")
    
    # Initialize with high memory
    init_hail("data_prep", driver_mem="16g", reference="GRCh38")

    # Paths
    ANCESTRY_PATH = config['inputs']['ancestry_pred']
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
    for _ in range(params['n_intervals']):
        chrom = random.choice(chroms)
        chrom_len = lengths[chrom]
        start = random.randint(1, chrom_len - params['window_size'])
        interval = hl.locus_interval(chrom, start, start + params['window_size'], reference_genome='GRCh38')
        hl_intervals.append(interval)

    interval_structs = [{'interval': i} for i in hl_intervals]
    interval_ht = hl.Table.parallelize(
        interval_structs, 
        schema=hl.tstruct(interval=hl.tinterval(hl.tlocus('GRCh38')))
    ).key_by('interval')

    # PASS 1: Count
    logger.info("PASS 1: Counting variants (Streaming)...")
    mt_query = get_filtered_mt(config, interval_ht, eur_ids)
    n_total = mt_query.count_rows()
    logger.info(f"Total variants found: {n_total}")

    target = params['target_variants']
    fraction = 1.0
    if n_total > target:
        fraction = target / n_total
        logger.info(f"Targeting {target} variants (fraction={fraction:.5f})")

    # PASS 2: Export
    logger.info("PASS 2: Sampling and Exporting (Streaming)...")
    mt_final = get_filtered_mt(config, interval_ht, eur_ids)
    
    if fraction < 1.0:
        mt_final = mt_final.sample_rows(fraction, seed=42)

    plink_output_base = OUTPUT_PATH.replace(".mt", "")
    logger.info(f"Exporting directly to PLINK: {plink_output_base}")
    
    hl.export_plink(
        mt_final,
        plink_output_base,
        ind_id=mt_final.s,
        fam_id=mt_final.s
    )
    
    logger.info("Job completed successfully.")

if __name__ == "__main__":
    main()
