#!/usr/bin/env python3
import os
import pandas as pd
import hail as hl
from utils import load_config, init_hail, setup_logger

logger = setup_logger("define_loci")

def main():
    config = load_config("config/config.json")
    # Initialize with reference genome to avoid warnings
    init_hail("define_loci", reference="GRCh38")
    
    # Paths
    GWAS_PATH = config['inputs']['phenotype_gwas']
    LOCI_OUTPUT = config['outputs']['loci_bed']
    
    # Params
    P_THRESH = config['params']['gwas_p_threshold']
    MHC_INTERVAL = config['params']['mhc_interval']
    FLANK = config['params']['locus_flank_window']

    # 1. Load GWAS
    logger.info(f"Loading GWAS results from {GWAS_PATH}...")
    gwas_ht = hl.read_table(GWAS_PATH)

    # 2. Filter Significant Hits
    logger.info(f"Filtering for P < {P_THRESH}...")
    sig_ht = gwas_ht.filter(gwas_ht.p_value < P_THRESH)
    
    # 3. Exclude MHC
    logger.info(f"Excluding MHC region ({MHC_INTERVAL})...")
    sig_ht = sig_ht.filter(~hl.parse_locus_interval(MHC_INTERVAL).contains(sig_ht.locus))
    
    n_sig = sig_ht.count()
    logger.info(f"Found {n_sig} significant SNPs outside MHC.")
    
    if n_sig == 0:
        logger.warning("No significant SNPs found. Exiting.")
        return

    # 4. Greedy Clumping (Distance-Based)
    logger.info("Collecting hits for clumping...")
    df_sig = sig_ht.select('p_value').to_pandas()
    df_sig['chrom'] = df_sig['locus'].apply(lambda x: x.contig)
    df_sig['pos'] = df_sig['locus'].apply(lambda x: x.position)
    df_sig = df_sig.sort_values('p_value') # Best P-value first
    
    defined_loci = []
    
    logger.info(f"Performing distance-based clumping (Window +/- {FLANK}bp)...")
    while not df_sig.empty:
        # Lead SNP is top of list
        lead = df_sig.iloc[0]
        
        chrom = lead['chrom']
        start = max(1, lead['pos'] - FLANK)
        end = lead['pos'] + FLANK
        
        defined_loci.append({
            'chrom': chrom,
            'start': start,
            'end': end,
            'lead_snp_pos': lead['pos'],
            'lead_snp_p': lead['p_value']
        })
        
        # Remove neighbors
        df_sig = df_sig[~((df_sig['chrom'] == chrom) & 
                          (df_sig['pos'] >= start) & 
                          (df_sig['pos'] <= end))]
    
    # 5. Save Results
    loci_df = pd.DataFrame(defined_loci)
    logger.info(f"Defined {len(loci_df)} independent loci.")
    
    # Ensure output dir exists
    os.makedirs(os.path.dirname(LOCI_OUTPUT), exist_ok=True)
    
    loci_df.to_csv(LOCI_OUTPUT, sep="\t", index=False)
    logger.info(f"Loci definitions saved to {LOCI_OUTPUT}")

if __name__ == "__main__":
    main()
