#!/usr/bin/env python3
import os
import subprocess
import pandas as pd
import hail as hl
from utils import load_config, init_hail, setup_logger

logger = setup_logger("define_loci")

def run_gcta_cojo(chrom, start, end, ma_path, bfile_path, out_prefix, gcta_bin):
    """Runs GCTA-COJO conditional analysis for a specific locus."""
    cmd = [
        gcta_bin,
        "--bfile", bfile_path,
        "--chr", chrom.replace("chr", ""), # GCTA often expects 1-22
        "--cojo-file", ma_path,
        "--cojo-slct", # Select independent signals
        "--cojo-p", "5e-8", # Threshold for genome-wide significance
        "--out", out_prefix
    ]
    
    # Add window restriction to GCTA to speed up loading
    # (Optional: extracting plink subset first might be safer if bfile is huge)
    
    try:
        logger.info(f"Running GCTA-COJO for {chrom}:{start}-{end}...")
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"GCTA-COJO failed for {out_prefix}: {e.stderr.decode()}")
        return False

def export_sumstats_for_gcta(ht, output_path):
    """Exports Hail Table to GCTA .ma format: SNP A1 A2 freq b se p N"""
    # GCTA expects: SNP A1 A2 freq b se p N
    # We assume 'ht' has these fields or we rename them.
    # Standardizing fields:
    # locus, alleles -> SNP (chr:pos:ref:alt), A1 (alt), A2 (ref)
    
    # Annotate necessary fields for export
    ht = ht.annotate(
        SNP = hl.str(ht.locus.contig) + ":" + hl.str(ht.locus.position) + ":" + ht.alleles[0] + ":" + ht.alleles[1],
        A1 = ht.alleles[1],
        A2 = ht.alleles[0],
        freq = ht.AF, # Assuming AF is present in GWAS results
        b = ht.beta,
        se = ht.standard_error,
        p = ht.p_value,
        N = ht.n_complete_samples # Or specific N field
    )
    
    ht.select('SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N').export(output_path)

def main():
    config = load_config("config/config.json")
    # Initialize with reference genome
    init_hail("define_loci", reference="GRCh38")
    
    # Paths
    GWAS_PATH = config['inputs']['phenotype_gwas']
    LOCI_OUTPUT = config['outputs']['loci_bed']
    COJO_DIR = f"{config['outputs']['base_dir']}/cojo_results"
    GCTA_BIN = os.path.abspath(config['tools']['gcta'])
    
    # LD Reference (From Phase 1)
    REF_BFILE = os.path.abspath(config['outputs']['data_dir'] + "/eur_common_snps_500k")
    
    # Params
    P_THRESH = config['params']['gwas_p_threshold']
    MHC_INTERVAL = config['params']['mhc_interval']
    FLANK = config['params']['locus_flank_window']

    # 1. Load GWAS
    logger.info(f"Loading GWAS results from {GWAS_PATH}...")
    gwas_ht = hl.read_table(GWAS_PATH)
    
    # Ensure N field exists (fix name if needed)
    if 'n_complete_samples' not in gwas_ht.row:
         # Fallback or check config for N field name
         logger.warning("Field 'n_complete_samples' not found. Checking for 'N'...")
         if 'N' in gwas_ht.row:
             gwas_ht = gwas_ht.annotate(n_complete_samples = gwas_ht.N)
         else:
             raise ValueError("GWAS HT must have sample size field (N or n_complete_samples)")

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

    # 4. Simple Clumping to define analysis windows
    logger.info("Defining windows for COJO...")
    df_sig = sig_ht.select('p_value').to_pandas()
    df_sig['chrom'] = df_sig['locus'].apply(lambda x: x.contig)
    df_sig['pos'] = df_sig['locus'].apply(lambda x: x.position)
    df_sig = df_sig.sort_values('p_value') 
    
    analysis_windows = []
    
    while not df_sig.empty:
        lead = df_sig.iloc[0]
        chrom = lead['chrom']
        # Define window around lead SNP
        start = max(1, lead['pos'] - FLANK)
        end = lead['pos'] + FLANK
        
        analysis_windows.append({
            'chrom': chrom,
            'start': start,
            'end': end,
            'lead_snp': lead
        })
        
        # Remove neighbors
        df_sig = df_sig[~((df_sig['chrom'] == chrom) & 
                          (df_sig['pos'] >= start) & 
                          (df_sig['pos'] <= end))]
    
    logger.info(f"Identified {len(analysis_windows)} candidate loci for COJO.")
    
    # 5. Run COJO per window
    os.makedirs(COJO_DIR, exist_ok=True)
    final_signals = []

    for i, window in enumerate(analysis_windows):
        chrom = window['chrom']
        start = window['start']
        end = window['end']
        locus_id = f"locus_{i}_{chrom}_{start}_{end}"
        
        logger.info(f"Processing {locus_id}...")
        
        # Filter GWAS to this window
        interval = hl.locus_interval(chrom, start, end, reference_genome='GRCh38')
        window_ht = gwas_ht.filter(interval.contains(gwas_ht.locus))
        
        # Export .ma file
        ma_file = f"{COJO_DIR}/{locus_id}.ma"
        export_sumstats_for_gcta(window_ht, ma_file)
        
        # Run GCTA
        cojo_out = f"{COJO_DIR}/{locus_id}"
        success = run_gcta_cojo(chrom, start, end, ma_file, REF_BFILE, cojo_out, GCTA_BIN)
        
        if success:
            # Read .jma.cojo file to get independent signals
            res_file = f"{cojo_out}.jma.cojo"
            if os.path.exists(res_file):
                try:
                    res_df = pd.read_csv(res_file, sep="\t")
                    logger.info(f"  -> Found {len(res_df)} independent signals.")
                    res_df['locus_id'] = locus_id
                    final_signals.append(res_df)
                except Exception as e:
                    logger.warning(f"  -> Failed to read COJO output: {e}")
            else:
                logger.warning("  -> No .jma.cojo file produced (possibly no variants in LD ref).")
    
    # 6. Aggregate and Save
    if final_signals:
        all_signals = pd.concat(final_signals, ignore_index=True)
        out_signals = f"{config['outputs']['base_dir']}/all_independent_signals.tsv"
        all_signals.to_csv(out_signals, sep="\t", index=False)
        logger.info(f"Saved all independent signals to {out_signals}")
        
        # Create final BED of loci based on COJO signals (or keep original windows)
        # For now, saving original windows + signal count
        loci_df = pd.DataFrame(analysis_windows)
        loci_df.to_csv(LOCI_OUTPUT, sep="\t", index=False)
        logger.info(f"Locus definitions saved to {LOCI_OUTPUT}")
    else:
        logger.warning("No independent signals found via COJO.")

if __name__ == "__main__":
    main()