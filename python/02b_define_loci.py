#!/usr/bin/env python3
"""Phase 2b: Lead SNP extraction, greedy clumping, and GCTA-COJO signal dissection.

Processes each chromosome independently using:
  - Pre-built per-chromosome .ma files (from 02a_export_gwas_ma.py)
  - Per-chromosome LD reference PLINK files (from 01b_qc_merge_background_snps.sh)

For each chromosome with significant GWAS hits:
  1. Greedy clumping to define candidate loci
  2. PLINK2 subset of LD ref to locus window (+500kb LD padding)
  3. GCTA-COJO --cojo-slct for stepwise conditional analysis
  4. Parse .jma.cojo output for independent signals

Outputs:
  - results/2-locus_definition/all_independent_signals.tsv
  - results/2-locus_definition/target_loci.bed
  - results/2-locus_definition/cojo/chrN/{locus_id}.jma.cojo  (per locus)

Usage:
    python python/02b_define_loci.py --config config/config.json
    python python/02b_define_loci.py --config config/config.json --force
    python python/02b_define_loci.py --config config/config.json --chrom chr6
"""

import argparse
import json
import logging
import os
import shutil
import subprocess
import sys
import tempfile
import time
from datetime import datetime
from typing import Dict, List, Optional, Tuple

import hail as hl
import pandas as pd

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

def setup_logger(name: str) -> logging.Logger:
    """Set up a logger writing to stdout with timestamps."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)],
    )
    return logging.getLogger(name)


logger = setup_logger("define_loci")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_FIELD_CANDIDATES: Dict[str, List[str]] = {
    'beta':   ['beta', 'effect_size', 'b'],
    'se':     ['standard_error', 'se', 'stderr'],
    'pval':   ['p_value', 'p_value_EUR', 'pval', 'p'],
    'af':     ['AF', 'AF_EUR', 'allele_frequency', 'freq'],
    'n':      ['n_complete_samples', 'N', 'n_samples', 'n'],
}

CHR_SIZES: Dict[str, int] = {
    'chr1': 248956422,  'chr2': 242193529,  'chr3': 198295559,
    'chr4': 190214555,  'chr5': 181538259,  'chr6': 170805979,
    'chr7': 159345973,  'chr8': 145138636,  'chr9': 138394717,
    'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
    'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189,
    'chr16': 90338345,  'chr17': 83257441,  'chr18': 80373285,
    'chr19': 58617616,  'chr20': 64444167,  'chr21': 46709983,
    'chr22': 50818468,
}

# Default GCTA-COJO parameters
DEFAULT_COJO_P = 5e-8
DEFAULT_COJO_WIND = 10000
DEFAULT_COJO_COLLINEAR = 0.9

# LD padding around each locus for PLINK subset
LD_PAD_BP = 500_000

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

def load_config(config_path: str) -> dict:
    """Load JSON config file."""
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Config file not found: {config_path}")
    with open(config_path, "r") as f:
        return json.load(f)

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Phase 2b: GCTA-COJO locus definition (per-chromosome)"
    )
    parser.add_argument(
        '--config', default='config/config.json',
        help='Path to config.json'
    )
    parser.add_argument(
        '--output-dir', default=None,
        help='Override base output directory for locus definition'
    )
    parser.add_argument(
        '--chrom', default=None,
        help='Process a single chromosome (e.g. chr6). Default: all with hits.'
    )
    parser.add_argument(
        '--force', action='store_true',
        help='Re-run COJO even if .jma.cojo output exists'
    )
    return parser.parse_args()

# ---------------------------------------------------------------------------
# Helper utilities
# ---------------------------------------------------------------------------

def _fmt_elapsed(seconds: float) -> str:
    """Format elapsed seconds as human-readable string."""
    if seconds < 60:
        return f"{seconds:.1f}s"
    m, s = divmod(int(seconds), 60)
    if m < 60:
        return f"{m}m {s:02d}s"
    h, m = divmod(m, 60)
    return f"{h}h {m:02d}m {s:02d}s"


def _resolve_field(ht: hl.Table, canonical: str, candidates: List[str]) -> str:
    """Return the first matching field name from candidates present in ht.row."""
    row_fields = set(ht.row)
    for name in candidates:
        if name in row_fields:
            logger.info(f"  Field '{canonical}' resolved to '{name}'")
            return name
    raise ValueError(
        f"Could not find field '{canonical}' in GWAS HT. "
        f"Tried: {candidates}. Available fields: {sorted(row_fields)}"
    )


def _run_command(
    cmd: List[str],
    log_path: str,
    label: str,
) -> Tuple[bool, str]:
    """Run a subprocess, tee stdout+stderr to log_path, return (success, stderr_tail)."""
    logger.info(f"  Running: {' '.join(cmd)}")
    try:
        with open(log_path, 'w') as fh:
            result = subprocess.run(
                cmd,
                stdout=fh,
                stderr=subprocess.STDOUT,
                check=True,
            )
        return True, ""
    except subprocess.CalledProcessError as exc:
        tail = ""
        try:
            with open(log_path) as fh:
                lines = fh.readlines()
            tail = "".join(lines[-20:])
        except Exception:
            pass
        logger.error(f"  {label} FAILED (exit {exc.returncode})\n{tail}")
        return False, tail


def _normalize_contig(contig: str) -> str:
    """Ensure contig has chr-prefix for GRCh38."""
    if not contig.startswith('chr'):
        return f"chr{contig}"
    return contig

# ---------------------------------------------------------------------------
# Step 1: Load and normalise GWAS HT
# ---------------------------------------------------------------------------

def load_gwas_ht(
    gwas_path: str,
    p_thresh: float,
    mhc_interval: str,
) -> Tuple[hl.Table, hl.Table, Dict[str, str]]:
    """Load GWAS HT, resolve field names, filter p < p_thresh, exclude MHC.

    Returns:
        gwas_ht    - full GWAS table (for reference)
        sig_ht     - significant-only table (p < p_thresh, no MHC)
        field_map  - mapping of canonical name -> actual HT field name
    """
    logger.info(f"Loading GWAS HT from {gwas_path} ...")
    gwas_ht = hl.read_table(gwas_path)

    logger.info("Resolving GWAS field names ...")
    field_map = {
        canonical: _resolve_field(gwas_ht, canonical, candidates)
        for canonical, candidates in _FIELD_CANDIDATES.items()
    }

    sample_contig = gwas_ht.locus.contig.collect()[0]
    has_chr_prefix = sample_contig.startswith('chr')
    logger.info(f"  Contig format in HT: '{sample_contig}' (chr-prefix: {has_chr_prefix})")

    if has_chr_prefix:
        mhc_str = mhc_interval if mhc_interval.startswith('chr') else f"chr{mhc_interval}"
    else:
        mhc_str = mhc_interval.replace('chr', '') if mhc_interval.startswith('chr') else mhc_interval

    p_field = field_map['pval']
    logger.info(f"Filtering GWAS: p_field='{p_field}', threshold={p_thresh}")
    sig_ht = gwas_ht.filter(gwas_ht[p_field] < p_thresh)

    logger.info(f"Excluding MHC region ({mhc_str}) ...")
    mhc_ivl = hl.parse_locus_interval(mhc_str, reference_genome='GRCh38')
    sig_ht = sig_ht.filter(~mhc_ivl.contains(sig_ht.locus))

    n_sig = sig_ht.count()
    logger.info(f"Significant SNPs outside MHC: {n_sig:,}")

    if n_sig == 0:
        raise RuntimeError(
            f"No significant SNPs found at p < {p_thresh} outside MHC. "
            "Check gwas_p_threshold and GWAS HT path in config."
        )

    return gwas_ht, sig_ht, field_map

# ---------------------------------------------------------------------------
# Step 2: Greedy clumping
# ---------------------------------------------------------------------------

def greedy_clump(
    sig_ht: hl.Table,
    p_field: str,
    flank: int,
    has_chr_prefix: bool,
) -> List[Dict]:
    """Greedy positional clumping of significant SNPs into non-overlapping windows.

    Returns list of dicts: {chrom, start, end, lead_snp_id, lead_pos, lead_p}
    """
    logger.info("Collecting significant SNPs for greedy clumping ...")

    df = sig_ht.select(
        p_value=sig_ht[p_field],
        contig=hl.str(sig_ht.locus.contig),
        position=sig_ht.locus.position,
    ).to_pandas()

    df = df.sort_values('p_value').reset_index(drop=True)
    df['snp_id'] = df['contig'] + ':' + df['position'].astype(str)

    windows = []
    consumed = [False] * len(df)

    for i, row in df.iterrows():
        if consumed[i]:
            continue
        chrom = row['contig']
        pos = int(row['position'])
        start = max(1, pos - flank)
        end = pos + flank

        windows.append({
            'chrom': _normalize_contig(chrom) if not chrom.startswith('chr') else chrom,
            'start': start,
            'end': end,
            'lead_snp_id': row['snp_id'],
            'lead_pos': pos,
            'lead_p': row['p_value'],
        })

        mask = (df['contig'] == chrom) & (df['position'] >= start) & (df['position'] <= end)
        consumed_idx = df[mask].index
        for idx in consumed_idx:
            consumed[idx] = True

    logger.info(f"Greedy clumping: {len(windows)} candidate loci from {len(df):,} significant SNPs")
    return windows

# ---------------------------------------------------------------------------
# Step 3: PLINK subset for locus (from per-chrom LD ref)
# ---------------------------------------------------------------------------

def extract_plink_subset(
    bfile: str,
    start: int,
    end: int,
    work_dir: str,
    locus_id: str,
    plink_threads: int = 4,
    plink_mem_mb: int = 8000,
) -> Optional[str]:
    """Extract a PLINK subset for the locus window (+LD padding) from a
    per-chromosome LD reference bfile.

    No --chr needed since bfile is already single-chromosome.

    Args:
        bfile: PLINK bfile prefix (per-chrom LD ref).
        start: Locus start position.
        end: Locus end position.
        work_dir: Temporary working directory.
        locus_id: Locus identifier for file naming.
        plink_threads: Number of threads for PLINK2.
        plink_mem_mb: Memory limit in MB for PLINK2.

    Returns:
        Path to subset bfile prefix, or None on failure.
    """
    ld_start = max(1, start - LD_PAD_BP)
    ld_end = end + LD_PAD_BP

    subset_prefix = os.path.join(work_dir, f"subset_{locus_id}")
    log_path = f"{subset_prefix}_plink.log"

    cmd = [
        'plink2',
        '--bfile', bfile,
        '--from-bp', str(ld_start),
        '--to-bp', str(ld_end),
        '--make-bed',
        '--out', subset_prefix,
        '--threads', str(plink_threads),
        '--memory', str(plink_mem_mb),
        '--no-psam-pheno',
    ]

    success, _ = _run_command(cmd, log_path, f"plink2 subset {locus_id}")
    if not success:
        return None

    if not os.path.exists(f"{subset_prefix}.bed"):
        logger.warning(f"  PLINK subset .bed not found: {subset_prefix}.bed")
        return None

    n_vars = sum(1 for _ in open(f"{subset_prefix}.bim"))
    logger.info(f"  PLINK subset: {n_vars:,} variants in LD window")
    return subset_prefix

# ---------------------------------------------------------------------------
# Step 4: Run GCTA-COJO
# ---------------------------------------------------------------------------

def run_gcta_cojo(
    gcta_bin: str,
    bfile: str,
    ma_path: str,
    out_prefix: str,
    p_thresh: float,
) -> bool:
    """Run GCTA-COJO --cojo-slct on a pre-extracted per-chrom PLINK subset.

    No --chr flag needed since bfile is already single-chromosome.

    Args:
        gcta_bin: Path to GCTA binary.
        bfile: PLINK bfile prefix for LD reference subset.
        ma_path: Path to full-chromosome .ma file.
        out_prefix: Output prefix for COJO results.
        p_thresh: P-value threshold for --cojo-p.

    Returns:
        True if COJO ran successfully.
    """
    log_path = f"{out_prefix}.gcta.log"

    cmd = [
        gcta_bin,
        '--bfile', bfile,
        '--cojo-file', ma_path,
        '--cojo-slct',
        '--cojo-p', str(p_thresh),
        '--cojo-wind', str(DEFAULT_COJO_WIND),
        '--cojo-collinear', str(DEFAULT_COJO_COLLINEAR),
        '--out', out_prefix,
    ]

    success, _ = _run_command(cmd, log_path, f"GCTA-COJO")
    return success

# ---------------------------------------------------------------------------
# Step 5: Parse COJO output
# ---------------------------------------------------------------------------

def parse_cojo_output(
    jma_path: str,
    locus_id: str,
    chrom: str,
) -> Optional[pd.DataFrame]:
    """Parse GCTA .jma.cojo file. Returns DataFrame with independent signals,
    or None if file missing/empty."""
    if not os.path.exists(jma_path):
        logger.warning(f"  .jma.cojo not found: {jma_path} (no independent signals selected)")
        return None

    try:
        df = pd.read_csv(jma_path, sep='\t')
    except Exception as exc:
        logger.warning(f"  Failed to read {jma_path}: {exc}")
        return None

    if df.empty:
        logger.warning(f"  .jma.cojo is empty for {locus_id}")
        return None

    df['locus_id'] = locus_id
    df['chrom'] = chrom
    logger.info(f"  COJO: {len(df)} independent signal(s) for {locus_id}")
    return df

# ---------------------------------------------------------------------------
# Step 6: Assemble BED from COJO signals
# ---------------------------------------------------------------------------

def build_loci_bed(
    all_signals_df: pd.DataFrame,
    flank: int,
    mhc_chrom: str,
    mhc_start: int,
    mhc_end: int,
) -> pd.DataFrame:
    """Build target_loci.bed from COJO independent signals.

    Each signal becomes a locus: signal_pos +/- flank, clipped to chr bounds.
    MHC signals are excluded. Columns: chrom, start, end, locus_id.
    """
    rows = []

    for _, sig in all_signals_df.iterrows():
        chrom = sig['chrom']
        pos = int(sig.get('bp', sig.get('pos', 0)))
        if pos == 0:
            logger.warning(
                f"  Could not determine position for signal in "
                f"{sig.get('locus_id')}; skipping"
            )
            continue

        start = max(1, pos - flank)
        end = pos + flank
        if chrom in CHR_SIZES:
            end = min(end, CHR_SIZES[chrom])

        if chrom == mhc_chrom and not (end < mhc_start or start > mhc_end):
            logger.info(f"  Excluding MHC-overlapping signal: {chrom}:{pos}")
            continue

        signal_id = f"{chrom}_{pos}_{sig.get('SNP', 'unknown')}"
        rows.append({
            'chrom':    chrom,
            'start':    start,
            'end':      end,
            'locus_id': signal_id,
        })

    return pd.DataFrame(rows, columns=['chrom', 'start', 'end', 'locus_id'])

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    """Run Phase 2b: lead SNP extraction + per-chromosome GCTA-COJO."""
    args = parse_args()

    run_start = time.time()
    logger.info("=" * 60)
    logger.info("PHASE 2b: LOCUS DEFINITION (GCTA-COJO, per-chromosome)")
    logger.info("=" * 60)
    logger.info(f"  Timestamp : {datetime.now().isoformat()}")
    logger.info(f"  Config    : {args.config}")

    config = load_config(args.config)

    # --- Paths ---
    locus_dir = args.output_dir or config['outputs'].get(
        'locus_def_dir', 'results/2-locus_definition'
    )
    cojo_dir = os.path.join(locus_dir, 'cojo')
    ma_dir = os.path.join(locus_dir, 'ma')
    ld_ref_dir = config['outputs'].get(
        'ld_reference', 'results/1-bg_snp/plink_qc'
    )

    os.makedirs(cojo_dir, exist_ok=True)

    loci_bed_path = os.path.join(locus_dir, 'target_loci.bed')
    signals_path = os.path.join(locus_dir, 'all_independent_signals.tsv')

    gwas_path = config['inputs']['phenotype_gwas']
    gcta_bin = os.path.abspath(config['tools']['gcta'])

    p_thresh = float(config['params']['gwas_p_threshold'])
    mhc_interval = config['params']['mhc_interval']
    flank = int(config['params']['locus_flank_window'])

    logger.info(f"  GWAS HT      : {gwas_path}")
    logger.info(f"  LD ref dir   : {ld_ref_dir}")
    logger.info(f"  .ma dir      : {ma_dir}")
    logger.info(f"  GCTA bin     : {gcta_bin}")
    logger.info(f"  P thresh     : {p_thresh}")
    logger.info(f"  COJO wind    : {DEFAULT_COJO_WIND} kb")
    logger.info(f"  MHC excl     : {mhc_interval}")
    logger.info(f"  Flank        : {flank:,} bp")
    logger.info(f"  Output       : {locus_dir}")
    logger.info(f"  Force        : {args.force}")
    logger.info("")

    # Validate dependencies
    if not os.path.exists(gcta_bin):
        raise FileNotFoundError(
            f"GCTA binary not found: {gcta_bin}\n"
            "Run python/install_gcta.py first."
        )

    # Initialize Hail (needed for loading GWAS HT)
    logger.info("Initializing Hail ...")
    hl.init(log='/tmp/hail_define_loci.log')
    hl.default_reference('GRCh38')
    logger.info("Hail initialized")
    logger.info("")

    try:
        # Step 1: Load GWAS and identify significant hits
        gwas_ht, sig_ht, field_map = load_gwas_ht(gwas_path, p_thresh, mhc_interval)

        # Step 2: Greedy clumping
        logger.info("")
        logger.info("Step 2: Greedy clumping ...")
        windows = greedy_clump(sig_ht, field_map['pval'], flank, True)
        logger.info(f"  -> {len(windows)} analysis windows")

        # Filter to single chromosome if requested
        if args.chrom:
            chrom_filter = args.chrom if args.chrom.startswith('chr') else f"chr{args.chrom}"
            windows = [w for w in windows if w['chrom'] == chrom_filter]
            logger.info(f"  -> {len(windows)} windows on {chrom_filter}")

        if not windows:
            logger.warning("No analysis windows to process. Exiting.")
            hl.stop()
            return

        # Group windows by chromosome
        chrom_windows: Dict[str, List[Dict]] = {}
        for win in windows:
            chrom_windows.setdefault(win['chrom'], []).append(win)

        logger.info(f"  Chromosomes with hits: {sorted(chrom_windows.keys())}")

        # Step 3-5: Per-chromosome COJO loop
        logger.info("")
        logger.info("Step 3: Running GCTA-COJO per chromosome ...")
        all_signals: List[pd.DataFrame] = []

        work_dir = tempfile.mkdtemp(prefix='cojo_work_')
        logger.info(f"  Working dir: {work_dir}")

        try:
            for chrom in sorted(chrom_windows.keys()):
                chr_windows = chrom_windows[chrom]
                chr_num = chrom.replace('chr', '')

                logger.info("")
                logger.info(f"--- {chrom} ({len(chr_windows)} loci) ---")

                # Validate per-chrom LD reference exists
                ld_ref_prefix = os.path.join(ld_ref_dir, f"{chrom}_background_qc")
                if not os.path.exists(f"{ld_ref_prefix}.bed"):
                    logger.error(
                        f"  LD reference not found: {ld_ref_prefix}.bed\n"
                        f"  Run bash/01b_qc_merge_background_snps.sh first."
                    )
                    continue

                # Validate per-chrom .ma file exists
                ma_path = os.path.join(ma_dir, f"{chrom}.ma")
                if not os.path.exists(ma_path):
                    logger.error(
                        f"  .ma file not found: {ma_path}\n"
                        f"  Run bash/02a_export_gwas_ma.sh first."
                    )
                    continue

                n_ld_vars = sum(1 for _ in open(f"{ld_ref_prefix}.bim"))
                n_ma_vars = sum(1 for _ in open(ma_path)) - 1
                logger.info(f"  LD ref: {n_ld_vars:,} variants")
                logger.info(f"  .ma   : {n_ma_vars:,} variants")

                # Create per-chrom cojo output directory
                chr_cojo_dir = os.path.join(cojo_dir, chrom)
                os.makedirs(chr_cojo_dir, exist_ok=True)

                for i, win in enumerate(chr_windows):
                    start = win['start']
                    end = win['end']
                    locus_id = f"locus_{chrom}_{start}_{end}"

                    logger.info(
                        f"  [{i+1}/{len(chr_windows)}] {locus_id} "
                        f"(lead: {win['lead_snp_id']}, p={win['lead_p']:.2e})"
                    )

                    # Resume: skip if COJO output already exists
                    jma_final = os.path.join(chr_cojo_dir, f"{locus_id}.jma.cojo")
                    if os.path.exists(jma_final) and not args.force:
                        logger.info(f"    SKIP (jma.cojo exists)")
                        sig_df = parse_cojo_output(jma_final, locus_id, chrom)
                        if sig_df is not None:
                            all_signals.append(sig_df)
                        continue

                    # 3a: Extract PLINK subset for locus window
                    subset_bfile = extract_plink_subset(
                        ld_ref_prefix, start, end, work_dir, locus_id,
                        plink_threads=4, plink_mem_mb=8000,
                    )
                    if subset_bfile is None:
                        logger.warning(
                            f"    PLINK subset failed for {locus_id}; skipping COJO"
                        )
                        continue

                    # 3b: Run GCTA-COJO
                    cojo_tmp_prefix = os.path.join(work_dir, locus_id)
                    success = run_gcta_cojo(
                        gcta_bin, subset_bfile, ma_path,
                        cojo_tmp_prefix, p_thresh,
                    )

                    # Copy COJO outputs to persistent cojo_dir
                    for ext in ['.jma.cojo', '.cma.cojo', '.ldr.cojo', '.gcta.log']:
                        src = f"{cojo_tmp_prefix}{ext}"
                        if os.path.exists(src):
                            shutil.copy2(
                                src, os.path.join(chr_cojo_dir, f"{locus_id}{ext}")
                            )

                    if not success:
                        logger.warning(f"    GCTA-COJO failed for {locus_id}")
                        continue

                    # 3c: Parse output
                    sig_df = parse_cojo_output(jma_final, locus_id, chrom)
                    if sig_df is not None:
                        all_signals.append(sig_df)

        finally:
            shutil.rmtree(work_dir, ignore_errors=True)
            logger.info(f"  Cleaned up working dir: {work_dir}")

        # Step 6: Assemble outputs
        logger.info("")
        logger.info("Step 6: Assembling outputs ...")

        if not all_signals:
            logger.warning(
                "No independent signals found via COJO. "
                "Check GWAS HT path, p-threshold, and LD reference."
            )
            hl.stop()
            return

        all_signals_df = pd.concat(all_signals, ignore_index=True)
        all_signals_df.to_csv(signals_path, sep='\t', index=False)
        logger.info(
            f"  Independent signals: {len(all_signals_df):,} -> {signals_path}"
        )

        # Parse MHC coords for BED exclusion
        mhc_parts = mhc_interval.replace('chr', '').split(':')
        mhc_chrom = _normalize_contig(mhc_parts[0])
        mhc_bp = mhc_parts[1].split('-') if len(mhc_parts) > 1 else ['0', '0']
        mhc_start_bp = int(mhc_bp[0])
        mhc_end_bp = int(mhc_bp[1]) if len(mhc_bp) > 1 else mhc_start_bp

        bed_df = build_loci_bed(
            all_signals_df, flank,
            mhc_chrom, mhc_start_bp, mhc_end_bp,
        )
        os.makedirs(os.path.dirname(loci_bed_path), exist_ok=True)
        bed_df.to_csv(loci_bed_path, sep='\t', index=False, header=True)
        logger.info(f"  Target loci BED: {len(bed_df):,} loci -> {loci_bed_path}")

    finally:
        try:
            hl.stop()
        except Exception:
            pass

    elapsed = time.time() - run_start
    logger.info("")
    logger.info("=" * 60)
    logger.info(f"Phase 2b complete in {_fmt_elapsed(elapsed)}")
    logger.info(f"  Signals : {signals_path}")
    logger.info(f"  BED     : {loci_bed_path}")
    logger.info(f"  COJO    : {cojo_dir}/")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
