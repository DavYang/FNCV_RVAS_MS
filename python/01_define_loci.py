#!/usr/bin/env python3
"""
Phase 2a: Lead SNP extraction, greedy clumping, and GCTA-COJO signal dissection.

Outputs:
  - results/2-locus_definition/all_independent_signals.tsv
  - results/2-locus_definition/target_loci.bed
  - results/2-locus_definition/cojo/{locus_id}.jma.cojo  (per-locus GCTA output)
"""

import argparse
import os
import shutil
import subprocess
import tempfile
import time
from datetime import datetime
from typing import Dict, List, Optional, Tuple

import hail as hl
import pandas as pd

from utils import load_config, setup_logger

logger = setup_logger("define_loci")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# AllxAll GWAS HT field name candidates (checked in order; first match wins)
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

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Phase 2a: GCTA-COJO locus definition")
    parser.add_argument('--config', default='config/config.json',
                        help='Path to config.json')
    parser.add_argument('--output-dir', default=None,
                        help='Override base output directory')
    return parser.parse_args()

# ---------------------------------------------------------------------------
# Helper utilities
# ---------------------------------------------------------------------------

def _fmt_elapsed(seconds: float) -> str:
    """Format elapsed seconds as human-readable string."""
    if seconds < 60:
        return f"{seconds:.1f}s"
    m, s = divmod(int(seconds), 60)
    return f"{m}m {s:02d}s"


def _resolve_field(ht: hl.Table, canonical: str, candidates: List[str]) -> str:
    """
    Return the first matching field name from candidates present in ht.row.
    Raises ValueError with a diagnostic message if none found.
    """
    row_fields = set(ht.row)
    for name in candidates:
        if name in row_fields:
            logger.info(f"  Field '{canonical}' resolved to '{name}'")
            return name
    raise ValueError(
        f"Could not find field '{canonical}' in GWAS HT. "
        f"Tried: {candidates}. Available fields: {sorted(row_fields)}"
    )


def _normalize_contig(contig: str) -> str:
    """Ensure contig has chr-prefix for GRCh38."""
    if not contig.startswith('chr'):
        return f"chr{contig}"
    return contig


def _run_command(
    cmd: List[str],
    log_path: str,
    label: str,
) -> Tuple[bool, str]:
    """
    Run a subprocess, tee stdout+stderr to log_path, and return (success, stderr_tail).
    """
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
        # Read last 20 lines of log for error context
        tail = ""
        try:
            with open(log_path) as fh:
                lines = fh.readlines()
            tail = "".join(lines[-20:])
        except Exception:
            pass
        logger.error(f"  {label} FAILED (exit {exc.returncode})\n{tail}")
        return False, tail

# ---------------------------------------------------------------------------
# Step 1: Load and normalise GWAS HT
# ---------------------------------------------------------------------------

def load_gwas_ht(
    gwas_path: str,
    p_thresh: float,
    mhc_interval: str,
) -> Tuple[hl.Table, hl.Table, Dict[str, str]]:
    """
    Load GWAS HT, resolve field names, filter p < p_thresh, exclude MHC.

    Returns:
        gwas_ht    - full GWAS table (used for per-window .ma export)
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

    # Detect contig format (chr-prefixed or bare)
    sample_contig = gwas_ht.locus.contig.collect()[0]
    has_chr_prefix = sample_contig.startswith('chr')
    logger.info(f"  Contig format in HT: '{sample_contig}' (chr-prefix: {has_chr_prefix})")

    # Normalize MHC interval to match HT contig format
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
    """
    Greedy positional clumping of significant SNPs into non-overlapping windows.

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
            'chrom': chrom,
            'start': start,
            'end': end,
            'lead_snp_id': row['snp_id'],
            'lead_pos': pos,
            'lead_p': row['p_value'],
        })

        # Mark all SNPs within this window as consumed
        mask = (df['contig'] == chrom) & (df['position'] >= start) & (df['position'] <= end)
        consumed_idx = df[mask].index
        for idx in consumed_idx:
            consumed[idx] = True

    logger.info(f"Greedy clumping: {len(windows)} candidate loci from {len(df):,} significant SNPs")
    return windows

# ---------------------------------------------------------------------------
# Step 3: Export .ma summary stats file (local)
# ---------------------------------------------------------------------------

def export_ma_file(
    gwas_ht: hl.Table,
    field_map: Dict[str, str],
    chrom: str,
    start: int,
    end: int,
    ma_path: str,
) -> bool:
    """
    Filter GWAS HT to window, write GCTA .ma format to a local file.

    GCTA .ma columns: SNP A1 A2 freq b se p N
    """
    logger.info(f"  Exporting .ma file for {chrom}:{start}-{end} -> {ma_path} ...")
    t0 = time.time()

    interval = hl.locus_interval(chrom, start, end, reference_genome='GRCh38',
                                  includes_end=True)
    win_ht = gwas_ht.filter(interval.contains(gwas_ht.locus))

    win_ht = win_ht.annotate(
        SNP=(hl.str(win_ht.locus.contig) + ':' +
             hl.str(win_ht.locus.position) + ':' +
             win_ht.alleles[0] + ':' + win_ht.alleles[1]),
        A1=win_ht.alleles[1],
        A2=win_ht.alleles[0],
        freq=win_ht[field_map['af']],
        b=win_ht[field_map['beta']],
        se=win_ht[field_map['se']],
        p=win_ht[field_map['pval']],
        N=hl.int32(win_ht[field_map['n']]),
    )

    # Export to a temp GCS path then localise, or export directly to local path
    # via pandas (avoids GCS write for ephemeral .ma files)
    df = win_ht.select('SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N').to_pandas()

    if df.empty:
        logger.warning(f"  No variants in window {chrom}:{start}-{end}; skipping.")
        return False

    df.to_csv(ma_path, sep='\t', index=False)
    logger.info(f"  .ma file: {len(df):,} variants ({_fmt_elapsed(time.time() - t0)})")
    return True

# ---------------------------------------------------------------------------
# Step 4: PLINK subset per locus
# ---------------------------------------------------------------------------

def extract_plink_subset(
    bfile: str,
    chrom: str,
    start: int,
    end: int,
    work_dir: str,
    plink_threads: int = 4,
    plink_mem_mb: int = 8000,
) -> Optional[str]:
    """
    Extract a PLINK subset for the locus window (±1 Mb for LD context) from the
    genome-wide LD reference bfile. Returns path to subset bfile prefix, or None
    on failure.
    """
    # Expand window by 500 kb on each side for LD context
    ld_start = max(1, start - 500_000)
    ld_end = end + 500_000
    chrom_num = chrom.replace('chr', '')

    subset_prefix = os.path.join(work_dir, f"subset_{chrom}_{start}_{end}")
    log_path = f"{subset_prefix}_plink.log"

    cmd = [
        'plink2',
        '--bfile', bfile,
        '--chr', chrom_num,
        '--from-bp', str(ld_start),
        '--to-bp', str(ld_end),
        '--make-bed',
        '--out', subset_prefix,
        '--threads', str(plink_threads),
        '--memory', str(plink_mem_mb),
        '--no-psam-pheno',
    ]

    success, _ = _run_command(cmd, log_path, f"plink2 subset {chrom}:{start}-{end}")
    if not success:
        return None

    if not os.path.exists(f"{subset_prefix}.bed"):
        logger.warning(f"  PLINK subset .bed not found: {subset_prefix}.bed")
        return None

    n_vars = sum(1 for _ in open(f"{subset_prefix}.bim"))
    logger.info(f"  PLINK subset: {n_vars:,} variants in LD window")
    return subset_prefix

# ---------------------------------------------------------------------------
# Step 5: Run GCTA-COJO
# ---------------------------------------------------------------------------

def run_gcta_cojo(
    gcta_bin: str,
    bfile: str,
    ma_path: str,
    out_prefix: str,
    chrom: str,
    p_thresh: float,
) -> bool:
    """
    Run GCTA-COJO --cojo-slct on a pre-extracted PLINK subset.
    Stdout/stderr written to {out_prefix}.gcta.log.
    """
    chrom_num = chrom.replace('chr', '')
    log_path = f"{out_prefix}.gcta.log"

    cmd = [
        gcta_bin,
        '--bfile', bfile,
        '--chr', chrom_num,
        '--cojo-file', ma_path,
        '--cojo-slct',
        '--cojo-p', str(p_thresh),
        '--out', out_prefix,
    ]

    success, _ = _run_command(cmd, log_path, f"GCTA-COJO {chrom}")
    return success

# ---------------------------------------------------------------------------
# Step 6: Parse COJO output
# ---------------------------------------------------------------------------

def parse_cojo_output(
    jma_path: str,
    locus_id: str,
    chrom: str,
) -> Optional[pd.DataFrame]:
    """
    Parse GCTA .jma.cojo file. Returns DataFrame with independent signals,
    or None if file missing/empty.
    """
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
# Step 7: Assemble BED from COJO signals
# ---------------------------------------------------------------------------

def build_loci_bed(
    all_signals_df: pd.DataFrame,
    flank: int,
    mhc_chrom: str,
    mhc_start: int,
    mhc_end: int,
) -> pd.DataFrame:
    """
    Build target_loci.bed from COJO independent signals.

    Each signal becomes a locus: signal_pos ± flank, clipped to chr bounds.
    MHC signals are excluded. Columns: chrom, start, end, locus_id.
    """
    rows = []

    for _, sig in all_signals_df.iterrows():
        chrom = sig['chrom']
        # GCTA .jma.cojo column for position is 'bp'
        pos = int(sig.get('bp', sig.get('pos', 0)))
        if pos == 0:
            logger.warning(f"  Could not determine position for signal in {sig.get('locus_id')}; skipping")
            continue

        start = max(1, pos - flank)
        end = pos + flank
        if chrom in CHR_SIZES:
            end = min(end, CHR_SIZES[chrom])

        # MHC exclusion on final signals
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
    """Run Phase 2a: lead SNP extraction + GCTA-COJO signal dissection."""
    args = parse_args()

    run_start = time.time()
    logger.info("=" * 60)
    logger.info("PHASE 2a: LOCUS DEFINITION (GCTA-COJO)")
    logger.info("=" * 60)
    logger.info(f"  Timestamp : {datetime.now().isoformat()}")
    logger.info(f"  Config    : {args.config}")

    config = load_config(args.config)

    # --- Paths ---
    project_dir = os.path.abspath(os.path.dirname(os.path.dirname(args.config)))
    base_dir = args.output_dir or os.path.join(project_dir, config['outputs']['base_dir'])
    locus_dir = os.path.join(base_dir, '2-locus_definition')
    cojo_dir = os.path.join(locus_dir, 'cojo')
    os.makedirs(cojo_dir, exist_ok=True)

    loci_bed_path = os.path.join(locus_dir, 'target_loci.bed')
    signals_path = os.path.join(locus_dir, 'all_independent_signals.tsv')

    gwas_path = config['inputs']['phenotype_gwas']
    ld_reference = os.path.abspath(
        os.path.join(project_dir, config['outputs']['ld_reference'])
    )
    gcta_bin = os.path.abspath(
        os.path.join(project_dir, config['tools']['gcta'])
    )

    p_thresh = float(config['params']['gwas_p_threshold'])
    mhc_interval = config['params']['mhc_interval']
    flank = int(config['params']['locus_flank_window'])

    logger.info(f"  GWAS HT   : {gwas_path}")
    logger.info(f"  LD ref    : {ld_reference}")
    logger.info(f"  GCTA bin  : {gcta_bin}")
    logger.info(f"  P thresh  : {p_thresh}")
    logger.info(f"  MHC excl  : {mhc_interval}")
    logger.info(f"  Flank     : {flank:,} bp")
    logger.info(f"  Output    : {locus_dir}")
    logger.info("")

    # Validate dependencies
    if not os.path.exists(f"{ld_reference}.bed"):
        raise FileNotFoundError(
            f"LD reference bfile not found: {ld_reference}.bed\n"
            "Ensure Phase 1 QC + merge has completed."
        )
    if not os.path.exists(gcta_bin):
        raise FileNotFoundError(
            f"GCTA binary not found: {gcta_bin}\n"
            "Run python/install_gcta.py first."
        )

    # Initialize Hail
    logger.info("Initializing Hail ...")
    hl.init(log='/tmp/hail_define_loci.log')
    hl.default_reference('GRCh38')
    logger.info("Hail initialized")
    logger.info("")

    try:
        # Step 1: Load GWAS
        gwas_ht, sig_ht, field_map = load_gwas_ht(gwas_path, p_thresh, mhc_interval)

        # Detect contig prefix for interval construction
        sample_contig = gwas_ht.locus.contig.collect()[0]
        has_chr_prefix = sample_contig.startswith('chr')

        # Step 2: Greedy clumping
        logger.info("")
        logger.info("Step 2: Greedy clumping ...")
        windows = greedy_clump(sig_ht, field_map['pval'], flank, has_chr_prefix)
        logger.info(f"  -> {len(windows)} analysis windows")

        # Step 3-5: Per-locus COJO loop
        logger.info("")
        logger.info("Step 3: Running GCTA-COJO per locus ...")
        all_signals: List[pd.DataFrame] = []

        # Temp dir for .ma files and PLINK subsets (local, ephemeral)
        work_dir = tempfile.mkdtemp(prefix='cojo_work_')
        logger.info(f"  Working dir: {work_dir}")

        try:
            for i, win in enumerate(windows):
                chrom = win['chrom']
                start = win['start']
                end = win['end']
                locus_id = f"locus_{i:04d}_{chrom}_{start}_{end}"

                logger.info(f"  [{i+1}/{len(windows)}] {locus_id}")

                # Resume: skip if COJO output already exists in cojo_dir
                jma_final = os.path.join(cojo_dir, f"{locus_id}.jma.cojo")
                if os.path.exists(jma_final):
                    logger.info(f"    SKIP (jma.cojo exists)")
                    sig_df = parse_cojo_output(jma_final, locus_id, chrom)
                    if sig_df is not None:
                        all_signals.append(sig_df)
                    continue

                # 3a: Export .ma file (local)
                ma_path = os.path.join(work_dir, f"{locus_id}.ma")
                ok = export_ma_file(gwas_ht, field_map, chrom, start, end, ma_path)
                if not ok:
                    continue

                # 3b: Extract PLINK subset
                subset_bfile = extract_plink_subset(
                    ld_reference, chrom, start, end, work_dir,
                    plink_threads=4,
                    plink_mem_mb=8000,
                )
                if subset_bfile is None:
                    logger.warning(f"    PLINK subset failed for {locus_id}; skipping COJO")
                    continue

                # 3c: Run GCTA-COJO (output to tmp, then copy to cojo_dir)
                cojo_tmp_prefix = os.path.join(work_dir, locus_id)
                success = run_gcta_cojo(
                    gcta_bin, subset_bfile, ma_path,
                    cojo_tmp_prefix, chrom, p_thresh,
                )

                # Copy COJO outputs to persistent cojo_dir regardless of success
                for ext in ['.jma.cojo', '.cma.cojo', '.gcta.log']:
                    src = f"{cojo_tmp_prefix}{ext}"
                    if os.path.exists(src):
                        shutil.copy2(src, os.path.join(cojo_dir, f"{locus_id}{ext}"))

                if not success:
                    logger.warning(f"    GCTA-COJO failed for {locus_id}")
                    continue

                # 3d: Parse output
                sig_df = parse_cojo_output(jma_final, locus_id, chrom)
                if sig_df is not None:
                    all_signals.append(sig_df)

        finally:
            shutil.rmtree(work_dir, ignore_errors=True)
            logger.info(f"  Cleaned up working dir: {work_dir}")

        # Step 4: Assemble outputs
        logger.info("")
        logger.info("Step 4: Assembling outputs ...")

        if not all_signals:
            logger.warning(
                "No independent signals found via COJO. "
                "Check GWAS HT path, p-threshold, and LD reference."
            )
            hl.stop()
            return

        all_signals_df = pd.concat(all_signals, ignore_index=True)
        all_signals_df.to_csv(signals_path, sep='\t', index=False)
        logger.info(f"  Independent signals: {len(all_signals_df):,} -> {signals_path}")

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
    logger.info(f"Phase 2a complete in {_fmt_elapsed(elapsed)}")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()