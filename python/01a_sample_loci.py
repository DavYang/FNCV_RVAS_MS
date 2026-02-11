#!/usr/bin/env python3
"""Sample loci from ACAF splitMT rows Table and export plink2-compatible files.

Phase 1 of the background SNP pipeline for Regenie Step 1. This script:
  1. Reads the ACAF splitMT rows Table (no genotype data).
  2. Bernoulli-samples ~500K loci proportionally across autosomes.
  3. Writes a Hail Table of sampled loci (provenance).
  4. Exports per-chromosome variant range files for plink2 --extract range.
  5. Exports a EUR sample keep-file for plink2 --keep.

Usage:
    python 01a_sample_loci.py --output-dir gs://bucket/path/to/output
"""

import argparse
import json
import os
import sys
import time

import hail as hl
import hailtop.fs as hfs
from datetime import datetime
from utils import load_config, setup_logger

logger = setup_logger("sample_loci")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
AUTOSOMES = [f"chr{i}" for i in range(1, 23)]

CHR_SIZES = {
    'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559,
    'chr4': 190214555, 'chr5': 181538259, 'chr6': 170805979,
    'chr7': 159345973, 'chr8': 145138636, 'chr9': 138394717,
    'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
    'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189,
    'chr16': 90338345, 'chr17': 83257441, 'chr18': 80373285,
    'chr19': 58617616, 'chr20': 64444167, 'chr21': 46709983,
    'chr22': 50818468,
}

ACAF_VARIANTS_PER_BP = 0.03826
SAMPLING_OVERSHOOT = 1.02


# ---------------------------------------------------------------------------
# Chromosome targets
# ---------------------------------------------------------------------------
def calculate_chromosome_targets(target_snps: int) -> dict:
    """Distribute target SNPs proportionally across autosomes (chr1-22).

    Args:
        target_snps: Total number of SNPs to sample genome-wide.

    Returns:
        dict mapping chromosome name to per-chromosome target count.
    """
    total_size = sum(CHR_SIZES.values())
    chr_targets = {}
    for chrom in AUTOSOMES:
        size = CHR_SIZES[chrom]
        chr_targets[chrom] = max(1, int(target_snps * (size / total_size)))

    current_total = sum(chr_targets.values())
    remaining = target_snps - current_total
    if remaining != 0:
        sorted_chroms = sorted(chr_targets, key=lambda c: chr_targets[c], reverse=True)
        for i, chrom in enumerate(sorted_chroms):
            if i >= abs(remaining):
                break
            chr_targets[chrom] += 1 if remaining > 0 else -1

    logger.info(f"Total target: {sum(chr_targets.values()):,} across {len(chr_targets)} autosomes")
    return chr_targets


# ---------------------------------------------------------------------------
# Pass 1: Sample loci from rows Table (no entry data)
# ---------------------------------------------------------------------------
def sample_loci(
    mt_path: str,
    output_ht_path: str,
    chr_targets: dict,
    config: dict,
) -> dict:
    """Sample loci proportionally across autosomes using rows Table only.

    This function never touches entry data. It reads the row keys of the ACAF
    splitMT, applies Bernoulli sampling per chromosome, and writes a small
    Hail Table of selected (locus, alleles) pairs.

    Args:
        mt_path: GCS path to ACAF splitMT.
        output_ht_path: GCS path to write the sampled loci Table.
        chr_targets: Dict of chromosome -> target SNP count.
        config: Loaded config dict.

    Returns:
        Summary dict with per-chromosome counts.
    """
    seed = config['sampling']['random_seed']
    avoid_hla = config['sampling'].get('avoid_hla', False)
    hla_region = config['sampling'].get('hla_region', {})

    logger.info("=" * 60)
    logger.info("PHASE 1: Sampling loci from rows Table (no entry data)")
    logger.info("=" * 60)
    pass1_start = time.time()

    logger.info(f"Reading rows Table from {mt_path} ...")
    mt = hl.read_matrix_table(mt_path)
    ht_rows = mt.rows().select()
    logger.info("Rows Table loaded (no entry data accessed)")

    per_chr_tables = []
    per_chr_counts = {}

    for chrom in AUTOSOMES:
        chr_target = chr_targets.get(chrom, 0)
        if chr_target == 0:
            continue

        chr_start = time.time()

        interval = hl.parse_locus_interval(
            f"{chrom}:1-{CHR_SIZES[chrom]}", reference_genome='GRCh38'
        )
        ht_chr = hl.filter_intervals(ht_rows, [interval])

        if avoid_hla and chrom == 'chr6' and hla_region:
            hla_chrom = f"chr{hla_region['chrom']}"
            hla_start = hla_region['start']
            hla_end = hla_region['end']
            hla_interval = hl.parse_locus_interval(
                f"{hla_chrom}:{hla_start}-{hla_end}", reference_genome='GRCh38'
            )
            ht_chr = ht_chr.filter(~hla_interval.contains(ht_chr.locus))
            logger.info(f"  {chrom}: Excluded HLA region "
                        f"{hla_chrom}:{hla_start}-{hla_end}")

        estimated_pool = int(CHR_SIZES[chrom] * ACAF_VARIANTS_PER_BP)
        fraction = min(1.0, (chr_target / estimated_pool) * SAMPLING_OVERSHOOT)

        chr_seed = seed + AUTOSOMES.index(chrom)
        ht_sampled = ht_chr.filter(hl.rand_unif(0.0, 1.0, seed=chr_seed) < fraction)

        per_chr_tables.append(ht_sampled)
        chr_time = time.time() - chr_start
        logger.info(f"  {chrom}: target={chr_target:,}, "
                    f"est_pool={estimated_pool:,}, "
                    f"fraction={fraction:.6f} ({chr_time:.1f}s)")

    global_target = sum(chr_targets.values())
    logger.info("Combining per-chromosome sampled loci ...")
    ht_combined = per_chr_tables[0]
    for ht in per_chr_tables[1:]:
        ht_combined = ht_combined.union(ht)

    logger.info(f"Writing sampled loci Table to {output_ht_path} ...")
    ht_combined.write(output_ht_path, overwrite=True)

    ht_final = hl.read_table(output_ht_path)
    total_count = ht_final.count()
    logger.info(f"Phase 1 complete: {total_count:,} loci sampled")

    chr_counts = ht_final.group_by(
        contig=ht_final.locus.contig
    ).aggregate(n=hl.agg.count())
    chr_counts_dict = {row.contig: row.n for row in chr_counts.collect()}

    for chrom in AUTOSOMES:
        n = chr_counts_dict.get(chrom, 0)
        target = chr_targets.get(chrom, 0)
        per_chr_counts[chrom] = {'sampled': n, 'target': target}
        if n > 0 or target > 0:
            logger.info(f"  {chrom}: {n:,} sampled (target: {target:,})")

    pass1_time = time.time() - pass1_start
    summary = {
        'phase': 1,
        'total_loci_sampled': total_count,
        'global_target': global_target,
        'per_chromosome': per_chr_counts,
        'output_ht_path': output_ht_path,
        'time_seconds': round(pass1_time, 1),
        'timestamp': datetime.now().isoformat(),
    }

    logger.info(f"Phase 1 done: {total_count:,} loci in {pass1_time:.1f}s")
    return summary


# ---------------------------------------------------------------------------
# Export plink2-compatible files
# ---------------------------------------------------------------------------
def export_plink2_inputs(
    loci_ht_path: str,
    output_dir: str,
    chr_targets: dict,
    config: dict,
) -> None:
    """Export variant range files and EUR keep-file for plink2.

    Reads the sampled loci Hail Table and writes:
      - Per-chromosome variant range files (chr, pos, pos) for --extract range
      - EUR sample keep-file (FID, IID) for --keep

    Args:
        loci_ht_path: GCS path to sampled loci Hail Table.
        output_dir: GCS base directory for output.
        chr_targets: Dict of chromosome -> target SNP count.
        config: Loaded config dict.
    """
    logger.info("=" * 60)
    logger.info("Exporting plink2-compatible input files")
    logger.info("=" * 60)

    plink2_dir = f"{output_dir}/plink2_inputs"

    # --- Export per-chromosome variant range files ---
    logger.info("Reading sampled loci Table ...")
    ht = hl.read_table(loci_ht_path)

    # Annotate with contig (both chr-prefixed and bare) and position for export
    ht = ht.annotate(
        contig_chr=ht.locus.contig,
        contig_bare=ht.locus.contig.replace('chr', ''),
        pos=hl.str(ht.locus.position),
    )

    for chrom in AUTOSOMES:
        if chr_targets.get(chrom, 0) == 0:
            continue

        interval = hl.parse_locus_interval(
            f"{chrom}:1-{CHR_SIZES[chrom]}", reference_genome='GRCh38'
        )
        ht_chr = hl.filter_intervals(ht, [interval])

        # plink2 --extract range format: chr pos_start pos_end (tab-delimited, no header)
        # Export with chr prefix (GRCh38 convention used by AoU PLINK files)
        range_path = f"{plink2_dir}/{chrom}_variants.range"
        ht_export = ht_chr.select(ht_chr.contig_chr, start=ht_chr.pos, end=ht_chr.pos)
        ht_export = ht_export.key_by()
        ht_export.export(range_path, header=False)

        # Also export without chr prefix (fallback if .bim uses bare numbers)
        range_path_bare = f"{plink2_dir}/{chrom}_variants_bare.range"
        ht_export_bare = ht_chr.select(ht_chr.contig_bare, start=ht_chr.pos, end=ht_chr.pos)
        ht_export_bare = ht_export_bare.key_by()
        ht_export_bare.export(range_path_bare, header=False)

        logger.info(f"  Exported {chrom} variant range files to {range_path}")

    # --- Export EUR sample keep-file ---
    logger.info("Exporting EUR sample keep-file ...")
    ancestry_path = config['inputs']['ancestry_pred']
    ancestry_ht = hl.import_table(
        ancestry_path,
        impute=True,
        types={'research_id': hl.tstr, 'ancestry_pred': hl.tstr},
    )
    eur_ht = ancestry_ht.filter(ancestry_ht.ancestry_pred == 'eur')
    # plink2 --keep format: FID IID (tab-delimited, no header)
    eur_keep = eur_ht.select(FID=eur_ht.research_id, IID=eur_ht.research_id)
    eur_keep = eur_keep.key_by()
    keep_path = f"{plink2_dir}/eur_samples.keep"
    eur_keep.export(keep_path, header=False)
    n_eur = eur_ht.count()
    logger.info(f"  Exported {n_eur:,} EUR samples to {keep_path}")

    logger.info("plink2 input files exported successfully")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Sample loci from ACAF splitMT and export plink2 input files"
    )
    parser.add_argument(
        '--output-dir', required=True,
        help='GCS output directory (set by bash wrapper)'
    )
    args = parser.parse_args()
    output_dir = args.output_dir

    config = load_config("config/config.json")

    hl.init(
        log='/tmp/hail_sample_loci.log',
        spark_conf={
            'spark.driver.memory': '12g',
            'spark.executor.memory': '8g',
            'spark.network.timeout': '600s',
            'spark.executor.heartbeatInterval': '120s',
            'spark.driver.maxResultSize': '4g',
        },
    )
    hl.default_reference('GRCh38')

    mt_path = config['inputs']['wgs_matrix_table']
    target_snps = config['sampling']['target_total_snps']
    test_mode = config['params'].get('test_mode', False)
    test_chromosome = config['params'].get('test_chromosome', None)

    logger.info(f"Source MT: {mt_path}")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Total target: {target_snps:,} across 22 autosomes")

    chr_targets = calculate_chromosome_targets(target_snps)

    if test_mode and test_chromosome:
        logger.info(f"=== TEST MODE: {test_chromosome} only ===")
        if test_chromosome not in chr_targets:
            logger.error(f"{test_chromosome} not in autosome targets")
            sys.exit(1)
        chr_targets = {test_chromosome: chr_targets[test_chromosome]}
        logger.info(f"Target for {test_chromosome}: {chr_targets[test_chromosome]:,}")
    else:
        logger.info("=== PRODUCTION MODE: chr1-22 ===")

    loci_ht_path = f"{output_dir}/sampled_loci.ht"

    try:
        # --- Phase 1: Sample loci ---
        if hfs.exists(f"{loci_ht_path}/_SUCCESS") or hfs.exists(f"{loci_ht_path}/metadata.json.gz"):
            logger.info(f"Sampled loci Table already exists at {loci_ht_path}, skipping sampling")
            pass1_summary = {'skipped': True, 'output_ht_path': loci_ht_path}
        else:
            pass1_summary = sample_loci(mt_path, loci_ht_path, chr_targets, config)

        # Write Phase 1 summary
        p1_summary_path = f"{output_dir}/pass1_summary.json"
        with hl.hadoop_open(p1_summary_path, 'w') as f:
            json.dump(pass1_summary, f, indent=2, default=str)
        logger.info(f"Phase 1 summary written to {p1_summary_path}")

        # --- Export plink2 input files ---
        export_plink2_inputs(loci_ht_path, output_dir, chr_targets, config)

        logger.info("=" * 60)
        logger.info("Phase 1 complete. Hail will shut down.")
        logger.info(f"Sampled loci: {loci_ht_path}")
        logger.info(f"plink2 inputs: {output_dir}/plink2_inputs/")
        logger.info("=" * 60)

    except Exception as e:
        logger.error(f"Pipeline failed: {e}")
        raise


if __name__ == "__main__":
    main()
