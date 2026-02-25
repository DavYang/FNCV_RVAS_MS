# 01 - Background SNP PLINK Pipeline

## Overview

This pipeline generates background common-SNP PLINK files (`.bed/.bim/.fam`) for use in
**Regenie Step 1** null model fitting. It samples a configurable number of common SNPs
(~1 million SNPs get downsampled to 500K) across autosomes (chr1-22) from the AoU v8 ACAF splitMT,
filters to EUR-ancestry samples, and exports per-chromosome PLINK binary files to GCS.
Proportional sampling based on chromosome size.
Results are copied and saved both ways, both on the current VM and GCS (virtual machine, google cloud storage).


### Spark/Hail Strategy

The pipeline uses **Hail** (backed by Apache Spark) running on the AoU Researcher Workbench
Dataproc cluster to process the AoU whole-genome MatrixTable. The key architectural decision
is a **per-chromosome serial execution model**:

- A **bash wrapper** (`bash/01_build_background_snps.sh`) loops over chr1-22 and spawns a
  **fresh Python/Hail/JVM process** for each chromosome.
- This avoids JVM state accumulation, metaspace leaks, and Spark DAG complexity that cause
  crashes in long-running multi-chromosome Hail sessions.
- Each chromosome invocation initializes Hail, processes one chromosome end-to-end, writes
  outputs to GCS, then shuts down the JVM cleanly before the next chromosome starts.

Each chromosome is processed in four steps:

1. **Count QC pool** (Step 0) -- Reads only the rows Table (no genotype data) and counts
   the number of variants passing `info.AF` and `variant_qc.call_rate` pre-filters for this
   chromosome. The count is cached to `pool_count.txt` to avoid recomputation on retries or
   re-invocations.
2. **Sample loci** (Step 1) -- Bernoulli-samples from the QC-passing pool at a fraction
   computed from the real pool size and a configurable `oversample_factor`. Uses rows Table
   only (no genotype data scanned).
3. **Extract and count** (Step 2a) -- Reads the full MatrixTable, joins to sampled loci and
   EUR samples, applies a EUR-specific genotype call-rate filter, checkpoints, and counts
   QC-passing variants. Does **not** export PLINK -- this avoids paying the expensive
   checkpoint+export cost on failed attempts.
4. **Export PLINK** (Step 2b) -- Called only once, after the variant count meets the target.
   Exports the checkpointed MatrixTable to PLINK binary format and cleans up the checkpoint.

If the count from Step 2a falls short of the target, the fraction is boosted and Steps 1-2a
are retried (up to `resample_max_attempts` times) before proceeding with available variants.

---

## Input Files

All input paths are defined in `config/config.json`:

| Input | Path | Description |
|-------|------|-------------|
| **ACAF splitMT** | `gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/acaf_threshold/splitMT/hail.mt` | AoU v8 WGS MatrixTable (multi-allelic sites split). Source of all variant and genotype data. |
| **Ancestry predictions** | `gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv` | TSV with `research_id` and `ancestry_pred` columns. Used to filter to EUR-ancestry samples (~234K). |

### Configuration (`config/config.json`)

Key sampling parameters:

```json
{
    "sampling": {
        "target_total_snps": 1000000,
        "random_seed": 42,
        "min_call_rate": 0.95,
        "hail_min_maf": 0.1,
        "hail_min_variant_call_rate": 0.99,
        "resample_max_attempts": 5,
        "oversample_factor": 2.5,
        "avoid_hla": true,
        "hla_region": {
            "chrom": "6",
            "start": 28510120,
            "end": 33480577
        }
    }
}
```

- **`target_total_snps`**: Total SNPs to sample genome-wide. Distributed proportionally by
  chromosome size (GRCh38 bp).
- **`random_seed`**: Base seed for reproducibility. Each chromosome gets `seed + chr_index`;
  retries add `(attempt - 1) * 100` to vary draws.
- **`min_call_rate`**: EUR-specific genotype call-rate threshold applied after extracting
  genotypes in Step 2a (default 0.95).
- **`hail_min_maf`**: Global `info.AF` threshold used as a proxy for EUR MAF in the rows-only
  pre-filter (default 0.1). Variants with pan-ancestry AF > 0.1 are very likely to have
  EUR MAF > 0.01-0.05, which is the downstream PLINK `--maf` target.
- **`hail_min_variant_call_rate`**: Pre-computed `variant_qc.call_rate` threshold in the
  rows-only pre-filter (default 0.99). Applied before any genotype data is read.
- **`resample_max_attempts`**: Maximum number of sample-then-count retry loops before
  proceeding with whatever variants are available (default 5).
- **`oversample_factor`**: Multiplier applied to the target when computing the initial
  Bernoulli sampling fraction: `fraction = min(1.0, target * oversample_factor / pool_size)`.
  Accounts for variants lost during the EUR call-rate filter (default 2.5).
- **`avoid_hla`**: If `true`, excludes the HLA/MHC region on chr6 from sampling.
- **`test_mode`** / **`test_chromosome`** (in `params`): If `test_mode` is `true`, only the
  specified chromosome is processed.

---

## Output Files

### GCS Output Directory Structure

```
gs://<WORKSPACE_BUCKET>/results/FNCV_RVAS_MS/<N>K_background_snps/
|
|-- shared/
|   |-- eur_samples.ht/          # Hail Table of EUR sample IDs (written once, reused)
|
|-- chr1/
|   |-- sampled_loci.ht/         # Hail Table of sampled locus keys for chr1
|   |-- chr1_background.bed      # PLINK binary genotype file
|   |-- chr1_background.bim      # PLINK variant info file
|   |-- chr1_background.fam      # PLINK sample info file
|   |-- pool_count.txt           # Cached QC-passing variant pool count (integer)
|   |-- summary.json             # Per-chromosome run metadata (status, counts, timing)
|
|-- chr2/
|   |-- sampled_loci.ht/
|   |-- chr2_background.bed
|   |-- chr2_background.bim
|   |-- chr2_background.fam
|   |-- pool_count.txt
|   |-- summary.json
|
|-- ...                          # chr3 through chr22 (same structure)
|
|-- tmp/                         # Temporary checkpoint files (cleaned up after each chrom)
|-- pipeline_summary.json        # Overall pipeline run summary
```

### Output File Descriptions

| File | Description |
|------|-------------|
| `chrN_background.bed/bim/fam` | PLINK binary fileset for chromosome N. Contains genotypes for ~N proportional SNPs across ~234K EUR samples. |
| `sampled_loci.ht/` | Hail Table storing the sampled variant loci (row keys only). Enables resume without re-sampling. |
| `pool_count.txt` | Plain-text file caching the QC-passing variant pool count for the chromosome. Reused across retries and re-invocations unless `--force-rerun` is passed. |
| `summary.json` | Per-chromosome JSON with status (`success`/`failed`), variant/sample counts, timing, pool size, oversample factor, number of attempts, and any error tracebacks. |
| `eur_samples.ht/` | Shared Hail Table of EUR sample IDs. Written once by the first chromosome, reused by all subsequent. |
| `pipeline_summary.json` | Overall run summary: succeeded/failed/skipped chromosome lists and total wall time. |

### Local Log Files

```
logs/
|-- 01_background_snps_<TIMESTAMP>.log           # Main wrapper log (all chromosomes)
|-- 01_background_snps_chr1_<TIMESTAMP>.log      # Per-chromosome Python/Spark output
|-- 01_background_snps_chr2_<TIMESTAMP>.log
|-- ...
|-- 01_background_snps.pid                       # PID file (removed on clean exit)
```

Hail JVM logs are written to `/tmp/hail_chrN.log` on the VM (not persisted to GCS).

---

## Instructions for Use

### 1. Configure the AoU Workbench Cloud Environment

Before running the pipeline, configure your cloud environment in the AoU Researcher Workbench:

1. Navigate to your workspace and open **Cloud analysis environment** settings.
2. Set the following:
   - **Compute type**: Dataproc Cluster
   - **CPUs**: 32, **RAM (GB)**: 208
   - **Workers**: 4, **Preemptible workers**: 40
   - **Worker CPUs**: 4, **Worker RAM (GB)**: 15, **Worker Disk (GB)**: 150
   - **Master Disk (GB)**: 500
3. **Set "Automatically pause after idle for" to 5 days** (or longer).
   - The default 30-minute idle timeout will kill the pipeline mid-run because background
     `nohup` processes are not detected as "active" by the Jupyter idle checker.
4. Click **Next** / **Update** to apply changes. This may restart the VM.

### 2. Edit Configuration

Edit `config/config.json` to set your desired parameters:

```bash
vi config/config.json
```

Key fields to review:
- `sampling.target_total_snps` -- Total SNPs to sample (e.g., 100000 or 500000)
- `sampling.avoid_hla` -- Set `true` for phenotypes affected by HLA (e.g., MS)
- `params.test_mode` -- Set `true` to test with a single chromosome first

### 3. Run the Pipeline

```bash
cd ~/FNCV_RVAS_MS
nohup bash bash/01_build_background_snps.sh > /dev/null 2>&1 &
```

### 4. Resume a Previous Run

The pipeline has built-in resume support. The Python script checks `summary.json` for each
chromosome and skips it if `status` is `"success"`. Failed or incomplete chromosomes are
automatically re-run. The cached `pool_count.txt` is also reused to avoid repeating the
pool-counting scan.

If the output directory from a previous run has a different name (e.g., from an older
date-stamped version), use `RESUME_OUTPUT_DIR`:

```bash
RESUME_OUTPUT_DIR="gs://<bucket>/results/FNCV_RVAS_MS/500K_background_snps_20260212" \
  nohup bash bash/01_build_background_snps.sh > /dev/null 2>&1 &
```

### 4b. Force Rerun

To ignore existing `summary.json` (even if `status=success`) and `pool_count.txt`, set
`FORCE_RERUN=true`. This re-counts the pool, re-samples, and re-exports PLINK for every
chromosome:

```bash
FORCE_RERUN=true RESUME_OUTPUT_DIR="gs://<bucket>/results/FNCV_RVAS_MS/500K_background_snps" \
  nohup bash bash/01_build_background_snps.sh > /dev/null 2>&1 &
```

The `FORCE_RERUN` env variable is passed as `--force-rerun` to the Python script.

### 5. Monitor Progress

```bash
# Main pipeline log (shows chromosome-level progress)
tail -f logs/01_background_snps_*.log

# Per-chromosome detailed log
tail -f logs/01_background_snps_chr6_*.log

# Hail/JVM log for the active chromosome
tail -f /tmp/hail_chr6.log

# Check if the wrapper process is still running
cat logs/01_background_snps.pid
ps -p $(cat logs/01_background_snps.pid) -o pid,etime,cmd
```

### 6. Validate Output PLINK Files

After the pipeline completes, verify the PLINK files are valid:

```bash
# Download a chromosome's files from GCS
gsutil -m cp gs://<bucket>/results/FNCV_RVAS_MS/100K_background_snps/chr21/chr21_background.* .

# Basic integrity check
plink --bfile chr21_background --freq --out chr21_check

# Full QC summary
plink --bfile chr21_background --missing --hardy --out chr21_qc
```

### 7. Stop the Pipeline

```bash
kill $(cat ~/FNCV_RVAS_MS/logs/01_background_snps.pid)
```

---

## Pipeline Scripts

| File | Description |
|------|-------------|
| `bash/01_build_background_snps.sh` | Bash wrapper. Loops chr1-22, computes per-chromosome targets, spawns one Python process per chromosome, passes `--force-rerun` when `FORCE_RERUN=true`, handles failure/summary. |
| `python/01_build_background_snps.py` | Python/Hail script. Processes a single chromosome end-to-end: counts QC pool, samples loci, extracts EUR GTs with call-rate filter, retries adaptively, exports PLINK only on success. Supports `--force-rerun` to override cached results. |
| `python/utils.py` | Shared utilities: `load_config()`, `setup_logger()`, `init_hail()`. |
| `config/config.json` | Central configuration for input paths, sampling parameters, and runtime options. |

---

## Performance Notes

- Each chromosome takes approximately **25-40 minutes** depending on size (benchmarked with
  500K target SNPs on the Dataproc configuration above).
- Full chr1-22 run completes in approximately **10-14 hours** (serial, with 60s pauses).
- The pipeline flushes OS page caches and sleeps 60 seconds between chromosomes to mitigate
  memory pressure from sequential JVM startups.
- Peak memory usage occurs during the checkpoint step (Step 2c), which materializes the
  filtered MatrixTable to GCS.

---

## Troubleshooting

| Symptom | Likely Cause | Fix |
|---------|-------------|-----|
| Pipeline stops after 6-7 chromosomes with no error in logs | AoU VM auto-paused due to idle timeout | Set "Automatically pause after idle for" to 5+ days in Cloud Environment settings |
| Resume starts from chr1 instead of where it stopped | Output directory name changed (e.g., date stamp) | Use `RESUME_OUTPUT_DIR` env var pointing to the existing output directory |
| `WORKSPACE_BUCKET environment variable not set` | VM restarted without AoU env vars | Re-open terminal from Jupyter; the env vars are set on container startup |
| `Hail init` fails or Spark errors | Dataproc cluster not running | Check cluster status in AoU workbench; restart if needed |
| chr6 has fewer SNPs than expected | HLA exclusion is enabled | Expected behavior when `avoid_hla: true`; the MHC region is excluded from sampling |
| Variant count below target after all attempts | Pool too small or call-rate filter too strict | Lower `hail_min_maf`, lower `min_call_rate`, or increase `oversample_factor` in config |
| Stale results after config change | Cached `pool_count.txt` or `summary.json` from prior run | Use `FORCE_RERUN=true` to ignore cached counts and re-process all chromosomes |

---

## Changelog

### 2026-02-25 -- Sampling Redesign (v2)

**Architecture changes:**

- **Real pool counting** (Step 0): Added `count_qc_pool()` to count the actual number of
  QC-passing variants per chromosome using the rows Table only (no genotype scan). The count
  is cached to `pool_count.txt` and reused across retries and re-invocations. Replaces the
  prior approach of estimating pool size from a genome-wide variants-per-bp constant.
- **Separated extract from export**: Split the old single `export_plink()` into two functions:
  `extract_and_count()` (Step 2a) checkpoints and counts variants without exporting PLINK,
  and `export_plink_final()` (Step 2b) exports PLINK only after the target is confirmed met.
  This avoids paying the expensive checkpoint+export cost on failed retry attempts.
- **Adaptive retry loop**: Fraction is computed as `min(1.0, target * oversample_factor / pool_size)`
  using the real pool count. On failure, the fraction is boosted proportionally and sampling is
  retried up to `resample_max_attempts` times.
- **`--force-rerun` flag**: Added CLI flag to override resume logic (`summary.json` with
  `status=success`) and cached pool counts. Passed from bash via `FORCE_RERUN=true` env var.
- **Resume logic centralized**: Removed the bash-level `gsutil stat` resume check. The Python
  script's own `summary.json` status check is now the single source of truth for skip/re-run
  decisions.

**Config additions:**

- `sampling.oversample_factor` (default 2.5): Controls initial oversampling to compensate for
  variants lost during EUR call-rate filtering.
- `sampling.min_call_rate` (default 0.95): EUR-specific genotype call-rate threshold.
- `sampling.hail_min_maf` (default 0.1): Global `info.AF` pre-filter (proxy for EUR MAF).
- `sampling.hail_min_variant_call_rate` (default 0.99): Pre-computed variant call-rate threshold.
- `sampling.resample_max_attempts` (default 5): Maximum retry attempts.

**Bash wrapper changes:**

- Added `FORCE_RERUN` env variable support; passes `--force-rerun` to Python.
- Removed bash-level resume skip (was checking `summary.json` via `gsutil stat`);
  resume is now handled entirely by the Python script.
