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

Each chromosome is processed in two steps:

1. **Sample loci** (Step 1) -- Bernoulli-samples from the QC-passing rows Table at a fraction
   computed from an empirical pool-size estimate and a configurable `oversample_factor`.
   Rows-only scan -- no genotype data is read. Pre-filters by `info.AF` (global AF > 0.1 as
   EUR MAF proxy) and `variant_qc.call_rate` (>= 0.99).
2. **Export PLINK** (Step 2) -- Reads the full MatrixTable, joins to sampled loci and EUR
   samples, applies a EUR-specific genotype call-rate filter, coalesces partitions, and
   exports PLINK **directly from the Spark DAG** (no intermediate checkpoint). The variant
   count is read from the `.bim` file after export.

With `oversample_factor >= 3`, the target is always met on the first attempt. Empirically,
the EUR call-rate filter drops <0.1% of pre-filtered variants, so no retry loop is needed.

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
        "oversample_factor": 3,
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
- **`oversample_factor`**: Multiplier applied to the target when computing the Bernoulli
  sampling fraction: `fraction = min(1.0, target * oversample_factor / estimated_pool)`.
  The estimated pool is derived from chromosome size (empirical ~1.75 QC-passing variants
  per kbp). With a value of 3, the target is always met on the first attempt (default 3).
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
|   |-- summary.json             # Per-chromosome run metadata (status, counts, timing)
|
|-- chr2/
|   |-- sampled_loci.ht/
|   |-- chr2_background.bed
|   |-- chr2_background.bim
|   |-- chr2_background.fam
|   |-- summary.json
|
|-- ...                          # chr3 through chr22 (same structure)
|
|-- pipeline_summary.json        # Overall pipeline run summary
```

### Output File Descriptions

| File | Description |
|------|-------------|
| `chrN_background.bed/bim/fam` | PLINK binary fileset for chromosome N. Contains genotypes for ~N proportional SNPs across ~234K EUR samples. |
| `sampled_loci.ht/` | Hail Table storing the sampled variant loci (row keys only). Enables resume without re-sampling. |
| `summary.json` | Per-chromosome JSON with status (`success`/`failed`), variant/sample counts, timing, oversample factor, fraction, and any error tracebacks. |
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
automatically re-run.

If the output directory from a previous run has a different name (e.g., from an older
date-stamped version), use `RESUME_OUTPUT_DIR`:

```bash
RESUME_OUTPUT_DIR="gs://<bucket>/results/FNCV_RVAS_MS/500K_background_snps_20260212" \
  nohup bash bash/01_build_background_snps.sh > /dev/null 2>&1 &
```

### 4b. Force Rerun

To ignore existing `summary.json` (even if `status=success`), set `FORCE_RERUN=true`.
This re-samples loci and re-exports PLINK for every chromosome:

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
| `python/01_build_background_snps.py` | Python/Hail script. Processes a single chromosome end-to-end: samples loci from rows Table, then builds a filtered DAG (EUR join + call-rate filter) and exports PLINK directly (no checkpoint). Supports `--force-rerun` to override cached results. |
| `python/utils.py` | Shared utilities: `load_config()`, `setup_logger()`, `init_hail()`. |
| `config/config.json` | Central configuration for input paths, sampling parameters, and runtime options. |

---

## Performance Notes

- Each chromosome takes approximately **30-35 minutes** depending on size (~17 min sampling +
  ~15 min export). Benchmarked with 1M target SNPs on the Dataproc configuration above.
- Full chr1-22 run completes in approximately **12-14 hours** (serial, with 60s pauses).
- The pipeline sleeps 60 seconds between chromosomes to allow OS memory reclamation from
  terminated JVM processes.
- Peak memory usage occurs during the PLINK export step, which evaluates the full
  filter+coalesce DAG and writes the binary genotype data to GCS.

---

## Troubleshooting

| Symptom | Likely Cause | Fix |
|---------|-------------|-----|
| Pipeline stops after 6-7 chromosomes with no error in logs | AoU VM auto-paused due to idle timeout | Set "Automatically pause after idle for" to 5+ days in Cloud Environment settings |
| Resume starts from chr1 instead of where it stopped | Output directory name changed (e.g., date stamp) | Use `RESUME_OUTPUT_DIR` env var pointing to the existing output directory |
| `WORKSPACE_BUCKET environment variable not set` | VM restarted without AoU env vars | Re-open terminal from Jupyter; the env vars are set on container startup |
| `Hail init` fails or Spark errors | Dataproc cluster not running | Check cluster status in AoU workbench; restart if needed |
| chr6 has fewer SNPs than expected | HLA exclusion is enabled | Expected behavior when `avoid_hla: true`; the MHC region is excluded from sampling |
| Variant count below target | `oversample_factor` too low or pre-filters too strict | Increase `oversample_factor` in config (3 is validated), or lower `hail_min_maf` / `min_call_rate` |
| Stale results after config change | `summary.json` from prior run shows `status=success` | Use `FORCE_RERUN=true` to re-process all chromosomes |

---

## Changelog

### 2026-02-26 -- Speed Optimization (v3)

**Architecture changes:**

- **Removed pool count step** (was Step 0): The `count_qc_pool()` function added ~2 min/chrom
  scanning the same rows Table that `sample_loci()` already scans. With `oversample_factor=3`,
  the target is always met on the first attempt, so knowing the exact pool size is unnecessary.
  The sampling fraction is now computed from an empirical estimate (~1.75 QC-passing variants
  per kbp of chromosome length).
- **Removed checkpoint**: The intermediate `checkpoint()` + `count_rows()` pattern in
  `extract_and_count()` added ~28 min/chrom to materialize the filtered MT to GCS. Since the
  retry loop is no longer needed (oversample_factor ensures first-attempt success), the full
  DAG (interval filter -> semi_join -> call-rate filter -> coalesce) is fed directly to
  `hl.export_plink()`. Variant count is read from the `.bim` file after export.
- **Removed retry loop**: With `oversample_factor=3`, chr9 yielded 144K variants against a
  48K target (3x overshoot) with zero variants lost to the EUR call-rate filter. The retry
  loop, fraction boosting, and `resample_max_attempts` config are no longer needed.
- **Merged extract_and_count + export_plink_final**: Replaced with a single `export_plink()`
  function that builds the DAG and exports in one pass.

**Config changes:**

- `sampling.oversample_factor`: Raised from 2 to 3 (empirically validated).
- `sampling.resample_max_attempts`: Removed (no longer used).

**Estimated speedup:** ~30 min saved per chromosome (from ~62 min to ~32 min). Full 22-chrom
run drops from ~23 hours to ~12 hours.

### 2026-02-25 -- Sampling Redesign (v2)

**Architecture changes:**

- **Real pool counting** (Step 0): Added `count_qc_pool()` to count the actual number of
  QC-passing variants per chromosome using the rows Table only (no genotype scan). The count
  was cached to `pool_count.txt` and reused across retries.
- **Separated extract from export**: Split `export_plink()` into `extract_and_count()` and
  `export_plink_final()` to avoid paying the export cost on failed retry attempts.
- **Adaptive retry loop**: Fraction boosted proportionally on failure, retried up to
  `resample_max_attempts` times.
- **`--force-rerun` flag**: Added CLI flag to override resume logic.
- **Resume logic centralized**: Python script's `summary.json` status check is the single
  source of truth.

**Config additions:**

- `sampling.oversample_factor`, `sampling.min_call_rate`, `sampling.hail_min_maf`,
  `sampling.hail_min_variant_call_rate`, `sampling.resample_max_attempts`.
