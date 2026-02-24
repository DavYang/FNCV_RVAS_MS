# 01 - Background SNP PLINK Pipeline

## Overview

This pipeline generates background common-SNP PLINK files (`.bed/.bim/.fam`) for use in
**Regenie Step 1** null model fitting. It samples a configurable number of common SNPs
(d~1 million SNPs get downsampled to 500K) across autosomes (chr1-22) from the AoU v8 ACAF splitMT,
filters to EUR-ancestry samples, and exports per-chromosome PLINK binary files to GCS.
Results are copied and saved both ways, both on the current VM and GCS (virtual machine, google cloud storage)


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

1. **Sample loci** -- Reads only the rows Table (no entry/genotype data) of the ACAF splitMT,
   filters to the target chromosome, applies Bernoulli sampling to select the target number
   of SNPs, and writes a sampled loci Hail Table.
2. **Export PLINK** -- Reads the full MatrixTable, filters to the sampled loci and EUR samples
   via `semi_join_rows`/`semi_join_cols`, checkpoints the filtered MT to break the Spark DAG,
   then calls `hl.export_plink()`.

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
        "target_total_snps": 100000,
        "random_seed": 42,
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
- **`random_seed`**: Base seed for reproducibility. Each chromosome gets `seed + chr_index`.
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
|-- tmp/                         # Temporary checkpoint files (cleaned up after each chrom)
|-- pipeline_summary.json        # Overall pipeline run summary
```

### Output File Descriptions

| File | Description |
|------|-------------|
| `chrN_background.bed/bim/fam` | PLINK binary fileset for chromosome N. Contains genotypes for ~N proportional SNPs across ~234K EUR samples. |
| `sampled_loci.ht/` | Hail Table storing the sampled variant loci (row keys only). Enables resume without re-sampling. |
| `summary.json` | Per-chromosome JSON with status (`success`/`failed`), variant/sample counts, timing, and any error tracebacks. |
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

The pipeline has built-in resume support. It checks for `summary.json` in each chromosome's
output directory and skips completed chromosomes.

If the output directory from a previous run has a different name (e.g., from an older
date-stamped version), use `RESUME_OUTPUT_DIR`:

```bash
RESUME_OUTPUT_DIR="gs://<bucket>/results/FNCV_RVAS_MS/500K_background_snps_20260212" \
  nohup bash bash/01_build_background_snps.sh > /dev/null 2>&1 &
```

Sub-step resume is also supported within each chromosome: if `sampled_loci.ht` or the `.bed`
file already exists, the corresponding step is skipped.

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
| `bash/01_build_background_snps.sh` | Bash wrapper. Loops chr1-22, computes per-chromosome targets, spawns one Python process per chromosome, handles resume/failure/summary. |
| `python/01_build_background_snps.py` | Python/Hail script. Processes a single chromosome: samples loci from the ACAF rows Table, filters MT to sampled loci + EUR samples, checkpoints, exports PLINK. |
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
