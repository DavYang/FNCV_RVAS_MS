# Enrichment Testing of AlleleStacker Candidate Variants for Blood-Trait GWAS Loci (EUR-only)

## Phenotype: Multiple Sclerosis

## Overview

This analysis evaluates whether candidate **functional non-coding variants (FNCVs)** are enriched for genetic loci associated with blood-related traits, compared to **non-functional non-coding variants (nFNCVs)**, using European-ancestry GWAS results from the All of Us (AoU) AllxAll resource.

The study tests if variants with specific methylation signatures (FNCVs) drive disease association at GWAS loci compared to matched background variants (nFNCVs). It begins with a single phenotype as a first pass (Multiple Sclerosis: NS_326.1) and will subsequently expand to multiple phenotypes.

---

## Software Stack

* **Hail:** Distributed data wrangling (extracting common SNPs from AoU MatrixTables).
* **PLINK (1.9/2.0):** Quality control, variant filtering, LD pruning, and fileset merging.
* **GCTA-COJO:** Conditional and joint analysis to isolate independent association signals.
* **SuSiE:** Bayesian fine-mapping to define 95% Credible Sets for target loci boundaries.
* **REGENIE:** Whole genome regression for the null model (Step 1) and SKAT-O aggregation testing (Step 2).

---

## Analysis Architecture

The pipeline is structured sequentially. The outputs of each phase directly dictate the search space, statistical background, or group masking for the subsequent phases.

### Phase 0: Phenotype & Covariate File Generation

**Goal:** Generate the strictly formatted phenotype and covariate text files required by REGENIE to correct for confounding factors like age, sex, and population stratification.

* **Methods:** * Query the AoU Curated Data Repository (CDR) to define NS_326.1 cases and controls.
* Extract primary demographic covariates (Age, Sex) and the top 10 Ancestry Principal Components (PCs) from the AoU ancestry pipeline.
* Format these into space-delimited text files with missing values encoded as `-9`.


* **Inputs:** AoU CDR (EHR/Survey data tables) + AoU Ancestry/PCA prediction tables.
* **Outputs:** `MS_phenotype.txt` and `MS_covariates.txt` (Columns must start with `FID` and `IID` to match PLINK `.fam` files).
* **Integration:** These files are mandatory inputs for Phase 1 to build the genetic background model and Phase 4 to align the rare variant genotypes with disease status and correct for confounding variables.

### Phase 1: Background (Null Model Construction)

**Goal:** Establish the genomic background model to correct for population structure and cryptic relatedness.

* **Status:** Common EUR background SNPs have been extracted via Hail, merged, and QC'ed into a single PLINK fileset saved locally, ready for bucket backup.
* **Methods:** Run REGENIE Step 1 (Whole Genome Regression) using the Leave-One-Chromosome-Out (LOCO) approach.
* **Inputs:** Merged PLINK files (`.bed`, `.bim`, `.fam`) + `MS_phenotype.txt` + `MS_covariates.txt`.
* **Outputs:** `MS_null_model_pred.list` and `.loco` files.
* **Integration:** These predictions act as the baseline null model for the Phase 4 burden tests, preventing false positives driven by ancestry or relatedness.

### Phase 2: Signal Dissection & Locus Definition

**Goal:** Define the exact genomic windows (GWAS Loci) where the enrichment test will occur by isolating independent association signals via conditional analysis.

* **Note:** For Multiple Sclerosis, it is standard practice to exclude the MHC region (Chr6: 25-35Mb) due to complex LD. This analysis excludes the MHC throughout.

Phase 2 is split into three sub-steps, each processed **per chromosome** for memory efficiency (recommended by GCTA authors):

#### Step 02: LD Reference QC (`02_qc_ld_reference.sh`)

Build per-chromosome LD reference panels from the raw PLINK files exported by Hail. Samples are already EUR-only (234K) from Phase 1 Hail export.

* **QC pipeline (PLINK2):**
  1. Variant QC: `--maf 0.01 --hwe 1e-6 --geno 0.05 --mind 0.05 --snps-only just-acgt --max-alleles 2`
  2. Remove duplicate variant IDs: `--rm-dup exclude-all`
* **Input:** `results/1-bg_snp/plink_no-qc/chrN_background.{bed,bim,fam}`
* **Output:** `results/2-locus_definition/ld_ref/chrN_ld_ref.{bed,bim,fam}` -> GCS

#### Step 02a: GWAS Summary Stats Export (`02a_export_gwas_ma.sh`)

Export full-chromosome `.ma` summary statistics from the AllxAll GWAS Hail Table. GCTA-COJO requires **all SNPs** per chromosome (not just per-window) to correctly estimate phenotypic variance.

* **Tool:** Hail (available at `/opt/conda` on AoU VM)
* **Input:** AllxAll GWAS HT (`gs://fc-aou-datasets-controlled/AllxAll/...`)
* **Output:** `results/2-locus_definition/ma/chrN.ma` (22 files) -> GCS
* `.ma` format: `SNP A1 A2 freq b se p N` (tab-separated)

#### Step 02b: GCTA-COJO Locus Definition (`02b_run_cojo.sh`)

Run conditional/joint analysis per chromosome to identify independent signals.

* **Methods:**
  1. **Lead SNP selection:** Identify significant hits (p < 5e-8) from GWAS HT, exclude MHC.
  2. **Greedy clumping:** Group hits into non-overlapping windows (lead SNP +/- 500kb flank).
  3. **Per-locus COJO:** For each locus, subset LD reference to window (+500kb LD padding), run `--cojo-slct` with `--cojo-p 5e-8`.
  4. **Assemble outputs:** Merge independent signals into BED file.
* *SuSiE Fine-Mapping*: Future step to narrow windows to 95% Credible Sets.
* **Inputs:** Per-chrom LD ref PLINK + per-chrom `.ma` files + GCTA binary
* **Outputs:**
  - `results/2-locus_definition/target_loci.bed` (chrom, start, end, locus_id)
  - `results/2-locus_definition/all_independent_signals.tsv`
  - `results/2-locus_definition/cojo/chrN/{locus}.jma.cojo`
* **Integration:** These coordinates restrict Phase 3 so that variant classification only occurs within strictly defined, disease-associated regions.

### Phase 3: Methylation-Based Variant Classification

**Goal:** Sort pre-generated rare variants within the defined GWAS loci into "Functional" (Query) and "Non-Functional" (Control) bins.

* **Selection Logic:**
* **FNCV (Functional Query):** Candidate FNCVs where variant alleles are exclusively observed as either hypermethylated/unmethylated AND are methylation outliers.
* **nFNCV (Matched Control):** Matched set of variants that are *not* methylation outliers.
* **Matching:** nFNCVs are selected to match FNCVs on MAF, GC content, and proximity to TSS.


* **Inputs:** FNCV/nFNCV variant lists + `target_loci.bed`.
* **Outputs:** `group_file.txt` (Regenie format `.set` and `.anno` files mapping Variant -> Set -> Locus).
* **Integration:** This file defines the specific "Masks" that REGENIE will evaluate in Phase 4.

### Phase 4: Statistical Testing & Enrichment

**Goal:** Test if the FNCV group carries more disease risk than the nFNCV group.

* **Methods:**
* **Regenie Step 2 (SKAT-O):** Run two specific masks per locus:
1. *Mask A*: Methylation-Defined Variants (FNCVs)
2. *Mask B*: Matched Background Variants (nFNCVs)


* **Enrichment Metric:** Compare the P-values `(-log10P_FNCV) - (-log10P_nFNCV)`.


* **Inputs:** Phase 1 Null Model (`.loco`), `MS_phenotype.txt`, `MS_covariates.txt`, WGS BGEN files (filtered to Loci), and Phase 3 Group Files.
* **Outputs:** Final enrichment delta and empirical significance calculated via label-shuffling permutation tests.


#### Order of Steps

| Step | Script | Description | Depends on |
|------|--------|-------------|------------|
| 00 | `00_run_phenotype_covariates.sh` | Phenotype + covariate files -> GCS | CDR access |
| 01a | `01a_build_background_snps.sh` | Per-chrom PLINK export via Hail -> GCS | Hail/Dataproc |
| 01b | `01b_qc_merge_background_snps.sh` | QC, LD-prune, merge, thin to 500K -> GCS | 01a |
| 01c | `01c_run_regenie_step1.sh` | REGENIE Step 1 null model (LOCO) -> GCS | 00, 01b |
| 02 | `02_qc_ld_reference.sh` | Per-chrom LD ref QC for GCTA-COJO -> GCS | 01a |
| 02a | `02a_export_gwas_ma.sh` | Per-chrom .ma summary stats via Hail -> GCS | CDR GWAS HT |
| 02b | `02b_run_cojo.sh` | GCTA-COJO per chromosome -> target_loci.bed -> GCS | 02, 02a |

Steps 02 and 02a can run in parallel (independent inputs). Step 02b depends on both.
Steps 01b-01c and 02-02b are independent branches and can run in parallel.

Each script uploads its outputs to GCS and can pull inputs from GCS if missing locally.

#### VM Crash Recovery

All intermediate and final outputs are backed up to GCS. After a VM crash:

```bash
# Restore phenotype/covariate files from GCS (no regeneration)
bash bash/00_run_phenotype_covariates.sh --restore

# 01b auto-downloads raw PLINK files from GCS if missing locally
nohup bash bash/01b_qc_merge_background_snps.sh > /dev/null 2>&1 &

# 01c auto-downloads phenotype + PLINK step1 files from GCS if missing
nohup bash bash/01c_run_regenie_step1.sh > /dev/null 2>&1 &

# 02_qc auto-downloads raw PLINK from GCS if missing
nohup bash bash/02_qc_ld_reference.sh > /dev/null 2>&1 &

# 02b auto-downloads LD ref + .ma files from GCS if missing
nohup bash bash/02b_run_cojo.sh > /dev/null 2>&1 &
```

All scripts support `--force` to re-run even if outputs already exist.
