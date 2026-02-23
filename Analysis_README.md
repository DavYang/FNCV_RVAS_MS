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

**Goal:** Define the exact genomic windows (GWAS Loci) where the enrichment test will occur by generating independent loci for testing, defined by 95% credible sets.

* **Note:** For Multiple Sclerosis, it is standard practice to exclude the MHC region (Chr6: 25-35Mb) from the general pipeline due to complex LD structures. This analysis will exclude this region.
* **Methods:**
* **Lead SNP Selection:** Identify significant hits (p < 5e-8) from `NS_326.1` summary statistics.
* **Locus Definition:** * *GCTA-COJO*: Used for conditional analysis to find secondary independent SNPs.
* *SuSiE Fine-Mapping*: Used to narrow the window to the 95% Credible Set.




* **Inputs:** AoU GWAS Summary Statistics + Per-chromosome PLINK files (as LD reference).
* **Outputs:** `target_loci.bed` (Chromosome, Start, End, Locus_ID).
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