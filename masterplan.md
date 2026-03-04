# Master Analysis Plan: AlleleStacker Enrichment Testing (EUR-only)

## Phenotype: Multiple Sclerosis (NS_326.1)

## 1. Overview

This analysis evaluates whether candidate **functional non-coding variants (FNCVs)** are enriched for genetic loci associated with Multiple Sclerosis (MS), compared to matched **non-functional non-coding variants (nFNCVs)**. The study leverages the *All of Us* (AoU) v7 European-ancestry cohort to perform a controlled Rare Variant Association Study (RVAS).

The central hypothesis is that rare variants exhibiting methylation outlier signatures (FNCVs) carry a higher burden of disease risk within established GWAS loci than strictly matched background variants (nFNCVs).

---

## 2. Software Stack & Toolchain

The pipeline utilizes a distributed cloud architecture (Hail, Google Cloud Services) integrating the following tools:

* **Hail (Python/Spark):** Used for genomic data wrangling, specifically extracting common background SNPs and and EUR sample genotypes from AoU MatrixTables (mt).
* **PLINK (1.9/2.0):** Used for Quality Control (QC), Linkage Disequilibrium (LD) pruning, and merging of chromosome-level filesets.
* **GCTA-COJO:** Performs Conditional and Joint Analysis to identify independent association signals within complex loci.
* **SuSiE (R):** Performs Bayesian Fine-Mapping to define the 95% Credible Sets that serve as the boundaries for the target loci.
* **REGENIE:** The primary statistical engine. Step 1 fits the whole-genome null model; Step 2 performs the SKAT-O aggregation tests.

---

## 3. Analysis Architecture

The pipeline is designed sequentially. The output of each phase forms the direct input or statistical foundation for the next.

### Phase 0: Phenotype & Covariate Engineering

**Goal:** Generate the strictly formatted phenotype and covariate files required by REGENIE to correct for confounding factors (Age, Sex, Ancestry).

* **Methodology:**
1. **Case/Control Definition:** Query the AoU Curated Data Repository (CDR) to define `NS_326.1` cases (via ICD codes/EHR) and clean controls.
2. **Covariate Extraction:** Join the phenotype data with the AoU Ancestry tables. Extract Age, Sex at Birth, and the top 10 Principal Components (PC1–PC10).
3. **Formatting:** Convert data into space-delimited text files. Missing values are encoded as `-9`. The first two columns must be `FID` and `IID` to match the genetic data.


* **Inputs:** AoU CDR Tables, AoU Ancestry Predictions.
* **Outputs:**
* `MS_phenotype.txt`: Contains binary case/control status.
* `MS_covariates.txt`: Contains Age, Sex, PC1..PC10.

* **Integration:** These files are inputs for **Phase 1** (to build the null model) and **Phase 4** (to align genotypes with clinical outcomes).

### Phase 1: Genomic Background (Null Model Construction)

**Goal:** Establish the genomic background model to correct for population structure and cryptic relatedness using common variants.

* **Status:** ~1,000,000 common EUR variants have been extracted via Hail, QC'ed, and merged into a single PLINK fileset (`merged_qc_common_snps`).
* **Methodology:**
1. **LD Pruning:** Prune the common variants to independent markers ($r^2 < 0.1$).
2. **MHC Exclusion:** Explicitly exclude the MHC region (Chr6: 25Mb–35Mb) to prevent model overfitting due to the massive MS signal in this region.
3. **REGENIE Step 1:** Run Whole Genome Regression using the Leave-One-Chromosome-Out (LOCO) scheme.


* **Inputs:**
* `merged_qc_common_snps.{bed,bim,fam}`
* `MS_phenotype.txt`
* `MS_covariates.txt`


* **Outputs:**
* `MS_null_model_pred.list`: List of prediction files.
* `*.loco` files: Chromosome-specific predictions.


* **Integration:** These LOCO predictions act as the baseline null model for **Phase 4**, ensuring that significant results are due to rare variant burden and not general relatedness.

### Phase 2: Signal Dissection & Locus Definition

**Goal:** Define the exact genomic windows (Target Loci) where the enrichment test will occur.

* **Methodology:**
1. **Lead SNP Selection:** Identify all significant hits ($p < 5 \times 10^{-8}$) from `NS_326.1` summary statistics. Threshold reduced to ($p < 5 \times 10^{-6}$) after reviw of the data. 
2. **Signal Separation (GCTA-COJO):** Use the 22 per-chromosome PLINK files as an LD reference to identify all independent secondary signals within GWAS peaks.
3. **Fine-Mapping (SuSiE):** For each independent signal, run SuSiE to identify the 95% Credible Set.
4. **Anchoring:** Define the "Locus Window" as the minimum/maximum position of the Credible Set, plus 1kb padding.
5. **MHC Exclusion:** Remove any loci falling within the Chr6 MHC region.


* **Inputs:**
* AoU GWAS Summary Statistics.
* Per-chromosome PLINK files (LD Reference, dense panel generated from re-QCing the sampled SNPs from the Hail mt).


* **Outputs:** `target_loci.bed` (Columns: Chrom, Start, End, Locus_ID).
* **Integration:** These coordinates restrict **Phase 3** so that variant classification only occurs within biologically relevant, disease-associated regions.

### Phase 3: Methylation-Based Variant Classification

**Goal:** Sort pre-generated rare variants within the defined GWAS loci into "Functional" (Query) and "Non-Functional" (Control) bins.

* **Methodology:**
1. **FNCV Selection:** Select rare variants ($MAF < 0.1\%$) within `target_loci.bed` that are AlleleStacker methylation outliers.
2. **nFNCV Selection (Matching):** Select rare variants within `target_loci.bed` that are *not* methylation outliers. These are strictly matched to FNCVs based on:
* Minor Allele Frequency (MAF).
* Local GC content.
* Distance to Transcription Start Site (TSS).


3. **Group File Generation:** Create the specific `.set` and `.anno` files required by REGENIE.


* **Inputs:**
* `target_loci.bed`
* List of FNCV/nFNCV Variant IDs.


* **Outputs:** `group_file.txt` (Maps Variant ID $\rightarrow$ Mask Name $\rightarrow$ Locus ID).
* **Integration:** This file defines the specific hypothesis tests (Masks) that REGENIE will evaluate in **Phase 4**.

### Phase 4: Statistical Testing & Enrichment

**Goal:** Quantify the disease burden difference between the Functional and Non-Functional groups.

* **Methodology:**
1. **Genotype Extraction:** Convert WGS VCFs (filtered to the Target Loci) into BGEN format.
2. **REGENIE Step 2:** Run the association test using **SKAT-O** (Sequence Kernel Association Test - Optimal). This is performed twice per locus: once for the FNCV mask, once for the nFNCV mask.
3. **Enrichment Calculation:** Calculate the delta statistic for each locus:

$$\Delta P = (-\log_{10} P_{\text{FNCV}}) - (-\log_{10} P_{\text{nFNCV}})$$


4. **Permutation Testing:** To assess significance, randomly shuffle the "Functional" and "Non-Functional" labels 1,000 times and re-run the test to generate an empirical P-value.


* **Inputs:**
* Phase 1 Null Model (`.loco`).
* `MS_phenotype.txt` and `MS_covariates.txt`.
* WGS BGEN files.
* Phase 3 Group Files.


* **Outputs:** Final Enrichment Report containing locus-specific $\Delta P$ and empirical significance.