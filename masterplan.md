# Master Analysis Plan: AlleleStacker Enrichment Testing (EUR-only)

## Phenotype: Multiple Sclerosis (NS_326.1)

## 1. Overview

This analysis evaluates whether candidate **functional non-coding variants (FNCVs)** are enriched for genetic loci associated with Multiple Sclerosis (MS), compared to matched **non-functional non-coding variants (nFNCVs)**. The study leverages the *All of Us* (AoU) v7 European-ancestry cohort to perform a controlled Rare Variant Association Study (RVAS).

The central hypothesis is that rare variants exhibiting methylation outlier signatures (FNCVs) carry a higher burden of disease risk within established GWAS loci than strictly matched background variants (nFNCVs).

---

## 2. Software Stack & Toolchain

The pipeline utilizes a distributed cloud architecture (Hail, Google Cloud Services) integrating the following tools:

* **Hail (Python/Spark):** Used for genomic data wrangling, specifically extracting common background SNPs and EUR sample genotypes from AoU MatrixTables (mt).
* **PLINK (1.9/2.0):** Used for Quality Control (QC), Linkage Disequilibrium (LD) pruning, and merging of chromosome-level filesets.
* **REGENIE:** The primary statistical engine. Step 1 fits the whole-genome null model; Step 2 performs the SKAT-O aggregation tests.
* **pyliftover (Python):** Lifts GWAS Catalog SNP positions from GRCh37 to GRCh38 for positional cross-referencing.

> **Benched (not currently used):** GCTA-COJO (conditional/joint analysis), SuSiE (Bayesian fine-mapping), PLINK LD clumping. These approaches were evaluated but deprioritized due to computational overhead and LD reference panel complexity.

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

### Phase 2: Locus Definition via GWAS Catalog Cross-Reference

**Goal:** Define the exact genomic windows (Target Loci) where the enrichment test will occur, anchored to known EUR MS GWAS loci from the literature.

* **Methodology:**

**Step A — Export AoU top SNPs (HPC/workbench):**
1. Export all AoU `NS_326.1` GWAS SNPs at $p < 5 \times 10^{-6}$ from per-chromosome `.ma` summary statistic files.
2. Output: `results/2-locus_definition/top_gwas_snps.tsv` (fields: SNP, chrom, pos, ref, alt, beta, se, p, freq).

**Step B — Parse GWAS Catalog (HPC):**
1. Load the full GWAS Catalog associations TSV (`gwas-catalog-download-associations-v1.0-full.tsv`, GRCh37 positions).
2. Filter to: MS trait (case-insensitive match on `DISEASE/TRAIT`), $p < 5 \times 10^{-8}$, autosomes, non-null positions.
3. **EUR-only study filter:** `INITIAL SAMPLE SIZE` must contain "European ancestry" and must **not** mention any non-EUR ancestry group (African, Asian, Japanese, Hispanic, admixed, Latino, multi, etc.). This excludes multi-ancestry GWAS entirely.
4. Liftover `CHR_POS` from GRCh37 to GRCh38 using `pyliftover` + UCSC `hg19ToHg38.over.chain.gz`.
5. Output: `gwas_catalog_ms_eur_hg38.tsv` → upload to GCS workspace bucket.

**Step C — Cross-reference & window generation (AoU workbench):**
1. For each AoU top SNP (GRCh38), check whether any EUR MS catalog SNP falls within $\pm 250$ kb (positional match, both now GRCh38).
2. **Validated SNPs** (catalog hit within ±250kb): generate $\pm 250$ kb window around the AoU SNP position.
3. **Unmatched SNPs** (no catalog hit): written to `novel_snps.tsv`, excluded from locus generation.
4. Merge overlapping windows (sort-and-merge on chromosome/start).
5. Exclude MHC region (Chr6: 25Mb–35Mb).
6. Assign `locus_id = chrN_<lead_snp_pos>` (lowest $p$ within merged window).

* **Inputs:**
    * AoU GWAS summary statistics (`.ma` files, per-chromosome).
    * GWAS Catalog full associations TSV (local HPC copy).

* **Outputs:**
    * `results/2-locus_definition/top_gwas_snps.tsv` — AoU top SNPs at $p < 5 \times 10^{-6}$.
    * `gwas_catalog_ms_eur_hg38.tsv` — filtered EUR MS catalog SNPs (GRCh38).
    * `results/2-locus_definition/gwas_catalog_validated_snps.tsv` — AoU SNPs matched to catalog.
    * `results/2-locus_definition/novel_snps.tsv` — AoU top SNPs with no catalog support.
    * `results/2-locus_definition/target_loci.bed` (Columns: Chrom, Start, End, Locus_ID).

* **Integration:** These coordinates restrict **Phase 3** so that variant classification only occurs within biologically relevant, literature-supported disease-associated regions.

> **Benched approaches:** GCTA-COJO signal dissection, SuSiE fine-mapping, and PLINK LD clumping were evaluated for locus definition but deprioritized due to computational overhead and LD reference panel complexity. Scripts `02b_run_cojo.sh` and `02b_define_loci.py` are retained but not currently used.

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