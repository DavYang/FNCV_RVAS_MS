# Analysis Memory: FNCV Enrichment in MS (EUR-only)

## Context
Evaluating enrichment of methylation-defined candidate functional non-coding variants (FNCVs) for association with Multiple Sclerosis (NS_326.1) in All of Us (AoU) v8 European samples.

## Current Progress

### Phase 0: Phenotype & Covariate Generation - COMPLETED
- **Goal:** strictly formatted files for REGENIE null model (Step 1) and association (Step 2).
- **Methods:**
  - Defined NS_326.1 cases and controls from AoU CDR.
  - Extracted Age, Sex, and top 10 Ancestry PCs.
  - Formatted as space-delimited text with `-9` for missing values.
- **Key Files:**
  - `MS_phenotype.txt`
  - `MS_covariates.txt`

### Phase 1: Background SNPs - COMPLETED
- **Background SNP Generation:** 
  - Extracted common SNPs from AoU v8 ACAF splitMT.
  - Applied call rate filter (>= 0.95) in Hail to address structural missingness in ACAF data.
  - Used 2.0x sampling overshoot to ensure target SNP counts post-QC.
  - Performed per-chromosome QC and merged into a genome-wide fileset.
- **QC Parameters:** 
  - `MAF >= 0.01`
  - `GENO <= 0.05`
  - `MAX-ALLELES 2`
  - Note: HWE and MIND filters were removed due to non-random missingness in ACAF data for EUR samples.
- **Key Files:**
  - Raw per-chrom PLINK: `results/1-bg_snp/plink_no-qc/`
  - QC'd per-chrom PLINK: `results/1-bg_snp/plink_qc/`
  - Merged genome-wide PLINK: `results/1-bg_snp/plink_qc/all_background_qc.{bed,bim,fam}`
  - Final post-merge QC PLINK: `results/1-bg_snp/plink_qc/all_background_final_qc.{bed,bim,fam}`

## Next Steps: Phase 2 (Signal Dissection & Locus Definition)
- **Objective:** Define genomic windows (GWAS Loci) for enrichment testing.
- **Tasks:**
  1. Identify lead SNPs (p < 5e-8) from NS_326.1 summary statistics.
  2. Use GCTA-COJO to find independent secondary signals.
  3. Use SuSiE fine-mapping to define 95% Credible Sets for locus boundaries.
  4. Exclude MHC region (Chr6: 25-35Mb).

## Critical Technical Notes
- **Missingness:** High sample-level missingness (~14% median) is expected/structural in AoU ACAF EUR data; do not use `--mind`.
- **HWE:** Do not use strict `--hwe` on ACAF-derived background SNPs as population-specific thresholds violate HWE assumptions.
- **Compute:** Phase 2 will require GCTA and SuSiE, using Dataproc for LD reference handling if necessary.
