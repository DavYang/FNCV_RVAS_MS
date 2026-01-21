
# Enrichment Testing of AlleleStacker Candidate Variants for Blood-Trait GWAS Loci (EUR-only)
## Phenotype: Multiple Sclerosis

## Overview

This analysis evaluates whether candidate **functional non-coding variants (FNCVs)** are enriched for genetic loci associated with blood-related traits, compared to **non-functional non-coding variants (nFNCVs)**, using European-ancestry GWAS results from the All of Us (AoU) AllxAll resource.

---

### 1. Analysis Architecture

The study tests if variants with specific methylation signatures (FNCVs) drive disease association at GWAS loci compared to matched background variants (nFNCVs), starting with:

- `Phenotype = NS_326.1`

#### **Phase 1: Background**

Establish the genomic background model to correct for relatedness, to be used for subsequent testing.

- **Status:** PLINK files (`.bed`, `.bim`, `.fam`) containing ~500K common SNPs are saved in `./null_model_data`.

#### **Phase 2: Signal Dissection & Locus Definition**
- Note: For this specific phenotype, it is standard practice to either exclude the MHC region from the general pipeline or handle it with a specialized HLA-specific tool. This analysis will exclude this region.

Define the exact genomic windows (GWAS Loci) where the enrichment test will occur, generating independent loci for testing, defined by 95% credible sets (may loosen this to increase loci width).

- **Methods:**
    - **Lead SNP Selection:** Identify significant hits (`p < 5e-8`) from `NS_326.1`
` summary statistics.
    - **Locus Definition:**  
        - *GCTA-COJO*: Used for conditional analysis to find secondary independent SNPs.  
        - *SuSiE Fine-Mapping*: Used to narrow the window to the 95% Credible Set.

- **Inputs:** AoU Summary Statistics  
- **Outputs:** A BED file or list of `Locus_ID` coordinates (Chromosome, Start, End)

---

### Phase 3: Methylation-Based Variant Classification

**Goal:** Sort variants within GWAS loci into "Functional" (Query) and "Non-Functional" (Control) bins.

- **Selection Logic:**
    - **FNCV (Functional Query):** Candidate FNCVs where variant alleles are exclusively observed as either hypermethylated/unmethylated AND are methylation outliers.
    - **nFNCV (Matched Control):** Matched set of variants that are *not* methylation outliers.
    - **Matching:** nFNCVs are selected to match FNCVs on MAF, GC content, and proximity to TSS.

- **Inputs:** FNCV/nFNCVs (Variant IDs), Defined GWAS Loci  
- **Outputs:** `group_file.txt` (Regenie format mapping Variant Set Locus)

---

### Phase 4: Statistical Testing & Enrichment

**Goal:** Test if the FNCV group carries more disease risk than the nFNCV group.

- **Methods:**
    - **Regenie Step 1:** Run the Null Model using Phase 1 output.
    - **Regenie Step 2 (SKAT-O):** Run two specific masks per locus:  
        1. *Mask A*: Methylation-Defined Variants (FNCVs)  
        2. *Mask B*: Matched Background Variants (nFNCVs)
    - **Enrichment Metric:** Compare the P-values (`(-log10P_FNCV) – (-log10P_nFNCV)`) or Beta estimates.

- **Inputs:** Regenie Null Model, WGS BGEN files (filtered to Loci)  
- **Outputs:** Enrichment delta (and empirical significance via permutation)

---

