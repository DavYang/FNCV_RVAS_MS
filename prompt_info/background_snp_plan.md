## Data Processing Flow for Sampling 1M Common SNPs (EUR Ancestry)

**1. Initialize Hail Environment**
- Start Hail session with appropriate cluster configuration
- Allocate sufficient memory for ~200K samples × 1M variants

**2. Load Ancestry Metadata**
- Read ancestry TSV file into Hail Table
- Key the table by sample ID to enable joining with genotype data

**3. Load Pre-filtered MatrixTable**
- Read ACAF threshold MT (AF > 0.01) directly from Google Cloud Storage
- No copying required - query in place

**4. Filter to EUR Ancestry Samples**
- Annotate MT columns with ancestry information from the TSV table
- Filter columns to retain only samples with EUR ancestry (~200K samples)
- This reduces the column dimension from ~400K to ~200K

**5. Random Sampling of Variants**
- Calculate sampling fraction: 1,000,000 / total_variant_count
- Use `mt.sample_rows()` to randomly select 1M SNPs
- This reduces the row dimension from ~99M to 1M variants

**6. Quality Check (Optional)**
- Verify final dimensions: ~200K samples × 1M variants
- Check allele frequency distribution to confirm common variants
- Verify no missingness issues

**7. Export Results**
- Write filtered MT to cloud storage, or
- Export to PLINK format if needed for downstream analysis, or
- Export to VCF/TSV depending on use case

**8. Clean Up**
- Stop Hail session to release cluster resources
- Verify output files are accessible