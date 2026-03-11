#!/usr/bin/env python3
"""DEPRECATED: GCTA-COJO locus definition has been replaced by a windowed approach.

The active Phase 2 pipeline is:
  02a_export_gwas_ma.py       -- export per-chromosome .ma files from AllxAll GWAS HT
  02c_export_top_gwas_snps.py -- filter to p < 5e-6, write top_gwas_snps.tsv
  02d_parse_gwas_catalog.py   -- parse EUR-only MS GWAS Catalog loci (runs on HPC)
  02e_define_loci_catalog.py  -- cross-reference AoU SNPs vs catalog, define +-250kb windows
"""

import sys

print("ERROR: 02b_define_loci.py is deprecated. GCTA-COJO is no longer used.")
print("Run: bash bash/02c_export_top_gwas_snps.sh")
print("Then: bash bash/02e_define_loci_catalog.sh")
sys.exit(1)
