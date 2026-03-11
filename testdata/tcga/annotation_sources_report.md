# Annotation Sources Report

Generated: 2026-03-11 15:31 UTC  
GENCODE transcripts: 254070  
Data: unified genomic index (SQLite)  
Workers: 4 (GOMAXPROCS)

## AlphaMissense Coverage

| Study | Variants | Missense | AM Hits | Coverage | likely_benign | ambiguous | likely_pathogenic |
|-------|----------|----------|---------|----------|---------------|-----------|-------------------|
| blca_tcga_gdc | 115850 | 75938 | 69937 | 92.1% | 41805 | 7662 | 20470 |
| brca_tcga_gdc | 89012 | 57703 | 53240 | 92.3% | 31416 | 5807 | 16017 |
| chol_tcga_gdc | 3764 | 2336 | 2151 | 92.1% | 1244 | 234 | 673 |
| coad_tcga_gdc | 244552 | 128655 | 118282 | 91.9% | 70342 | 12727 | 35213 |
| gbm_tcga_gdc | 54870 | 32684 | 30155 | 92.3% | 18286 | 3259 | 8610 |
| luad_tcga_gdc | 190868 | 120249 | 111597 | 92.8% | 65087 | 12682 | 33828 |
| skcm_tcga_gdc | 353450 | 204372 | 185786 | 90.9% | 115483 | 20014 | 50289 |
| **Total** | **1052366** | **621937** | **571148** | **91.8%** | **343663** | **62385** | **165100** |

## ClinVar Coverage

| Study | Variants | ClinVar Hits | Hit Rate | Pathogenic | Likely_pathogenic | Uncertain | Benign | Likely_benign | Other |
|-------|----------|--------------|----------|------------|-------------------|-----------|--------|---------------|-------|
| blca_tcga_gdc | 115850 | 8211 | 7.09% | 605 | 153 | 4999 | 47 | 1649 | 758 |
| brca_tcga_gdc | 89012 | 8025 | 9.02% | 915 | 167 | 4559 | 56 | 1478 | 850 |
| chol_tcga_gdc | 3764 | 395 | 10.49% | 39 | 9 | 230 | 2 | 70 | 45 |
| coad_tcga_gdc | 244552 | 32906 | 13.46% | 2415 | 776 | 19501 | 224 | 6643 | 3347 |
| gbm_tcga_gdc | 54870 | 8298 | 15.12% | 549 | 121 | 4984 | 59 | 1680 | 905 |
| luad_tcga_gdc | 190868 | 9334 | 4.89% | 622 | 311 | 5635 | 54 | 1859 | 853 |
| skcm_tcga_gdc | 353450 | 28341 | 8.02% | 1023 | 498 | 17041 | 148 | 7319 | 2312 |
| **Total** | **1052366** | **95510** | **9.08%** | **6168** | **2035** | **56949** | **590** | **20698** | **9070** |

## SIGNAL Coverage

*Note: SIGNAL data uses GRCh37 coordinates. TCGA GDC MAFs use GRCh38, so few hits are expected.*

| Study | Variants | SIGNAL Hits | Hit Rate |
|-------|----------|-------------|----------|
| blca_tcga_gdc | 115850 | 0 | 0.00% |
| brca_tcga_gdc | 89012 | 0 | 0.00% |
| chol_tcga_gdc | 3764 | 0 | 0.00% |
| coad_tcga_gdc | 244552 | 0 | 0.00% |
| gbm_tcga_gdc | 54870 | 0 | 0.00% |
| luad_tcga_gdc | 190868 | 0 | 0.00% |
| skcm_tcga_gdc | 353450 | 0 | 0.00% |
| **Total** | **1052366** | **0** | **0.00%** |

## Genomic Index Lookup Performance

| Study | Variants | Base Time | Index Lookup Time | Overhead | Lookups/sec |
|-------|----------|-----------|-------------------|----------|-------------|
| blca_tcga_gdc | 115850 | 21.017s | 52.384s | 249.2% | 2212 |
| brca_tcga_gdc | 89012 | 17.437s | 24.928s | 143.0% | 3571 |
| chol_tcga_gdc | 3764 | 725ms | 921ms | 127.1% | 4086 |
| coad_tcga_gdc | 244552 | 33.095s | 58.798s | 177.7% | 4159 |
| gbm_tcga_gdc | 54870 | 9.415s | 14.882s | 158.1% | 3687 |
| luad_tcga_gdc | 190868 | 29.932s | 51.475s | 172.0% | 3708 |
| skcm_tcga_gdc | 353450 | 46.805s | 1m37.97s | 209.3% | 3608 |
| **Total** | **1052366** | **2m38.426s** | **5m1.357s** | **190.2%** | **3492** |
