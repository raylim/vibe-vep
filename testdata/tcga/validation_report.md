# GRCh38 Validation Report

Generated: 2026-03-11 16:49 UTC  
Assembly: GRCh38  
GENCODE transcripts: 254070 (loaded from GTF/FASTA in 31.732s)  
Workers: 4 (GOMAXPROCS)

## Match Rates

| Study | Variants | Conseq Match | Conseq Mismatch | Conseq Rate | HGVSp Match | HGVSp Mismatch | HGVSp Rate | HGVSc Match | HGVSc Mismatch | HGVSc Rate |
|-------|----------|-------------|-----------------|-------------|-------------|----------------|------------|-------------|----------------|------------|
| blca_tcga_gdc | 115850 | 115654 | 1 | 99.8% | 110031 | 0 | 95.0% | 114214 | 0 | 98.6% |
| brca_tcga_gdc | 89012 | 88854 | 1 | 99.8% | 84213 | 0 | 94.6% | 87820 | 0 | 98.7% |
| chol_tcga_gdc | 3764 | 3758 | 0 | 99.8% | 3485 | 0 | 92.6% | 3696 | 0 | 98.2% |
| coad_tcga_gdc | 244552 | 244128 | 3 | 99.8% | 229301 | 0 | 93.8% | 241086 | 0 | 98.6% |
| gbm_tcga_gdc | 54870 | 54754 | 0 | 99.8% | 51933 | 0 | 94.6% | 54160 | 0 | 98.7% |
| luad_tcga_gdc | 190868 | 190466 | 2 | 99.8% | 181823 | 0 | 95.3% | 188319 | 0 | 98.7% |
| skcm_tcga_gdc | 353450 | 352616 | 3 | 99.8% | 337041 | 0 | 95.4% | 349039 | 0 | 98.8% |
| **Total** | **1052366** | **1050230** | **10** | **99.8%** | **997827** | **0** | **94.8%** | **1038334** | **0** | **98.7%** |

## Consequence Category Breakdown

| Study | delins_normalized | gene_model_change | match | mismatch | no_cds_data | position_shift | transcript_model_change | upstream_reclassified |
|-------|------|------|------|------|------|------|------|------|
| blca_tcga_gdc | 1 | 2 | 114755 | 1 | 74 | 30 | 88 | 899 |
| brca_tcga_gdc | 2 | 0 | 88310 | 1 | 78 | 27 | 50 | 544 |
| chol_tcga_gdc | 0 | 0 | 3727 | 0 | 4 | 1 | 1 | 31 |
| coad_tcga_gdc | 8 | 3 | 242446 | 3 | 193 | 64 | 153 | 1682 |
| gbm_tcga_gdc | 0 | 0 | 54383 | 0 | 59 | 25 | 32 | 371 |
| luad_tcga_gdc | 1 | 2 | 189115 | 2 | 219 | 77 | 101 | 1351 |
| skcm_tcga_gdc | 5 | 5 | 350048 | 3 | 400 | 214 | 207 | 2568 |
| **Total** | **17** | **12** | **1042784** | **10** | **1027** | **438** | **632** | **7446** |

## HGVSp Category Breakdown

| Study | both_empty | delins_normalized | dup_vs_ins | fuzzy_fs | maf_empty | maf_nonstandard | match | no_cds_data | position_shift | splice_no_protein | splice_vs_predicted | splice_vs_syn | transcript_model_change | vep_empty |
|-------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|
| blca_tcga_gdc | 473 | 7 | 1 | 645 | 745 | 2608 | 110031 | 74 | 456 | 57 | 27 | 612 | 59 | 55 |
| brca_tcga_gdc | 329 | 19 | 4 | 1041 | 454 | 1565 | 84213 | 79 | 616 | 94 | 45 | 476 | 27 | 50 |
| chol_tcga_gdc | 24 | 2 | 0 | 77 | 23 | 84 | 3485 | 4 | 31 | 8 | 1 | 21 | 2 | 2 |
| coad_tcga_gdc | 2284 | 8 | 2 | 3024 | 1339 | 5442 | 229301 | 196 | 1133 | 43 | 34 | 1373 | 94 | 279 |
| gbm_tcga_gdc | 192 | 0 | 2 | 582 | 310 | 1114 | 51933 | 59 | 273 | 14 | 11 | 326 | 25 | 29 |
| luad_tcga_gdc | 711 | 10 | 0 | 1059 | 1090 | 4115 | 181823 | 223 | 594 | 53 | 21 | 1005 | 65 | 99 |
| skcm_tcga_gdc | 1373 | 17 | 0 | 493 | 2050 | 7932 | 337041 | 401 | 1037 | 52 | 18 | 2676 | 123 | 237 |
| **Total** | **5386** | **63** | **9** | **6921** | **6011** | **22860** | **997827** | **1036** | **4140** | **321** | **157** | **6489** | **395** | **751** |

## HGVSc Category Breakdown

| Study | both_empty | delins_normalized | dup_vs_ins | maf_empty | match | position_shift | transcript_model_change | vep_empty |
|-------|------|------|------|------|------|------|------|------|
| blca_tcga_gdc | 10 | 81 | 25 | 892 | 114214 | 334 | 292 | 2 |
| brca_tcga_gdc | 8 | 106 | 34 | 538 | 87820 | 360 | 146 | 0 |
| chol_tcga_gdc | 0 | 7 | 3 | 31 | 3696 | 19 | 8 | 0 |
| coad_tcga_gdc | 25 | 90 | 175 | 1653 | 241086 | 1097 | 422 | 4 |
| gbm_tcga_gdc | 3 | 9 | 22 | 369 | 54160 | 212 | 95 | 0 |
| luad_tcga_gdc | 25 | 207 | 45 | 1329 | 188319 | 509 | 431 | 3 |
| skcm_tcga_gdc | 47 | 152 | 26 | 2525 | 349039 | 1037 | 619 | 5 |
| **Total** | **118** | **652** | **330** | **7337** | **1038334** | **3568** | **2013** | **14** |

## Cancer Gene Mismatches

No mismatches across all 1207 cancer genes tested.

## Performance

Transcript load: 31.732s from GTF/FASTA

| Study | Variants | Sequential | Seq v/s | Parallel | Par v/s | Speedup |
|-------|----------|-----------|---------|----------|---------|--------|
| blca_tcga_gdc | 115850 | 8.026s | 14434 | 4.631s | 25016 | 1.73x |
| brca_tcga_gdc | 89012 | 7.477s | 11904 | 3.302s | 26956 | 2.26x |
| chol_tcga_gdc | 3764 | 210ms | 17919 | 159ms | 23743 | 1.32x |
| coad_tcga_gdc | 244552 | 14.694s | 16643 | 10.195s | 23987 | 1.44x |
| gbm_tcga_gdc | 54870 | 3.354s | 16360 | 1.565s | 35063 | 2.14x |
| luad_tcga_gdc | 190868 | 12.118s | 15750 | 7.728s | 24699 | 1.57x |
| skcm_tcga_gdc | 353450 | 23.092s | 15306 | 15.852s | 22297 | 1.46x |
| **Total** | **1052366** | **1m8.972s** | **15258** | **43.432s** | **24230** | **1.59x** |
