# GRCh38 Validation Report

Generated: 2026-03-11 15:32 UTC  
Assembly: GRCh38  
GENCODE transcripts: 507365 (loaded from gob cache in 2.705s)  
Workers: 160 (GOMAXPROCS)

## Match Rates

| Study | Variants | Conseq Match | Conseq Mismatch | Conseq Rate | HGVSp Match | HGVSp Mismatch | HGVSp Rate | HGVSc Match | HGVSc Mismatch | HGVSc Rate |
|-------|----------|-------------|-----------------|-------------|-------------|----------------|------------|-------------|----------------|------------|
| blca_tcga_gdc | 116684 | 116598 | 1 | 99.9% | 110950 | 0 | 95.1% | 115140 | 0 | 98.7% |
| brca_tcga_gdc | 89012 | 88945 | 1 | 99.9% | 84478 | 0 | 94.9% | 87772 | 0 | 98.6% |
| chol_tcga_gdc | 3764 | 3761 | 0 | 99.9% | 3500 | 0 | 93.0% | 3694 | 0 | 98.1% |
| coad_tcga_gdc | 244552 | 244350 | 3 | 99.9% | 229450 | 0 | 93.8% | 240951 | 0 | 98.5% |
| gbm_tcga_gdc | 54870 | 54810 | 0 | 99.9% | 51949 | 0 | 94.7% | 54140 | 0 | 98.7% |
| luad_tcga_gdc | 190868 | 190705 | 1 | 99.9% | 181808 | 0 | 95.3% | 188382 | 0 | 98.7% |
| skcm_tcga_gdc | 353450 | 353027 | 3 | 99.9% | 336999 | 0 | 95.3% | 348942 | 0 | 98.7% |
| **Total** | **1053200** | **1052196** | **9** | **99.9%** | **999134** | **0** | **94.9%** | **1039021** | **0** | **98.7%** |

## Consequence Category Breakdown

| Study | delins_normalized | gene_model_change | match | mismatch | no_cds_data | position_shift | transcript_model_change | upstream_reclassified |
|-------|------|------|------|------|------|------|------|------|
| blca_tcga_gdc | 1 | 0 | 115692 | 1 | 2 | 46 | 36 | 906 |
| brca_tcga_gdc | 2 | 0 | 88401 | 1 | 1 | 38 | 25 | 544 |
| chol_tcga_gdc | 0 | 0 | 3730 | 0 | 0 | 1 | 2 | 31 |
| coad_tcga_gdc | 8 | 3 | 242677 | 3 | 2 | 107 | 79 | 1673 |
| gbm_tcga_gdc | 0 | 0 | 54439 | 0 | 0 | 37 | 23 | 371 |
| luad_tcga_gdc | 1 | 1 | 189355 | 1 | 0 | 113 | 47 | 1350 |
| skcm_tcga_gdc | 5 | 2 | 350460 | 3 | 1 | 287 | 125 | 2567 |
| **Total** | **17** | **6** | **1044754** | **9** | **6** | **629** | **337** | **7442** |

## HGVSp Category Breakdown

| Study | both_empty | delins_normalized | dup_vs_ins | fuzzy_fs | maf_empty | maf_nonstandard | match | no_cds_data | position_shift | splice_no_protein | splice_vs_predicted | splice_vs_syn | transcript_model_change | vep_empty |
|-------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|
| blca_tcga_gdc | 485 | 7 | 1 | 481 | 739 | 2672 | 110950 | 2 | 566 | 57 | 28 | 616 | 15 | 65 |
| brca_tcga_gdc | 337 | 18 | 3 | 732 | 446 | 1583 | 84478 | 2 | 730 | 94 | 45 | 475 | 10 | 59 |
| chol_tcga_gdc | 24 | 2 | 0 | 62 | 23 | 85 | 3500 | 0 | 33 | 8 | 1 | 21 | 1 | 4 |
| coad_tcga_gdc | 2306 | 8 | 2 | 2711 | 1317 | 5497 | 229450 | 4 | 1461 | 46 | 33 | 1376 | 43 | 298 |
| gbm_tcga_gdc | 202 | 0 | 2 | 544 | 300 | 1124 | 51949 | 0 | 349 | 14 | 11 | 326 | 15 | 34 |
| luad_tcga_gdc | 731 | 10 | 0 | 1036 | 1070 | 4162 | 181808 | 4 | 837 | 54 | 21 | 1002 | 20 | 113 |
| skcm_tcga_gdc | 1384 | 17 | 0 | 403 | 2039 | 7984 | 336999 | 2 | 1548 | 53 | 21 | 2675 | 72 | 253 |
| **Total** | **5469** | **62** | **8** | **5969** | **5934** | **23107** | **999134** | **14** | **5524** | **326** | **160** | **6491** | **176** | **826** |

## HGVSc Category Breakdown

| Study | both_empty | delins_normalized | dup_vs_ins | maf_empty | match | position_shift | transcript_model_change | vep_empty |
|-------|------|------|------|------|------|------|------|------|
| blca_tcga_gdc | 9 | 82 | 25 | 899 | 115140 | 442 | 87 | 0 |
| brca_tcga_gdc | 8 | 106 | 32 | 538 | 87772 | 477 | 79 | 0 |
| chol_tcga_gdc | 0 | 7 | 3 | 31 | 3694 | 21 | 8 | 0 |
| coad_tcga_gdc | 20 | 90 | 171 | 1658 | 240951 | 1431 | 227 | 4 |
| gbm_tcga_gdc | 3 | 9 | 21 | 369 | 54140 | 288 | 40 | 0 |
| luad_tcga_gdc | 22 | 207 | 43 | 1332 | 188382 | 750 | 131 | 1 |
| skcm_tcga_gdc | 45 | 152 | 26 | 2527 | 348942 | 1566 | 190 | 2 |
| **Total** | **107** | **653** | **321** | **7354** | **1039021** | **4975** | **762** | **7** |

## Cancer Gene Mismatches

No mismatches across all 1207 cancer genes tested.

## Performance

Transcript load: 2.705s from gob cache

| Study | Variants | Sequential | Seq v/s | Parallel | Par v/s | Speedup |
|-------|----------|-----------|---------|----------|---------|--------|
| blca_tcga_gdc | 116684 | 5.689s | 20511 | 617ms | 189032 | 9.22x |
| brca_tcga_gdc | 89012 | 4.381s | 20316 | 495ms | 179822 | 8.85x |
| chol_tcga_gdc | 3764 | 189ms | 19923 | 23ms | 160193 | 8.04x |
| coad_tcga_gdc | 244552 | 11.38s | 21489 | 1.338s | 182795 | 8.51x |
| gbm_tcga_gdc | 54870 | 2.565s | 21393 | 283ms | 194062 | 9.07x |
| luad_tcga_gdc | 190868 | 8.667s | 22022 | 1.255s | 152057 | 6.90x |
| skcm_tcga_gdc | 353450 | 15.33s | 23057 | 1.675s | 211044 | 9.15x |
| **Total** | **1053200** | **48.201s** | **21850** | **5.686s** | **185215** | **8.48x** |
