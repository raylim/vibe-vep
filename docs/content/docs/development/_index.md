---
title: Development
description: Setting up a development environment and contributing to vibe-vep.
weight: 4
---

## Prerequisites

- Go 1.24+ with `CGO_ENABLED=1` (required for DuckDB and SQLite)
- Git

## Building

```bash
git clone https://github.com/inodb/vibe-vep.git
cd vibe-vep
CGO_ENABLED=1 go build -o vibe-vep ./cmd/vibe-vep
```

## Running Tests

```bash
# Run all tests (fast mode, ~1s)
go test ./... -short -count=1

# Run benchmarks
go test ./internal/annotate/ -bench . -benchmem

# Run full validation benchmark (~30min)
go test ./internal/output/ -run TestValidationBenchmark -v -count=1 -timeout 30m

# Run annotation sources benchmark
go test ./internal/output/ -run TestAnnotationSourcesBenchmark -v -count=1 -timeout 30m
```

## Test Data

Download TCGA test data for validation (~1.6GB):

```bash
make download-testdata
```

This provides 7 TCGA GDC studies (1,052,366 total variants):

| Study | Variants |
|-------|----------|
| CHOL | 3,800 |
| GBM | 55,000 |
| BRCA | 89,000 |
| BLCA | 116,000 |
| LUAD | 191,000 |
| COAD | 245,000 |
| SKCM | 353,000 |

## Project Structure

```
cmd/vibe-vep/       CLI entry point (annotate, download commands)
internal/
  annotate/         Consequence prediction (PredictConsequence, Annotator)
  cache/            Transcript cache (GENCODE GTF/FASTA loader)
  duckdb/           DuckDB cache for transcripts and variant results
  genomicindex/     Unified SQLite index for annotation source lookups (AM, ClinVar, SIGNAL)
  maf/              MAF file parser
  output/           Output formatting and validation comparison
  vcf/              VCF file parser
testdata/
  cache/            Test transcript data (JSON)
  tcga/             TCGA GDC MAF files for validation
```

## Roadmap

- [ ] **Feature parity for MAF annotation** — Match the annotation capabilities of the [genome-nexus-annotation-pipeline](https://github.com/genome-nexus/genome-nexus-annotation-pipeline)
  - [x] Consequence prediction (99.8% concordance)
  - [x] HGVSp/HGVSc notation
  - [x] Full MAF output format
  - [x] Cancer gene annotations (OncoKB)
  - [x] VCF to MAF conversion
  - [x] `--pick` / `--most-severe` annotation filtering
- [ ] **Additional annotation sources**
  - [x] AlphaMissense pathogenicity scores
  - [x] ClinVar clinical significance
  - [x] Cancer Hotspots
  - [x] SIGNAL germline frequencies
  - [ ] SIFT predictions
  - [ ] PolyPhen-2 predictions
  - [ ] gnomAD allele frequencies
- [ ] **Re-annotate datahub GDC studies**
- [ ] **Replace genome-nexus-annotation-pipeline for datahub**

## License

MIT License
