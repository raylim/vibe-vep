---
title: Output Formats
description: VCF and MAF output format details.
weight: 4
aliases:
  - /docs/output-formats/
---

## VCF Output (default for VCF input)

When the input is a VCF file, vibe-vep outputs VCF format by default. Original VCF header and data lines are preserved, and a `CSQ` INFO field is appended with consequence annotations for all overlapping transcripts. The CSQ format follows the VEP convention — pipe-delimited fields per transcript, comma-separated between transcripts:

```
CSQ=Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|CANONICAL
```

Canonical transcripts are flagged with `CANONICAL=YES`.

## MAF Output (default for MAF input)

When the input is a MAF file, vibe-vep outputs MAF format by default. All original columns are preserved exactly as-is, and `vibe.*` namespaced columns are appended with fresh predictions:

| Column | Description |
|--------|-------------|
| `vibe.hugo_symbol` | Gene symbol |
| `vibe.consequence` | SO consequence term |
| `vibe.variant_classification` | MAF variant classification |
| `vibe.transcript_id` | Ensembl transcript ID |
| `vibe.hgvsc` | HGVS coding DNA notation |
| `vibe.hgvsp` | HGVS protein notation (3-letter) |
| `vibe.hgvsp_short` | HGVS protein notation (1-letter) |

When annotation sources are configured, additional `vibe.{source}.{column}` columns are appended (e.g., `vibe.oncokb.gene_type`, `vibe.alphamissense.score`).

Use `vibe-vep version --maf-columns` to see the full column mapping.

## Performance

- **End-to-end throughput**: ~14,000 variants/sec parallel (4 workers), ~5,000 single-threaded
- **Cache loading**: ~25 seconds to load 254k GENCODE v46 transcripts
- **Memory**: Proportional to transcript count (~254k for GENCODE v46)
- **Total benchmark**: 1,052,366 variants across 7 TCGA studies
