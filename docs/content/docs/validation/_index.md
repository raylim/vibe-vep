---
title: Validation
description: Validation methodology and results against TCGA datasets.
weight: 3
---

## How Validation Works

The `vibe-vep compare` command compares existing MAF annotations against fresh predictions. For each MAF variant, the tool:

1. Parses the variant and its existing annotation (consequence, gene, transcript ID)
2. Re-annotates the variant against all overlapping GENCODE transcripts
3. Selects the best matching annotation for comparison (see below)
4. Normalizes consequence terms and compares

## Transcript Selection for Validation

When comparing against a MAF entry, the tool selects the best VEP annotation using this priority:

1. **Exact transcript ID match** — if the MAF specifies a transcript ID (e.g., ENST00000311936), use the annotation for that transcript. However, if the transcript's biotype has changed between GENCODE versions (e.g., was `protein_coding` but is now `retained_intron`) and the MAF has a coding consequence, the match is skipped.

2. **Same gene, best transcript** — among annotations for the same gene (by Hugo symbol), prefer:
   - Canonical transcripts (marked by GENCODE/Ensembl)
   - Protein-coding biotypes over non-coding
   - Higher-impact consequences (HIGH > MODERATE > LOW > MODIFIER)

3. **Any transcript, best match** — if no gene match, use the same ranking across all annotations.

## Consequence Normalization

MAF files may use different consequence terms than SO standard. The validation normalizes both sides before comparison:

- **MAF to SO mapping**: `Missense_Mutation` to `missense_variant`, `Silent` to `synonymous_variant`, etc.
- **Modifier stripping**: Drops terms like `non_coding_transcript_variant`, `NMD_transcript_variant`, `coding_sequence_variant`, `start_retained_variant`, `stop_retained_variant` that are secondary modifiers when a higher-impact term is present
- **Splice normalization**: Maps `splice_donor_region_variant` and `splice_donor_5th_base_variant` to `splice_region_variant`; drops `splice_region_variant` when a primary consequence is present
- **Impact-based stripping**: Drops `intron_variant` when splice donor/acceptor is present; drops UTR terms when a HIGH-impact term is present; drops `stop_gained`/`stop_lost` when co-occurring with `frameshift_variant`
- **Inframe grouping**: Normalizes `protein_altering_variant`, `inframe_deletion`, and `inframe_insertion` to a common term
- **Upstream/downstream tolerance**: MAF upstream/downstream calls are always accepted as matching, since different canonical transcript sets produce different transcript boundaries
- **Sorting**: Comma-separated terms are sorted alphabetically for consistent comparison

## GRCh38 Results (TCGA)

Tested against 7 TCGA GDC studies from [cBioPortal/datahub](https://github.com/cBioPortal/datahub):

{{< validation-report assembly="grch38" >}}

## GRCh37 Results (MSK-IMPACT)

{{< validation-report assembly="grch37" >}}

## Reproducing Validation

To download test data and regenerate the reports:

```bash
# Download all test data
make download-testdata

# Run GRCh38 validation
go test ./internal/output/ -run TestValidationBenchmark$ -v -count=1 -timeout 30m

# Run GRCh37 validation
go test ./internal/output/ -run TestValidationBenchmarkGRCh37 -v -count=1 -timeout 30m
```

Full markdown reports are also available:
- [GRCh38 report](https://github.com/inodb/vibe-vep/blob/main/testdata/tcga/validation_report.md)
- [GRCh37 report](https://github.com/inodb/vibe-vep/blob/main/testdata/grch37/validation_report.md)
