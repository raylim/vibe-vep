---
title: Validation
description: Validation methodology and results against TCGA datasets.
weight: 5
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

## Validation Results

Tested against 7 TCGA GDC studies from [cBioPortal/datahub](https://github.com/cBioPortal/datahub) (1,052,366 total variants):

| Column | Match | Mismatch | Rate |
|--------|-------|----------|------|
| Consequence | 1,050,079 | 46 | 99.8% |
| HGVSp | 997,292 | 282 | 94.8% |
| HGVSc | 1,037,748 | 278 | 98.6% |

Mismatches that are due to GENCODE version differences (not algorithm bugs) are reclassified into separate categories:
- **transcript_model_change** (501): transcript biotype changed between versions (e.g. protein_coding to retained_intron)
- **gene_model_change** (4): gene boundary differences (coding vs intergenic)
- **position_shift** (695 consequence, 4,733 HGVSp, 5,753 HGVSc): CDS coordinate changes between GENCODE versions

The validation also tracks mismatches in [OncoKB cancer genes](https://www.oncokb.org/cancerGenes) — only 9 cancer genes have any mismatches (1 each) across 1M+ variants.

See the full [validation report](https://github.com/inodb/vibe-vep/blob/main/testdata/tcga/validation_report.md) for per-study breakdowns and category details.

### Reproducing Validation

To download the TCGA test data and regenerate the report:

```bash
make download-testdata
go test ./internal/output/ -run TestValidationBenchmark -v -count=1
```

### Remaining Mismatches

The 46 consequence mismatches are primarily:
- CDS sequence differences between GENCODE versions (synonymous vs missense)
- Transcript structure differences (exon/intron boundary changes)
- Complex multi-region indels with ambiguous consequence priority
