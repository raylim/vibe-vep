---
title: Consequence Prediction
description: How vibe-vep predicts variant consequences using Sequence Ontology terms.
weight: 1
aliases:
  - /docs/consequence-prediction/
---

## Supported Consequences

The tool predicts variant consequences using Sequence Ontology (SO) terms, ordered by impact:

| Impact | Consequences |
|--------|-------------|
| HIGH | stop_gained, frameshift_variant, stop_lost, start_lost, splice_donor_variant, splice_acceptor_variant |
| MODERATE | missense_variant, inframe_insertion, inframe_deletion |
| LOW | synonymous_variant, splice_region_variant, stop_retained_variant |
| MODIFIER | intron_variant, 5_prime_UTR_variant, 3_prime_UTR_variant, upstream/downstream_gene_variant, non_coding_transcript_exon_variant, mature_miRNA_variant |

## How It Works

For each variant, the tool determines which transcripts overlap the variant position, then classifies the effect on each transcript:

1. **Upstream/Downstream**: Variant outside transcript boundaries
2. **Intronic**: Variant between exons. Checks for splice site overlap:
   - **Splice donor/acceptor** (HIGH): within +/-1-2bp of exon boundary on intron side (donor at 5' end of intron, acceptor at 3' end, strand-aware)
   - **Splice region** (LOW): within 3bp exon side or 3-8bp intron side of splice junction
3. **UTR**: Variant in 5' or 3' untranslated region
4. **Coding**: Variant in CDS — calculates codon/amino acid change to determine missense, synonymous, stop_gained, etc.

For **indels**, the tool checks the entire deletion span (not just the start position) for:
- Splice site overlap across the full `[pos, pos+len(ref)-1]` range
- Start codon deletion (for deletions spanning from UTR into CDS)
- Stop codon overlap (frameshift at stop codon produces frameshift_variant,stop_lost)
- In-frame insertions that create stop codons

## Transcript Biotype Handling

The tool treats transcripts as protein-coding if they have defined CDS coordinates (`CDSStart > 0 && CDSEnd > 0`), which covers:
- `protein_coding`
- `nonsense_mediated_decay` (appends `NMD_transcript_variant` modifier)
- `IG_*_gene` / `TR_*_gene` (immunoglobulin/T-cell receptor segments)
- `protein_coding_LoF`, `non_stop_decay`

For `miRNA` biotype transcripts, exonic variants are classified as `mature_miRNA_variant` instead of `non_coding_transcript_exon_variant`.
