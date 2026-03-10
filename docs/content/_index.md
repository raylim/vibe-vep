---
title: vibe-vep
---

{{< blocks/cover title="vibe-vep" image_anchor="top" height="med" >}}
<a class="btn btn-lg btn-primary me-3 mb-4" href="docs/how-to-use/getting-started/">
Get Started <i class="fas fa-arrow-alt-circle-right ms-2"></i>
</a>
<a class="btn btn-lg btn-secondary me-3 mb-4" href="https://github.com/inodb/vibe-vep">
GitHub <i class="fab fa-github ms-2 "></i>
</a>
<p class="lead mt-4">A lightweight variant effect predictor for cBioPortal, Genome Nexus, and OncoKB</p>
{{< /blocks/cover >}}

{{% blocks/lead color="primary" %}}

**vibe-vep** is a single Go binary with no Perl dependencies that predicts variant consequences using GENCODE annotations (~95MB download vs 17GB VEP cache).

It incorporates transcript prioritization similar to Genome Nexus — selecting a single gene and protein change by prioritizing coding transcripts and highest-impact consequences — while also providing effect predictions for all overlapping transcripts.

Achieves **99.8% concordance** with GDC/VEP annotations across 1M+ TCGA variants.

{{% /blocks/lead %}}

{{% blocks/section color="dark" type="row" %}}

{{% blocks/feature icon="fa-bolt" title="Fast" %}}
~14,000 variants/sec parallel (4 workers), ~5,000 single-threaded.
Cache loading in ~25 seconds for 254k GENCODE v46 transcripts.
{{% /blocks/feature %}}

{{% blocks/feature icon="fa-check-circle" title="Validated" %}}
Tested against 7 TCGA GDC studies (1,052,366 variants).
99.8% consequence match rate, 0 cancer gene mismatches.
{{% /blocks/feature %}}

{{% blocks/feature icon="fa-database" title="Multiple Data Sources" %}}
AlphaMissense, ClinVar, Cancer Hotspots, SIGNAL, and OncoKB annotation sources.
Extensible architecture for adding new sources.
{{% /blocks/feature %}}

{{% /blocks/section %}}

{{% blocks/section %}}

## Key Features

- **VCF/MAF Parsing**: Annotate variants from VCF or MAF files
- **GENCODE Annotations**: Uses GENCODE GTF/FASTA (~95MB download vs 17GB VEP cache)
- **Consequence Prediction**: Classifies variants using Sequence Ontology terms
- **Validation Mode**: Compare existing MAF annotations against predictions
- **VCF to MAF Conversion**: Convert VCF files to MAF format with annotations
- **Annotation Sources**: AlphaMissense, ClinVar, Cancer Hotspots, SIGNAL, OncoKB

{{% /blocks/section %}}
