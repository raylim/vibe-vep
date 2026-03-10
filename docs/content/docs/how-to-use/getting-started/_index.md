---
title: Getting Started
description: Install vibe-vep and annotate your first variants.
weight: 1
aliases:
  - /docs/getting-started/
---

## Installation

```bash
go install github.com/inodb/vibe-vep/cmd/vibe-vep@latest
```

Or build from source:

```bash
git clone https://github.com/inodb/vibe-vep.git
cd vibe-vep
go build -o vibe-vep ./cmd/vibe-vep
```

## Quick Start

### 1. Download GENCODE Annotations (one-time setup)

```bash
# For GRCh38 (default)
vibe-vep download --assembly GRCh38

# For GRCh37
vibe-vep download --assembly GRCh37
```

This downloads ~95MB of annotation data to `~/.vibe-vep/`.

### 2. Annotate Variants

```bash
# Annotate a VCF file
vibe-vep annotate vcf input.vcf

# Annotate a MAF file
vibe-vep annotate maf input.maf

# Output to file
vibe-vep annotate vcf -o output.vcf input.vcf

# Convert VCF to MAF
vibe-vep convert vcf2maf input.vcf
```

### 3. Validate MAF Annotations

Compare existing MAF annotations against VEP predictions:

```bash
vibe-vep compare data_mutations.txt
```

## Usage

```
vibe-vep - Variant Effect Predictor

Commands:
  annotate    Annotate variants (vcf, maf, or variant subcommands)
  compare     Compare MAF annotations against predictions
  config      Manage configuration (show/set/get)
  convert     Convert between formats (vcf2maf)
  download    Download GENCODE annotation files
  prepare     Build transcript cache for fast startup
  version     Show version and data source information

Annotate Options:
  --assembly      Genome assembly: GRCh37 or GRCh38 (default: GRCh38)
  -o, --output    Output file (default: stdout)
  --canonical     Only report canonical transcript annotations
  --pick          One annotation per variant (best transcript)
  --most-severe   One annotation per variant (highest impact)
  --save-results  Save annotation results to DuckDB for later lookup
  --no-cache      Skip transcript cache, always load from GTF/FASTA
  --clear-cache   Clear and rebuild transcript and variant caches

Download Options:
  --assembly      Genome assembly: GRCh37 or GRCh38 (default: GRCh38)
  --output        Output directory (default: ~/.vibe-vep/)
```

## Examples

### Annotating a VCF file

The repository includes a small [example VCF](https://github.com/inodb/vibe-vep/blob/main/examples/example.vcf) containing the KRAS G12C (p.Gly12Cys) hotspot mutation:

```bash
vibe-vep annotate vcf examples/example.vcf
```

VCF input produces VCF output by default — the original VCF lines are preserved and a `CSQ` INFO field is added with consequence annotations for all overlapping transcripts (canonical flagged with `YES`):

```
##fileformat=VCFv4.2
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from vibe-vep. Format: Allele|Consequence|IMPACT|SYMBOL|...">
#CHROM  POS       ID           REF  ALT  QUAL  FILTER  INFO
12      25245351  rs121913529  C    A    100   PASS    DP=50;CSQ=A|missense_variant|MODERATE|KRAS|...|YES,...
```

### Annotating a MAF file

The repository also includes a small [example MAF](https://github.com/inodb/vibe-vep/blob/main/examples/example.maf) with the same KRAS G12C variant:

```bash
vibe-vep annotate maf examples/example.maf
```

When the input is MAF, the output defaults to MAF format — all original columns are preserved and annotation columns are updated with fresh VEP predictions:

```
Hugo_Symbol  Entrez_Gene_Id  Center  ...  Consequence       Variant_Classification  ...  HGVSc    HGVSp        HGVSp_Short  Transcript_ID
KRAS         3845            .       ...  missense_variant  Missense_Mutation       ...  c.34G>T  p.Gly12Cys   p.G12C       ENST00000311936
```

### Other examples

```bash
# Validate TCGA MAF file
vibe-vep compare data_mutations.txt

# Show all validation results (not just mismatches)
vibe-vep compare --all data_mutations.txt

# Annotate with GRCh37
vibe-vep annotate vcf --assembly GRCh37 sample.vcf

# Pick one annotation per variant (best transcript)
vibe-vep annotate vcf --pick input.vcf

# Convert VCF to MAF format
vibe-vep convert vcf2maf input.vcf -o output.maf
```

## Configuration

vibe-vep reads configuration from `~/.vibe-vep.yaml` (or `--config` flag). Environment variables with `VIBE_VEP_` prefix also work (e.g. `VIBE_VEP_ONCOKB_CANCER_GENE_LIST`).

```yaml
# OncoKB cancer gene list — adds Gene_Type column to MAF output
oncokb:
  cancer-gene-list: /path/to/cancerGeneList.tsv

# Optional annotation sources
annotations:
  alphamissense: true   # AlphaMissense pathogenicity scores (CC BY 4.0)
  clinvar: true         # ClinVar clinical significance
  hotspots: /path/to/hotspots_v2_and_3d.txt  # Cancer Hotspots (path to TSV)
  signal: true          # SIGNAL germline frequencies (GRCh37 only)
```

The `cancerGeneList.tsv` file can be downloaded from [OncoKB](https://www.oncokb.org/cancerGenes) or is included in the repository.

**Enabling annotation sources:**

```bash
# AlphaMissense (GRCh38): download + prepare + enable
vibe-vep config set annotations.alphamissense true
vibe-vep download  # fetches ~643 MB
vibe-vep prepare   # builds SQLite index

# ClinVar (GRCh38): download + enable
vibe-vep config set annotations.clinvar true
vibe-vep download  # fetches ~182 MB clinvar.vcf.gz

# Hotspots: point to TSV file
vibe-vep config set annotations.hotspots /path/to/hotspots_v2_and_3d.txt

# SIGNAL (GRCh37 only): enable
vibe-vep config set annotations.signal true
```

Use `vibe-vep version` to see which sources are loaded and `vibe-vep version --maf-columns` for the full column mapping.
