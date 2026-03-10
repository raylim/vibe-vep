#!/usr/bin/env bash
# prepare_clinvar.sh — Download and filter ClinVar data for independent benchmarking.
#
# Downloads:
#   ~/.vibe-vep/clinvar/variant_summary.txt.gz   (ClinVar variant summary, ~435 MB)
#   ~/.vibe-vep/clinvar/MANE.GRCh38.v1.5.summary.txt.gz  (MANE Select mapping, ~1 MB)
#
# Produces:
#   testdata/clinvar/clinvar_patho.vcf.gz  — VCF of pathogenic SNVs and small indels with CLNP/CLNTX INFO
#
# Usage:
#   bash scripts/prepare_clinvar.sh [--force]
#
# The output VCF is used by run_snpeff_clinvar.sh and run_vep_clinvar.sh.
# The benchmark test reads variant_summary.txt.gz directly.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

CLINVAR_DIR="${HOME}/.vibe-vep/clinvar"
OUTDIR="${REPO_ROOT}/testdata/clinvar"
FORCE=false

for arg in "$@"; do
  case "$arg" in
    --force) FORCE=true ;;
  esac
done

mkdir -p "$CLINVAR_DIR" "$OUTDIR"

SUMMARY_URL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
MANE_URL="https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/MANE.GRCh38.v1.5.summary.txt.gz"

SUMMARY_PATH="${CLINVAR_DIR}/variant_summary.txt.gz"
MANE_PATH="${CLINVAR_DIR}/MANE.GRCh38.v1.5.summary.txt.gz"

# ── Download variant_summary.txt.gz ──────────────────────────────────────────
if [[ "$FORCE" == "true" || ! -f "$SUMMARY_PATH" ]]; then
  echo "[prepare_clinvar] Downloading variant_summary.txt.gz (~435 MB) ..."
  curl -fsSL -o "$SUMMARY_PATH" "$SUMMARY_URL"
  echo "[prepare_clinvar] Downloaded to $SUMMARY_PATH"
else
  echo "[prepare_clinvar] Using cached $SUMMARY_PATH"
fi

# ── Download MANE summary ─────────────────────────────────────────────────────
if [[ "$FORCE" == "true" || ! -f "$MANE_PATH" ]]; then
  echo "[prepare_clinvar] Downloading MANE summary (~1 MB) ..."
  curl -fsSL -o "$MANE_PATH" "$MANE_URL"
  echo "[prepare_clinvar] Downloaded to $MANE_PATH"
else
  echo "[prepare_clinvar] Using cached $MANE_PATH"
fi

OUT_VCF="${OUTDIR}/clinvar_patho.vcf.gz"

if [[ "$FORCE" == "true" || ! -f "$OUT_VCF" ]]; then
  echo "[prepare_clinvar] Filtering variant_summary → $OUT_VCF ..."

  # Column indices (1-based for awk):
  #  2=Type  3=Name  5=GeneSymbol  7=ClinicalSignificance  25=ReviewStatus
  # 17=Assembly 19=Chromosome 32=PositionVCF 33=ReferenceAlleleVCF 34=AlternateAlleleVCF
  #
  # Filters applied:
  #   Assembly == GRCh38
  #   Type in (single nucleotide variant, Deletion, Insertion, Indel)
  #   ClinicalSignificance contains "athogenic" (catches Pathogenic, Likely_pathogenic)
  #     but NOT "not Pathogenic" or "Conflicting"
  #   Name contains "(p." (has protein-level HGVS)
  #   PositionVCF is a positive integer
  #   Ref and Alt are non-empty nucleotide sequences (no "N/A" or symbolic alleles)
  #
  # Output VCF INFO fields:
  #   CLNP=p.Arg273Cys   (URL-encoded protein change)
  #   CLNTX=NM_000546    (unversioned RefSeq transcript)
  #   CLNG=TP53          (gene symbol)
  #   CLNSIG=Pathogenic  (clinical significance)

  zcat "$SUMMARY_PATH" | awk -F'\t' '
NR==1 {
  # Find column indices from header
  for (i=1; i<=NF; i++) {
    h = $i
    sub(/^#/,"",h)
    col[h] = i
  }
  # Verify required columns exist
  if (!col["Assembly"] || !col["Type"] || !col["Name"] ||
      !col["ClinicalSignificance"] || !col["Chromosome"] ||
      !col["PositionVCF"] || !col["ReferenceAlleleVCF"] || !col["AlternateAlleleVCF"]) {
    print "ERROR: missing required column in variant_summary header" > "/dev/stderr"
    exit 1
  }
  next
}
{
  assembly  = $col["Assembly"]
  vartype   = $col["Type"]
  name      = $col["Name"]
  sig       = $col["ClinicalSignificance"]
  chrom     = $col["Chromosome"]
  posvcf    = $col["PositionVCF"]
  ref       = $col["ReferenceAlleleVCF"]
  alt       = $col["AlternateAlleleVCF"]
  gene      = $col["GeneSymbol"]

  # Filters
  if (assembly != "GRCh38") next
  if (vartype != "single nucleotide variant" &&
      vartype != "Deletion" && vartype != "Insertion" && vartype != "Indel") next
  if (sig !~ /athogenic/) next
  if (sig ~ /not Pathogenic|Conflicting/) next
  if (name !~ /\(p\./) next
  if (posvcf+0 <= 0) next
  if (ref ~ /^[Nn][Aa]$/ || alt ~ /^[Nn][Aa]$/) next
  if (ref !~ /^[ACGTNacgtn]+$/ || alt !~ /^[ACGTNacgtn]+$/) next

  # Extract protein change: everything between "(p." and ")"
  if (match(name, /\(p\.([^)]+)\)/, arr)) {
    protein = "p." arr[1]
  } else next
  if (protein == "p.?" || protein ~ /\?$/) next

  # Extract unversioned RefSeq transcript
  tx = name
  if (match(tx, /^([A-Z]{2}_[0-9]+)/, arr)) {
    tx = arr[1]
  } else {
    tx = "."
  }

  # URL-encode protein (replace spaces with %20, commas with %2C)
  gsub(/ /, "%20", protein)
  gsub(/,/, "%2C", protein)

  # Emit VCF line
  info = "CLNP=" protein ";CLNTX=" tx ";CLNG=" gene ";CLNSIG=" sig
  print chrom "\t" posvcf "\t.\t" ref "\t" alt "\t.\t.\t" info
}
' | sort -k1,1V -k2,2n | gzip > "$OUT_VCF"

  COUNT=$(zcat "$OUT_VCF" | wc -l)
  echo "[prepare_clinvar] Generated $OUT_VCF with $COUNT variants"
else
  echo "[prepare_clinvar] Using cached $OUT_VCF"
fi

echo "[prepare_clinvar] Done."
echo ""
echo "Next steps:"
echo "  bash scripts/run_snpeff_clinvar.sh"
echo "  # (optionally) bash scripts/run_vep_clinvar.sh"
echo "  CC=~/miniforge3/envs/gcc/bin/x86_64-conda-linux-gnu-cc CGO_ENABLED=1 \\"
echo "    go test ./internal/output/ -run TestClinVarBenchmark -v -timeout 60m"
