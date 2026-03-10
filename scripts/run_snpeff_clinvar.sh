#!/usr/bin/env bash
# run_snpeff_clinvar.sh — Annotate ClinVar pathogenic VCF with snpEff.
#
# Input:  testdata/clinvar/clinvar_patho.vcf.gz  (from prepare_clinvar.sh)
# Output: testdata/clinvar/snpeff.vcf.gz
#
# Usage:
#   bash scripts/run_snpeff_clinvar.sh [--snpeff-jar PATH] [--db GENOME] [--config FILE] [--force]
#
# Options:
#   --snpeff-jar PATH   Path to snpEff.jar  [default: ~/snpEff/snpEff.jar or $SNPEFF_JAR]
#   --db GENOME         snpEff genome database  [default: GRCh38.115]
#   --config FILE       Path to custom snpEff.config (e.g. /tmp/snpEff_custom.config)
#   --force             Re-annotate even if output exists
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

CLINVAR_DIR="$REPO_ROOT/testdata/clinvar"
INPUT_VCF="$CLINVAR_DIR/clinvar_patho.vcf.gz"
OUTPUT_VCF="$CLINVAR_DIR/snpeff.vcf.gz"

SNPEFF_DB="${SNPEFF_DB:-GRCh38.115}"
SNPEFF_DATA_DIR="${SNPEFF_DATA_DIR:-}"
SNPEFF_CONFIG="${SNPEFF_CONFIG:-}"
FORCE=0

# Locate snpEff.jar
if [[ -z "${SNPEFF_JAR:-}" ]]; then
  for candidate in \
    "$HOME/snpEff/snpEff.jar" \
    "$(dirname "$0")/snpEff.jar" \
    "/opt/snpEff/snpEff.jar" \
    "/usr/local/share/snpEff/snpEff.jar"; do
    if [[ -f "$candidate" ]]; then
      SNPEFF_JAR="$candidate"
      break
    fi
  done
fi

for arg in "$@"; do
  case "$arg" in
    --snpeff-jar) shift; SNPEFF_JAR="$1" ;;
    --db)         shift; SNPEFF_DB="$1" ;;
    --config)     shift; SNPEFF_CONFIG="$1" ;;
    --force)      FORCE=1 ;;
  esac
done

if [[ -z "${SNPEFF_JAR:-}" || ! -f "$SNPEFF_JAR" ]]; then
  echo "ERROR: snpEff.jar not found. Set SNPEFF_JAR or pass --snpeff-jar." >&2
  exit 1
fi

if [[ ! -f "$INPUT_VCF" ]]; then
  echo "ERROR: $INPUT_VCF not found. Run scripts/prepare_clinvar.sh first." >&2
  exit 1
fi

if [[ -f "$OUTPUT_VCF" && "$FORCE" -eq 0 ]]; then
  echo "[clinvar-snpeff] $OUTPUT_VCF already exists (use --force to re-annotate)"
  exit 0
fi

echo "[clinvar-snpeff] Annotating $(zcat "$INPUT_VCF" | wc -l) variants with $SNPEFF_DB ..."
echo "[clinvar-snpeff] Using snpEff: $SNPEFF_JAR"

config_arg=()
if [[ -n "${SNPEFF_CONFIG:-}" ]]; then
  config_arg=(-config "$SNPEFF_CONFIG")
fi

data_dir_arg=()
if [[ -n "${SNPEFF_DATA_DIR:-}" ]]; then
  data_dir_arg=(-dataDir "$SNPEFF_DATA_DIR")
fi

# snpEff cannot read /dev/stdin; write input to a temp file.
TMP_VCF=$(mktemp /tmp/clinvar_input_XXXXXX.vcf)
trap 'rm -f "$TMP_VCF"' EXIT

{
  printf '##fileformat=VCFv4.1\n'
  printf '##source=prepare_clinvar.sh\n'
  printf '##INFO=<ID=CLNP,Number=1,Type=String,Description="ClinVar protein HGVS">\n'
  printf '##INFO=<ID=CLNTX,Number=1,Type=String,Description="ClinVar RefSeq transcript">\n'
  printf '##INFO=<ID=CLNG,Number=1,Type=String,Description="ClinVar gene symbol">\n'
  printf '##INFO=<ID=CLNSIG,Number=1,Type=String,Description="ClinVar clinical significance">\n'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
  zcat "$INPUT_VCF"
} > "$TMP_VCF"

java -Xmx12g -jar "$SNPEFF_JAR" ann \
    "${config_arg[@]}" \
    "${data_dir_arg[@]}" \
    -v \
    -noStats \
    -noLog \
    -nodownload \
    "$SNPEFF_DB" \
    "$TMP_VCF" \
  | gzip -c > "$OUTPUT_VCF"

COUNT=$(zcat "$OUTPUT_VCF" | grep -v '^#' | wc -l)
echo "[clinvar-snpeff] Done → $OUTPUT_VCF ($COUNT variants annotated)"
