#!/usr/bin/env bash
# run_annovar_clinvar.sh — Annotate ClinVar pathogenic VCF with ANNOVAR.
#
# Input:  testdata/clinvar/clinvar_patho.vcf.gz  (from prepare_clinvar.sh)
# Output: testdata/clinvar/annovar.txt.gz
#
# Usage:
#   bash scripts/run_annovar_clinvar.sh [--annovar-dir DIR] [--db BUILDVER] [--force]
#
# Options:
#   --annovar-dir DIR   Path to ANNOVAR directory containing table_annovar.pl  [default: ~/annovar]
#   --db BUILDVER       Genome build version  [default: hg38]
#   --force             Re-annotate even if output exists
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

CLINVAR_DIR="$REPO_ROOT/testdata/clinvar"
INPUT_VCF="$CLINVAR_DIR/clinvar_patho.vcf.gz"
OUTPUT_TXT="$CLINVAR_DIR/annovar.txt.gz"

ANNOVAR_DIR="${ANNOVAR_DIR:-$HOME/annovar}"
BUILDVER="${BUILDVER:-hg38}"
PROTOCOL="refGeneWithVer"
FORCE=0

for arg in "$@"; do
  case "$arg" in
    --annovar-dir) shift; ANNOVAR_DIR="$1" ;;
    --db)          shift; BUILDVER="$1" ;;
    --force)       FORCE=1 ;;
  esac
done

TABLE_ANNOVAR="$ANNOVAR_DIR/table_annovar.pl"
HUMANDB="$ANNOVAR_DIR/humandb"

if [[ ! -f "$TABLE_ANNOVAR" ]]; then
  echo "ERROR: table_annovar.pl not found at $TABLE_ANNOVAR. Set ANNOVAR_DIR or pass --annovar-dir." >&2
  exit 1
fi

if [[ ! -f "$INPUT_VCF" ]]; then
  echo "ERROR: $INPUT_VCF not found. Run scripts/prepare_clinvar.sh first." >&2
  exit 1
fi

if [[ ! -f "$HUMANDB/${BUILDVER}_${PROTOCOL}.txt" ]]; then
  echo "ERROR: Database $HUMANDB/${BUILDVER}_${PROTOCOL}.txt not found." >&2
  exit 1
fi

if [[ -f "$OUTPUT_TXT" && "$FORCE" -eq 0 ]]; then
  echo "[clinvar-annovar] $OUTPUT_TXT already exists (use --force to re-annotate)"
  exit 0
fi

N=$(zcat "$INPUT_VCF" | wc -l)
echo "[clinvar-annovar] Annotating $N variants with ANNOVAR $BUILDVER $PROTOCOL ..."
echo "[clinvar-annovar] Using ANNOVAR: $TABLE_ANNOVAR"

# ANNOVAR requires a proper VCF header for --vcfinput.
TMP_VCF=$(mktemp /tmp/clinvar_annovar_XXXXXX.vcf)
TMP_OUT=$(mktemp -d /tmp/annovar_out_XXXXXX)
trap 'rm -f "$TMP_VCF"; rm -rf "$TMP_OUT"' EXIT

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

ANNO_START=$SECONDS
perl "$TABLE_ANNOVAR" "$TMP_VCF" "$HUMANDB" \
  -buildver "$BUILDVER" \
  -out "$TMP_OUT/out" \
  -remove \
  -protocol "$PROTOCOL" \
  -operation g \
  -nastring . \
  -vcfinput \
  -polish \
  2>&1
ELAPSED=$(( SECONDS - ANNO_START ))

MULTIANNO_TXT="$TMP_OUT/out.${BUILDVER}_multianno.txt"
if [[ ! -f "$MULTIANNO_TXT" ]]; then
  echo "ERROR: Expected output not found: $MULTIANNO_TXT" >&2
  exit 1
fi

gzip -c "$MULTIANNO_TXT" > "$OUTPUT_TXT"
COUNT=$(zcat "$OUTPUT_TXT" | grep -v "^Chr" | wc -l)
echo "$ELAPSED" > "${OUTPUT_TXT%.txt.gz}.elapsed"
echo "[clinvar-annovar] Done → $OUTPUT_TXT ($COUNT variants annotated, ${ELAPSED}s)"
