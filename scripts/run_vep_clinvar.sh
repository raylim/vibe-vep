#!/usr/bin/env bash
# run_vep_clinvar.sh — Annotate ClinVar pathogenic VCF with Ensembl VEP.
#
# Input:  testdata/clinvar/clinvar_patho.vcf.gz  (from prepare_clinvar.sh)
# Output: testdata/clinvar/vep.vcf.gz
#
# Usage:
#   bash scripts/run_vep_clinvar.sh [--vep-cmd CMD] [--cache-dir DIR] [--threads N] [--force]
#
# Options:
#   --vep-cmd CMD   Path or command for vep  [default: vep or $VEP_CMD]
#   --cache-dir DIR VEP cache directory      [default: ~/.vep or $VEP_CACHE_DIR]
#   --assembly ASM  Genome assembly          [default: GRCh38]
#   --threads N     Parallel forks for VEP  [default: 4]
#   --force         Re-annotate even if output exists
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

CLINVAR_DIR="$REPO_ROOT/testdata/clinvar"
INPUT_VCF="$CLINVAR_DIR/clinvar_patho.vcf.gz"
OUTPUT_VCF="$CLINVAR_DIR/vep.vcf.gz"

MAMBA="${MAMBA:-$(command -v mamba 2>/dev/null || command -v conda 2>/dev/null || echo "")}"
MAMBA_ENV="${MAMBA_ENV:-ensembl-vep}"
VEP_CMD="${VEP_CMD:-}"
VEP_CACHE_DIR="${VEP_CACHE_DIR:-$HOME/.vep}"
ASSEMBLY="${ASSEMBLY:-GRCh38}"
THREADS="${THREADS:-4}"
FORCE=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --vep-cmd)   VEP_CMD="$2"; shift 2 ;;
    --mamba-env) MAMBA_ENV="$2"; shift 2 ;;
    --cache-dir) VEP_CACHE_DIR="$2"; shift 2 ;;
    --assembly)  ASSEMBLY="$2"; shift 2 ;;
    --threads)   THREADS="$2"; shift 2 ;;
    --force)     FORCE=1; shift ;;
    -h|--help)
      sed -n '2,/^set/p' "$0" | grep '^#' | sed 's/^# \{0,1\}//'
      exit 0 ;;
    *) echo "Unknown option: $1" >&2; exit 1 ;;
  esac
done

# Resolve vep command.
VEP_RUN=()
if [[ -n "$VEP_CMD" ]]; then
  VEP_RUN=("$VEP_CMD")
elif command -v vep &>/dev/null; then
  VEP_RUN=(vep)
elif [[ -n "$MAMBA" ]]; then
  VEP_RUN=("$MAMBA" run -n "$MAMBA_ENV" vep)
else
  echo "ERROR: vep not found. Set --vep-cmd or install in a conda env named '$MAMBA_ENV'." >&2
  exit 1
fi

if [[ ! -f "$INPUT_VCF" ]]; then
  echo "ERROR: $INPUT_VCF not found. Run scripts/prepare_clinvar.sh first." >&2
  exit 1
fi

if [[ -f "$OUTPUT_VCF" && "$FORCE" -eq 0 ]]; then
  echo "[clinvar-vep] $OUTPUT_VCF already exists (use --force to re-annotate)"
  exit 0
fi

# Remove existing output so VEP doesn't skip due to its own file-exists check.
[[ -f "$OUTPUT_VCF" ]] && rm -f "$OUTPUT_VCF"

COUNT=$(zcat "$INPUT_VCF" | wc -l)
echo "[clinvar-vep] Annotating $COUNT variants with Ensembl VEP (assembly=$ASSEMBLY, forks=$THREADS) ..."
echo "[clinvar-vep] VEP command: ${VEP_RUN[*]}"

# Build a proper VCF with header from the headerless input.
TMP_VCF=$(mktemp /tmp/clinvar_vep_input_XXXXXX.vcf)
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

ANNO_START=$SECONDS
"${VEP_RUN[@]}" \
  --input_file "$TMP_VCF" \
  --output_file "$OUTPUT_VCF" \
  --vcf \
  --compress_output bgzip \
  --offline \
  --cache \
  --dir_cache "$VEP_CACHE_DIR" \
  --assembly "$ASSEMBLY" \
  --species homo_sapiens \
  --fork "$THREADS" \
  --no_stats \
  --hgvs \
  --transcript_version \
  --canonical \
  --biotype \
  --numbers \
  --quiet \
  2>&1 | grep -v "^#" | tail -5 || true
ELAPSED=$(( SECONDS - ANNO_START ))

if [[ -f "$OUTPUT_VCF" ]]; then
  ANNOTATED=$(zcat "$OUTPUT_VCF" | grep -v '^#' | wc -l)
  echo "$ELAPSED" > "${OUTPUT_VCF%.vcf.gz}.elapsed"
  echo "[clinvar-vep] Done → $OUTPUT_VCF ($ANNOTATED variants annotated, ${ELAPSED}s)"
else
  echo "[clinvar-vep] ERROR: output VCF not created" >&2
  exit 1
fi
