#!/usr/bin/env bash
# run_vep.sh — Convert TCGA MAF files to VCF and annotate with Ensembl VEP.
#
# Output VCFs are written to testdata/tcga/vep/<study>.vcf.gz.
# Run this script once; then use the benchmark test:
#
#   go test ./internal/output/ -run TestSnpEffBenchmark -v -count=1
#
# Requirements:
#   - Ensembl VEP (vep command, or set --vep-cmd / VEP_CMD)
#   - VEP cache for GRCh38 in ~/.vep (install with INSTALL.pl --AUTO cf --CACHEDIR ~/.vep)
#
# Usage:
#   scripts/run_vep.sh [--vep-cmd CMD] [--cache-dir DIR] [--threads N] [--studies A,B,...]
#
# Options:
#   --vep-cmd CMD     Path or command for vep  [default: vep or $VEP_CMD]
#   --cache-dir DIR   VEP cache directory       [default: ~/.vep or $VEP_CACHE_DIR]
#   --assembly ASM    Genome assembly           [default: GRCh38]
#   --threads N       Parallel forks for VEP   [default: 4]
#   --studies A,B,... Comma-separated studies   [default: all TCGA studies]
#   --force           Re-annotate even if output exists
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
TCGA_DIR="$REPO_ROOT/testdata/tcga"
VEP_DIR="$TCGA_DIR/vep"

MAMBA="${MAMBA:-$(command -v mamba || command -v conda || echo "")}"
MAMBA_ENV="${MAMBA_ENV:-ensembl-vep}"
VEP_CMD="${VEP_CMD:-}"
VEP_CACHE_DIR="${VEP_CACHE_DIR:-$HOME/.vep}"
ASSEMBLY="${ASSEMBLY:-GRCh38}"
THREADS="${THREADS:-4}"
FORCE=0

ALL_STUDIES=(
  blca_tcga_gdc
  brca_tcga_gdc
  chol_tcga_gdc
  coad_tcga_gdc
  gbm_tcga_gdc
  luad_tcga_gdc
  skcm_tcga_gdc
)
STUDIES=("${ALL_STUDIES[@]}")

# Parse arguments.
while [[ $# -gt 0 ]]; do
  case "$1" in
    --vep-cmd)    VEP_CMD="$2"; shift 2 ;;
    --mamba-env)  MAMBA_ENV="$2"; shift 2 ;;
    --cache-dir)  VEP_CACHE_DIR="$2"; shift 2 ;;
    --assembly)   ASSEMBLY="$2"; shift 2 ;;
    --threads)    THREADS="$2"; shift 2 ;;
    --force)      FORCE=1; shift ;;
    --studies)    IFS=',' read -ra STUDIES <<< "$2"; shift 2 ;;
    -h|--help)
      sed -n '2,/^set/p' "$0" | grep '^#' | sed 's/^# \{0,1\}//'
      exit 0 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

# Build the VEP invocation array.
# Prefer: explicit VEP_CMD > mamba run -n MAMBA_ENV vep > bare vep.
declare -a VEP_RUN
if [[ -n "$VEP_CMD" ]]; then
  VEP_RUN=("$VEP_CMD")
elif [[ -n "$MAMBA" && -n "$MAMBA_ENV" ]]; then
  VEP_RUN=("$MAMBA" "run" "-n" "$MAMBA_ENV" "vep")
else
  VEP_RUN=("vep")
fi

# Validate VEP is runnable.
if ! "${VEP_RUN[@]}" --help &>/dev/null; then
  echo "ERROR: VEP not runnable (tried: ${VEP_RUN[*]})" >&2
  echo "       Install: mamba create -n ensembl-vep -c bioconda ensembl-vep" >&2
  exit 1
fi

if [[ ! -d "$VEP_CACHE_DIR" ]]; then
  echo "ERROR: VEP cache directory not found: $VEP_CACHE_DIR" >&2
  echo "       Download cache with INSTALL.pl:" >&2
  echo "         conda run -n ensembl-vep vep_install -a cf -s homo_sapiens -y GRCh38 --CACHEDIR ~/.vep" >&2
  exit 1
fi

if ! command -v python3 &>/dev/null; then
  echo "ERROR: python3 not found (required for MAF → VCF conversion)." >&2
  exit 1
fi

VEP_VERSION=$("${VEP_RUN[@]}" --help 2>&1 | grep -oP 'VEP version \K[0-9]+' | head -1 || echo "unknown")

echo "VEP command : ${VEP_RUN[*]} (version $VEP_VERSION)"
echo "Cache dir   : $VEP_CACHE_DIR"
echo "Assembly    : $ASSEMBLY"
echo "Output dir  : $VEP_DIR"
echo "Threads     : $THREADS"
echo ""

mkdir -p "$VEP_DIR"

# Python helper: converts a TCGA GDC MAF file to a minimal VCF suitable for VEP.
# (Same MAF→VCF logic as run_snpeff.sh.)
MAF2VCF_PY='
import sys, gzip, re

def open_maybe_gz(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)

def maf_to_vcf(maf_path, vcf_path):
    seen = set()
    records = []

    with open_maybe_gz(maf_path) as fh:
        header = None
        col = {}
        for line in fh:
            line = line.rstrip("\n\r")
            if line.startswith("#"):
                continue
            if header is None:
                header = line.split("\t")
                col = {h: i for i, h in enumerate(header)}
                for req in ("Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2"):
                    if req not in col:
                        sys.exit(f"Missing required MAF column: {req}")
                continue
            f = line.split("\t")
            if len(f) <= max(col["Chromosome"], col["Start_Position"],
                             col["Reference_Allele"], col["Tumor_Seq_Allele2"]):
                continue
            chrom = f[col["Chromosome"]]
            if not chrom.startswith("chr"):
                chrom = "chr" + chrom
            try:
                start = int(f[col["Start_Position"]])
            except ValueError:
                continue
            ref_allele = f[col["Reference_Allele"]]
            alt_allele = f[col["Tumor_Seq_Allele2"]]

            if ref_allele == "-":
                pos = start
                vcf_ref = "N"
                vcf_alt = "N" + alt_allele
            elif alt_allele == "-":
                pos = start - 1
                vcf_ref = "N" + ref_allele
                vcf_alt = "N"
            else:
                pos = start
                vcf_ref = ref_allele
                vcf_alt = alt_allele

            key = (chrom, pos, vcf_ref, vcf_alt)
            if key in seen:
                continue
            seen.add(key)
            records.append((chrom, pos, vcf_ref, vcf_alt))

    records.sort(key=lambda r: (chrom_sort_key(r[0]), r[1]))

    with gzip.open(vcf_path, "wt") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write(f"##source=run_vep.sh (MAF->VCF)\n")
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom, pos, ref, alt in records:
            out.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\n")

    print(f"  Wrote {len(records)} unique variants to {vcf_path}", file=sys.stderr)

def chrom_sort_key(c):
    c = c.lstrip("chr")
    try:
        return (0, int(c), "")
    except ValueError:
        return (1, 0, c)

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("maf")
    ap.add_argument("vcf")
    args = ap.parse_args()
    maf_to_vcf(args.maf, args.vcf)
'

MAF2VCF_SCRIPT=$(mktemp /tmp/maf2vcf_vep_XXXXXX.py)
echo "$MAF2VCF_PY" > "$MAF2VCF_SCRIPT"
trap 'rm -f "$MAF2VCF_SCRIPT"' EXIT

run_study() {
  local study="$1"
  local maf_file="$TCGA_DIR/${study}_data_mutations.txt"
  local out_vcf="$VEP_DIR/${study}.vcf.gz"

  if [[ ! -f "$maf_file" ]]; then
    echo "[$study] Skipping: MAF not found ($maf_file)" >&2
    echo "         Run: make download-testdata" >&2
    return
  fi

  if [[ -f "$out_vcf" && "$FORCE" -eq 0 ]]; then
    echo "[$study] Already exists, skipping (use --force to re-run)"
    return
  fi

  echo "[$study] Converting MAF → VCF..."
  local input_vcf
  input_vcf=$(mktemp /tmp/vep_input_XXXXXX.vcf.gz)
  python3 "$MAF2VCF_SCRIPT" "$maf_file" "$input_vcf"

  echo "[$study] Running VEP..."
  "${VEP_RUN[@]}" \
    --input_file "$input_vcf" \
    --output_file "$out_vcf" \
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

  rm -f "$input_vcf"

  if [[ -f "$out_vcf" ]]; then
    local size
    size=$(du -sh "$out_vcf" | cut -f1)
    echo "[$study] Done → $out_vcf ($size)"
  else
    echo "[$study] ERROR: output VCF not created" >&2
    return 1
  fi
}

export -f run_study
export TCGA_DIR VEP_DIR VEP_CACHE_DIR ASSEMBLY MAF2VCF_SCRIPT FORCE

for study in "${STUDIES[@]}"; do
  run_study "$study"
done

echo ""
echo "All studies processed. Run the benchmark:"
echo "  go test ./internal/output/ -run TestSnpEffBenchmark -v -count=1"
