#!/usr/bin/env bash
# run_snpeff.sh — Convert TCGA MAF files to VCF and annotate with snpEff.
#
# Output VCFs are written to testdata/tcga/snpeff/<study>.vcf.gz.
# Run this script once; then use the benchmark test:
#
#   go test ./internal/output/ -run TestSnpEffBenchmark -v -count=1
#
# Requirements:
#   - Java (java command on PATH)
#   - snpEff JAR (set SNPEFF_JAR or place snpEff.jar next to this script)
#   - snpEff GRCh38 database (downloaded automatically on first run)
#
# Usage:
#   scripts/run_snpeff.sh [--snpeff-jar PATH] [--db GENOME] [--data-dir DIR] [--config FILE] [--threads N]
#
# Options:
#   --snpeff-jar PATH   Path to snpEff.jar  [default: ~/snpEff/snpEff.jar or $SNPEFF_JAR]
#   --db GENOME         snpEff genome database  [default: GRCh38.99]
#   --data-dir DIR      snpEff data directory containing genome databases  [default: snpEff built-in]
#   --config FILE       Path to custom snpEff.config  [default: snpEff built-in]
#   --threads N         Parallel jobs  [default: 4]
#   --studies A,B,...   Comma-separated study names  [default: all TCGA studies]
#   --force             Re-annotate even if output exists
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
TCGA_DIR="$REPO_ROOT/testdata/tcga"
SNPEFF_DIR="$TCGA_DIR/snpeff"

# Default snpEff genome database (GRCh38 Ensembl 115, aligned with GENCODE v49 and VEP v115).
SNPEFF_DB="${SNPEFF_DB:-GRCh38.115}"
SNPEFF_DATA_DIR="${SNPEFF_DATA_DIR:-}"
SNPEFF_CONFIG="${SNPEFF_CONFIG:-}"
THREADS="${THREADS:-4}"
FORCE=0

# Locate snpEff.jar.
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
    --snpeff-jar) SNPEFF_JAR="$2"; shift 2 ;;
    --db)         SNPEFF_DB="$2"; shift 2 ;;
    --data-dir)   SNPEFF_DATA_DIR="$2"; shift 2 ;;
    --config)     SNPEFF_CONFIG="$2"; shift 2 ;;
    --threads)    THREADS="$2"; shift 2 ;;
    --force)      FORCE=1; shift ;;
    --studies)    IFS=',' read -ra STUDIES <<< "$2"; shift 2 ;;
    -h|--help)
      sed -n '2,/^set/p' "$0" | grep '^#' | sed 's/^# \{0,1\}//'
      exit 0 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

# Validate requirements.
if ! command -v java &>/dev/null; then
  echo "ERROR: java not found. Install Java 8+ and ensure it is on PATH." >&2
  exit 1
fi

if [[ -z "${SNPEFF_JAR:-}" || ! -f "$SNPEFF_JAR" ]]; then
  echo "ERROR: snpEff.jar not found." >&2
  echo "       Set SNPEFF_JAR or pass --snpeff-jar PATH." >&2
  echo "       Download: https://pcingola.github.io/SnpEff/#download" >&2
  exit 1
fi

if ! command -v python3 &>/dev/null; then
  echo "ERROR: python3 not found (required for MAF → VCF conversion)." >&2
  exit 1
fi

echo "snpEff JAR : $SNPEFF_JAR"
echo "Database   : $SNPEFF_DB"
echo "Data dir   : ${SNPEFF_DATA_DIR:-(default)}"
echo "Config     : ${SNPEFF_CONFIG:-(default)}"
echo "Output dir : $SNPEFF_DIR"
echo "Threads    : $THREADS"
echo ""

mkdir -p "$SNPEFF_DIR"

# Python helper: converts a TCGA GDC MAF file to a minimal VCF suitable for snpEff.
#
# Handles:
#   SNVs       — ref/alt are single bases, direct conversion.
#   Deletions  — MAF ref is non-dash, alt is "-".  VCF convention: prepend "N"
#                anchor and set pos = Start_Position - 1.  Without a reference
#                genome the anchor is "N"; snpEff normalises these correctly.
#   Insertions — MAF ref is "-", alt is the inserted sequence.  VCF convention:
#                pos = Start_Position, ref = "N", alt = "N" + inserted_seq.
#
# Note: snpEff handles "N" anchors when the reference FASTA is not provided.
#       For exact VCF compliance, pass a reference FASTA via --ref-fasta.
MAF2VCF_PY='
import sys, gzip, re

def open_maybe_gz(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)

def maf_to_vcf(maf_path, vcf_path, ref_fasta=None):
    """Convert a TCGA GDC MAF file to a sorted, unique VCF for snpEff."""
    # Load reference FASTA index for anchor bases (optional).
    fasta_seq = {}
    if ref_fasta:
        try:
            fasta_seq = load_fasta(ref_fasta)
        except Exception as e:
            print(f"Warning: could not load reference FASTA: {e}", file=sys.stderr)

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
                # Required columns.
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

            # Convert to VCF representation.
            if ref_allele == "-":
                # Insertion: MAF start = anchor position.
                pos = start
                anchor = get_anchor(fasta_seq, chrom, pos)
                vcf_ref = anchor
                vcf_alt = anchor + alt_allele
            elif alt_allele == "-":
                # Deletion: MAF start = first deleted base; VCF pos = start-1.
                pos = start - 1
                anchor = get_anchor(fasta_seq, chrom, pos)
                vcf_ref = anchor + ref_allele
                vcf_alt = anchor
            else:
                # SNV or complex.
                pos = start
                vcf_ref = ref_allele
                vcf_alt = alt_allele

            key = (chrom, pos, vcf_ref, vcf_alt)
            if key in seen:
                continue
            seen.add(key)
            records.append((chrom, pos, vcf_ref, vcf_alt))

    # Sort by chromosome (natural) then position.
    records.sort(key=lambda r: (chrom_sort_key(r[0]), r[1]))

    with gzip.open(vcf_path, "wt") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write(f"##source=run_snpeff.sh (MAF→VCF)\n")
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom, pos, ref, alt in records:
            out.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\n")

    print(f"  Wrote {len(records)} unique variants to {vcf_path}", file=sys.stderr)

def get_anchor(fasta_seq, chrom, pos):
    """Return the reference base at pos (1-based), or N if unavailable."""
    seq = fasta_seq.get(chrom) or fasta_seq.get(chrom.lstrip("chr"))
    if seq and 0 < pos <= len(seq):
        return seq[pos - 1].upper()
    return "N"

def chrom_sort_key(c):
    c = c.lstrip("chr")
    try:
        return (0, int(c), "")
    except ValueError:
        return (1, 0, c)

def load_fasta(path):
    """Load a FASTA file into a dict of chrom → sequence string."""
    seqs = {}
    current = None
    parts = []
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if current is not None:
                    seqs[current] = "".join(parts)
                current = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
    if current is not None:
        seqs[current] = "".join(parts)
    return seqs

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("maf")
    ap.add_argument("vcf")
    ap.add_argument("--ref-fasta", default=None)
    args = ap.parse_args()
    maf_to_vcf(args.maf, args.vcf, args.ref_fasta)
'

# Write the Python helper to a temp file so we can call it per-study.
MAF2VCF_SCRIPT=$(mktemp /tmp/maf2vcf_XXXXXX.py)
echo "$MAF2VCF_PY" > "$MAF2VCF_SCRIPT"
trap 'rm -f "$MAF2VCF_SCRIPT"' EXIT

# Process studies (optionally in parallel).
run_study() {
  local study="$1"
  local maf_file="$TCGA_DIR/${study}_data_mutations.txt"
  local out_vcf="$SNPEFF_DIR/${study}.vcf.gz"

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
  local input_vcf="$SNPEFF_DIR/${study}_input.vcf.gz"
  python3 "$MAF2VCF_SCRIPT" "$maf_file" "$input_vcf"

  echo "[$study] Running snpEff..."
  local config_arg=()
  if [[ -n "${SNPEFF_CONFIG:-}" ]]; then
    config_arg=(-config "$SNPEFF_CONFIG")
  fi
  local data_dir_arg=()
  if [[ -n "${SNPEFF_DATA_DIR:-}" ]]; then
    data_dir_arg=(-dataDir "$SNPEFF_DATA_DIR")
  fi
  java -Xmx12g -jar "$SNPEFF_JAR" ann \
    "${config_arg[@]}" \
    "${data_dir_arg[@]}" \
    -v \
    -noStats \
    -noLog \
    -nodownload \
    "$SNPEFF_DB" \
    "$input_vcf" \
    | gzip -c > "$out_vcf"

  rm -f "$input_vcf"
  echo "[$study] Done → $out_vcf"
}

export -f run_study
export TCGA_DIR SNPEFF_DIR SNPEFF_JAR SNPEFF_DB SNPEFF_DATA_DIR SNPEFF_CONFIG MAF2VCF_SCRIPT FORCE

if command -v parallel &>/dev/null && [[ "$THREADS" -gt 1 ]]; then
  printf '%s\n' "${STUDIES[@]}" | parallel -j "$THREADS" run_study {}
else
  for study in "${STUDIES[@]}"; do
    run_study "$study"
  done
fi

echo ""
echo "All studies processed. Run the benchmark:"
echo "  go test ./internal/output/ -run TestSnpEffBenchmark -v -count=1"
