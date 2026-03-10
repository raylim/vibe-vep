#!/usr/bin/env bash
#
# Download GRCh37 MAF files from cBioPortal/datahub for validation testing.
# Files are saved to testdata/grch37/.
#
set -euo pipefail

STUDIES=(
  msk_impact_50k_2026
)

BASE_URL="https://media.githubusercontent.com/media/cBioPortal/datahub/master/public"
DEST_DIR="$(cd "$(dirname "$0")/.." && pwd)/testdata/grch37"

mkdir -p "$DEST_DIR"

for study in "${STUDIES[@]}"; do
  dest="$DEST_DIR/${study}_data_mutations.txt"
  if [[ -f "$dest" ]]; then
    echo "Skipping $study (already exists)"
    continue
  fi
  url="$BASE_URL/$study/data_mutations.txt"
  echo "Downloading $study..."
  curl -fL --progress-bar -o "$dest" "$url"
done

echo "Done. Files saved to $DEST_DIR"
