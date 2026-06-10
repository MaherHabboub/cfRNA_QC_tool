#!/bin/bash

set -euo pipefail

CONFIG="${1:-}"

if [[ -z "$CONFIG" ]]; then
    echo "Usage: bash GTF_to_BED12.sh path/to/config.sh"
    exit 1
fi

if [[ ! -f "$CONFIG" ]]; then
    echo "ERROR: Config file not found: $CONFIG"
    exit 1
fi

source "$CONFIG"

# -----------------------------
# Required config variables
# -----------------------------
: "${GTF:?ERROR: GTF not set in config}"
: "${OUTDIR:?ERROR: OUTDIR not set in config}"

module purge
module load Kent_tools/479-GCC-13.3.0

ANNOT_DIR="${OUTDIR}/annotation"
mkdir -p "$ANNOT_DIR"

# Use GTF filename as prefix, minus .gtf/.gff
BASENAME="$(basename "$GTF")"
PREFIX="${BASENAME%.gtf}"
PREFIX="${PREFIX%.gff}"

GENEPRED="${ANNOT_DIR}/${PREFIX}.genePred"
BED12_UNSORTED="${ANNOT_DIR}/${PREFIX}.bed12.unsorted.bed"
BED12="${ANNOT_DIR}/${PREFIX}.bed12.bed"

echo "Using gtfToGenePred: $(command -v gtfToGenePred)"
echo "Using genePredToBed: $(command -v genePredToBed)"
echo "Input GTF: $GTF"
echo "Output BED12: $BED12"
echo

[[ -f "$GTF" ]] || {
    echo "ERROR: GTF not found: $GTF"
    exit 1
}

# Convert GTF -> genePred (keeps CDS info)
gtfToGenePred -genePredExt -allErrors "$GTF" "$GENEPRED"

# Convert genePred -> BED12
genePredToBed "$GENEPRED" "$BED12_UNSORTED"

# Sort BED
sort -k1,1 -k2,2n "$BED12_UNSORTED" > "$BED12"

# Quick validation: BED12 must be 12 columns and thickStart/thickEnd numeric
echo "Validation (first line):"
awk 'BEGIN{FS="\t"} NR==1{print "NF="NF, "thickStart(col7)="$7, "thickEnd(col8)="$8}' "$BED12"

BAD_LINE=$(awk 'BEGIN{FS="\t"} ($7 !~ /^[0-9]+$/ || $8 !~ /^[0-9]+$/){print; exit}' "$BED12" || true)

if [[ -n "${BAD_LINE}" ]]; then
  echo "ERROR: Found non-numeric thickStart/thickEnd in BED12 (this will break RSeQC)."
  echo "First offending line:"
  echo "$BAD_LINE"
  exit 2
fi

# Write BED12 path to a small helper file for the master script / user
echo "$BED12" > "${ANNOT_DIR}/BED12.path.txt"

echo
echo "Done."
echo "Created: $BED12"
echo "Path saved to: ${ANNOT_DIR}/BED12.path.txt"