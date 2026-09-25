#!/bin/bash

set -euo pipefail

CONFIG="${1:-}"
TARGET_SAMPLE="${2:-ALL}"

if [[ -z "$CONFIG" ]]; then
    echo "Usage: bash Read_Distribution.sh path/to/config.sh [sample_id]"
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
: "${SAMPLESHEET:?ERROR: SAMPLESHEET not set in config}"
: "${OUTDIR:?ERROR: OUTDIR not set in config}"

# -----------------------------
# Paths
# -----------------------------
BED12_PATH_FILE="${OUTDIR}/annotation/BED12.path.txt"
RESULT_DIR="${OUTDIR}/read_distribution"

# -----------------------------
# Software environment
# -----------------------------
module purge
module load RSeQC/5.0.1-foss-2023a
# Match RSeQC's foss/2023a compiler toolchain.
module load SAMtools/1.18-GCC-12.3.0

command -v samtools >/dev/null 2>&1 || {
    echo "ERROR: samtools is unavailable after loading SAMtools." >&2
    exit 1
}

# -----------------------------
# Validate prerequisites
# -----------------------------
if [[ ! -s "$BED12_PATH_FILE" ]]; then
    echo "ERROR: Generated BED12 path file is missing or empty: $BED12_PATH_FILE" >&2
    echo "Run GTF_to_BED12.sh successfully before this module." >&2
    exit 1
fi

BED12="$(<"$BED12_PATH_FILE")"

if [[ -z "$BED12" || ! -f "$BED12" ]]; then
    echo "ERROR: Generated BED12 file is missing: ${BED12:-<empty path>}" >&2
    echo "Path was read from: $BED12_PATH_FILE" >&2
    exit 1
fi

mkdir -p "$RESULT_DIR"

echo "Running read distribution QC..."
echo "Target sample: $TARGET_SAMPLE"
echo "Generated BED12: $BED12"

tail -n +2 "$SAMPLESHEET" | while IFS=$'\t' read -r SAMPLE FASTQ1 FASTQ2 BAM STARLOG SJTAB LAYOUT CONDITION TRANSCRIPTOME_BAM
do
    if [[ "$TARGET_SAMPLE" != "ALL" && "$SAMPLE" != "$TARGET_SAMPLE" ]]; then
        continue
    fi


    echo "------------------------------------"
    echo "Processing: $SAMPLE"

    if [[ ! -f "$BAM" ]]; then
        echo "WARNING: BAM not found for $SAMPLE, skipping"
        continue
    fi

    samtools quickcheck -v "$BAM" || {
        echo "ERROR: BAM failed samtools quickcheck for $SAMPLE: $BAM" >&2
        exit 1
    }

    if ! samtools view -H "$BAM" | awk '$1 == "@SQ" {found=1} END {exit !found}'; then
        echo "ERROR: BAM header contains no @SQ reference-sequence records for $SAMPLE: $BAM" >&2
        exit 1
    fi

    SAMPLE_OUTDIR="${RESULT_DIR}/${SAMPLE}"
    mkdir -p "$SAMPLE_OUTDIR"

    OUTTXT="${SAMPLE_OUTDIR}/${SAMPLE}.read_distribution.txt"
    SORTED_INPUT=""
    INPUT_BAM="$BAM"
    sort_order="$(samtools view -H "$BAM" | awk '$1 == "@HD" {for (i = 1; i <= NF; i++) if ($i ~ /^SO:/) {print substr($i, 4); exit}}')"

    if [[ "$sort_order" != "coordinate" ]]; then
        SORTED_INPUT="${SAMPLE_OUTDIR}/${SAMPLE}.coordinate_sorted.tmp.bam"
        echo "Input BAM is not declared coordinate-sorted (SO=${sort_order:-unspecified}); creating a temporary coordinate-sorted BAM."
        samtools sort -o "$SORTED_INPUT" "$BAM"
        INPUT_BAM="$SORTED_INPUT"
    fi

    read_distribution.py \
      -i "$INPUT_BAM" \
      -r "$BED12" \
      > "$OUTTXT"

    [[ -z "$SORTED_INPUT" ]] || rm -f "$SORTED_INPUT"

    echo "Done: $SAMPLE"
    echo "Wrote: $OUTTXT"

done

echo "Read distribution QC complete."
