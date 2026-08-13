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

tail -n +2 "$SAMPLESHEET" | while IFS=$'\t' read -r SAMPLE FASTQ1 FASTQ2 BAM STARLOG SJTAB LAYOUT CONDITION
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

    SAMPLE_OUTDIR="${RESULT_DIR}/${SAMPLE}"
    mkdir -p "$SAMPLE_OUTDIR"

    OUTTXT="${SAMPLE_OUTDIR}/${SAMPLE}.read_distribution.txt"

    read_distribution.py \
      -i "$BAM" \
      -r "$BED12" \
      > "$OUTTXT"

    echo "Done: $SAMPLE"
    echo "Wrote: $OUTTXT"

done

echo "Read distribution QC complete."
