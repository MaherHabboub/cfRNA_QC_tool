#!/bin/bash

set -euo pipefail

CONFIG="${1:-}"

if [[ -z "$CONFIG" ]]; then
    echo "Usage: bash QC_read_dist.sh path/to/config.sh"
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
: "${BED12:?ERROR: BED12 not set in config}"

module purge
module load RSeQC/5.0.1-foss-2023a

RESULT_DIR="${OUTDIR}/read_distribution"
mkdir -p "$RESULT_DIR"

echo "Running read distribution QC..."

tail -n +2 "$SAMPLESHEET" | while IFS=$'\t' read -r SAMPLE FASTQ1 FASTQ2 BAM STARLOG SJTAB LAYOUT CONDITION
do

    echo "------------------------------------"
    echo "Processing: $SAMPLE"

    if [[ ! -f "$BAM" ]]; then
        echo "WARNING: BAM not found for $SAMPLE, skipping"
        continue
    fi

    if [[ ! -f "$BED12" ]]; then
        echo "ERROR: BED12 not found: $BED12"
        exit 1
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