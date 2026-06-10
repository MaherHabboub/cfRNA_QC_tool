#!/bin/bash

set -euo pipefail

CONFIG="${1:-}"

if [[ -z "$CONFIG" ]]; then
    echo "Usage: bash QC_strandedness.sh path/to/config.sh"
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
: "${EXON_BED:?ERROR: EXON_BED not set in config}"

module purge
module load RSeQC/5.0.1-foss-2023a

RESULT_DIR="${OUTDIR}/strandedness"
mkdir -p "$RESULT_DIR"

echo "Running strandedness QC..."

if [[ ! -f "$EXON_BED" ]]; then
    echo "ERROR: EXON_BED not found: $EXON_BED"
    exit 1
fi

tail -n +2 "$SAMPLESHEET" | while IFS=$'\t' read -r SAMPLE FASTQ1 FASTQ2 BAM STARLOG SJTAB LAYOUT CONDITION
do

    echo "------------------------------------"
    echo "Processing: $SAMPLE"

    if [[ ! -f "$BAM" ]]; then
        echo "WARNING: BAM not found for $SAMPLE, skipping"
        continue
    fi

    SAMPLE_OUTDIR="${RESULT_DIR}/${SAMPLE}"
    mkdir -p "$SAMPLE_OUTDIR"

    FULL_OUT="${SAMPLE_OUTDIR}/${SAMPLE}_RSeQC_output_all.txt"
    SHORT_OUT="${SAMPLE_OUTDIR}/${SAMPLE}_RSeQC_output.txt"

    infer_experiment.py \
      -r "$EXON_BED" \
      -i "$BAM" \
      > "$FULL_OUT"

    # Same extraction logic as original:
    # extract the "+-" line like the lab script does
    out=$(grep "+-" "$FULL_OUT" | cut -d":" -f2 | tr -d ' ' || true)

    if [[ -z "$out" ]]; then
        echo "WARNING: Could not extract '+-' strandedness value for $SAMPLE"
        out="NA"
    fi

    echo -e "${SAMPLE}\t${out}" > "$SHORT_OUT"

    echo "Done: $SAMPLE"
    echo "Wrote: $FULL_OUT"
    echo "Wrote: $SHORT_OUT"

done

echo "Strandedness QC complete."