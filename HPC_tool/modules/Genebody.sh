#!/bin/bash

set -euo pipefail

CONFIG="${1:-}"
TARGET_SAMPLE="${2:-ALL}"

if [[ -z "$CONFIG" ]]; then
    echo "Usage: bash Genebody.sh path/to/config.sh [sample_id]"
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
module load picard/3.0.0-Java-17
module load RSeQC/5.0.1-foss-2023a

RESULT_DIR="${OUTDIR}/gene_body_coverage"
mkdir -p "$RESULT_DIR"

echo "Running gene body coverage QC..."
echo "Target sample: $TARGET_SAMPLE"

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

    PREFIX="${RESULT_DIR}/${SAMPLE}"

    echo "[1/2] Building BAM index..."

    # Same logic as original
    java -jar "$EBROOTPICARD/picard.jar" BuildBamIndex \
      I="$BAM" \
      O="${BAM}.bai"

    # Same check as original
    test -f "${BAM}.bai"

    echo "[2/2] Running RSeQC geneBody_coverage..."

    # Same tolerant behavior as original
    set +e

    geneBody_coverage.py \
      -r "$BED12" \
      -i "$BAM" \
      -o "$PREFIX"

    rc=$?

    set -e

    # Same logic as original
    if [[ $rc -ne 0 && ! -f "${PREFIX}.geneBodyCoverage.txt" ]]; then
        echo "WARNING: geneBody_coverage failed for $SAMPLE"
        continue
    fi

    echo "Done: $SAMPLE"

done

echo "Gene body coverage QC complete."