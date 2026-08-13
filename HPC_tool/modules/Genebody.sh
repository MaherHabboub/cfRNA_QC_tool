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

# -----------------------------
# Paths
# -----------------------------
BED12_PATH_FILE="${OUTDIR}/annotation/BED12.path.txt"
RESULT_DIR="${OUTDIR}/gene_body_coverage"

# -----------------------------
# Software environment
# -----------------------------
module purge
module load picard/3.0.0-Java-17
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

echo "Running gene body coverage QC..."
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

    PREFIX="${RESULT_DIR}/${SAMPLE}"

    echo "[1/2] Building BAM index..."

    java -jar "$EBROOTPICARD/picard.jar" BuildBamIndex \
      I="$BAM" \
      O="${BAM}.bai"

    test -f "${BAM}.bai"

    echo "[2/2] Running RSeQC geneBody_coverage..."

    set +e

    geneBody_coverage.py \
      -r "$BED12" \
      -i "$BAM" \
      -o "$PREFIX"

    rc=$?

    set -e

    if [[ $rc -ne 0 && ! -f "${PREFIX}.geneBodyCoverage.txt" ]]; then
        echo "WARNING: geneBody_coverage failed for $SAMPLE"
        continue
    fi

    echo "Done: $SAMPLE"

done

echo "Gene body coverage QC complete."
