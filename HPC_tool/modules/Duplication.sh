#!/bin/bash

set -euo pipefail

CONFIG="${1:-}"

if [[ -z "$CONFIG" ]]; then
    echo "Usage: bash 03_duplication_qc.sh path/to/config.sh"
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

module purge
module load picard/3.0.0-Java-17

RESULT_DIR="${OUTDIR}/duplication"
mkdir -p "$RESULT_DIR"

echo "Running duplication QC..."

# Skip header
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

    METRICS="${SAMPLE_OUTDIR}/${SAMPLE}.markdup.metrics.txt"
    TMP_BAM="${SAMPLE_OUTDIR}/${SAMPLE}.markdup.tmp.bam"
    SUMMARY="${SAMPLE_OUTDIR}/${SAMPLE}.duplication_summary.tsv"

    echo "Running Picard MarkDuplicates..."

    # Same logic as original
    java -jar "$EBROOTPICARD/picard.jar" MarkDuplicates \
      I="$BAM" \
      O="$TMP_BAM" \
      M="$METRICS" \
      ASSUME_SORTED=true \
      VALIDATION_STRINGENCY=SILENT \
      REMOVE_DUPLICATES=false \
      CREATE_INDEX=false

    # Same cleanup logic
    rm -f "$TMP_BAM"

    # Same metric extraction logic
    PCT=$(awk '
        BEGIN{FS="\t"}
        $1=="LIBRARY"{hdr=1; next}
        hdr && $0!=""{print $9; exit}
    ' "$METRICS")

    echo -e "sample\tpercent_duplication" > "$SUMMARY"
    echo -e "${SAMPLE}\t${PCT}" >> "$SUMMARY"

    echo "Done: $SAMPLE"
    echo "Duplication: $PCT"

done

echo "Duplication QC complete."