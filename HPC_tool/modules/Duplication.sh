#!/bin/bash

set -euo pipefail

CONFIG="${1:-}"
TARGET_SAMPLE="${2:-ALL}"

if [[ -z "$CONFIG" ]]; then
    echo "Usage: bash Duplication.sh path/to/config.sh [sample_id]"
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

DOWNSAMPLE_ENABLED="${DOWNSAMPLE_ENABLED:-yes}"

case "$DOWNSAMPLE_ENABLED" in
    yes|no) ;;
    *)
        echo "ERROR: DOWNSAMPLE_ENABLED must be 'yes' or 'no': $DOWNSAMPLE_ENABLED" >&2
        exit 1
        ;;
esac

# -----------------------------
# Paths
# -----------------------------
RESULT_DIR="${OUTDIR}/duplication"
DOWNSAMPLE_MANIFEST="${OUTDIR}/downsampled_bams/downsampling_manifest.tsv"

# -----------------------------
# Software environment
# -----------------------------
module purge
module load picard/3.0.0-Java-17

# -----------------------------
# Validate downsampling prerequisite
# -----------------------------
if [[ "$DOWNSAMPLE_ENABLED" == "yes" ]]; then
    if [[ ! -s "$DOWNSAMPLE_MANIFEST" ]]; then
        echo "ERROR: Downsampling is enabled, but its manifest is missing or empty: $DOWNSAMPLE_MANIFEST" >&2
        echo "Run Downsample.sh successfully before this module." >&2
        exit 1
    fi

    expected_header=$'sample_id\toriginal_bam\tselected_bam\toriginal_alignments\tretained_alignments\trequested_fraction\tobserved_fraction\tseed\tstatus'
    if [[ "$(head -n 1 "$DOWNSAMPLE_MANIFEST")" != "$expected_header" ]]; then
        echo "ERROR: Downsampling manifest has an unexpected header: $DOWNSAMPLE_MANIFEST" >&2
        exit 1
    fi
fi

mkdir -p "$RESULT_DIR"

echo "Running duplication QC..."
echo "Target sample: $TARGET_SAMPLE"

tail -n +2 "$SAMPLESHEET" | while IFS=$'\t' read -r SAMPLE FASTQ1 FASTQ2 BAM STARLOG SJTAB LAYOUT CONDITION
do
    if [[ "$TARGET_SAMPLE" != "ALL" && "$SAMPLE" != "$TARGET_SAMPLE" ]]; then
        continue
    fi


    echo "------------------------------------"
    echo "Processing: $SAMPLE"

    QC_BAM="$BAM"

    if [[ "$DOWNSAMPLE_ENABLED" == "yes" ]]; then
        manifest_match_count="$(awk -F '\t' -v sample="$SAMPLE" 'NR > 1 && $1 == sample {count++} END {print count+0}' "$DOWNSAMPLE_MANIFEST")"
        if [[ "$manifest_match_count" != "1" ]]; then
            echo "ERROR: Expected one downsampling manifest record for $SAMPLE; found $manifest_match_count." >&2
            exit 1
        fi

        QC_BAM="$(awk -F '\t' -v sample="$SAMPLE" 'NR > 1 && $1 == sample {print $3; exit}' "$DOWNSAMPLE_MANIFEST")"
        if [[ -z "$QC_BAM" || ! -f "$QC_BAM" ]]; then
            echo "ERROR: Selected BAM is missing for $SAMPLE: ${QC_BAM:-<empty path>}" >&2
            exit 1
        fi
    elif [[ ! -f "$QC_BAM" ]]; then
        echo "WARNING: BAM not found for $SAMPLE, skipping"
        continue
    fi

    SAMPLE_OUTDIR="${RESULT_DIR}/${SAMPLE}"
    mkdir -p "$SAMPLE_OUTDIR"

    METRICS="${SAMPLE_OUTDIR}/${SAMPLE}.markdup.metrics.txt"
    TMP_BAM="${SAMPLE_OUTDIR}/${SAMPLE}.markdup.tmp.bam"
    SUMMARY="${SAMPLE_OUTDIR}/${SAMPLE}.duplication_summary.tsv"

    echo "Running Picard MarkDuplicates..."
    echo "Input BAM: $QC_BAM"

    java -jar "$EBROOTPICARD/picard.jar" MarkDuplicates \
      I="$QC_BAM" \
      O="$TMP_BAM" \
      M="$METRICS" \
      ASSUME_SORTED=true \
      VALIDATION_STRINGENCY=SILENT \
      REMOVE_DUPLICATES=false \
      CREATE_INDEX=false

    rm -f "$TMP_BAM"

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
