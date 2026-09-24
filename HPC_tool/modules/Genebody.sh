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
BED12_PATH_FILE="${OUTDIR}/annotation/BED12.path.txt"
RESULT_DIR="${OUTDIR}/gene_body_coverage"
DOWNSAMPLE_MANIFEST="${OUTDIR}/downsampled_bams/downsampling_manifest.tsv"

# -----------------------------
# Software environment
# -----------------------------
module purge
module load picard/3.0.0-Java-17
module load SAMtools
module load RSeQC/5.0.1-foss-2023a

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

echo "Running gene body coverage QC..."
echo "Target sample: $TARGET_SAMPLE"
echo "Generated BED12: $BED12"

tail -n +2 "$SAMPLESHEET" | while IFS=$'\t' read -r SAMPLE FASTQ1 FASTQ2 BAM STARLOG SJTAB LAYOUT CONDITION TRANSCRIPTOME_BAM
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

    samtools quickcheck -v "$QC_BAM" || {
        echo "ERROR: BAM failed samtools quickcheck for $SAMPLE: $QC_BAM" >&2
        exit 1
    }

    if ! samtools view -H "$QC_BAM" | awk '$1 == "@SQ" {found=1} END {exit !found}'; then
        echo "ERROR: BAM header contains no @SQ reference-sequence records for $SAMPLE: $QC_BAM" >&2
        exit 1
    fi

    PREFIX="${RESULT_DIR}/${SAMPLE}"
    SORTED_INPUT=""
    INPUT_BAM="$QC_BAM"
    sort_order="$(samtools view -H "$QC_BAM" | awk '$1 == "@HD" {for (i = 1; i <= NF; i++) if ($i ~ /^SO:/) {print substr($i, 4); exit}}')"

    if [[ "$sort_order" != "coordinate" ]]; then
        SORTED_INPUT="${RESULT_DIR}/${SAMPLE}.coordinate_sorted.tmp.bam"
        echo "Input BAM is not declared coordinate-sorted (SO=${sort_order:-unspecified}); creating a temporary coordinate-sorted BAM."
        java -jar "$EBROOTPICARD/picard.jar" SortSam \
          I="$QC_BAM" \
          O="$SORTED_INPUT" \
          SORT_ORDER=coordinate \
          VALIDATION_STRINGENCY=SILENT
        INPUT_BAM="$SORTED_INPUT"
    fi

    echo "[1/2] Building BAM index..."
    echo "Input BAM: $INPUT_BAM"

    java -jar "$EBROOTPICARD/picard.jar" BuildBamIndex \
      I="$INPUT_BAM" \
      O="${INPUT_BAM}.bai"

    test -f "${INPUT_BAM}.bai"

    echo "[2/2] Running RSeQC geneBody_coverage..."

    set +e

    geneBody_coverage.py \
      -r "$BED12" \
      -i "$INPUT_BAM" \
      -o "$PREFIX"

    rc=$?

    set -e

    if [[ $rc -ne 0 && ! -f "${PREFIX}.geneBodyCoverage.txt" ]]; then
        echo "WARNING: geneBody_coverage failed for $SAMPLE"
        [[ -z "$SORTED_INPUT" ]] || rm -f "$SORTED_INPUT" "${SORTED_INPUT}.bai"
        continue
    fi

    [[ -z "$SORTED_INPUT" ]] || rm -f "$SORTED_INPUT" "${SORTED_INPUT}.bai"

    echo "Done: $SAMPLE"

done

echo "Gene body coverage QC complete."
