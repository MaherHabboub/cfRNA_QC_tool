#!/bin/bash

set -euo pipefail

CONFIG="${1:-}"

if [[ -z "$CONFIG" || ! -f "$CONFIG" ]]; then
    echo "Usage: bash Downsample.sh path/to/config_downsampling.sh" >&2
    exit 1
fi

source "$CONFIG"

: "${SOURCE_SAMPLESHEET:?ERROR: SOURCE_SAMPLESHEET not set in config}"
: "${SAMPLESHEET:?ERROR: SAMPLESHEET not set in config}"
: "${OUTDIR:?ERROR: OUTDIR not set in config}"
: "${SOURCE_BED12_PATH_FILE:?ERROR: SOURCE_BED12_PATH_FILE not set in config}"
: "${DOWNSAMPLE_SEED:?ERROR: DOWNSAMPLE_SEED not set in config}"

DOWNSAMPLE_THREADS="${DOWNSAMPLE_THREADS:-4}"
EXPECTED_SOURCE_SAMPLES="${EXPECTED_SOURCE_SAMPLES:-2}"

if (( ${#DOWNSAMPLE_LABELS[@]} != ${#DOWNSAMPLE_FRACTIONS[@]} )); then
    echo "ERROR: DOWNSAMPLE_LABELS and DOWNSAMPLE_FRACTIONS differ in length." >&2
    exit 1
fi

if [[ ! -s "$SOURCE_SAMPLESHEET" ]]; then
    echo "ERROR: Source samplesheet missing or empty: $SOURCE_SAMPLESHEET" >&2
    exit 1
fi

source_count=$(tail -n +2 "$SOURCE_SAMPLESHEET" | awk -F '\t' 'NF && $1 != "" {n++} END {print n+0}')
if [[ "$source_count" -ne "$EXPECTED_SOURCE_SAMPLES" ]]; then
    echo "ERROR: Expected $EXPECTED_SOURCE_SAMPLES source samples, found $source_count." >&2
    exit 1
fi

module purge
module load "${SAMTOOLS_MODULE:-SAMtools}"

BAM_DIR="${OUTDIR}/downsampled_bams"
INPUT_DIR="$(dirname "$SAMPLESHEET")"
BENCHMARK_DIR="${OUTDIR}/benchmark"
ANNOTATION_DIR="${OUTDIR}/annotation"
MANIFEST="${BENCHMARK_DIR}/downsampling_manifest.tsv"

mkdir -p "$BAM_DIR" "$INPUT_DIR" "$BENCHMARK_DIR" "$ANNOTATION_DIR"

if [[ ! -s "$SOURCE_BED12_PATH_FILE" ]]; then
    echo "ERROR: BED12 path file missing or empty: $SOURCE_BED12_PATH_FILE" >&2
    exit 1
fi

BED12="$(<"$SOURCE_BED12_PATH_FILE")"
if [[ -z "$BED12" || ! -f "$BED12" ]]; then
    echo "ERROR: BED12 referenced by $SOURCE_BED12_PATH_FILE does not exist: ${BED12:-<empty>}" >&2
    exit 1
fi
printf '%s\n' "$BED12" > "${ANNOTATION_DIR}/BED12.path.txt"

source_header=$(head -n 1 "$SOURCE_SAMPLESHEET")
printf '%s\n' "$source_header" > "$SAMPLESHEET"
printf 'source_sample\ttest_sample\tlabel\trequested_fraction\tsource_bam\ttest_bam\toriginal_alignments\tretained_alignments\tobserved_fraction\tdownsampling_seconds\tbam_bytes\n' > "$MANIFEST"

echo "Creating benchmark BAMs with samtools seed $DOWNSAMPLE_SEED"

while IFS=$'\t' read -r SAMPLE FASTQ1 FASTQ2 BAM STARLOG SJTAB LAYOUT CONDITION
do
    [[ -z "${SAMPLE:-}" ]] && continue

    if [[ ! -f "$BAM" ]]; then
        echo "ERROR: BAM not found for $SAMPLE: $BAM" >&2
        exit 1
    fi

    samtools quickcheck -v "$BAM"
    original_count=$(samtools view -@ "$DOWNSAMPLE_THREADS" -c "$BAM")

    for i in "${!DOWNSAMPLE_LABELS[@]}"
    do
        label="${DOWNSAMPLE_LABELS[$i]}"
        fraction="${DOWNSAMPLE_FRACTIONS[$i]}"
        test_sample="${SAMPLE}__${label}"
        test_bam="${BAM_DIR}/${test_sample}.bam"
        start_epoch=$(date +%s)

        echo "[$SAMPLE] $label ($fraction)"

        if [[ "$fraction" == "1" || "$fraction" == "1.0" || "$fraction" == "1.00" ]]; then
            if [[ -e "$test_bam" && ! -L "$test_bam" ]]; then
                echo "ERROR: Refusing to replace non-symlink ds100 BAM: $test_bam" >&2
                exit 1
            fi
            ln -sfn "$BAM" "$test_bam"
        else
            fraction_digits="${fraction#0.}"
            tmp_bam="${test_bam}.tmp.${SLURM_JOB_ID:-$$}"
            rm -f "$tmp_bam"
            samtools view \
                -@ "$DOWNSAMPLE_THREADS" \
                -b \
                -s "${DOWNSAMPLE_SEED}.${fraction_digits}" \
                -o "$tmp_bam" \
                "$BAM"
            samtools quickcheck -v "$tmp_bam"
            mv -f "$tmp_bam" "$test_bam"
        fi

        retained_count=$(samtools view -@ "$DOWNSAMPLE_THREADS" -c "$test_bam")
        end_epoch=$(date +%s)
        elapsed=$((end_epoch - start_epoch))
        bam_bytes=$(stat -Lc '%s' "$test_bam")
        observed_fraction=$(awk -v n="$retained_count" -v d="$original_count" 'BEGIN {if (d>0) printf "%.8f", n/d; else print "NA"}')

        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$SAMPLE" "$test_sample" "$label" "$fraction" "$BAM" "$test_bam" \
            "$original_count" "$retained_count" "$observed_fraction" "$elapsed" "$bam_bytes" \
            >> "$MANIFEST"

        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$test_sample" "$FASTQ1" "$FASTQ2" "$test_bam" "$STARLOG" "$SJTAB" "$LAYOUT" "$CONDITION" \
            >> "$SAMPLESHEET"
    done
done < <(tail -n +2 "$SOURCE_SAMPLESHEET")

echo "Downsampling complete."
echo "Generated samplesheet: $SAMPLESHEET"
echo "Manifest: $MANIFEST"
