#!/bin/bash

set -euo pipefail

CONFIG="${1:-}"

if [[ -z "$CONFIG" ]]; then
    echo "Usage: bash Downsample.sh path/to/config.sh" >&2
    exit 1
fi

if [[ ! -f "$CONFIG" ]]; then
    echo "ERROR: Config file not found: $CONFIG" >&2
    exit 1
fi

source "$CONFIG"

# -----------------------------
# Required config variables
# -----------------------------
: "${SAMPLESHEET:?ERROR: SAMPLESHEET not set in config}"
: "${OUTDIR:?ERROR: OUTDIR not set in config}"

DOWNSAMPLE_ENABLED="${DOWNSAMPLE_ENABLED:-yes}"
DOWNSAMPLE_TARGET_ALIGNMENTS="${DOWNSAMPLE_TARGET_ALIGNMENTS:-1000000}"
DOWNSAMPLE_SEED="${DOWNSAMPLE_SEED:-42}"
DOWNSAMPLE_THREADS="${DOWNSAMPLE_THREADS:-4}"

case "$DOWNSAMPLE_ENABLED" in
    yes|no) ;;
    *)
        echo "ERROR: DOWNSAMPLE_ENABLED must be 'yes' or 'no': $DOWNSAMPLE_ENABLED" >&2
        exit 1
        ;;
esac

for value_name in DOWNSAMPLE_TARGET_ALIGNMENTS DOWNSAMPLE_SEED DOWNSAMPLE_THREADS; do
    value="${!value_name}"
    if [[ ! "$value" =~ ^[0-9]+$ ]]; then
        echo "ERROR: $value_name must be a non-negative integer: $value" >&2
        exit 1
    fi
done

if (( DOWNSAMPLE_TARGET_ALIGNMENTS < 1 || DOWNSAMPLE_THREADS < 1 )); then
    echo "ERROR: DOWNSAMPLE_TARGET_ALIGNMENTS and DOWNSAMPLE_THREADS must be at least 1." >&2
    exit 1
fi

if [[ "$DOWNSAMPLE_ENABLED" != "yes" ]]; then
    echo "Downsampling is disabled by DOWNSAMPLE_ENABLED=no; no BAMs were created."
    exit 0
fi

# -----------------------------
# Paths
# -----------------------------
BAM_DIR="${OUTDIR}/downsampled_bams"
MANIFEST="${BAM_DIR}/downsampling_manifest.tsv"
MANIFEST_TMP="${MANIFEST}.tmp.$$"

# -----------------------------
# Software environment
# -----------------------------
module purge
module load SAMtools

command -v samtools >/dev/null 2>&1 || {
    echo "ERROR: samtools is unavailable after loading SAMtools." >&2
    exit 1
}

if [[ ! -s "$SAMPLESHEET" ]]; then
    echo "ERROR: Samplesheet is missing or empty: $SAMPLESHEET" >&2
    exit 1
fi

duplicate_samples="$(
    tail -n +2 "$SAMPLESHEET" | awk -F '\t' 'NF && $1 != "" {count[$1]++} END {for (sample in count) if (count[sample] > 1) print sample}'
)"
if [[ -n "$duplicate_samples" ]]; then
    echo "ERROR: Samplesheet contains duplicate sample IDs:" >&2
    echo "$duplicate_samples" >&2
    exit 1
fi

mkdir -p "$BAM_DIR"
trap 'rm -f "$MANIFEST_TMP"' EXIT

HEADER='sample_id\toriginal_bam\tselected_bam\toriginal_alignments\tretained_alignments\trequested_fraction\tobserved_fraction\tseed\tstatus'
printf '%b\n' "$HEADER" > "$MANIFEST_TMP"

echo "Running BAM downsampling..."
echo "Target alignments: $DOWNSAMPLE_TARGET_ALIGNMENTS"
echo "Seed: $DOWNSAMPLE_SEED"
echo "Threads: $DOWNSAMPLE_THREADS"

while IFS=$'\t' read -r SAMPLE FASTQ1 FASTQ2 BAM STARLOG SJTAB LAYOUT CONDITION
do
    [[ -z "${SAMPLE:-}" ]] && continue

    if [[ ! -f "$BAM" ]]; then
        echo "ERROR: BAM not found for $SAMPLE: $BAM" >&2
        exit 1
    fi

    samtools quickcheck -v "$BAM"
    original_count="$(samtools view -@ "$DOWNSAMPLE_THREADS" -c "$BAM")"

    selected_bam="$BAM"
    retained_count="$original_count"
    requested_fraction="1.000000000000"
    observed_fraction="1.000000000000"
    status="not_downsampled"

    if (( original_count > DOWNSAMPLE_TARGET_ALIGNMENTS )); then
        requested_fraction="$(
            awk -v target="$DOWNSAMPLE_TARGET_ALIGNMENTS" -v original="$original_count" \
                'BEGIN {printf "%.12f", target / original}'
        )"
        fraction_digits="${requested_fraction#0.}"
        selected_bam="${BAM_DIR}/${SAMPLE}.downsampled.bam"
        tmp_bam="${selected_bam}.tmp.${SLURM_JOB_ID:-$$}"

        rm -f "$tmp_bam"
        samtools view \
            -@ "$DOWNSAMPLE_THREADS" \
            -b \
            -s "${DOWNSAMPLE_SEED}.${fraction_digits}" \
            -o "$tmp_bam" \
            "$BAM"
        samtools quickcheck -v "$tmp_bam"

        retained_count="$(samtools view -@ "$DOWNSAMPLE_THREADS" -c "$tmp_bam")"
        observed_fraction="$(
            awk -v retained="$retained_count" -v original="$original_count" \
                'BEGIN {printf "%.12f", retained / original}'
        )"

        mv -f "$tmp_bam" "$selected_bam"
        rm -f "${selected_bam}.bai" "${selected_bam%.bam}.bai"
        status="downsampled"
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$SAMPLE" "$BAM" "$selected_bam" "$original_count" "$retained_count" \
        "$requested_fraction" "$observed_fraction" "$DOWNSAMPLE_SEED" "$status" \
        >> "$MANIFEST_TMP"

    echo "$SAMPLE: $status ($retained_count of $original_count alignments)"
done < <(tail -n +2 "$SAMPLESHEET")

mv -f "$MANIFEST_TMP" "$MANIFEST"

echo "BAM downsampling complete."
echo "Manifest: $MANIFEST"
