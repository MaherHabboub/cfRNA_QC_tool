#!/bin/bash

set -euo pipefail

CONFIG="${1:-}"
TARGET_SAMPLE="${2:-ALL}"

if [[ -z "$CONFIG" ]]; then
    echo "Usage: bash Splice_Junction.sh path/to/config.sh [sample_id]"
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
MAPQ_MIN=30
RESULT_DIR="${OUTDIR}/splice_junctions"

# -----------------------------
# Software environment
# -----------------------------
module purge
module load SAMtools
module load Anaconda3/2024.06-1

for cmd in samtools python; do
    command -v "$cmd" >/dev/null 2>&1 || {
        echo "ERROR: Required command not available after module loading: $cmd" >&2
        exit 1
    }
done

mkdir -p "$RESULT_DIR"

HEADER='sample\tcondition\ttotal_unique_mapped_reads\tspliced_reads\tfraction_spliced\tmapq_min'

echo "Running splice junction QC..."
echo "Sample sheet: $SAMPLESHEET"
echo "Output: $RESULT_DIR"
echo "MAPQ cutoff: $MAPQ_MIN"

N_SAMPLES=0

# Expected samplesheet positions: sample ID (1), BAM (4), STAR log (5),
# SJ.out.tab (6), and condition (8). Optional columns may follow condition.
while IFS= read -r SAMPLE_LINE; do
    IFS=$'\t' read -r -a FIELDS <<< "$SAMPLE_LINE"

    if (( ${#FIELDS[@]} < 8 )); then
        echo "ERROR: Expected at least 8 tab-separated columns, found ${#FIELDS[@]}" >&2
        echo "Problematic row: $SAMPLE_LINE" >&2
        exit 1
    fi

    SAMPLE="${FIELDS[0]}"
    BAM="${FIELDS[3]}"
    STARLOG="${FIELDS[4]}"
    SJTAB="${FIELDS[5]}"
    CONDITION="${FIELDS[7]}"

    [[ -z "${SAMPLE:-}" ]] && continue
    [[ "$TARGET_SAMPLE" == "ALL" || "$SAMPLE" == "$TARGET_SAMPLE" ]] || continue

    SAMPLE="${SAMPLE//$'\r'/}"
    BAM="${BAM//$'\r'/}"
    STARLOG="${STARLOG//$'\r'/}"
    SJTAB="${SJTAB//$'\r'/}"
    CONDITION="${CONDITION//$'\r'/}"

    if [[ -z "$CONDITION" ]]; then
        echo "ERROR: Condition is empty for sample: $SAMPLE" >&2
        exit 1
    fi

    SAMPLE_OUTDIR="${RESULT_DIR}/${SAMPLE}"
    mkdir -p "$SAMPLE_OUTDIR"
    rm -f "${SAMPLE_OUTDIR}/.complete" \
        "${SAMPLE_OUTDIR}/${SAMPLE}.splice_junction_summary.tsv" \
        "${SAMPLE_OUTDIR}/${SAMPLE}.Log.final.out"

    echo "------------------------------------"
    echo "Processing: $SAMPLE"

    # Preserve the existing STAR junction-table summary.
    if [[ -f "$SJTAB" ]]; then
        JUNCTION_TSV="${SAMPLE_OUTDIR}/${SAMPLE}.splice_junction_summary.tsv"

        awk -v sample="$SAMPLE" 'BEGIN {
            total=0; annotated=0; novel=0; uniq=0; multi=0
        }
        {
            total++
            if ($6 == 1) annotated++; else novel++
            uniq += $7
            multi += $8
        }
        END {
            fraction_annotated=(total ? annotated/total : 0)
            fraction_novel=(total ? novel/total : 0)
            print "sample\ttotal_junctions\tannotated_junctions\tnovel_junctions\tfraction_annotated\tfraction_novel\tsum_unique_support\tsum_multi_support"
            printf "%s\t%d\t%d\t%d\t%.6f\t%.6f\t%d\t%d\n", sample, total, annotated, novel, fraction_annotated, fraction_novel, uniq, multi
        }' "$SJTAB" > "$JUNCTION_TSV"

        echo "Junction summary: $JUNCTION_TSV"
    else
        echo "WARNING: SJ.out.tab not found; junction summary skipped"
    fi

    if [[ -f "$STARLOG" ]]; then
        cp -f "$STARLOG" "${SAMPLE_OUTDIR}/${SAMPLE}.Log.final.out"
    else
        echo "WARNING: STAR Log.final.out not found; log copy skipped"
    fi

    # A CIGAR containing N crosses one or more splice junctions. Restrict to
    # primary, mapped, non-duplicate, QC-passing alignments with MAPQ >= 30.
    if [[ ! -f "$BAM" ]]; then
        echo "WARNING: BAM not found; read-fraction calculation skipped"
        continue
    fi

    samtools quickcheck -v "$BAM" || {
        echo "ERROR: BAM failed samtools quickcheck: $BAM" >&2
        exit 1
    }

    if ! samtools view -H "$BAM" | awk '$1 == "@SQ" {found=1} END {exit !found}'; then
        echo "ERROR: BAM header contains no @SQ reference-sequence records for $SAMPLE: $BAM" >&2
        exit 1
    fi

    COUNTS="$({
        samtools view -F 3844 -q "$MAPQ_MIN" "$BAM" |
        awk 'BEGIN {total=0; nonsplice=0; splice=0}
             {total++; if ($6 ~ /[0-9]+N/) splice++; else nonsplice++}
             END {printf "%d\t%d\t%d", total, nonsplice, splice}'
    })"

    IFS=$'\t' read -r TOTAL_UNIQUE NONSPLICE_READS SPLICE_READS <<< "$COUNTS"

    if (( TOTAL_UNIQUE == 0 )); then
        FRACTION_SPLICED="0.000000"
        echo "WARNING: No reads passed the primary/unique filters"
    else
        FRACTION_SPLICED="$(
            awk -v nonsplice="$NONSPLICE_READS" -v total="$TOTAL_UNIQUE" \
                'BEGIN {printf "%.6f", 1 - (nonsplice / total)}'
        )"
    fi

    SAMPLE_TSV="${SAMPLE_OUTDIR}/${SAMPLE}.splice_read_fraction.tsv"
    {
        printf '%b\n' "$HEADER"
        printf '%s\t%s\t%d\t%d\t%s\t%d\n' \
            "$SAMPLE" "$CONDITION" "$TOTAL_UNIQUE" "$SPLICE_READS" \
            "$FRACTION_SPLICED" "$MAPQ_MIN"
    } > "$SAMPLE_TSV"

    printf '%s\n' "${HPC_RUN_ID:-manual}" > "${SAMPLE_OUTDIR}/.complete.tmp.$"
    mv -f "${SAMPLE_OUTDIR}/.complete.tmp.$" "${SAMPLE_OUTDIR}/.complete"
    ((N_SAMPLES+=1))

    echo "Unique mapped reads: $TOTAL_UNIQUE"
    echo "Condition: $CONDITION"
    echo "Spliced reads: $SPLICE_READS ($FRACTION_SPLICED)"
done < <(tail -n +2 "$SAMPLESHEET")

(( N_SAMPLES > 0 )) || {
    echo "ERROR: No sample produced a spliced read fraction" >&2
    exit 1
}

echo "Splice junction QC complete for $N_SAMPLES sample(s). Cohort summaries are created during reporting."
