#!/bin/bash

set -euo pipefail

CONFIG="${1:-}"
TARGET_SAMPLE="${2:-ALL}"

if [[ -z "$CONFIG" ]]; then
    echo "Usage: bash Star_Alignment.sh path/to/config.sh [sample_id]" >&2
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
case "${STAR_ENABLED:-yes}" in
    no) echo "STAR alignment is disabled."; exit 0 ;;
    yes) ;;
    *) echo "ERROR: STAR_ENABLED must be 'yes' or 'no'" >&2; exit 1 ;;
esac
: "${SAMPLESHEET:?ERROR: SAMPLESHEET not set}"
: "${OUTDIR:?ERROR: OUTDIR not set}"
: "${STAR_INDEX:?ERROR: STAR_INDEX not set}"
: "${GTF:?ERROR: GTF not set}"
STAR_THREADS="${SLURM_CPUS_PER_TASK:-8}"
[[ "$STAR_THREADS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: STAR_THREADS must be positive" >&2; exit 1; }
[[ -d "$STAR_INDEX" && -r "$GTF" ]] || { echo "ERROR: STAR_INDEX or GTF is missing" >&2; exit 1; }

# -----------------------------
# Paths
# -----------------------------
STAR_OUTDIR="${STAR_OUTDIR:-${OUTDIR}/star}"

# -----------------------------
# Software environment
# -----------------------------
module purge
module load env/software/doduo
module load STAR/2.7.11b-GCC-13.2.0

command -v STAR >/dev/null || { echo "ERROR: STAR is unavailable" >&2; exit 1; }

# -----------------------------
# Read samples
# -----------------------------
# Header-based reading also supports standalone calls with FASTQ-only sheets.
# Capture first so a parser failure cannot be hidden by process substitution.
ROWS="$(python3 - "$SAMPLESHEET" "$TARGET_SAMPLE" <<'PY'
import csv
from pathlib import Path
import re
import sys

sheet = Path(sys.argv[1]).resolve()
target = sys.argv[2]
selected, seen = [], set()
try:
    with sheet.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if not any(row.values()):
                continue
            sample = (row.get("sample_id") or "").strip()
            if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", sample) or sample == "ALL" or sample in seen:
                raise ValueError(f"Invalid or duplicate sample_id: {sample!r}")
            seen.add(sample)
            if target != "ALL" and sample != target:
                continue
            layout = (row.get("layout") or "").strip()
            condition = (row.get("condition") or "").strip()
            if layout not in ("PE", "SE") or condition in ("", "NA", "."):
                raise ValueError(f"Valid layout and condition are required for {sample}")
            reads = []
            for key in ("fastq_r1", "fastq_r2"):
                value = (row.get(key) or "").strip()
                reads.append("NA" if value in ("", "NA", ".") else str((sheet.parent / value).resolve()))
            selected.append([sample, *reads, layout, condition])
    if not selected:
        raise ValueError(f"No samples found for {target!r}")
    csv.writer(sys.stdout, delimiter="\t", lineterminator="\n").writerows(selected)
except (OSError, ValueError) as error:
    sys.exit(f"ERROR: {error}")
PY
)"

# -----------------------------
# Process samples
# -----------------------------
while IFS=$'\t' read -r SAMPLE R1 R2 LAYOUT CONDITION; do
    [[ -r "$R1" && "$R1" != "NA" ]] || { echo "ERROR: FASTQ not found: $R1" >&2; exit 1; }
    reads=("$R1")
    if [[ "$LAYOUT" == "PE" ]]; then
        [[ -r "$R2" && "$R2" != "NA" ]] || { echo "ERROR: FASTQ not found: $R2" >&2; exit 1; }
        if [[ "$R1" == *.gz && "$R2" != *.gz || "$R1" != *.gz && "$R2" == *.gz ]]; then
            echo "ERROR: STAR mates must use the same compression for $SAMPLE" >&2
            exit 1
        fi
        reads+=("$R2")
    fi
    read_command=()
    [[ "$R1" != *.gz ]] || read_command=(--readFilesCommand zcat)
    sample_dir="${STAR_OUTDIR}/${SAMPLE}"
    BAM="${sample_dir}/${SAMPLE}.Aligned.sortedByCoord.out.bam"
    mkdir -p "$sample_dir"
    echo "Running STAR for $SAMPLE ($LAYOUT, $CONDITION)"
    STAR \
        --runThreadN "$STAR_THREADS" \
        --genomeDir "$STAR_INDEX" \
        --sjdbGTFfile "$GTF" \
        --readFilesIn "${reads[@]}" \
        ${read_command[@]+"${read_command[@]}"} \
        --outFileNamePrefix "${sample_dir}/${SAMPLE}." \
        --outSAMtype BAM SortedByCoordinate \
        --outReadsUnmapped Fastx \
        --twopassMode Basic \
        --outMultimapperOrder Random \
        --outSAMmultNmax -1 \
        --outFilterMultimapNmax 10 \
        --outSAMprimaryFlag AllBestScore \
        --outFilterScoreMinOverLread 0.66 \
        --outFilterMatchNminOverLread 0.66 \
        --outFilterMatchNmin 20
    echo "Done: $SAMPLE; BAM: $BAM"
done <<< "$ROWS"
