#!/bin/bash

set -euo pipefail

CONFIG="${1:-}"
TARGET_SAMPLE="${2:-ALL}"

if [[ -z "$CONFIG" ]]; then
    echo "Usage: bash HTSeq_Counts.sh path/to/config.sh [sample_id]" >&2
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
case "${HTSEQ_ENABLED:-yes}" in
    no) echo "HTSeq counting is disabled."; exit 0 ;;
    yes) ;;
    *) echo "ERROR: HTSEQ_ENABLED must be 'yes' or 'no'" >&2; exit 1 ;;
esac
: "${SAMPLESHEET:?ERROR: SAMPLESHEET not set}"
: "${OUTDIR:?ERROR: OUTDIR not set}"
: "${GTF:?ERROR: GTF not set}"
HTSEQ_STRANDED="${HTSEQ_STRANDED:-yes}"
HTSEQ_THREADS="${SLURM_CPUS_PER_TASK:-1}"
HPC_RUN_ID="${HPC_RUN_ID:-manual}"
case "$HTSEQ_STRANDED" in
    yes|no|reverse) ;;
    *) echo "ERROR: HTSEQ_STRANDED must be yes, no, or reverse" >&2; exit 1 ;;
esac
[[ "$HTSEQ_THREADS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: HTSEQ_THREADS must be positive" >&2; exit 1; }
[[ -r "$GTF" ]] || { echo "ERROR: GTF not found: $GTF" >&2; exit 1; }

# -----------------------------
# Paths
# -----------------------------
HTSEQ_OUTDIR="${HTSEQ_OUTDIR:-${OUTDIR}/htseq}"

# -----------------------------
# Software environment
# -----------------------------
module purge
module load env/software/doduo
module load HTSeq/2.0.7-foss-2023a
# Match HTSeq's foss/2023a compiler toolchain.
module load SAMtools/1.18-GCC-12.3.0

command -v htseq-count >/dev/null || { echo "ERROR: htseq-count is unavailable" >&2; exit 1; }
command -v samtools >/dev/null || { echo "ERROR: samtools is unavailable" >&2; exit 1; }

# -----------------------------
# Read samples
# -----------------------------
# Read the named columns so standalone and resolved samplesheets both work.
ROWS="$(python3 - "$SAMPLESHEET" "$TARGET_SAMPLE" "${STAR_ENABLED:-yes}" \
    "${STAR_OUTDIR:-${OUTDIR}/star}" <<'PY'
import csv
from pathlib import Path
import re
import sys

sheet = Path(sys.argv[1]).resolve()
target, star_enabled, star_outdir = sys.argv[2:]
selected, seen = [], set()
try:
    if star_enabled not in ("yes", "no"):
        raise ValueError("STAR_ENABLED must be 'yes' or 'no'")
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
            if star_enabled == "yes":
                bam = Path(star_outdir).resolve() / sample / f"{sample}.Aligned.sortedByCoord.out.bam"
            else:
                value = (row.get("bam") or "").strip()
                if value in ("", "NA", "."):
                    raise ValueError(f"bam is required for {sample}")
                bam = (sheet.parent / value).resolve()
            selected.append([sample, str(bam)])
    if not selected:
        raise ValueError(f"No samples found for {target!r}")
    csv.writer(sys.stdout, delimiter="\t", lineterminator="\n").writerows(selected)
except (OSError, ValueError) as error:
    sys.exit(f"ERROR: {error}")
PY
)"

# -----------------------------
# Helper functions
# -----------------------------
# Each submitted job receives exactly one sample. Manual ALL mode still counts
# separately; combining samples belongs exclusively to Aggregate.sh.
process_sample() (
    local sample="$1" bam="$2"
    local sample_dir="${HTSEQ_OUTDIR}/${sample}"
    mkdir -p "$sample_dir"
    local outfile="${sample_dir}/${sample}_htseq_counts.txt"
    local marker="${outfile}.complete"
    rm -f "$marker"
    local workdir
    workdir="$(mktemp -d "${sample_dir}/.htseq.XXXXXX")"
    trap 'rm -rf "$workdir"' EXIT

    [[ -f "$bam" ]] || { echo "ERROR: BAM not found for $sample: $bam" >&2; exit 1; }
    samtools quickcheck -v "$bam"
    samtools view -H "$bam" > "${workdir}/header.sam"
    awk '$1 == "@SQ" {found=1} END {exit !found}' "${workdir}/header.sam" || {
        echo "ERROR: BAM has no @SQ records for $sample" >&2; exit 1;
    }
    if ! awk '$1 == "@HD" {for (i=2; i<=NF; i++) if ($i == "SO:coordinate") found=1} END {exit !found}' "${workdir}/header.sam"; then
        samtools sort -@ "$HTSEQ_THREADS" -o "${workdir}/sorted.bam" "$bam"
        bam="${workdir}/sorted.bam"
        samtools quickcheck -v "$bam"
    fi

    echo "Running HTSeq for $sample; full BAM: $2; strandedness: $HTSEQ_STRANDED"
    htseq-count \
        --format bam \
        --order pos \
        --nonunique none \
        --stranded "$HTSEQ_STRANDED" \
        "$bam" "$GTF" > "${workdir}/counts.txt"
    [[ -s "${workdir}/counts.txt" ]] || { echo "ERROR: Empty HTSeq counts for $sample" >&2; exit 1; }
    mv -f "${workdir}/counts.txt" "$outfile"
    printf '%s\n' "$HPC_RUN_ID" > "${workdir}/complete"
    mv -f "${workdir}/complete" "$marker"
    echo "Wrote: $outfile"
)

# -----------------------------
# Process samples
# -----------------------------
while IFS=$'\t' read -r SAMPLE BAM; do
    process_sample "$SAMPLE" "$BAM"
done <<< "$ROWS"
