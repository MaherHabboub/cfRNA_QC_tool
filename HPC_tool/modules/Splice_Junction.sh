#!/bin/bash

set -euo pipefail

CONFIG="${1:-}"

if [[ -z "$CONFIG" ]]; then
    echo "Usage: bash Splice_Junction.sh path/to/config.sh"
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
COMBINED_TSV="${RESULT_DIR}/splice_read_fractions.tsv"
CONDITION_SUMMARY_TSV="${RESULT_DIR}/splice_read_fraction_cohort_summary.tsv"
PLOT_PNG="${RESULT_DIR}/splice_read_fractions.png"
PLOT_PDF="${RESULT_DIR}/splice_read_fractions.pdf"

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
TMP_TSV="${COMBINED_TSV}.tmp.$$"
trap 'rm -f "$TMP_TSV"' EXIT

HEADER='sample\tcondition\ttotal_unique_mapped_reads\tspliced_reads\tfraction_spliced\tmapq_min'
printf '%b\n' "$HEADER" > "$TMP_TSV"

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

    tail -n 1 "$SAMPLE_TSV" >> "$TMP_TSV"
    ((N_SAMPLES+=1))

    echo "Unique mapped reads: $TOTAL_UNIQUE"
    echo "Condition: $CONDITION"
    echo "Spliced reads: $SPLICE_READS ($FRACTION_SPLICED)"
done < <(tail -n +2 "$SAMPLESHEET")

(( N_SAMPLES > 0 )) || {
    echo "ERROR: No sample produced a spliced read fraction" >&2
    exit 1
}

mv -f "$TMP_TSV" "$COMBINED_TSV"

# Plot one spliced-read fraction distribution per condition. Individual
# points are samples; the annotation above each box is the unweighted mean.
python - "$COMBINED_TSV" "$CONDITION_SUMMARY_TSV" "$PLOT_PNG" "$PLOT_PDF" <<'PYTHON'
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

table, summary_file, png_file, pdf_file = sys.argv[1:]
df = pd.read_csv(table, sep="\t")

if df.empty:
    raise SystemExit("Combined splice-read fraction table is empty")

valid = df["total_unique_mapped_reads"] > 0
if not valid.all():
    excluded = ", ".join(df.loc[~valid, "sample"].astype(str))
    print(f"WARNING: Excluding samples with zero qualifying reads: {excluded}", file=sys.stderr)
    df = df.loc[valid].copy()

if df.empty:
    raise SystemExit("No samples with qualifying reads are available for plotting")

if df["condition"].isna().any() or df["condition"].astype(str).str.strip().eq("").any():
    raise SystemExit("At least one sample has an empty condition")

df["condition"] = df["condition"].astype(str)
conditions = list(dict.fromkeys(df["condition"]))
groups = [
    df.loc[df["condition"] == condition, "fraction_spliced"].to_numpy(dtype=float)
    for condition in conditions
]
means = np.array([values.mean() for values in groups])

summary = (
    df.groupby("condition", sort=False)["fraction_spliced"]
      .agg(n_samples="size", mean_fraction="mean", median_fraction="median",
           standard_deviation="std", minimum_fraction="min", maximum_fraction="max")
      .reset_index()
)
summary.to_csv(summary_file, sep="\t", index=False, float_format="%.6f")

rng = np.random.default_rng(42)
x_positions = np.arange(1, len(conditions) + 1)
cmap = plt.get_cmap("tab10")
colors = [cmap(i % 10) for i in range(len(conditions))]

fig_width = max(6.6, 1.35 * len(conditions) + 2.0)
fig, ax = plt.subplots(figsize=(fig_width, 5.7))
ax.set_facecolor("#EBEBEB")

boxplot = ax.boxplot(
    groups,
    positions=x_positions,
    widths=0.55,
    patch_artist=True,
    showfliers=False,
    medianprops={"color": "#000000", "linewidth": 1.5},
    whiskerprops={"color": "#4D4D4D", "linewidth": 1.1},
    capprops={"color": "#4D4D4D", "linewidth": 1.1},
    boxprops={"edgecolor": "#4D4D4D", "linewidth": 1.1},
)

for box, color in zip(boxplot["boxes"], colors):
    box.set_facecolor(color)
    box.set_alpha(0.55)

for x, values, color in zip(x_positions, groups, colors):
    jitter = rng.uniform(-0.075, 0.075, size=len(values))
    ax.scatter(
        np.full(len(values), x) + jitter,
        values,
        s=34,
        color=color,
        alpha=0.95,
        edgecolors="#333333",
        linewidths=0.4,
        zorder=3,
    )

label_heights = [min(1.075, values.max() + 0.055) for values in groups]
for x, mean, label_y in zip(x_positions, means, label_heights):
    ax.text(x, label_y, f"{mean:.4f}", ha="center", va="bottom", fontsize=15, fontweight="bold")

ax.set_xlim(0.4, len(conditions) + 0.6)
ax.set_ylim(0.0, 1.12)
ax.set_xticks(x_positions, conditions)
ax.set_yticks(np.arange(0, 1.01, 0.20))
ax.set_yticks(np.arange(0, 1.01, 0.05), minor=True)
ax.set_ylabel("Fraction of uniquely mapped reads crossing splice junctions")
ax.set_xlabel("Condition")
ax.spines[["top", "right"]].set_visible(False)
for side in ["left", "bottom"]:
    ax.spines[side].set_visible(True)
    ax.spines[side].set_color("#000000")
    ax.spines[side].set_linewidth(1.1)

ax.tick_params(axis="both", which="major", color="#000000", width=1.0)
ax.tick_params(axis="y", which="minor", length=0)
ax.grid(axis="y", which="major", color="#A8A8A8", linewidth=0.8)
ax.grid(axis="y", which="minor", color="#CACACA", linewidth=0.45)
ax.set_axisbelow(True)

fig.tight_layout()
fig.savefig(png_file, dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(pdf_file, bbox_inches="tight", facecolor="white")
plt.close(fig)
PYTHON

echo "Splice junction QC complete."
echo "Samples plotted: $N_SAMPLES"
echo "Combined table: $COMBINED_TSV"
echo "Condition summary: $CONDITION_SUMMARY_TSV"
echo "PNG plot: $PLOT_PNG"
echo "PDF plot: $PLOT_PDF"
