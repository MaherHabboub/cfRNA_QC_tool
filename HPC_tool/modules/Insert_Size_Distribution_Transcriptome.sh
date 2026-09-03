#!/bin/bash

set -euo pipefail

CONFIG="${1:-}"
TARGET_SAMPLE="${2:-ALL}"

if [[ -z "$CONFIG" ]]; then
    echo "Usage: bash Insert_Size_Distribution_Transcriptome.sh path/to/config.sh [sample_id]"
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
RESULT_DIR="${OUTDIR}/insert_size_distribution"
THREADS="${SLURM_CPUS_PER_TASK:-4}"

# -----------------------------
# Software environment
# -----------------------------
module purge
module load SAMtools/1.19.2-GCC-13.2.0
module load Anaconda3/2024.06-1

if ! command -v samtools >/dev/null 2>&1; then
    echo "ERROR: samtools is not available after loading the software environment"
    exit 1
fi

if ! command -v python >/dev/null 2>&1; then
    echo "ERROR: python is not available after loading the software environment"
    exit 1
fi

mkdir -p "$RESULT_DIR"

echo "============================================================"
echo "Transcriptome insert size distribution QC"
echo "Collapse concordant transcript placements"
echo "Exclude fragment-length ambiguous read pairs"
echo "Output directory: $RESULT_DIR"
echo "Target sample: $TARGET_SAMPLE"
echo "============================================================"
echo

tail -n +2 "$SAMPLESHEET" | while IFS=$'\t' read -r SAMPLE FASTQ1 FASTQ2 GENOMIC_BAM STARLOG SJTAB LAYOUT CONDITION TRANSCRIPTOME_BAM
do

    if [[ "$TARGET_SAMPLE" != "ALL" && "$SAMPLE" != "$TARGET_SAMPLE" ]]; then
        continue
    fi

    echo "============================================================"
    echo "Processing sample: $SAMPLE"
    echo "============================================================"

    if [[ "$LAYOUT" != "PE" ]]; then
        echo "Skipping $SAMPLE (insert size distribution QC requires paired-end data)"
        continue
    fi

    TRANSCRIPTOME_BAM="${TRANSCRIPTOME_BAM//$'\r'/}"

    if [[ -z "$TRANSCRIPTOME_BAM" || "$TRANSCRIPTOME_BAM" == "NA" || "$TRANSCRIPTOME_BAM" == "." ]]; then
        echo "ERROR: No transcriptome_bam is specified for $SAMPLE in $SAMPLESHEET"
        exit 1
    fi

    BAM="$TRANSCRIPTOME_BAM"

    if [[ ! -f "$BAM" ]]; then
        echo "ERROR: Transcriptome BAM not found for $SAMPLE: $BAM"
        exit 1
    fi

    SAMPLE_OUTDIR="${RESULT_DIR}/${SAMPLE}"
    mkdir -p "$SAMPLE_OUTDIR"

    QBAM="${SAMPLE_OUTDIR}/${SAMPLE}.qnamesort.tmp.bam"

    HISTO="${SAMPLE_OUTDIR}/${SAMPLE}.insert_size_distribution_histogram.tsv"
    CLASSIFICATION="${SAMPLE_OUTDIR}/${SAMPLE}.insert_size_distribution_classification.tsv"
    AMBIG_EXAMPLES="${SAMPLE_OUTDIR}/${SAMPLE}.insert_size_distribution_ambiguous_examples.tsv"

    SUMMARY="${SAMPLE_OUTDIR}/${SAMPLE}.insert_size_distribution_summary.tsv"
    PLOT="${SAMPLE_OUTDIR}/${SAMPLE}.insert_size_distribution_hist.png"

    # ========================================================
    # 1. Query-name sort
    # ========================================================

    echo
    echo "[1/3] Query-name sorting transcriptome BAM..."
    echo "Input BAM:"
    echo "$BAM"

    samtools sort \
        -n \
        -@ "$THREADS" \
        -m 6G \
        -T "${SAMPLE_OUTDIR}/${SAMPLE}.sorttmp" \
        -o "$QBAM" \
        "$BAM"

    # ========================================================
    # 2. Collapse transcript placements by original read pair
    #
    # Example A:
    #
    #   ENST1 -> 60 bp
    #   ENST2 -> 60 bp
    #   ENST3 -> 60 bp
    #
    #   unique lengths = {60}
    #   => count ONE 60-bp fragment
    #
    # Example B:
    #
    #   ENST1 -> 168 bp
    #   ENST2 -> 168 bp
    #   ENST3 -> 240 bp
    #
    #   unique lengths = {168,240}
    #   => ambiguous
    #   => exclude
    #
    # STAR transcript-coordinate TLEN is SAM column 9.
    # abs(TLEN) is used because mates have opposite signs.
    # ========================================================

    echo
    echo "[2/3] Collapsing transcript placements by original read pair..."

    python - \
        "$QBAM" \
        "$HISTO" \
        "$CLASSIFICATION" \
        "$AMBIG_EXAMPLES" \
        "$THREADS" <<'PY'

import sys
import subprocess
from collections import Counter

qbam = sys.argv[1]
histo_out = sys.argv[2]
classification_out = sys.argv[3]
ambiguous_examples_out = sys.argv[4]
threads = sys.argv[5]

histogram = Counter()

total_read_pairs = 0
accepted = 0
ambiguous = 0
no_valid_length = 0

current_read = None
current_lengths = set()
current_records = 0

# Only save a small set of examples for inspection
MAX_AMBIGUOUS_EXAMPLES = 100
ambiguous_examples = []


def finish_read(read_name, lengths, n_records):
    global total_read_pairs
    global accepted
    global ambiguous
    global no_valid_length

    if read_name is None:
        return

    total_read_pairs += 1

    # All transcript placements agree
    if len(lengths) == 1:

        fragment_length = next(iter(lengths))

        histogram[fragment_length] += 1
        accepted += 1

    # Transcript placements disagree on fragment length
    elif len(lengths) > 1:

        ambiguous += 1

        if len(ambiguous_examples) < MAX_AMBIGUOUS_EXAMPLES:
            ambiguous_examples.append(
                (
                    read_name,
                    n_records,
                    len(lengths),
                    ",".join(str(x) for x in sorted(lengths))
                )
            )

    # No useful paired TLEN found
    else:

        no_valid_length += 1


# ============================================================
# Stream the query-name-sorted BAM through samtools
# ============================================================

cmd = [
    "samtools",
    "view",
    "-@",
    threads,
    qbam
]

proc = subprocess.Popen(
    cmd,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE,
    text=True
)

for line in proc.stdout:

    if not line.strip():
        continue

    fields = line.rstrip("\n").split("\t")

    if len(fields) < 9:
        continue

    read_name = fields[0]

    try:
        flag = int(fields[1])
        tlen = int(fields[8])
    except ValueError:
        continue

    # --------------------------------------------------------
    # New QNAME means previous physical fragment is complete
    # --------------------------------------------------------

    if current_read is None:
        current_read = read_name

    elif read_name != current_read:

        finish_read(
            current_read,
            current_lengths,
            current_records
        )

        current_read = read_name
        current_lengths = set()
        current_records = 0

    current_records += 1

    # --------------------------------------------------------
    # Ignore unmapped read or unmapped mate
    #
    # 0x4 = read unmapped
    # 0x8 = mate unmapped
    # --------------------------------------------------------

    if flag & 0x4:
        continue

    if flag & 0x8:
        continue

    # --------------------------------------------------------
    # TLEN = 0 gives no usable paired-fragment estimate
    # --------------------------------------------------------

    if tlen == 0:
        continue

    current_lengths.add(abs(tlen))


# Final read name
finish_read(
    current_read,
    current_lengths,
    current_records
)

proc.stdout.close()

stderr_output = proc.stderr.read()
proc.stderr.close()

return_code = proc.wait()

if return_code != 0:
    raise RuntimeError(
        "samtools view failed with exit code "
        f"{return_code}\n{stderr_output}"
    )


# ============================================================
# Write collapsed histogram
# ============================================================

with open(histo_out, "w") as out:

    for length in sorted(histogram):
        out.write(
            f"{length}\t{histogram[length]}\n"
        )


# ============================================================
# Write classification summary
# ============================================================

with open(classification_out, "w") as out:

    out.write("metric\tvalue\n")

    out.write(
        f"total_read_pairs\t{total_read_pairs}\n"
    )

    out.write(
        f"accepted_same_length\t{accepted}\n"
    )

    out.write(
        f"ambiguous_different_lengths\t{ambiguous}\n"
    )

    out.write(
        f"no_valid_fragment_length\t{no_valid_length}\n"
    )

    if total_read_pairs > 0:

        out.write(
            f"accepted_fraction\t"
            f"{accepted / total_read_pairs:.8f}\n"
        )

        out.write(
            f"ambiguous_fraction\t"
            f"{ambiguous / total_read_pairs:.8f}\n"
        )

        out.write(
            f"no_valid_length_fraction\t"
            f"{no_valid_length / total_read_pairs:.8f}\n"
        )


# ============================================================
# Save examples of ambiguous fragments
# ============================================================

with open(ambiguous_examples_out, "w") as out:

    out.write(
        "read_name\t"
        "bam_records\t"
        "number_of_unique_lengths\t"
        "observed_fragment_lengths\n"
    )

    for row in ambiguous_examples:

        out.write(
            "\t".join(map(str, row)) + "\n"
        )


print("Read-pair classification complete")
print(f"Total read pairs: {total_read_pairs}")
print(f"Accepted same-length pairs: {accepted}")
print(f"Ambiguous different-length pairs: {ambiguous}")
print(f"No valid fragment length: {no_valid_length}")

PY

    # Query-name BAM no longer needed
    rm -f "$QBAM"

    # ========================================================
    # Check histogram
    # ========================================================

    if [[ ! -s "$HISTO" ]]; then
        echo "ERROR: No accepted fragments were produced for $SAMPLE"
        exit 1
    fi

    # ========================================================
    # 3. Summary + plot
    # ========================================================

    echo
    echo "[3/3] Generating insert-size summary and plot..."

    python - \
        "$SAMPLE" \
        "$HISTO" \
        "$CLASSIFICATION" \
        "$SUMMARY" \
        "$PLOT" <<'PY'

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

sample = sys.argv[1]
histo = sys.argv[2]
classification_file = sys.argv[3]
summary_out = sys.argv[4]
plot_out = sys.argv[5]

# ============================================================
# Fragment histogram
# ============================================================

df = pd.read_csv(
    histo,
    sep="\t",
    header=None,
    names=["fragment_length", "count"]
)

df["fragment_length"] = pd.to_numeric(
    df["fragment_length"],
    errors="coerce"
)

df["count"] = pd.to_numeric(
    df["count"],
    errors="coerce"
).fillna(0)

df = (
    df
    .dropna(subset=["fragment_length"])
    .sort_values("fragment_length")
)

total = int(df["count"].sum())

if total == 0:
    raise RuntimeError(
        f"No accepted fragments for {sample}"
    )

# ============================================================
# Classification statistics
# ============================================================

class_df = pd.read_csv(
    classification_file,
    sep="\t"
)

classification = dict(
    zip(
        class_df["metric"],
        class_df["value"]
    )
)

# ============================================================
# Helpers
# ============================================================

def frac(lo, hi):

    mask = (
        (df["fragment_length"] >= lo) &
        (df["fragment_length"] <= hi)
    )

    return float(
        df.loc[mask, "count"].sum() / total
    )


def frac_gt(x):

    return float(
        df.loc[
            df["fragment_length"] > x,
            "count"
        ].sum() / total
    )


def weighted_percentile(values, weights, percentile):

    values = np.asarray(values)
    weights = np.asarray(weights)

    order = np.argsort(values)

    values = values[order]
    weights = weights[order]

    cumulative = np.cumsum(weights)

    cutoff = percentile / 100.0 * cumulative[-1]

    index = np.searchsorted(
        cumulative,
        cutoff,
        side="left"
    )

    return float(values[index])


# ============================================================
# Shared insert-size metrics
# ============================================================

frac_20_120 = frac(20, 120)
frac_150_180 = frac(150, 180)
frac_300_1000 = frac(300, 1000)

peak = df.loc[
    (df["fragment_length"] >= 165) &
    (df["fragment_length"] <= 170),
    "count"
].mean()

flank = pd.concat([
    df.loc[
        (df["fragment_length"] >= 130) &
        (df["fragment_length"] <= 150),
        "count"
    ],

    df.loc[
        (df["fragment_length"] >= 180) &
        (df["fragment_length"] <= 220),
        "count"
    ]

]).mean()

peak_enrich = (
    float(peak / flank)
    if (
        not np.isnan(peak) and
        not np.isnan(flank) and
        flank > 0
    )
    else np.nan
)


# ============================================================
# Additional distribution metrics
# ============================================================

frac_le_500 = frac(1, 500)
frac_500_1000 = frac(501, 1000)
frac_1000_5000 = frac(1001, 5000)
frac_gt_5000 = frac_gt(5000)

lengths = df["fragment_length"].to_numpy()
counts = df["count"].to_numpy()

median = weighted_percentile(
    lengths,
    counts,
    50
)

p90 = weighted_percentile(
    lengths,
    counts,
    90
)

p95 = weighted_percentile(
    lengths,
    counts,
    95
)

p99 = weighted_percentile(
    lengths,
    counts,
    99
)


# ============================================================
# Final summary
# ============================================================

out = pd.DataFrame([{

    "sample": sample,

    "coordinate_system":
        "transcriptome_collapsed",

    # Number entering the actual insert-size distribution
    "total_fragments":
        total,

    # Read-pair classification
    "total_read_pairs_seen":
        int(float(classification["total_read_pairs"])),

    "accepted_same_length":
        int(float(classification["accepted_same_length"])),

    "ambiguous_different_lengths":
        int(float(classification["ambiguous_different_lengths"])),

    "no_valid_fragment_length":
        int(float(classification["no_valid_fragment_length"])),

    "accepted_fraction":
        float(classification["accepted_fraction"]),

    "ambiguous_fraction":
        float(classification["ambiguous_fraction"]),

    "no_valid_length_fraction":
        float(classification["no_valid_length_fraction"]),

    # Shared insert-size metrics
    "fraction_20_120":
        frac_20_120,

    "fraction_150_180":
        frac_150_180,

    "fraction_300_1000":
        frac_300_1000,

    "peak167_enrichment":
        peak_enrich,

    # Additional distribution metrics
    "fraction_le_500":
        frac_le_500,

    "fraction_500_1000":
        frac_500_1000,

    "fraction_1000_5000":
        frac_1000_5000,

    "fraction_gt_5000":
        frac_gt_5000,

    "median_fragment_length":
        median,

    "p90_fragment_length":
        p90,

    "p95_fragment_length":
        p95,

    "p99_fragment_length":
        p99,

    "max_fragment_length":
        float(df["fragment_length"].max())

}])

out.to_csv(
    summary_out,
    sep="\t",
    index=False
)


# ============================================================
# Plot
# ============================================================

df_plot = df[
    (df["fragment_length"] >= 0) &
    (df["fragment_length"] <= 500)
].copy()

plt.figure(
    figsize=(7, 4)
)

plt.plot(
    df_plot["fragment_length"],
    df_plot["count"]
)

plt.axvspan(
    150,
    180,
    alpha=0.2
)

plt.xlabel(
    "Transcriptome insert size (bp)"
)

plt.ylabel(
    "Unique physical fragments"
)

plt.title(
    f"Transcriptome Insert Size Distribution: {sample}"
)

plt.tight_layout()

plt.savefig(
    plot_out,
    dpi=200
)

plt.close()

PY

    echo
    echo "Finished: $SAMPLE"
    echo
    echo "Histogram:"
    echo "$HISTO"
    echo
    echo "Classification:"
    echo "$CLASSIFICATION"
    echo
    echo "Ambiguous examples:"
    echo "$AMBIG_EXAMPLES"
    echo
    echo "Summary:"
    echo "$SUMMARY"
    echo
    echo "Plot:"
    echo "$PLOT"
    echo
    echo "Finished: $(date)"
    echo

done

echo "============================================================"
echo "All transcriptome insert-size analyses complete."
echo "Finished: $(date)"
echo "Results:"
echo "$RESULT_DIR"
echo "============================================================"
