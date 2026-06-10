#!/bin/bash

set -euo pipefail

CONFIG="${1:-}"

if [[ -z "$CONFIG" ]]; then
    echo "Usage: bash 05_fragment_size_qc.sh path/to/config.sh"
    exit 1
fi

if [[ ! -f "$CONFIG" ]]; then
    echo "ERROR: Config file not found: $CONFIG"
    exit 1
fi

source "$CONFIG"

: "${SAMPLESHEET:?ERROR: SAMPLESHEET not set in config}"
: "${OUTDIR:?ERROR: OUTDIR not set in config}"

module purge
module load picard/3.0.0-Java-17
module load BEDTools/2.31.1-GCC-13.2.0
module load Anaconda3/2024.06-1

RESULT_DIR="${OUTDIR}/fragment_size"
mkdir -p "$RESULT_DIR"

echo "Running fragment size QC..."

tail -n +2 "$SAMPLESHEET" | while IFS=$'\t' read -r SAMPLE FASTQ1 FASTQ2 BAM STARLOG SJTAB LAYOUT CONDITION
do

    echo "------------------------------------"
    echo "Processing: $SAMPLE"

    if [[ "$LAYOUT" != "PE" ]]; then
        echo "Skipping $SAMPLE (fragment size QC requires paired-end data)"
        continue
    fi

    if [[ ! -f "$BAM" ]]; then
        echo "WARNING: BAM not found for $SAMPLE, skipping"
        continue
    fi

    SAMPLE_OUTDIR="${RESULT_DIR}/${SAMPLE}"
    mkdir -p "$SAMPLE_OUTDIR"

    QBAM="${SAMPLE_OUTDIR}/${SAMPLE}.qnamesort.tmp.bam"

    HISTO="${SAMPLE_OUTDIR}/${SAMPLE}.insert_size_histogram.tsv"
    SUMMARY="${SAMPLE_OUTDIR}/${SAMPLE}.fragment_size_summary.tsv"
    PLOT="${SAMPLE_OUTDIR}/${SAMPLE}.fragment_size_hist.png"

    echo "[1/3] Queryname sort..."

    # Same logic as original
    java -jar "$EBROOTPICARD/picard.jar" SortSam \
      I="$BAM" \
      O="$QBAM" \
      SORT_ORDER=queryname \
      VALIDATION_STRINGENCY=SILENT

    echo "[2/3] BEDPE fragment extraction..."

    # Same logic as original
    bedtools bamtobed -bedpe -i "$QBAM" \
      | awk 'BEGIN{OFS="\t"}
             $1==$4 {
               s=$2; if($5<s) s=$5;
               e=$3; if($6>e) e=$6;
               len=e-s;
               if(len>0) c[len]++
             }
             END{
               for(k in c) print k,c[k]
             }' \
      | sort -k1,1n > "$HISTO"

    rm -f "$QBAM"

    echo "[3/3] Summarizing..."

    python - "$SAMPLE" "$HISTO" "$SUMMARY" "$PLOT" <<'PY'
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

sample = sys.argv[1]
histo = sys.argv[2]
summary_out = sys.argv[3]
plot_out = sys.argv[4]

# Same logic as original
df = pd.read_csv(
    histo,
    sep="\t",
    header=None,
    names=["insert_size","count"]
)

df["insert_size"] = pd.to_numeric(
    df["insert_size"],
    errors="coerce"
)

df["count"] = pd.to_numeric(
    df["count"],
    errors="coerce"
).fillna(0.0)

df = df.dropna(
    subset=["insert_size"]
).sort_values(
    "insert_size"
)

total = float(
    df["count"].sum()
)

if total == 0:
    total = 1.0

def frac(lo, hi):
    m = (
        (df["insert_size"] >= lo) &
        (df["insert_size"] <= hi)
    )

    return float(
        df.loc[m, "count"].sum() / total
    )

# Same biological windows
frac_20_120 = frac(20, 120)
frac_150_180 = frac(150, 180)
frac_300_1000 = frac(300, 1000)

# Same peak enrichment logic
peak = df.loc[
    (df["insert_size"] >= 165) &
    (df["insert_size"] <= 170),
    "count"
].mean()

flank = pd.concat([
    df.loc[
        (df["insert_size"] >= 130) &
        (df["insert_size"] <= 150),
        "count"
    ],
    df.loc[
        (df["insert_size"] >= 180) &
        (df["insert_size"] <= 220),
        "count"
    ]
]).mean()

peak_enrich = (
    float(peak / flank)
    if (not np.isnan(flank) and flank > 0)
    else np.nan
)

out = pd.DataFrame([{
    "sample": sample,
    "total_fragments": total,
    "fraction_20_120": frac_20_120,
    "fraction_150_180": frac_150_180,
    "fraction_300_1000": frac_300_1000,
    "peak167_enrichment": peak_enrich
}])

out.to_csv(
    summary_out,
    sep="\t",
    index=False
)

# Same plotting logic
df_plot = df[
    (df["insert_size"] >= 0) &
    (df["insert_size"] <= 500)
].copy()

plt.figure(figsize=(7,4))

plt.plot(
    df_plot["insert_size"],
    df_plot["count"]
)

plt.axvspan(
    150,
    180,
    alpha=0.2
)

plt.xlabel("Fragment length (bp)")
plt.ylabel("Count")
plt.title(f"Fragment size distribution: {sample}")

plt.tight_layout()

plt.savefig(
    plot_out,
    dpi=200
)
PY

    echo "Done: $SAMPLE"

done

echo "Fragment size QC complete."