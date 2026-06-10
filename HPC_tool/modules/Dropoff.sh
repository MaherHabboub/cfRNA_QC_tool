#!/bin/bash

set -euo pipefail

CONFIG="${1:-}"
TARGET_SAMPLE="${2:-ALL}"

if [[ -z "$CONFIG" ]]; then
    echo "Usage: bash Dropoff.sh path/to/config.sh [sample_id]"
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

BINS="${OUTDIR}/annotation/exon_intron_bins/exon_intron_bins.bed"

[[ -f "$BINS" ]] || {
    echo "ERROR: Bin file not found: $BINS"
    exit 1
}

module purge
module load BEDTools/2.31.1-GCC-13.2.0
module load Anaconda3/2024.06-1

# Try loading SAMtools. If the exact version is unavailable, try default SAMtools.
module load SAMtools/1.19.2-GCC-13.2.0 2>/dev/null || module load SAMtools 2>/dev/null || true

if ! command -v samtools >/dev/null 2>&1; then
    echo "ERROR: samtools not available. Please load a SAMtools module in Dropoff.sh."
    exit 1
fi

RESULT_DIR="${OUTDIR}/dropoff"
mkdir -p "$RESULT_DIR"

echo "Running exon-intron dropoff QC"
echo "Target sample: $TARGET_SAMPLE"
echo "Bins file: $BINS"

# Skip header, read samplesheet
tail -n +2 "$SAMPLESHEET" | while IFS=$'\t' read -r SAMPLE FASTQ1 FASTQ2 BAM STARLOG SJTAB LAYOUT CONDITION
do

    if [[ "$TARGET_SAMPLE" != "ALL" && "$SAMPLE" != "$TARGET_SAMPLE" ]]; then
        continue
    fi

    echo "------------------------------------"
    echo "Processing: $SAMPLE"

    if [[ ! -f "$BAM" ]]; then
        echo "WARNING: BAM not found for $SAMPLE, skipping"
        continue
    fi

    SAMPLE_OUTDIR="${RESULT_DIR}/${SAMPLE}"
    mkdir -p "$SAMPLE_OUTDIR"

    RAW="${SAMPLE_OUTDIR}/${SAMPLE}.bin_coverage.tsv"
    PROFILE="${SAMPLE_OUTDIR}/${SAMPLE}.dropoff_profile.tsv"
    PLOT="${SAMPLE_OUTDIR}/${SAMPLE}.dropoff_profile.png"

    GENOME_FILE="${SAMPLE_OUTDIR}/${SAMPLE}.genome_from_bam.tsv"
    FILTERED_BINS="${SAMPLE_OUTDIR}/${SAMPLE}.exon_intron_bins.filtered.bed"
    SORTED_BINS="${SAMPLE_OUTDIR}/${SAMPLE}.exon_intron_bins.sorted.bed"

    echo "Creating genome file from BAM header..."

    samtools view -H "$BAM" \
      | awk '$1=="@SQ"{
          sn="";
          ln="";
          for(i=1;i<=NF;i++){
              if($i ~ /^SN:/){sn=substr($i,4)}
              if($i ~ /^LN:/){ln=substr($i,4)}
          }
          if(sn!="" && ln!=""){print sn"\t"ln}
        }' > "$GENOME_FILE"

    if [[ ! -s "$GENOME_FILE" ]]; then
        echo "ERROR: Genome file could not be created from BAM header for $SAMPLE"
        exit 1
    fi

    echo "Filtering bins to chromosomes present in BAM header..."

    awk 'BEGIN{FS=OFS="\t"}
         NR==FNR {keep[$1]=1; next}
         ($1 in keep)' "$GENOME_FILE" "$BINS" > "$FILTERED_BINS"

    if [[ ! -s "$FILTERED_BINS" ]]; then
        echo "ERROR: Filtered bins file is empty for $SAMPLE"
        exit 1
    fi

    echo "Sorting bins according to BAM genome order..."

    bedtools sort \
      -i "$FILTERED_BINS" \
      -g "$GENOME_FILE" \
      > "$SORTED_BINS"

    if [[ ! -s "$SORTED_BINS" ]]; then
        echo "ERROR: Sorted bins file is empty for $SAMPLE"
        exit 1
    fi

    echo "Running bedtools coverage with -sorted and genome file..."

    bedtools coverage \
      -a "$SORTED_BINS" \
      -b "$BAM" \
      -split \
      -sorted \
      -g "$GENOME_FILE" \
      -counts > "$RAW"

    echo "Generating profile..."

    python - "$SAMPLE" "$RAW" "$PROFILE" "$PLOT" <<'PY'
import pandas as pd
import matplotlib.pyplot as plt
import re
import sys

sample = sys.argv[1]
raw = sys.argv[2]
out_tsv = sys.argv[3]
out_png = sys.argv[4]

df = pd.read_csv(raw, sep="\t", header=None)

df.columns = [
    "chrom",
    "start",
    "end",
    "bin_label",
    "boundary_id",
    "strand",
    "count"
]

mids = []

for lab in df["bin_label"]:
    m = re.match(r'^(exon|intron)_(-?\d+)_(-?\d+)$', lab)

    if not m:
        mids.append(None)
        continue

    a = int(m.group(2))
    b = int(m.group(3))

    mids.append((a + b) / 2.0)

df["dist_mid"] = mids

profile = df.groupby(
    ["bin_label", "dist_mid"],
    as_index=False
)["count"].mean()

profile = profile.sort_values("dist_mid")

exon_mean = profile.loc[
    profile["dist_mid"] < 0,
    "count"
].mean()

if exon_mean == 0 or pd.isna(exon_mean):
    exon_mean = 1.0

profile["norm_count"] = profile["count"] / exon_mean

profile.to_csv(
    out_tsv,
    sep="\t",
    index=False
)

plt.figure(figsize=(7,4))

plt.plot(
    profile["dist_mid"],
    profile["norm_count"],
    marker="o"
)

plt.axvline(
    0,
    linestyle="--"
)

# Force y-axis to start at 0 and end at the maximum value, with slight padding
ymax = profile["norm_count"].max()

if pd.isna(ymax) or ymax <= 0:
    ymax = 1.0

plt.ylim(0, ymax * 1.05)

plt.xlabel(
    "Distance from exon–intron boundary (bp)\n"
    "(negative=exon side, positive=intron side)"
)

plt.ylabel(
    "Normalized mean read count (exon mean = 1.0)"
)

plt.title(
    f"Exon–Intron Drop-off: {sample}"
)

plt.tight_layout()

plt.savefig(
    out_png,
    dpi=200
)
PY

    echo "Done: $SAMPLE"

done

echo "Exon-intron dropoff QC complete."