#!/bin/bash

set -euo pipefail

CONFIG="${1:-}"

if [[ -z "$CONFIG" ]]; then
    echo "Usage: bash Make_Dropoff_Bins.sh path/to/config.sh"
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
: "${GTF:?ERROR: GTF is not set in config}"
: "${OUTDIR:?ERROR: OUTDIR is not set in config}"

# -----------------------------
# Output paths
# -----------------------------
BIN_OUTDIR="${OUTDIR}/annotation/exon_intron_bins"
mkdir -p "$BIN_OUTDIR"

BINS_BED="${BIN_OUTDIR}/exon_intron_bins.bed"
BINS_RAW="${BIN_OUTDIR}/exon_intron_bins.raw_transcript_level.bed"

echo "Creating exon-intron boundary bins"
echo "GTF: $GTF"
echo "Raw transcript-level BED: $BINS_RAW"
echo "Final deduplicated BED: $BINS_BED"

[[ -f "$GTF" ]] || { echo "ERROR: GTF not found: $GTF"; exit 1; }

module purge
module load Anaconda3/2024.06-1

python - "$GTF" "$BINS_RAW" "$BINS_BED" <<'PY'
import re
import sys
from collections import defaultdict

gtf = sys.argv[1]
raw_bed = sys.argv[2]
out_bed = sys.argv[3]

FLANK = 200
BIN = 50

tx_re = re.compile(r'transcript_id "([^"]+)"')

exons = defaultdict(list)

with open(gtf, "r") as f:
    for line in f:
        if not line or line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            continue

        chrom, source, feature, start, end, score, strand, frame, attrs = fields

        if feature != "exon":
            continue

        m = tx_re.search(attrs)
        if not m:
            continue

        tx = m.group(1)
        s0 = int(start) - 1
        e0 = int(end)

        exons[(tx, chrom, strand)].append((s0, e0))


def make_bins_for_boundary(chrom, strand, boundary, exon_interval, intron_interval, boundary_id):
    ex_s, ex_e = exon_interval
    in_s, in_e = intron_interval

    rows = []

    for d0 in range(-FLANK, 0, BIN):
        d1 = d0 + BIN
        label = f"exon_{d0}_{d1}"

        if strand == "+":
            s = boundary + d0
            e = boundary + d1
        else:
            s = boundary - d1
            e = boundary - d0

        s = max(s, ex_s)
        e = min(e, ex_e)

        if s < e:
            rows.append((chrom, s, e, label, boundary_id, strand))

    for d0 in range(0, FLANK, BIN):
        d1 = d0 + BIN
        label = f"intron_{d0}_{d1}"

        if strand == "+":
            s = boundary + d0
            e = boundary + d1
        else:
            s = boundary - d1
            e = boundary - d0

        s = max(s, in_s)
        e = min(e, in_e)

        if s < e:
            rows.append((chrom, s, e, label, boundary_id, strand))

    return rows


out_rows = []
boundary_counter = 0

for (tx, chrom, strand), ex_list in exons.items():
    ex_list.sort()

    for (s1, e1), (s2, e2) in zip(ex_list, ex_list[1:]):
        intron_s = e1
        intron_e = s2

        if intron_e <= intron_s:
            continue

        if strand == "+":
            boundary = e1
            exon_interval = (s1, e1)
            intron_interval = (intron_s, intron_e)
        else:
            boundary = s2
            exon_interval = (s2, e2)
            intron_interval = (intron_s, intron_e)

        boundary_id = f"{tx}|B{boundary_counter}"
        boundary_counter += 1

        out_rows.extend(
            make_bins_for_boundary(
                chrom=chrom,
                strand=strand,
                boundary=boundary,
                exon_interval=exon_interval,
                intron_interval=intron_interval,
                boundary_id=boundary_id,
            )
        )

# -----------------------------
# Write raw transcript-level BED
# -----------------------------
with open(raw_bed, "w") as out:
    for chrom, s, e, label, bid, strand in out_rows:
        out.write(f"{chrom}\t{s}\t{e}\t{label}\t{bid}\t{strand}\n")

# -----------------------------
# Deduplicate bins
# -----------------------------
# The original raw rows are transcript-level.
# For QC, we deduplicate by genomic bin identity:
# chrom, start, end, bin_label, strand.
#
# This prevents genes with many transcript isoforms from producing many
# duplicate or redundant bins and making bedtools coverage explode.
# -----------------------------
seen = set()
dedup_rows = []

for chrom, s, e, label, bid, strand in out_rows:
    key = (chrom, s, e, label, strand)

    if key in seen:
        continue

    seen.add(key)

    # Make a stable artificial boundary/bin ID.
    # Downstream code only needs this column to exist.
    dedup_id = f"DEDUP_{len(dedup_rows)}"

    dedup_rows.append((chrom, s, e, label, dedup_id, strand))

# Sort final rows by chrom/start/end, same as before.
dedup_rows.sort(key=lambda x: (x[0], x[1], x[2], x[3], x[5]))

with open(out_bed, "w") as out:
    for chrom, s, e, label, bid, strand in dedup_rows:
        out.write(f"{chrom}\t{s}\t{e}\t{label}\t{bid}\t{strand}\n")

print(f"Wrote raw transcript-level bins: {raw_bed}")
print(f"Wrote final deduplicated bins: {out_bed}")
print(f"Total transcript boundaries: {boundary_counter}")
print(f"Total raw bins: {len(out_rows)}")
print(f"Total deduplicated bins: {len(dedup_rows)}")
print(f"Reduction factor: {len(out_rows) / len(dedup_rows):.2f}x" if dedup_rows else "Reduction factor: NA")
PY

# Final sort for bedtools safety
sort -k1,1 -k2,2n "$BINS_BED" -o "$BINS_BED"

echo
echo "Done."
echo "Raw transcript-level bins:"
ls -lh "$BINS_RAW"
wc -l "$BINS_RAW"
echo
echo "Final deduplicated bins:"
ls -lh "$BINS_BED"
wc -l "$BINS_BED"