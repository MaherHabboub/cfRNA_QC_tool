#!/bin/bash

set -euo pipefail

CONFIG="${1:-}"

if [[ -z "$CONFIG" ]]; then
    echo "Usage: bash QC_map.sh path/to/config.sh"
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

RESULT_DIR="${OUTDIR}/mapping"
mkdir -p "$RESULT_DIR"

echo "Running mapping QC..."

# Same parsing helper as original
get_val() {
  local logfile="$1"
  local key="$2"

  awk -F'\\|' -v k="$key" '
    {
      left=$1
      gsub(/^[ \t]+|[ \t]+$/, "", left)

      if (left == k) {

        val=$2
        gsub(/^[ \t]+|[ \t]+$/, "", val)

        print val
        exit
      }
    }
  ' "$logfile"
}

tail -n +2 "$SAMPLESHEET" | while IFS=$'\t' read -r SAMPLE FASTQ1 FASTQ2 BAM STARLOG SJTAB LAYOUT CONDITION
do

    echo "------------------------------------"
    echo "Processing: $SAMPLE"

    if [[ ! -f "$STARLOG" ]]; then
        echo "WARNING: STAR Log.final.out not found for $SAMPLE, skipping"
        continue
    fi

    SAMPLE_OUTDIR="${RESULT_DIR}/${SAMPLE}"
    mkdir -p "$SAMPLE_OUTDIR"

    LOCAL_LOG="${SAMPLE_OUTDIR}/${SAMPLE}.Log.final.out"

    # Same copy behavior
    cp -f "$STARLOG" "$LOCAL_LOG"

    # Same metrics extraction
    input_reads="$(get_val "$STARLOG" 'Number of input reads')"

    uniq_n="$(get_val "$STARLOG" 'Uniquely mapped reads number')"

    uniq_pct="$(get_val "$STARLOG" 'Uniquely mapped reads %')"

    multi_pct="$(get_val "$STARLOG" '% of reads mapped to multiple loci')"

    too_many_pct="$(get_val "$STARLOG" '% of reads mapped to too many loci')"

    unmap_mismatch_pct="$(get_val "$STARLOG" '% of reads unmapped: too many mismatches')"

    unmap_short_pct="$(get_val "$STARLOG" '% of reads unmapped: too short')"

    unmap_other_pct="$(get_val "$STARLOG" '% of reads unmapped: other')"

    # Same mapped % calculation
    mapped_pct="$(python3 - <<PY
def p(x):
    return float(x.replace("%","").strip()) if x else 0.0

u = p("${uniq_pct}")
m = p("${multi_pct}")
t = p("${too_many_pct}")

print(f"{u+m+t:.2f}%")
PY
)"

    # Same unmapped % calculation
    unmapped_pct="$(python3 - <<PY
def p(x):
    return float(x.replace("%","").strip()) if x else 0.0

mm = p("${unmap_mismatch_pct}")
sh = p("${unmap_short_pct}")
ot = p("${unmap_other_pct}")

print(f"{mm+sh+ot:.2f}%")
PY
)"

    OUTTSV="${SAMPLE_OUTDIR}/${SAMPLE}.mapping_summary.tsv"

    {
      echo -e "sample\tinput_reads\tuniquely_mapped_reads\tuniquely_mapped_pct\tmapped_pct\tmulti_mapped_pct\ttoo_many_loci_pct\tunmapped_pct\tunmapped_too_many_mismatches_pct\tunmapped_too_short_pct\tunmapped_other_pct"

      echo -e "${SAMPLE}\t${input_reads}\t${uniq_n}\t${uniq_pct}\t${mapped_pct}\t${multi_pct}\t${too_many_pct}\t${unmapped_pct}\t${unmap_mismatch_pct}\t${unmap_short_pct}\t${unmap_other_pct}"

    } > "$OUTTSV"

    echo "Done: $SAMPLE"

done

echo "Mapping QC complete."