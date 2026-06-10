#!/bin/bash

set -euo pipefail

CONFIG="${1:-}"

if [[ -z "$CONFIG" ]]; then
    echo "Usage: bash QC_fastqc.sh path/to/config.sh"
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

module purge
module load FastQC/0.11.9-Java-11

RESULT_DIR="${OUTDIR}/fastqc/raw"
mkdir -p "$RESULT_DIR"

echo "Running FastQC..."

# -----------------------------
# Helper functions
# -----------------------------

get_data_from_zip() {
  local zip="$1"
  unzip -p "$zip" '*fastqc/fastqc_data.txt'
}

fastqc_zip_name_from_fastq() {
  local fq="$1"
  local base

  base="$(basename "$fq")"

  base="${base%.gz}"
  base="${base%.fastq}"
  base="${base%.fq}"

  echo "${base}_fastqc.zip"
}

parse_fastqc_stream() {
  awk '
    BEGIN{
      FS="\t";
      in_q=0; in_gc=0; in_ad=0;
      n=0; minmean=""; maxad=0; have_ad=0;
      best_gc="NA"; best_c=-1;
    }

    /^>>Per base sequence quality/ {in_q=1; next}
    in_q && /^>>END_MODULE/ {in_q=0}

    /^>>Per sequence GC content/ {in_gc=1; next}
    in_gc && /^>>END_MODULE/ {in_gc=0}

    /^>>Adapter Content/ {in_ad=1; next}
    in_ad && /^>>END_MODULE/ {in_ad=0}

    in_q && $1 ~ /^#/ {next}
    in_q && NF>=2 {
      n++;
      means[n]=$2+0;
      if (minmean=="" || ($2+0)<minmean) minmean=$2+0;
    }

    in_gc && $1 ~ /^#/ {next}
    in_gc && NF>=2 {
      gc=$1+0;
      c=$2+0;
      if (c>best_c){
        best_c=c;
        best_gc=gc;
      }
    }

    in_ad && $1 ~ /^#/ {next}
    in_ad && NF>=2 {
      for (i=2; i<=NF; i++){
        v=$i+0;
        if (!have_ad || v>maxad) maxad=v;
        have_ad=1;
      }
    }

    END{
      if (n==0){
        mean_last10="NA";
        minmean_out="NA";
      } else {
        start=(n>10 ? n-9 : 1);
        sum=0;
        cnt=0;

        for (i=start; i<=n; i++){
          sum+=means[i];
          cnt++;
        }

        mean_last10=sum/cnt;
        minmean_out=minmean;
      }

      maxad_out = (have_ad ? maxad : "NA");

      if (mean_last10=="NA" || minmean_out=="NA") {
        printf("%s\t%s\t%s\t%s\n", mean_last10, minmean_out, best_gc, maxad_out);
      } else {
        printf("%.6f\t%.6f\t%s\t%s\n", mean_last10, minmean_out, best_gc, maxad_out);
      }
    }
  '
}

# -----------------------------
# Process samples
# -----------------------------
tail -n +2 "$SAMPLESHEET" | while IFS=$'\t' read -r SAMPLE FASTQ1 FASTQ2 BAM STARLOG SJTAB LAYOUT CONDITION
do
    echo "------------------------------------"
    echo "Processing: $SAMPLE"
    echo "Layout: $LAYOUT"
    echo "FASTQ1: $FASTQ1"
    echo "FASTQ2: $FASTQ2"

    SAMPLE_OUTDIR="${RESULT_DIR}/${SAMPLE}"
    mkdir -p "$SAMPLE_OUTDIR"

    METRICS_TSV="${SAMPLE_OUTDIR}/${SAMPLE}.fastqc_parsed_metrics.tsv"

    if [[ ! -f "$FASTQ1" ]]; then
        echo "WARNING: FASTQ1 not found for $SAMPLE, skipping"
        continue
    fi

    # -----------------------------
    # Run FastQC
    # -----------------------------
    if [[ "$LAYOUT" == "PE" ]]; then

        if [[ "$FASTQ2" == "NA" || -z "$FASTQ2" ]]; then
            echo "WARNING: FASTQ2 is NA/empty for PE sample $SAMPLE, skipping"
            continue
        fi

        if [[ ! -f "$FASTQ2" ]]; then
            echo "WARNING: FASTQ2 missing for PE sample $SAMPLE: $FASTQ2"
            continue
        fi

        fastqc \
          -t "${THREADS:-4}" \
          -o "$SAMPLE_OUTDIR" \
          "$FASTQ1" "$FASTQ2"

    elif [[ "$LAYOUT" == "SE" ]]; then

        fastqc \
          -t "${THREADS:-4}" \
          -o "$SAMPLE_OUTDIR" \
          "$FASTQ1"

    else
        echo "WARNING: Unknown layout for $SAMPLE: $LAYOUT"
        echo "Expected layout to be SE or PE. Skipping."
        continue
    fi

    # -----------------------------
    # Parse results
    # -----------------------------
    {
        echo -e "sample\tread\tmean_q_last10bp\tmin_mean_q_anypos\tgc_peak_percent\tmax_adapter_percent"

        ZIP1="${SAMPLE_OUTDIR}/$(fastqc_zip_name_from_fastq "$FASTQ1")"

        if [[ ! -f "$ZIP1" ]]; then
            echo "WARNING: FastQC zip not found for R1: $ZIP1" >&2
        else
            m1=$(get_data_from_zip "$ZIP1" | parse_fastqc_stream)
            echo -e "${SAMPLE}\tR1\t${m1}"
        fi

        if [[ "$LAYOUT" == "PE" ]]; then
            ZIP2="${SAMPLE_OUTDIR}/$(fastqc_zip_name_from_fastq "$FASTQ2")"

            if [[ ! -f "$ZIP2" ]]; then
                echo "WARNING: FastQC zip not found for R2: $ZIP2" >&2
            else
                m2=$(get_data_from_zip "$ZIP2" | parse_fastqc_stream)
                echo -e "${SAMPLE}\tR2\t${m2}"
            fi
        fi

    } > "$METRICS_TSV"

    echo "Wrote: $METRICS_TSV"
    echo "Done: $SAMPLE"

done

echo "FastQC complete."