#!/bin/bash

set -euo pipefail

CONFIG="${1:-}"

if [[ -z "$CONFIG" ]]; then
    echo "Usage: bash QC_splice_junction.sh path/to/config.sh"
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

RESULT_DIR="${OUTDIR}/splice_junctions"
mkdir -p "$RESULT_DIR"

echo "Running splice junction QC..."

tail -n +2 "$SAMPLESHEET" | while IFS=$'\t' read -r SAMPLE FASTQ1 FASTQ2 BAM STARLOG SJTAB LAYOUT CONDITION
do

    echo "------------------------------------"
    echo "Processing: $SAMPLE"

    if [[ ! -f "$SJTAB" ]]; then
        echo "WARNING: SJ.out.tab not found for $SAMPLE, skipping"
        continue
    fi

    SAMPLE_OUTDIR="${RESULT_DIR}/${SAMPLE}"
    mkdir -p "$SAMPLE_OUTDIR"

    OUTTSV="${SAMPLE_OUTDIR}/${SAMPLE}.splice_junction_summary.tsv"

    # Keep a copy of STAR final log for QC records, if available
    if [[ -f "$STARLOG" ]]; then
        cp -f "$STARLOG" "$SAMPLE_OUTDIR/${SAMPLE}.Log.final.out"
    else
        echo "WARNING: STAR Log.final.out not found for $SAMPLE, continuing without copying"
    fi

    # Same splice junction summarization logic as original
    awk -v sample="$SAMPLE" 'BEGIN{
      total=0; annotated=0; novel=0; uniq=0; multi=0;
    }
    {
      total++;
      if ($6==1) annotated++; else novel++;
      uniq += $7;
      multi += $8;
    }
    END{
      frac_annot = (total>0 ? annotated/total : 0);
      frac_novel = (total>0 ? novel/total : 0);
      frac_annot_uniq = ((uniq+multi)>0 ? annotated/total : frac_annot); # keep simple
      printf("sample\ttotal_junctions\tannotated_junctions\tnovel_junctions\tfraction_annotated\tfraction_novel\tsum_unique_support\tsum_multi_support\n");
      printf("%s\t%d\t%d\t%d\t%.6f\t%.6f\t%d\t%d\n", sample,total,annotated,novel,frac_annot,frac_novel,uniq,multi);
    }' "$SJTAB" > "$OUTTSV"

    echo "Done: $SAMPLE"
    echo "Wrote: $OUTTSV"

done

echo "Splice junction QC complete."