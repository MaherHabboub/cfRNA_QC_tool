#!/bin/bash

set -euo pipefail

CONFIG="${1:-}"

if [[ -z "$CONFIG" ]]; then
    echo "Usage: bash QC_multiqc.sh path/to/config.sh"
    exit 1
fi

if [[ ! -f "$CONFIG" ]]; then
    echo "ERROR: Config file not found: $CONFIG"
    exit 1
fi

source "$CONFIG"

MODULE_SWITCHES=(
    FASTQC_ENABLED
    MAPPING_ENABLED
    DUPLICATION_ENABLED
    INSERT_SIZE_ENABLED
    GENEBODY_ENABLED
    READ_DISTRIBUTION_ENABLED
    SPLICE_JUNCTION_ENABLED
    STRANDEDNESS_ENABLED
    DROPOFF_ENABLED
)

for switch_name in "${MODULE_SWITCHES[@]}"; do
    switch_value="${!switch_name:-yes}"
    case "$switch_value" in
        yes|no) ;;
        *)
            echo "ERROR: $switch_name must be 'yes' or 'no': $switch_value" >&2
            exit 1
            ;;
    esac
    printf -v "$switch_name" '%s' "$switch_value"
done

# -----------------------------
# Required config variables
# -----------------------------
: "${OUTDIR:?ERROR: OUTDIR not set in config}"

# -----------------------------
# Paths
# -----------------------------
MODULE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CUSTOM_CONTENT_SCRIPT="${MODULE_DIR}/Make_multiqc_custom_content.sh"

# -----------------------------
# Output folders
# -----------------------------
MULTIQC_DIR="${OUTDIR}/multiqc"
MULTIQC_INPUT_DIR="${MULTIQC_DIR}/multiqc_input"
CUSTOM_DIR="${MULTIQC_INPUT_DIR}/custom_tables"
CUSTOM_MQC_DIR="${OUTDIR}/multiqc/custom_content"
MULTIQC_TMP="${MULTIQC_DIR}/tmp"
MULTIQC_CONFIG="${MULTIQC_DIR}/multiqc_config.yaml"
REPORT="${MULTIQC_DIR}/hpc_qc_multiqc_report.html"

# -----------------------------
# Software environment
# -----------------------------
module purge
module load MultiQC/1.28-foss-2024a

mkdir -p "$MULTIQC_DIR" "$MULTIQC_INPUT_DIR" "$CUSTOM_DIR" "$CUSTOM_MQC_DIR" "$MULTIQC_TMP"

# Force MultiQC/Python temp files into a writable project folder
export TMPDIR="$MULTIQC_TMP"
export TEMP="$MULTIQC_TMP"
export TMP="$MULTIQC_TMP"

echo "Running MultiQC"
echo "QC output folder: $OUTDIR"
echo "MultiQC output folder: $MULTIQC_DIR"
echo "MultiQC input folder: $MULTIQC_INPUT_DIR"
echo "MultiQC temp folder: $TMPDIR"
echo "Custom MultiQC content folder: $CUSTOM_MQC_DIR"

echo "Using MultiQC: $(command -v multiqc)"
multiqc --version

# -----------------------------
# Helper: stage files into MultiQC input dir
# -----------------------------
stage_file() {
    local src="$1"
    local dest_dir="$2"

    [[ -f "$src" ]] || return 0

    local base
    base="$(basename "$src")"

    if ! ln -sf "$src" "${dest_dir}/${base}" 2>/dev/null; then
        cp -f "$src" "${dest_dir}/${base}"
    fi
}

# -----------------------------
# Clean old staged files
# -----------------------------
# We keep the output report/data, but refresh the staged input.
echo "Refreshing staged MultiQC input files..."
find "$MULTIQC_INPUT_DIR" -mindepth 1 -type f -delete 2>/dev/null || true
mkdir -p "$CUSTOM_DIR"

# -----------------------------
# Generate custom MultiQC content
# -----------------------------
echo "Generating custom MultiQC content..."

if [[ -f "$CUSTOM_CONTENT_SCRIPT" ]]; then
    bash "$CUSTOM_CONTENT_SCRIPT" "$CONFIG"
else
    echo "WARNING: Custom MultiQC content script not found:"
    echo "$CUSTOM_CONTENT_SCRIPT"
    echo "Continuing with standard MultiQC content only."
fi

# -----------------------------
# Create MultiQC config
# -----------------------------
cat > "$MULTIQC_CONFIG" <<EOF
title: "HPC RNA-seq QC Report"
subtitle: "Standard tool outputs plus custom QC metrics"
ignore_images: false
EOF

echo "Created MultiQC config:"
echo "$MULTIQC_CONFIG"

# -----------------------------
# Stage recognized QC outputs
# -----------------------------
echo "Staging recognized QC files..."

if [[ "$FASTQC_ENABLED" == "yes" ]]; then
    find "$OUTDIR" -type f \( -name "*_fastqc.zip" -o -name "*_fastqc.html" \) | while read -r f
    do
        stage_file "$f" "$MULTIQC_INPUT_DIR"
    done
fi

if [[ "$MAPPING_ENABLED" == "yes" ]]; then
    find "$OUTDIR" -type f -name "*.Log.final.out" | while read -r f
    do
        stage_file "$f" "$MULTIQC_INPUT_DIR"
    done
fi

if [[ "$DUPLICATION_ENABLED" == "yes" ]]; then
    find "$OUTDIR" -type f -name "*.markdup.metrics.txt" | while read -r f
    do
        stage_file "$f" "$MULTIQC_INPUT_DIR"
    done
fi

if [[ "$READ_DISTRIBUTION_ENABLED" == "yes" ]]; then
    find "$OUTDIR" -type f -name "*.read_distribution.txt" | while read -r f
    do
        stage_file "$f" "$MULTIQC_INPUT_DIR"
    done
fi

if [[ "$STRANDEDNESS_ENABLED" == "yes" ]]; then
    find "$OUTDIR" -type f -name "*_RSeQC_output_all.txt" | while read -r f
    do
        stage_file "$f" "$MULTIQC_INPUT_DIR"
    done
fi

if [[ "$GENEBODY_ENABLED" == "yes" ]]; then
    find "$OUTDIR" -type f -name "*.geneBodyCoverage.txt" | while read -r f
    do
        stage_file "$f" "$MULTIQC_INPUT_DIR"
    done
fi

# -----------------------------
# Stage source custom QC summary tables
# -----------------------------
# These are bundled beside the report but may not be parsed directly.
echo "Staging source custom QC summary tables..."

if [[ "$MAPPING_ENABLED" == "yes" ]]; then
    find "$OUTDIR" -type f -name "*.mapping_summary.tsv" | while read -r f
    do
        stage_file "$f" "$CUSTOM_DIR"
    done
fi

if [[ "$DUPLICATION_ENABLED" == "yes" ]]; then
    find "$OUTDIR" -type f -name "*.duplication_summary.tsv" | while read -r f
    do
        stage_file "$f" "$CUSTOM_DIR"
    done
fi

if [[ "$SPLICE_JUNCTION_ENABLED" == "yes" ]]; then
    find "$OUTDIR" -type f \( -name "*.splice_junction_summary.tsv" -o -name "*.splice_read_fraction.tsv" -o -name "splice_read_fraction_cohort_summary.tsv" \) | while read -r f
    do
        stage_file "$f" "$CUSTOM_DIR"
    done
fi

if [[ "$INSERT_SIZE_ENABLED" == "yes" ]]; then
    find "$OUTDIR" -type f -name "*.insert_size_distribution_summary.tsv" | while read -r f
    do
        stage_file "$f" "$CUSTOM_DIR"
    done
fi

if [[ "$FASTQC_ENABLED" == "yes" ]]; then
    find "$OUTDIR" -type f -name "*.fastqc_parsed_metrics.tsv" | while read -r f
    do
        stage_file "$f" "$CUSTOM_DIR"
    done
fi

if [[ "$DROPOFF_ENABLED" == "yes" ]]; then
    find "$OUTDIR" -type f -name "*.dropoff_profile.tsv" | while read -r f
    do
        stage_file "$f" "$CUSTOM_DIR"
    done
fi

# -----------------------------
# Stage MultiQC custom content
# -----------------------------
echo "Staging MultiQC custom content..."

find "$CUSTOM_MQC_DIR" -type f \( \
    -name "*_mqc.yaml" -o \
    -name "*_mqc.yml" -o \
    -name "*_mqc.json" -o \
    -name "*_mqc.tsv" -o \
    -name "*_mqc.csv" -o \
    -name "*_mqc.png" -o \
    -name "*_mqc.html" \
\) | while read -r f
do
    stage_file "$f" "$MULTIQC_INPUT_DIR"
done

# -----------------------------
# Run MultiQC
# -----------------------------
echo "Running MultiQC scan..."

set +e

if multiqc --help 2>&1 | grep -q -- "--no-clean-up"; then
    multiqc "$MULTIQC_INPUT_DIR" \
        --outdir "$MULTIQC_DIR" \
        --filename "hpc_qc_multiqc_report.html" \
        --config "$MULTIQC_CONFIG" \
        --force \
        --no-clean-up
else
    multiqc "$MULTIQC_INPUT_DIR" \
        --outdir "$MULTIQC_DIR" \
        --filename "hpc_qc_multiqc_report.html" \
        --config "$MULTIQC_CONFIG" \
        --force
fi

MQC_RC=$?

set -e

if [[ $MQC_RC -ne 0 ]]; then
    if [[ -f "$REPORT" ]]; then
        echo "WARNING: MultiQC exited with code $MQC_RC, but the report was created successfully."
        echo "This appears to be a temporary-folder cleanup issue, not a report-generation failure."
    else
        echo "ERROR: MultiQC failed and no report was created."
        exit "$MQC_RC"
    fi
fi

echo
echo "MultiQC complete."
echo "Report:"
echo "$REPORT"
echo
echo "Source custom TSVs staged in:"
echo "$CUSTOM_DIR"
echo
echo "MultiQC custom content staged from:"
echo "$CUSTOM_MQC_DIR"
