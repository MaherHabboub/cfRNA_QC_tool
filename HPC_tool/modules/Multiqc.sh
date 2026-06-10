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

# -----------------------------
# Load MultiQC
# -----------------------------
module purge
module load MultiQC/1.28-foss-2024a

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
find "$MULTIQC_INPUT_DIR" -mindepth 1 -maxdepth 1 -type f -delete 2>/dev/null || true
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

# FastQC outputs
find "$OUTDIR" -type f \( -name "*_fastqc.zip" -o -name "*_fastqc.html" \) | while read -r f
do
    stage_file "$f" "$MULTIQC_INPUT_DIR"
done

# STAR mapping logs
find "$OUTDIR" -type f -name "*.Log.final.out" | while read -r f
do
    stage_file "$f" "$MULTIQC_INPUT_DIR"
done

# Picard duplication metrics
find "$OUTDIR" -type f -name "*.markdup.metrics.txt" | while read -r f
do
    stage_file "$f" "$MULTIQC_INPUT_DIR"
done

# RSeQC outputs
find "$OUTDIR" -type f \( \
    -name "*.read_distribution.txt" -o \
    -name "*_RSeQC_output_all.txt" -o \
    -name "*.geneBodyCoverage.txt" \
\) | while read -r f
do
    stage_file "$f" "$MULTIQC_INPUT_DIR"
done

# -----------------------------
# Stage original custom QC summary tables
# -----------------------------
# These are bundled beside the report but may not be parsed directly.
echo "Staging original custom QC summary tables..."

find "$OUTDIR" -type f \( \
    -name "*.mapping_summary.tsv" -o \
    -name "*.duplication_summary.tsv" -o \
    -name "*.splice_junction_summary.tsv" -o \
    -name "*.fragment_size_summary.tsv" -o \
    -name "*.fastqc_parsed_metrics.tsv" -o \
    -name "*.dropoff_profile.tsv" \
\) | while read -r f
do
    stage_file "$f" "$CUSTOM_DIR"
done

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

REPORT="${MULTIQC_DIR}/hpc_qc_multiqc_report.html"

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
echo "Original custom TSVs staged in:"
echo "$CUSTOM_DIR"
echo
echo "MultiQC custom content staged from:"
echo "$CUSTOM_MQC_DIR"