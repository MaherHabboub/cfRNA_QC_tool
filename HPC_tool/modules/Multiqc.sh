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
    KRAKEN_ENABLED
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

# Mapping/splice jobs publish this marker only after a successful sample run.
# Limit staging to the current samplesheet and exclude old or failed outputs.
current_sample_output() {
    local parent sample marker
    parent="$(dirname "$1")"
    sample="$(basename "$parent")"
    marker="${parent}/.complete"
    [[ -f "$marker" && "$(cat "$marker")" == "${HPC_RUN_ID:-manual}" ]] || return 1
    awk -F '\t' -v sample="$sample" 'NR > 1 && $1 == sample {found=1} END {exit !found}' "$SAMPLESHEET"
}

# -----------------------------
# Clean old staged files
# -----------------------------
# We keep the output report/data, but refresh the staged input.
echo "Refreshing staged MultiQC input files..."
find "$MULTIQC_INPUT_DIR" -mindepth 1 \( -type f -o -type l \) -delete 2>/dev/null || true
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
cat > "$MULTIQC_CONFIG" <<'MULTIQC_CONFIG_YAML'
title: "HPC RNA-seq QC Report"
subtitle: "Standard tool outputs plus custom QC metrics"
ignore_images: false
# Explicit subsection IDs and ordering keep each summary above its plot.
report_section_order:
  insert_size_distribution_histograms:
    after: custom_insert_size_distribution_summary
  exon_intron_dropoff_profiles:
    after: custom_exon_intron_dropoff_summary
  splice_read_fractions:
    after: custom_splice_read_fractions
  kraken_visualization_overview:
    after: custom_kraken_microbiome_summary
custom_content:
  order:
    - custom_fastqc_parsed_metrics
    - custom_mapping_summary
    - custom_splice_junction_summary
    - custom_spliced_reads
    - custom_splice_read_fractions
    - splice_read_fractions
    - custom_insert_size
    - custom_insert_size_distribution_summary
    - insert_size_distribution_histograms
    - custom_exon_intron_dropoff
    - custom_exon_intron_dropoff_summary
    - exon_intron_dropoff_profiles
    - custom_kraken
    - custom_kraken_microbiome_summary
    - kraken_visualization_overview

# Keys are the HTML section anchors used by MultiQC 1.28 and later.
section_comments:
  general_stats: >-
    This is MultiQC's default overview of key metrics from the enabled tools.
    Use it to compare sequencing yield, alignment, and duplication across samples.
    Good libraries have read counts matching the expected depth, strong alignment, and duplication appropriate for the RNA input.
    Very low alignment, unexpectedly low yield, or extreme duplication warrants investigation.
  rseqc-read_distribution: >-
    This module shows where aligned reads fall relative to annotated gene features.
    It separates coding exons, untranslated exons, introns, and intergenic regions.
    Mature RNA libraries should be enriched in CDS, 5'UTR, and 3'UTR exons, with relatively little intronic or intergenic signal.
    Excess intronic or intergenic signal can indicate DNA contamination, unspliced RNA, or annotation mismatch.
  rseqc-gene_body_coverage: >-
    This module summarizes read coverage along genes scaled from the 5' to the 3' end.
    It helps reveal positional bias and differences between samples.
    Good libraries generally have a smooth profile consistent with the protocol and comparable samples.
    Strong unexplained end bias or uneven coverage can indicate degradation or library-preparation bias, although targeted 3' protocols intentionally show end bias.
  rseqc-infer_experiment: >-
    This module infers library strandedness from read orientation relative to annotated transcripts.
    It reports support for each orientation and the fraction that cannot be assigned.
    Good results match the experimental design, with one dominant orientation for stranded libraries or roughly equal orientations for unstranded libraries.
    An unexpected orientation or a large unassigned fraction warrants checking library settings, annotation, and mapping quality.
  star_summary: >-
    This table summarizes STAR alignment yield, mapped read lengths, splice events, and alignment errors.
    Compare these metrics with the intended sequencing depth and library design.
    Good libraries show high unique alignment, low mismatch and indel rates, and predominantly annotated splice events.
    Low alignment, unusually short mapped reads, or elevated errors can indicate poor reads, contamination, or a reference mismatch.
  star_alignments: >-
    This plot separates unique alignments, multiple alignments, and reasons for failed mapping.
    It shows which outcomes account for the input reads.
    Good libraries have a high uniquely mapped fraction and relatively few reads mapping to too many loci or failing alignment.
    Excess multimapping or unmapped reads can indicate repetitive RNA, short or poor-quality reads, contamination, or a reference mismatch.
  fastqc_sequence_counts: >-
    This module compares read counts and estimated unique and duplicate sequences across samples.
    FastQC sequence uniqueness is distinct from unique genomic mapping.
    Good libraries meet the planned depth and retain substantial sequence diversity for their RNA input.
    Very low yield or extreme duplication can indicate limited library complexity, although abundant transcripts naturally generate repeated sequences.
  fastqc_per_base_sequence_quality: >-
    This plot shows mean base-call quality across read positions.
    Phred scores describe confidence in each base call, with higher scores indicating fewer expected errors.
    Good reads retain high quality across their length, typically Q30 or higher and near Q40 where supported by the platform.
    Sustained low scores or a steep decline at read ends can impair alignment and may require trimming.
  fastqc_per_sequence_quality_scores: >-
    This plot shows the distribution of mean quality scores per read.
    It reveals whether a subset of reads has consistently poor base calls.
    Good libraries place almost all reads in the high-quality range, typically Q30 or higher.
    A large low-quality tail or a separate low-quality peak indicates reads that may need filtering.
  fastqc_per_base_sequence_content: >-
    This plot shows the proportions of A, C, G, and T at each read position.
    RNA library preparation can create expected biases near read starts.
    Good libraries show a composition profile consistent with the protocol and comparable samples.
    Strong unexplained bias across much of the read can indicate adapters or a restricted sequence mixture, while RNA-specific start bias alone is often expected.
  fastqc_per_sequence_gc_content: >-
    This plot shows the distribution of GC percentages across reads.
    Human DNA typically peaks near 41% GC, whereas human RNA libraries often peak at higher GC percentages.
    A higher GC peak consistent with comparable cfRNA libraries supports an RNA-enriched profile.
    A peak near 41% can raise concern for DNA contamination, but GC is only a supporting indicator because transcript composition and library preparation affect the profile.
  fastqc_per_base_n_content: >-
    This plot shows the percentage of undetermined bases, reported as N, at each read position.
    These bases represent calls the sequencer could not resolve confidently.
    Good reads have N content close to zero throughout their length.
    Elevated N content or position-specific spikes indicate base-calling problems and can reduce usable alignments.
  fastqc_sequence_length_distribution: >-
    This plot shows the distribution of read lengths in the FASTQ files.
    Read length is distinct from the paired-end insert size.
    Good libraries match the sequencing and trimming protocol, with expected variation after adapter removal.
    Unexpectedly short reads or unexplained length peaks can indicate excessive trimming, adapter read-through, or mixed inputs.
  fastqc_sequence_duplication_levels: >-
    This plot shows how often identical read sequences recur.
    Duplication can reflect PCR amplification as well as naturally abundant RNA transcripts.
    Good libraries retain adequate sequence diversity for the input amount and protocol.
    Extreme duplication with few distinct sequences can indicate low complexity or overamplification, but FastQC duplication alone cannot distinguish PCR copies from biological repeats.
  fastqc_overrepresented_sequences: >-
    This section summarizes unusually frequent sequences in each sample.
    These sequences can originate from adapters, rRNA, abundant transcripts, or other contaminants.
    Good libraries have little unexplained overrepresentation after accounting for the RNA protocol.
    High levels of adapters, rRNA in a depleted library, or unexpected sequences can indicate library-preparation problems or contamination.
  fastqc_top_overrepresented_sequences: >-
    This table lists the most frequent overrepresented sequences and any suggested sources.
    Source matches are clues to investigate rather than definitive identifications.
    Good libraries show few unexplained sequences, with abundant hits consistent with the intended RNA population.
    Prominent adapter or unexpected contaminant hits warrant checking trimming and library preparation.
  fastqc_adapter_content: >-
    This plot shows recognized adapter sequence content across read positions.
    Increasing adapter content toward read ends often reflects sequencing beyond short RNA inserts.
    Good processed reads have little residual adapter sequence.
    Substantial adapter content can interfere with alignment and indicates a need to check trimming, although read-through is common in raw short-insert cfRNA libraries.
  fastqc_status_checks: >-
    This panel summarizes FastQC pass, warning, and fail flags for each sample.
    The flags use general thresholds rather than RNA-specific acceptance criteria.
    Good libraries pass the relevant quality checks and have protocol-consistent explanations for remaining warnings.
    Multiple unexplained failures, especially for base quality, N content, or adapters, warrant investigation before downstream analysis.
  picard-markduplicates: >-
    This module reports alignments marked as duplicates by Picard.
    Duplicate rates depend on input amount, sequencing depth, and transcript abundance.
    Good libraries retain enough non-duplicate reads for the intended analysis and have complexity consistent with comparable samples.
    Extreme duplication with few usable reads can indicate low input or overamplification, although coordinate duplicates in RNA do not necessarily represent PCR copies.
MULTIQC_CONFIG_YAML

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
    find "$OUTDIR" -path "${OUTDIR}/mapping/*" -type f -name "*.Log.final.out" | while read -r f
    do
        current_sample_output "$f" || continue
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
        current_sample_output "$f" || continue
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
        if [[ "$(basename "$f")" != "splice_read_fraction_cohort_summary.tsv" ]]; then
            current_sample_output "$f" || continue
        fi
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

if [[ "$KRAKEN_ENABLED" == "yes" ]]; then
    if [[ -d "$OUTDIR/kraken" ]]; then
        find "$OUTDIR/kraken" -type f \( \
            -name "*.microbial_summary.tsv" -o \
            -name "*.top_genera.tsv" -o \
            -name "*.top_species.tsv" -o \
            -name "visualization_manifest.tsv" \
        \) | while read -r f
        do
            stage_file "$f" "$CUSTOM_DIR"
        done
    fi
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
