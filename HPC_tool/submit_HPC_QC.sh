#!/bin/bash

set -euo pipefail

CONFIG="${1:-}"

if [[ -z "$CONFIG" ]]; then
    echo "Usage: bash submit_HPC_QC.sh path/to/config.sh"
    exit 1
fi

if [[ ! -f "$CONFIG" ]]; then
    echo "ERROR: Config file not found: $CONFIG"
    exit 1
fi

# Convert config path to absolute path
CONFIG="$(readlink -f "$CONFIG")"

TOOL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULE_DIR="${TOOL_DIR}/modules"
WRAPPER_DIR="${TOOL_DIR}/wrappers"
LOG_DIR="${TOOL_DIR}/logs"

mkdir -p "$WRAPPER_DIR" "$LOG_DIR"

source "$CONFIG"

# Slurm resource settings belong to the submitter. Edit these allocations here;
# modules take their thread counts from SLURM_CPUS_PER_TASK.
FASTQC_THREADS=2
DOWNSAMPLE_THREADS=4
STAR_THREADS=8
STAR_MEM="60G"
STAR_TIME="36:00:00"
HTSEQ_THREADS=1
HTSEQ_MEM="60G"
HTSEQ_TIME="36:00:00"

: "${OUTDIR:?ERROR: OUTDIR not set in config}"
: "${SAMPLESHEET:?ERROR: SAMPLESHEET not set in config}"
DOWNSAMPLE_ENABLED="${DOWNSAMPLE_ENABLED:-yes}"
DOWNSAMPLE_TARGET_ALIGNMENTS="${DOWNSAMPLE_TARGET_ALIGNMENTS:-1000000}"
DOWNSAMPLE_SEED="${DOWNSAMPLE_SEED:-42}"
FASTQC_ENABLED="${FASTQC_ENABLED:-yes}"
MAPPING_ENABLED="${MAPPING_ENABLED:-yes}"
DUPLICATION_ENABLED="${DUPLICATION_ENABLED:-yes}"
INSERT_SIZE_ENABLED="${INSERT_SIZE_ENABLED:-yes}"
GENEBODY_ENABLED="${GENEBODY_ENABLED:-yes}"
READ_DISTRIBUTION_ENABLED="${READ_DISTRIBUTION_ENABLED:-yes}"
SPLICE_JUNCTION_ENABLED="${SPLICE_JUNCTION_ENABLED:-yes}"
STRANDEDNESS_ENABLED="${STRANDEDNESS_ENABLED:-yes}"
DROPOFF_ENABLED="${DROPOFF_ENABLED:-yes}"
KRAKEN_ENABLED="${KRAKEN_ENABLED:-yes}"
STAR_ENABLED="${STAR_ENABLED:-yes}"
HTSEQ_ENABLED="${HTSEQ_ENABLED:-yes}"
STAR_OUTDIR="${STAR_OUTDIR:-${OUTDIR}/star}"
HTSEQ_OUTDIR="${HTSEQ_OUTDIR:-${OUTDIR}/htseq}"
HTSEQ_STRANDED="${HTSEQ_STRANDED:-yes}"

if [[ ! "$FASTQC_THREADS" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: FASTQC_THREADS must be a positive integer: $FASTQC_THREADS" >&2
    exit 1
fi

case "$DOWNSAMPLE_ENABLED" in
    yes|no) ;;
    *)
        echo "ERROR: DOWNSAMPLE_ENABLED must be 'yes' or 'no': $DOWNSAMPLE_ENABLED" >&2
        exit 1
        ;;
esac

MODULE_SWITCHES=(
    STAR_ENABLED
    HTSEQ_ENABLED
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
    switch_value="${!switch_name}"
    case "$switch_value" in
        yes|no) ;;
        *)
            echo "ERROR: $switch_name must be 'yes' or 'no': $switch_value" >&2
            exit 1
            ;;
    esac
done

for prefix in STAR HTSEQ; do
    switch_name="${prefix}_ENABLED"
    [[ "${!switch_name}" == "yes" ]] || continue
    value_name="${prefix}_THREADS"
    [[ "${!value_name}" =~ ^[1-9][0-9]*$ ]] || {
        echo "ERROR: $value_name must be a positive integer" >&2; exit 1;
    }
    value_name="${prefix}_MEM"
    [[ "${!value_name}" =~ ^[1-9][0-9]*[KMGTkmgt]?$ ]] || {
        echo "ERROR: $value_name must be a positive Slurm memory value (e.g. 60G)" >&2; exit 1;
    }
    value_name="${prefix}_TIME"
    [[ "${!value_name}" =~ ^([0-9]+-)?[0-9]+:[0-5][0-9]:[0-5][0-9]$ ]] || {
        echo "ERROR: $value_name must use [days-]hours:mm:ss" >&2; exit 1;
    }
done
if [[ "$HTSEQ_ENABLED" == "yes" ]]; then
    case "$HTSEQ_STRANDED" in
        yes|no|reverse) ;;
        *) echo "ERROR: HTSEQ_STRANDED must be yes, no, or reverse" >&2; exit 1 ;;
    esac
fi

for value_name in DOWNSAMPLE_TARGET_ALIGNMENTS DOWNSAMPLE_SEED DOWNSAMPLE_THREADS; do
    value="${!value_name}"
    if [[ ! "$value" =~ ^[0-9]+$ ]]; then
        echo "ERROR: $value_name must be a non-negative integer: $value" >&2
        exit 1
    fi
done

if (( DOWNSAMPLE_TARGET_ALIGNMENTS < 1 || DOWNSAMPLE_THREADS < 1 )); then
    echo "ERROR: DOWNSAMPLE_TARGET_ALIGNMENTS and DOWNSAMPLE_THREADS must be at least 1." >&2
    exit 1
fi

if [[ ! -f "$SAMPLESHEET" ]]; then
    echo "ERROR: Samplesheet not found: $SAMPLESHEET" >&2
    exit 1
fi

# -----------------------------
# Optional HPC cluster setup
# -----------------------------
# Keep cluster selection in the user config so the workflow can be used on a
# different Slurm cluster without editing this submission script. Leave both
# variables empty when the cluster has already been selected externally.
if [[ -n "${CLUSTER_MODULE:-}" || -n "${CLUSTER_ENV_MODULE:-}" ]]; then
    if ! command -v module >/dev/null 2>&1; then
        echo "ERROR: A cluster module was set in the config, but the 'module' command is unavailable." >&2
        exit 1
    fi

    if [[ -n "${CLUSTER_MODULE:-}" ]]; then
        module swap "$CLUSTER_MODULE"
    fi

    if [[ -n "${CLUSTER_ENV_MODULE:-}" ]]; then
        module load "$CLUSTER_ENV_MODULE"
    fi
fi

mkdir -p "${OUTDIR}/logs"

# Resolve relative config paths from the submission working directory before
# Slurm changes directories. Sample paths are relative to the samplesheet.
for path_name in OUTDIR SAMPLESHEET STAR_OUTDIR HTSEQ_OUTDIR GTF EXON_BED STAR_INDEX; do
    path_value="${!path_name:-}"
    if [[ -n "$path_value" && "$path_value" != /* ]]; then
        printf -v "$path_name" '%s/%s' "$PWD" "$path_value"
    fi
done
RUN_DIR="$(mktemp -d "${OUTDIR}/logs/run.XXXXXXXX")"
HPC_RUN_ID="${RUN_DIR##*/}"
RESOLVED_SAMPLESHEET="${RUN_DIR}/samplesheet.tsv"
(
    export SAMPLESHEET STAR_OUTDIR GTF EXON_BED STAR_INDEX DOWNSAMPLE_ENABLED
    export "${MODULE_SWITCHES[@]}"
    # Keep samplesheet preparation here, like the embedded parsers in modules.
    python3 - "$RESOLVED_SAMPLESHEET" <<'PY'
import csv
import os
from pathlib import Path
import re
import sys

COLUMNS = [
    "sample_id", "fastq_r1", "fastq_r2", "bam", "star_log", "sj_tab",
    "layout", "condition", "transcriptome_bam",
]
PATH_COLUMNS = ["fastq_r1", "fastq_r2", "bam", "star_log", "sj_tab", "transcriptome_bam"]


def enabled(name):
    return os.environ.get(name, "yes") == "yes"


def require_file(value, label):
    if value == "NA" or not Path(value).is_file() or not os.access(value, os.R_OK):
        raise ValueError(f"{label} not found or unreadable: {value}")


def read_samples():
    sheet = Path(os.environ["SAMPLESHEET"]).resolve()
    samples, seen = [], set()
    with sheet.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError("Samplesheet header is missing or contains duplicate columns")
        for column in ("sample_id", "layout", "condition"):
            if column not in reader.fieldnames:
                raise ValueError(f"Samplesheet missing required column: {column}")
        for raw in reader:
            if not any(raw.values()):
                continue
            if None in raw:
                raise ValueError(f"Too many fields on samplesheet line {reader.line_num}")
            row = {key: (raw.get(key) or "").strip() for key in COLUMNS}
            sample = row["sample_id"]
            if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", sample) or sample == "ALL":
                raise ValueError(f"Invalid sample_id: {sample!r}; use letters, digits, _, . or - (not ALL)")
            if sample in seen:
                raise ValueError(f"Duplicate sample_id: {sample}")
            seen.add(sample)
            if row["layout"] not in ("PE", "SE"):
                raise ValueError(f"layout must be PE or SE for {sample}")
            if row["condition"] in ("", "NA", "."):
                raise ValueError(f"condition is required for {sample}")
            if any("\n" in value or "\r" in value or "\t" in value for value in row.values()):
                raise ValueError(f"Embedded tabs/newlines are not supported for {sample}")
            for key in PATH_COLUMNS:
                value = row[key]
                # Relative sample paths are relative to the input samplesheet.
                row[key] = "NA" if value in ("", "NA", ".") else str((sheet.parent / value).resolve())
            if enabled("STAR_ENABLED"):
                prefix = Path(os.environ["STAR_OUTDIR"]).resolve() / sample / f"{sample}."
                row["bam"] = f"{prefix}Aligned.sortedByCoord.out.bam"
                row["star_log"] = f"{prefix}Log.final.out"
                row["sj_tab"] = f"{prefix}SJ.out.tab"
            samples.append(row)
    if not samples:
        raise ValueError("Samplesheet contains no samples")
    return samples


def validate_inputs(samples):
    # Annotation jobs always run, regardless of the analysis switches.
    require_file(os.environ.get("GTF", "NA"), "GTF")
    if enabled("STRANDEDNESS_ENABLED"):
        require_file(os.environ.get("EXON_BED", "NA"), "EXON_BED")
    if enabled("STAR_ENABLED"):
        index = os.environ.get("STAR_INDEX", "")
        if not index or not Path(index).is_dir():
            raise ValueError(f"STAR_INDEX directory not found: {index}")
    bam_needed = any(enabled(name) for name in (
        "HTSEQ_ENABLED", "DOWNSAMPLE_ENABLED", "DUPLICATION_ENABLED",
        "INSERT_SIZE_ENABLED", "GENEBODY_ENABLED", "READ_DISTRIBUTION_ENABLED",
        "SPLICE_JUNCTION_ENABLED", "STRANDEDNESS_ENABLED", "DROPOFF_ENABLED",
    ))
    for row in samples:
        sample = row["sample_id"]
        if enabled("STAR_ENABLED") or enabled("FASTQC_ENABLED"):
            require_file(row["fastq_r1"], f"fastq_r1 for {sample}")
            if row["layout"] == "PE":
                require_file(row["fastq_r2"], f"fastq_r2 for {sample}")
                if enabled("STAR_ENABLED") and row["fastq_r1"].endswith(".gz") != row["fastq_r2"].endswith(".gz"):
                    raise ValueError(f"STAR mates must use the same compression for {sample}")
        if row["transcriptome_bam"] != "NA":
            require_file(row["transcriptome_bam"], f"transcriptome_bam for {sample}")
        if not enabled("STAR_ENABLED"):
            if bam_needed:
                require_file(row["bam"], f"bam for {sample}")
            if enabled("MAPPING_ENABLED") or enabled("KRAKEN_ENABLED"):
                require_file(row["star_log"], f"star_log for {sample} (disable mapping/Kraken for non-STAR BAMs)")
            if enabled("SPLICE_JUNCTION_ENABLED"):
                for key in ("sj_tab", "star_log"):
                    if row[key] != "NA":
                        require_file(row[key], f"{key} for {sample}")
            if enabled("KRAKEN_ENABLED"):
                prefix = Path(row["star_log"]).parent / f"{sample}.Unmapped.out.mate"
                require_file(f"{prefix}1", f"STAR unmapped mate1 for {sample}")
                if row["layout"] == "PE":
                    require_file(f"{prefix}2", f"STAR unmapped mate2 for {sample}")


try:
    samples = read_samples()
    validate_inputs(samples)
    with open(sys.argv[1], "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(samples)
except (OSError, ValueError, KeyError) as error:
    sys.exit(f"ERROR: {error}")
PY
)
ORIGINAL_CONFIG="$CONFIG"
CONFIG="${RUN_DIR}/config.sh"
# Snapshot the config and pin resolved inputs so every job sees the same sheet.
cat "$ORIGINAL_CONFIG" > "$CONFIG"
printf '\n# Resolved inputs for this submission.\n' >> "$CONFIG"
SAMPLESHEET="$RESOLVED_SAMPLESHEET"
for value_name in SAMPLESHEET OUTDIR STAR_OUTDIR HTSEQ_OUTDIR GTF EXON_BED STAR_INDEX HPC_RUN_ID \
    HTSEQ_STRANDED DOWNSAMPLE_ENABLED "${MODULE_SWITCHES[@]}"; do
    printf '%s=%q\n' "$value_name" "${!value_name:-}" >> "$CONFIG"
done
# Keep wrappers from previous submissions available for inspection and reruns.
WRAPPER_DIR="${WRAPPER_DIR}/${HPC_RUN_ID}"
mkdir -p "$WRAPPER_DIR"

SUBMIT_LOG="${OUTDIR}/logs/submitted_jobs.tsv"

echo -e "step\tjob_name\tjob_id\twrapper_script\tdependency\textra_args" > "$SUBMIT_LOG"

echo "============================================================"
echo "Submitting HPC QC workflow"
echo "TOOL_DIR: $TOOL_DIR"
echo "MODULE_DIR: $MODULE_DIR"
echo "WRAPPER_DIR: $WRAPPER_DIR"
echo "CONFIG: $CONFIG"
echo "SAMPLESHEET: $SAMPLESHEET"
echo "OUTDIR: $OUTDIR"
echo "DOWNSAMPLE_ENABLED: $DOWNSAMPLE_ENABLED"
if [[ "$DOWNSAMPLE_ENABLED" == "yes" ]]; then
    echo "DOWNSAMPLE_TARGET_ALIGNMENTS: $DOWNSAMPLE_TARGET_ALIGNMENTS"
fi
for switch_name in "${MODULE_SWITCHES[@]}"; do
    echo "$switch_name: ${!switch_name}"
done
echo "SUBMIT_LOG: $SUBMIT_LOG"
echo "============================================================"
echo

submit_step() {
    local step_id="$1"
    local job_name="$2"
    local module_script="$3"
    local time="$4"
    local mem="$5"
    local cpus="$6"
    local dependency_type="${7:-}"
    local dependency_jobs="${8:-}"
    local extra_args="${9:-}"

    local wrapper="${WRAPPER_DIR}/${step_id}_${job_name}.sh"
    local dependency_arg=""

    if [[ -n "$dependency_type" && -n "$dependency_jobs" ]]; then
        dependency_arg="--dependency=${dependency_type}:${dependency_jobs}"
    fi

    if [[ ! -f "${MODULE_DIR}/${module_script}" ]]; then
        echo "ERROR: Module script not found: ${MODULE_DIR}/${module_script}" >&2
        exit 1
    fi

    cat > "$wrapper" <<EOF
#!/bin/bash

#SBATCH -J ${job_name}
#SBATCH -D ${TOOL_DIR}
#SBATCH -t ${time}
#SBATCH --mem=${mem}
#SBATCH --cpus-per-task=${cpus}
#SBATCH --mail-type=FAIL
#SBATCH -o ${LOG_DIR}/${step_id}_${job_name}_%j.out
#SBATCH -e ${LOG_DIR}/${step_id}_${job_name}_%j.err

set -euo pipefail

echo "============================================================"
echo "Running step: ${step_id}_${job_name}"
echo "Module: ${MODULE_DIR}/${module_script}"
echo "Config: ${CONFIG}"
echo "Extra args: ${extra_args}"
echo "Started: \$(date)"
echo "Host: \$(hostname)"
echo "============================================================"

bash "${MODULE_DIR}/${module_script}" "${CONFIG}" ${extra_args}

echo "============================================================"
echo "Finished step: ${step_id}_${job_name}"
echo "Finished: \$(date)"
echo "============================================================"
EOF

    chmod +x "$wrapper"

    echo "Submitting ${step_id}_${job_name}" >&2
    echo "Wrapper: $wrapper" >&2

    if [[ -n "$dependency_arg" ]]; then
        echo "Dependency: $dependency_arg" >&2
        sbatch_output="$(sbatch --kill-on-invalid-dep=yes "$dependency_arg" "$wrapper")"
    else
        echo "Dependency: none" >&2
        sbatch_output="$(sbatch "$wrapper")"
    fi

    echo "$sbatch_output" >&2

    job_id="$(echo "$sbatch_output" | awk '{print $4}')"

    echo -e "${step_id}\t${job_name}\t${job_id}\t${wrapper}\t${dependency_arg:-none}\t${extra_args:-none}" >> "$SUBMIT_LOG"

    # Important: only print job ID to stdout, so command substitution works cleanly
    echo "$job_id"
}

# ============================================================
# Independent annotation preparation
# ============================================================
gtf_job=$(submit_step "01" "gtf_to_bed12" "GTF_to_BED12.sh" "01:00:00" "16G" "1")
bins_job=$(submit_step "02" "make_dropoff_bins" "Make_Dropoff_Bins.sh" "01:00:00" "16G" "1")

# Every submitted preparation/analysis job belongs in the reporting barrier.
qc_jobs=("$gtf_job" "$bins_job")

# ============================================================
# Per-sample processing: wait only for actual input producers.
# An empty dependency list means the inputs already exist.
# ============================================================
while IFS=$'\t' read -r SAMPLE FASTQ1 FASTQ2 BAM STARLOG SJTAB LAYOUT CONDITION TRANSCRIPTOME_BAM
do
    [[ -z "${SAMPLE:-}" ]] && continue
    star_job=""
    downsample_job=""

    if [[ "$STAR_ENABLED" == "yes" ]]; then
        star_job=$(submit_step "00a" "star_${SAMPLE}" "Star_Alignment.sh" \
            "$STAR_TIME" "$STAR_MEM" "$STAR_THREADS" "" "" "$SAMPLE")
        qc_jobs+=("$star_job")
    fi

    if [[ "$FASTQC_ENABLED" == "yes" ]]; then
        fastqc_job=$(submit_step "03" "fastqc_${SAMPLE}" "Fastqc.sh" \
            "01:00:00" "3G" "$FASTQC_THREADS" "" "" "$SAMPLE")
        qc_jobs+=("$fastqc_job")
    fi

    if [[ "$DOWNSAMPLE_ENABLED" == "yes" ]]; then
        downsample_job=$(submit_step "00" "downsample_${SAMPLE}" "Downsample.sh" \
            "04:00:00" "16G" "$DOWNSAMPLE_THREADS" "afterok" "$star_job" "$SAMPLE")
        qc_jobs+=("$downsample_job")
    fi

    if [[ "$HTSEQ_ENABLED" == "yes" ]]; then
        htseq_job=$(submit_step "00b" "htseq_${SAMPLE}" "HTSeq_Counts.sh" \
            "$HTSEQ_TIME" "$HTSEQ_MEM" "$HTSEQ_THREADS" "afterok" "$star_job" "$SAMPLE")
        qc_jobs+=("$htseq_job")
    fi

    if [[ "$MAPPING_ENABLED" == "yes" ]]; then
        map_job=$(submit_step "04" "mapping_${SAMPLE}" "Map.sh" \
            "02:00:00" "16G" "1" "afterok" "$star_job" "$SAMPLE")
        qc_jobs+=("$map_job")
    fi

    # Downsampling already depends on this sample's STAR job when needed.
    duplication_dependency="${downsample_job:-$star_job}"
    genebody_dependency="${gtf_job}${duplication_dependency:+:${duplication_dependency}}"
    read_distribution_dependency="${gtf_job}${star_job:+:${star_job}}"
    dropoff_dependency="${bins_job}${star_job:+:${star_job}}"

    if [[ "$DUPLICATION_ENABLED" == "yes" ]]; then
        dup_job=$(submit_step "05" "duplication_${SAMPLE}" "Duplication.sh" \
            "2:00:00" "40G" "2" "afterok" "$duplication_dependency" "$SAMPLE")
        qc_jobs+=("$dup_job")
    fi

    TRANSCRIPTOME_BAM="${TRANSCRIPTOME_BAM//$'\r'/}"
    if [[ -n "$TRANSCRIPTOME_BAM" && "$TRANSCRIPTOME_BAM" != "NA" && "$TRANSCRIPTOME_BAM" != "." ]]; then
        insert_size_module="Insert_Size_Distribution_Transcriptome.sh"
        insert_size_time="24:00:00"
        insert_size_memory="40G"
        insert_size_cpus="4"
        insert_size_dependency=""
    else
        insert_size_module="Insert_Size_Distribution_Genomic.sh"
        insert_size_time="2:00:00"
        insert_size_memory="20G"
        insert_size_cpus="2"
        insert_size_dependency="$star_job"
    fi
    if [[ "$INSERT_SIZE_ENABLED" == "yes" ]]; then
        insert_size_job=$(submit_step "06" "insert_size_distribution_${SAMPLE}" "$insert_size_module" \
            "$insert_size_time" "$insert_size_memory" "$insert_size_cpus" \
            "afterok" "$insert_size_dependency" "$SAMPLE")
        qc_jobs+=("$insert_size_job")
    fi

    if [[ "$GENEBODY_ENABLED" == "yes" ]]; then
        genebody_job=$(submit_step "07" "genebody_${SAMPLE}" "Genebody.sh" \
            "6:00:00" "8G" "2" "afterok" "$genebody_dependency" "$SAMPLE")
        qc_jobs+=("$genebody_job")
    fi

    if [[ "$READ_DISTRIBUTION_ENABLED" == "yes" ]]; then
        read_dist_job=$(submit_step "08" "read_distribution_${SAMPLE}" "Read_Distribution.sh" \
            "02:00:00" "8G" "1" "afterok" "$read_distribution_dependency" "$SAMPLE")
        qc_jobs+=("$read_dist_job")
    fi

    if [[ "$SPLICE_JUNCTION_ENABLED" == "yes" ]]; then
        splice_job=$(submit_step "09" "splice_junction_${SAMPLE}" "Splice_Junction.sh" \
            "04:00:00" "16G" "1" "afterok" "$star_job" "$SAMPLE")
        qc_jobs+=("$splice_job")
    fi

    if [[ "$STRANDEDNESS_ENABLED" == "yes" ]]; then
        strand_job=$(submit_step "10" "strandedness_${SAMPLE}" "Strandedness.sh" \
            "02:00:00" "4G" "1" "afterok" "$star_job" "$SAMPLE")
        qc_jobs+=("$strand_job")
    fi

    if [[ "$DROPOFF_ENABLED" == "yes" ]]; then
        dropoff_job=$(submit_step "11" "dropoff_${SAMPLE}" "Dropoff.sh" \
            "3:00:00" "6G" "2" "afterok" "$dropoff_dependency" "$SAMPLE")
        qc_jobs+=("$dropoff_job")
    fi

    if [[ "$KRAKEN_ENABLED" == "yes" ]]; then
        kraken_job=$(submit_step "12" "kraken_${SAMPLE}" "Kraken.sh" \
            "04:00:00" "90G" "4" "afterok" "$star_job" "$SAMPLE")
        qc_jobs+=("$kraken_job")
    fi
done < <(tail -n +2 "$SAMPLESHEET")

# Join QC jobs with colon for Slurm dependency
qc_dep="$(IFS=:; echo "${qc_jobs[*]}")"

# ============================================================
# 3. Reporting and aggregation
# Use afterany so these still run even if one QC module fails.
# ============================================================

multiqc_job=$(
submit_step \
    "13" \
    "multiqc" \
    "Multiqc.sh" \
    "02:00:00" \
    "16G" \
    "2" \
    "afterany" \
    "$qc_dep"
)

aggregate_job=$(
submit_step \
    "14" \
    "aggregate" \
    "Aggregate.sh" \
    "01:00:00" \
    "16G" \
    "1" \
    "afterany" \
    "$multiqc_job"
)

echo
echo "============================================================"
echo "Submitted HPC QC workflow."
echo "Submitted jobs log:"
echo "$SUBMIT_LOG"
echo
echo "Final aggregate job:"
echo "$aggregate_job"
echo "Runtime config for manual module reruns: $CONFIG"
echo "============================================================"
