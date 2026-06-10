#!/bin/bash

set -euo pipefail

CONFIG="${1:-}"

if [[ -z "$CONFIG" ]]; then
    echo "Usage: bash submit_dropoff_multiqc_aggregate.sh path/to/config.sh"
    exit 1
fi

if [[ ! -f "$CONFIG" ]]; then
    echo "ERROR: Config file not found: $CONFIG"
    exit 1
fi

# -----------------------------
# Select UGent HPC cluster
# -----------------------------
module swap cluster/doduo || true
module load env/software/doduo || true

# Convert config path to absolute path
CONFIG="$(readlink -f "$CONFIG")"

TOOL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULE_DIR="${TOOL_DIR}/modules"
WRAPPER_DIR="${TOOL_DIR}/wrappers"
LOG_DIR="${TOOL_DIR}/logs"

mkdir -p "$WRAPPER_DIR" "$LOG_DIR"

source "$CONFIG"

: "${OUTDIR:?ERROR: OUTDIR not set in config}"
: "${SAMPLESHEET:?ERROR: SAMPLESHEET not set in config}"

mkdir -p "${OUTDIR}/logs"

SUBMIT_LOG="${OUTDIR}/logs/submitted_dropoff_multiqc_aggregate_jobs.tsv"

echo -e "step\tjob_name\tjob_id\twrapper_script\tdependency\textra_args" > "$SUBMIT_LOG"

echo "============================================================"
echo "Submitting per-sample Dropoff + MultiQC + Aggregate workflow"
echo "TOOL_DIR: $TOOL_DIR"
echo "MODULE_DIR: $MODULE_DIR"
echo "WRAPPER_DIR: $WRAPPER_DIR"
echo "CONFIG: $CONFIG"
echo "SAMPLESHEET: $SAMPLESHEET"
echo "OUTDIR: $OUTDIR"
echo "SUBMIT_LOG: $SUBMIT_LOG"
echo "============================================================"
echo

# -----------------------------
# Sanity checks
# -----------------------------

BINS="${OUTDIR}/annotation/exon_intron_bins/exon_intron_bins.bed"

if [[ ! -f "$BINS" ]]; then
    echo "ERROR: Dropoff bins file not found:"
    echo "$BINS"
    echo
    echo "You need to run the full workflow first, or at least:"
    echo "  GTF_to_BED12.sh"
    echo "  Make_Dropoff_Bins.sh"
    exit 1
fi

echo "Found dropoff bins:"
echo "$BINS"
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
        sbatch_output="$(sbatch "$dependency_arg" "$wrapper")"
    else
        echo "Dependency: none" >&2
        sbatch_output="$(sbatch "$wrapper")"
    fi

    echo "$sbatch_output" >&2

    job_id="$(echo "$sbatch_output" | awk '{print $4}')"

    echo -e "${step_id}\t${job_name}\t${job_id}\t${wrapper}\t${dependency_arg:-none}\t${extra_args:-none}" >> "$SUBMIT_LOG"

    echo "$job_id"
}

# ============================================================
# 1. Run exon-intron dropoff per sample
# ============================================================

dropoff_jobs=()

while IFS=$'\t' read -r SAMPLE FASTQ1 FASTQ2 BAM STARLOG SJTAB LAYOUT CONDITION
do
    [[ -z "${SAMPLE:-}" ]] && continue

    dropoff_job=$(
    submit_step \
        "11" \
        "dropoff_${SAMPLE}" \
        "Dropoff.sh" \
        "12:00:00" \
        "80G" \
        "2" \
        "" \
        "" \
        "$SAMPLE"
    )

    dropoff_jobs+=("$dropoff_job")

done < <(tail -n +2 "$SAMPLESHEET")

dropoff_dep="$(IFS=:; echo "${dropoff_jobs[*]}")"

# ============================================================
# 2. Run MultiQC after all dropoff jobs finish
# Use afterany so MultiQC still runs even if one dropoff job fails.
# ============================================================

multiqc_job=$(
submit_step \
    "12" \
    "multiqc" \
    "Multiqc.sh" \
    "02:00:00" \
    "16G" \
    "2" \
    "afterany" \
    "$dropoff_dep"
)

# ============================================================
# 3. Run aggregation after MultiQC finishes
# ============================================================

aggregate_job=$(
submit_step \
    "13" \
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
echo "Submitted per-sample Dropoff + MultiQC + Aggregate workflow."
echo "Submitted jobs log:"
echo "$SUBMIT_LOG"
echo
echo "Dropoff jobs:"
printf '%s\n' "${dropoff_jobs[@]}"
echo
echo "MultiQC job:"
echo "$multiqc_job"
echo
echo "Final aggregate job:"
echo "$aggregate_job"
echo "============================================================"