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
        sbatch_output="$(sbatch "$dependency_arg" "$wrapper")"
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
# 1. Annotation jobs
# ============================================================

gtf_job=$(
submit_step \
    "01" \
    "gtf_to_bed12" \
    "GTF_to_BED12.sh" \
    "01:00:00" \
    "16G" \
    "1"
)

bins_job=$(
submit_step \
    "02" \
    "make_dropoff_bins" \
    "Make_Dropoff_Bins.sh" \
    "01:00:00" \
    "16G" \
    "1" \
    "afterok" \
    "$gtf_job"
)

# ============================================================
# 2. QC module jobs
# These run as separate jobs after annotation setup.
# ============================================================

qc_jobs=()

fastqc_job=$(
submit_step \
    "03" \
    "fastqc" \
    "Fastqc.sh" \
    "02:00:00" \
    "16G" \
    "4" \
    "afterok" \
    "$bins_job"
)
qc_jobs+=("$fastqc_job")

map_job=$(
submit_step \
    "04" \
    "mapping" \
    "Map.sh" \
    "02:00:00" \
    "16G" \
    "1" \
    "afterok" \
    "$bins_job"
)
qc_jobs+=("$map_job")

dup_job=$(
submit_step \
    "05" \
    "duplication" \
    "Duplication.sh" \
    "12:00:00" \
    "60G" \
    "2" \
    "afterok" \
    "$bins_job"
)
qc_jobs+=("$dup_job")

frag_job=$(
submit_step \
    "06" \
    "fragmentation" \
    "Fragmentation.sh" \
    "12:00:00" \
    "30G" \
    "2" \
    "afterok" \
    "$bins_job"
)
qc_jobs+=("$frag_job")

# ============================================================
# 2b. Gene body coverage
# One job per sample, because this module can be slow.
# ============================================================

echo "Submitting one Genebody job per sample..." >&2

while IFS=$'\t' read -r SAMPLE FASTQ1 FASTQ2 BAM STARLOG SJTAB LAYOUT CONDITION
do
    [[ -z "${SAMPLE:-}" ]] && continue

    genebody_job=$(
    submit_step \
        "07" \
        "genebody_${SAMPLE}" \
        "Genebody.sh" \
        "24:00:00" \
        "20G" \
        "2" \
        "afterok" \
        "$bins_job" \
        "$SAMPLE"
    )

    qc_jobs+=("$genebody_job")

done < <(tail -n +2 "$SAMPLESHEET")

read_dist_job=$(
submit_step \
    "08" \
    "read_distribution" \
    "Read_Distribution.sh" \
    "02:00:00" \
    "16G" \
    "2" \
    "afterok" \
    "$bins_job"
)
qc_jobs+=("$read_dist_job")

splice_job=$(
submit_step \
    "09" \
    "splice_junction" \
    "Splice_Junction.sh" \
    "02:00:00" \
    "16G" \
    "1" \
    "afterok" \
    "$bins_job"
)
qc_jobs+=("$splice_job")

strand_job=$(
submit_step \
    "10" \
    "strandedness" \
    "Strandedness.sh" \
    "02:00:00" \
    "16G" \
    "2" \
    "afterok" \
    "$bins_job"
)
qc_jobs+=("$strand_job")

# ============================================================
# 2c. Exon-intron dropoff
# One job per sample, because this module can be slow and memory-intensive.
# ============================================================

echo "Submitting one Dropoff job per sample..." >&2

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
        "afterok" \
        "$bins_job" \
        "$SAMPLE"
    )

    qc_jobs+=("$dropoff_job")

done < <(tail -n +2 "$SAMPLESHEET")

# Join QC jobs with colon for Slurm dependency
qc_dep="$(IFS=:; echo "${qc_jobs[*]}")"

# ============================================================
# 3. Reporting and aggregation
# Use afterany so these still run even if one QC module fails.
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
    "$qc_dep"
)

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
echo "Submitted HPC QC workflow."
echo "Submitted jobs log:"
echo "$SUBMIT_LOG"
echo
echo "Final aggregate job:"
echo "$aggregate_job"
echo "============================================================"