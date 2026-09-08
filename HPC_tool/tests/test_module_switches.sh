#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOL_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
SUBMITTER="${TOOL_DIR}/submit_HPC_QC.sh"
TEST_TMP="$(mktemp -d)"
FAKE_BIN="${TEST_TMP}/bin"
CALL_LOG="${TEST_TMP}/sbatch_calls.tsv"
JOB_COUNTER="${TEST_TMP}/job_counter"
SYSTEM_PATH="$PATH"
OUTDIR_PATH="${TEST_TMP}/out"

cleanup() {
    rm -rf "$TEST_TMP"
}
trap cleanup EXIT

mkdir -p "$FAKE_BIN"
printf '%s\n' 1000 > "$JOB_COUNTER"

printf '%s\n' \
    '#!/bin/bash' \
    'set -euo pipefail' \
    'printf "%s\\n" "$*" >> "$TEST_SBATCH_CALL_LOG"' \
    'job_id="$(cat "$TEST_SBATCH_JOB_COUNTER")"' \
    'printf "%s\\n" "$((job_id + 1))" > "$TEST_SBATCH_JOB_COUNTER"' \
    'echo "Submitted batch job $job_id"' \
    > "${FAKE_BIN}/sbatch"
chmod +x "${FAKE_BIN}/sbatch"

printf '%s\n' \
    '#!/bin/bash' \
    'exit 0' \
    > "${FAKE_BIN}/module"
chmod +x "${FAKE_BIN}/module"

printf '%s\n' \
    '#!/bin/bash' \
    'while IFS= read -r _; do :; done' \
    'exit 0' \
    > "${FAKE_BIN}/python"
chmod +x "${FAKE_BIN}/python"

printf '%s\n' \
    '#!/bin/bash' \
    'set -euo pipefail' \
    'if [[ "${1:-}" == "--help" ]]; then echo "--no-clean-up"; exit 0; fi' \
    'if [[ "${1:-}" == "--version" ]]; then echo "fake MultiQC"; exit 0; fi' \
    'outdir=""' \
    'filename="hpc_qc_multiqc_report.html"' \
    'while [[ $# -gt 0 ]]; do' \
    '    case "$1" in' \
    '        --outdir) outdir="$2"; shift 2 ;;' \
    '        --filename) filename="$2"; shift 2 ;;' \
    '        *) shift ;;' \
    '    esac' \
    'done' \
    'mkdir -p "$outdir"' \
    ': > "${outdir}/${filename}"' \
    > "${FAKE_BIN}/multiqc"
chmod +x "${FAKE_BIN}/multiqc"

write_samplesheet() {
    printf '%s\n' \
        $'sample_id\tfastq_r1\tfastq_r2\tbam\tstar_log\tsj_tab\tlayout\tcondition' \
        $'sample_01\tr1.fastq.gz\tr2.fastq.gz\tsample.bam\tsample.Log.final.out\tsample.SJ.out.tab\tPE\tCONTROL' \
        > "${TEST_TMP}/samples.tsv"
}

write_config() {
    local config_path="$1"
    shift

    printf '%s\n' \
        'CLUSTER_MODULE=""' \
        'CLUSTER_ENV_MODULE=""' \
        'FASTQC_THREADS=2' \
        'DOWNSAMPLE_ENABLED="no"' \
        'DOWNSAMPLE_TARGET_ALIGNMENTS=1000000' \
        'DOWNSAMPLE_SEED=42' \
        'DOWNSAMPLE_THREADS=4' \
        "SAMPLESHEET=\"${TEST_TMP}/samples.tsv\"" \
        "OUTDIR=\"${OUTDIR_PATH}\"" \
        "$@" \
        > "$config_path"
}

assert_contains() {
    local pattern="$1"
    local file="$2"
    grep -q -- "$pattern" "$file" || {
        echo "Expected '$pattern' in $file" >&2
        exit 1
    }
}

assert_not_contains() {
    local pattern="$1"
    local file="$2"
    if grep -q -- "$pattern" "$file"; then
        echo "Did not expect '$pattern' in $file" >&2
        exit 1
    fi
}

run_submitter() {
    PATH="${FAKE_BIN}:$PATH" \
    TEST_SBATCH_CALL_LOG="$CALL_LOG" \
    TEST_SBATCH_JOB_COUNTER="$JOB_COUNTER" \
    bash "$SUBMITTER" "$1" >/dev/null
}

write_samplesheet
OUTDIR_PATH="${TEST_TMP}/out"
write_config "${TEST_TMP}/default_config.sh"
run_submitter "${TEST_TMP}/default_config.sh"
assert_contains 'fastqc_sample_01' "$CALL_LOG"
assert_contains 'mapping' "$CALL_LOG"
assert_contains 'duplication_sample_01' "$CALL_LOG"
assert_contains 'insert_size_distribution_sample_01' "$CALL_LOG"
assert_contains 'genebody_sample_01' "$CALL_LOG"
assert_contains 'read_distribution_sample_01' "$CALL_LOG"
assert_contains 'splice_junction' "$CALL_LOG"
assert_contains 'strandedness_sample_01' "$CALL_LOG"
assert_contains 'dropoff_sample_01' "$CALL_LOG"

: > "$CALL_LOG"
OUTDIR_PATH="${TEST_TMP}/out"
write_config "${TEST_TMP}/disabled_config.sh" \
    'FASTQC_ENABLED="no"' \
    'MAPPING_ENABLED="no"' \
    'DUPLICATION_ENABLED="no"' \
    'INSERT_SIZE_ENABLED="no"' \
    'GENEBODY_ENABLED="no"' \
    'READ_DISTRIBUTION_ENABLED="no"' \
    'SPLICE_JUNCTION_ENABLED="no"' \
    'STRANDEDNESS_ENABLED="no"' \
    'DROPOFF_ENABLED="no"'
run_submitter "${TEST_TMP}/disabled_config.sh"
assert_contains 'gtf_to_bed12' "$CALL_LOG"
assert_contains 'make_dropoff_bins' "$CALL_LOG"
assert_contains 'multiqc' "$CALL_LOG"
assert_contains 'aggregate' "$CALL_LOG"
assert_contains '--dependency=afterany:' "$CALL_LOG"
assert_not_contains 'fastqc_sample_01' "$CALL_LOG"
assert_not_contains 'mapping' "$CALL_LOG"
assert_not_contains 'duplication_sample_01' "$CALL_LOG"
assert_not_contains 'insert_size_distribution_sample_01' "$CALL_LOG"
assert_not_contains 'genebody_sample_01' "$CALL_LOG"
assert_not_contains 'read_distribution_sample_01' "$CALL_LOG"
assert_not_contains 'splice_junction' "$CALL_LOG"
assert_not_contains 'strandedness_sample_01' "$CALL_LOG"
assert_not_contains 'dropoff_sample_01' "$CALL_LOG"

OUTDIR_PATH="${TEST_TMP}/out"
write_config "${TEST_TMP}/invalid_config.sh" 'FASTQC_ENABLED="maybe"'
if PATH="${FAKE_BIN}:$PATH" bash "$SUBMITTER" "${TEST_TMP}/invalid_config.sh" >"${TEST_TMP}/invalid.out" 2>&1; then
    echo "Expected invalid config to fail" >&2
    exit 1
fi
assert_contains "FASTQC_ENABLED must be 'yes' or 'no'" "${TEST_TMP}/invalid.out"

# A stale mapping file must not be staged when mapping is disabled, while an
# enabled FastQC file remains available to MultiQC.
OUTDIR_PATH="${TEST_TMP}/report_out"
mkdir -p \
    "${OUTDIR_PATH}/fastqc/raw/sample_01" \
    "${OUTDIR_PATH}/mapping/sample_01" \
    "${OUTDIR_PATH}/multiqc/multiqc_input/custom_tables"
printf '%s\n' 'FastQC fixture' > "${OUTDIR_PATH}/fastqc/raw/sample_01/sample_01_fastqc.html"
printf '%s\n' 'stale mapping fixture' > "${OUTDIR_PATH}/mapping/sample_01/sample_01.Log.final.out"
printf '%s\n' 'stale staged fixture' > "${OUTDIR_PATH}/multiqc/multiqc_input/custom_tables/stale.mapping_summary.tsv"
write_config "${TEST_TMP}/report_config.sh" \
    'FASTQC_ENABLED="yes"' \
    'MAPPING_ENABLED="no"' \
    'DUPLICATION_ENABLED="no"' \
    'INSERT_SIZE_ENABLED="no"' \
    'GENEBODY_ENABLED="no"' \
    'READ_DISTRIBUTION_ENABLED="no"' \
    'SPLICE_JUNCTION_ENABLED="no"' \
    'STRANDEDNESS_ENABLED="no"' \
    'DROPOFF_ENABLED="no"'
PATH="${FAKE_BIN}:$SYSTEM_PATH" bash "${TOOL_DIR}/modules/Multiqc.sh" "${TEST_TMP}/report_config.sh" >/dev/null
[[ -f "${OUTDIR_PATH}/multiqc/multiqc_input/sample_01_fastqc.html" ]] || {
    echo "Expected enabled FastQC fixture to be staged" >&2
    exit 1
}
[[ ! -e "${OUTDIR_PATH}/multiqc/multiqc_input/sample_01.Log.final.out" ]] || {
    echo "Disabled mapping fixture was staged" >&2
    exit 1
}
[[ ! -e "${OUTDIR_PATH}/multiqc/multiqc_input/custom_tables/stale.mapping_summary.tsv" ]] || {
    echo "Stale staged fixture was not removed" >&2
    exit 1
}

# Run the aggregate fixture assertions wherever its Python dependencies are
# available (normally in the HPC Anaconda module). This environment does not
# install them, so the check is intentionally skipped when unavailable.
if PATH="$SYSTEM_PATH" command -v python >/dev/null 2>&1 && PATH="$SYSTEM_PATH" python -c 'import numpy, pandas' >/dev/null 2>&1; then
    OUTDIR_PATH="${TEST_TMP}/aggregate_out"
    mkdir -p "${OUTDIR_PATH}/mapping/sample_01"
    printf '%s\n' \
        $'sample\tinput_reads\tuniquely_mapped_pct\tmapped_pct\tunmapped_pct\tunmapped_too_short_pct\tmulti_mapped_pct' \
        $'sample_01\t100\t90\t95\t5\t1\t2' \
        > "${OUTDIR_PATH}/mapping/sample_01/sample_01.mapping_summary.tsv"
    write_config "${TEST_TMP}/aggregate_config.sh" \
        'FASTQC_ENABLED="no"' \
        'MAPPING_ENABLED="no"' \
        'DUPLICATION_ENABLED="no"' \
        'INSERT_SIZE_ENABLED="no"' \
        'GENEBODY_ENABLED="no"' \
        'READ_DISTRIBUTION_ENABLED="no"' \
        'SPLICE_JUNCTION_ENABLED="no"' \
        'STRANDEDNESS_ENABLED="no"' \
        'DROPOFF_ENABLED="no"'
    PATH="$SYSTEM_PATH:${FAKE_BIN}" bash "${TOOL_DIR}/modules/Aggregate.sh" "${TEST_TMP}/aggregate_config.sh" >/dev/null
    aggregate_value="$(awk -F '\t' 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == "uniquely_mapped_pct") c = i } NR == 2 { print $c }' "${OUTDIR_PATH}/summary/hpc_qc_summary.tsv")"
    [[ "$aggregate_value" == "NA" ]] || {
        echo "Disabled mapping metric should be NA, found: $aggregate_value" >&2
        exit 1
    }
else
    echo "Skipping aggregate fixture: Python pandas/numpy unavailable."
fi

echo "HPC module switch tests passed."
