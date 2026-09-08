#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOCAL_TOOL_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
REPO_DIR="$(cd "${LOCAL_TOOL_DIR}/.." && pwd)"
RUNNER="${LOCAL_TOOL_DIR}/scripts/run_all_local_qc.R"
REPORT_SCRIPT="${LOCAL_TOOL_DIR}/scripts/08_generate_html_report.R"
REPORT_TEMPLATE="${LOCAL_TOOL_DIR}/scripts/08_report_template.Rmd"
REAL_RSCRIPT="$(command -v Rscript)"
TEST_TMP="$(mktemp -d)"
FAKE_BIN="${TEST_TMP}/bin"
CALL_LOG="${TEST_TMP}/child_rscript_calls.tsv"

cleanup() {
    rm -rf "$TEST_TMP"
}
trap cleanup EXIT

mkdir -p "$FAKE_BIN"

printf '%s\n' \
    '#!/bin/bash' \
    'set -euo pipefail' \
    'printf "%s\\n" "$*" >> "$LOCAL_QC_CHILD_CALL_LOG"' \
    'exit 0' \
    > "${FAKE_BIN}/Rscript"
chmod +x "${FAKE_BIN}/Rscript"

write_config() {
    local path="$1"
    local out_root="$2"
    shift 2

    printf '%s\n' \
        'COUNTS <- "/test/counts.tsv"' \
        "OUT_ROOT <- \"${out_root}\"" \
        'METADATA <- ""' \
        'HPC_SUMMARY <- ""' \
        'BIOTYPE_ENABLED <- TRUE' \
        'PCA_SSGSEA_ENABLED <- TRUE' \
        'SEX_INFERENCE_ENABLED <- TRUE' \
        'TOP_N <- 1000L' \
        'MAX_PCS <- 5L' \
        'TOP_SCATTER <- 10L' \
        'REPORT_TITLE <- "Test report"' \
        'STOP_ON_FAIL <- FALSE' \
        "$@" \
        > "$path"
}

assert_contains() {
    local pattern="$1"
    local path="$2"
    grep -q -- "$pattern" "$path" || {
        echo "Expected '$pattern' in $path" >&2
        exit 1
    }
}

assert_not_contains() {
    local pattern="$1"
    local path="$2"
    if grep -q -- "$pattern" "$path"; then
        echo "Did not expect '$pattern' in $path" >&2
        exit 1
    fi
}

# Default selection runs every local analysis module except HPC correlation,
# because no HPC_SUMMARY path was configured.
DEFAULT_OUT="${TEST_TMP}/default_out"
write_config "${TEST_TMP}/default_config.R" "$DEFAULT_OUT"
PATH="${FAKE_BIN}:$PATH" LOCAL_QC_CHILD_CALL_LOG="$CALL_LOG" \
    "$REAL_RSCRIPT" "$RUNNER" --config "${TEST_TMP}/default_config.R" >/dev/null
assert_contains '01_normalize_log.R' "$CALL_LOG"
assert_contains '02_biotype_distribution.R' "$CALL_LOG"
assert_contains '03_pca_top_variable_genes.R' "$CALL_LOG"
assert_contains '04_ssgsea_score.R' "$CALL_LOG"
assert_contains '05_sex_inference_XIST_vs_Ypanel.R' "$CALL_LOG"
assert_contains '08_generate_html_report.R' "$CALL_LOG"
assert_not_contains '07_hpc_metric_pc_correlation.R' "$CALL_LOG"
assert_contains $'07_hpc_metric_pc_correlation\t.*\tSKIPPED\t' "${DEFAULT_OUT}/run_logs/local_qc_run_status.tsv"

# Disabling modules prevents child invocation and records explicit reasons.
: > "$CALL_LOG"
DISABLED_OUT="${TEST_TMP}/disabled_out"
write_config "${TEST_TMP}/disabled_config.R" "$DISABLED_OUT" \
    'BIOTYPE_ENABLED <- FALSE' \
    'PCA_SSGSEA_ENABLED <- FALSE' \
    'SEX_INFERENCE_ENABLED <- FALSE' \
    'HPC_SUMMARY <- "/test/hpc_summary.tsv"'
PATH="${FAKE_BIN}:$PATH" LOCAL_QC_CHILD_CALL_LOG="$CALL_LOG" \
    "$REAL_RSCRIPT" "$RUNNER" --config "${TEST_TMP}/disabled_config.R" >/dev/null
assert_contains '01_normalize_log.R' "$CALL_LOG"
assert_contains '08_generate_html_report.R' "$CALL_LOG"
assert_not_contains '02_biotype_distribution.R' "$CALL_LOG"
assert_not_contains '03_pca_top_variable_genes.R' "$CALL_LOG"
assert_not_contains '04_ssgsea_score.R' "$CALL_LOG"
assert_not_contains '05_sex_inference_XIST_vs_Ypanel.R' "$CALL_LOG"
assert_not_contains '07_hpc_metric_pc_correlation.R' "$CALL_LOG"
assert_contains 'Disabled by BIOTYPE_ENABLED.' "${DISABLED_OUT}/run_logs/local_qc_run_status.tsv"
assert_contains 'Disabled by PCA_SSGSEA_ENABLED.' "${DISABLED_OUT}/run_logs/local_qc_run_status.tsv"
assert_contains 'Disabled by SEX_INFERENCE_ENABLED.' "${DISABLED_OUT}/run_logs/local_qc_run_status.tsv"
assert_contains 'PCA_SSGSEA_ENABLED is FALSE; HPC metric correlation requires PCA scores.' "${DISABLED_OUT}/run_logs/local_qc_run_status.tsv"

# Config validation fails before any child analysis script is scheduled.
: > "$CALL_LOG"
INVALID_OUT="${TEST_TMP}/invalid_out"
write_config "${TEST_TMP}/invalid_config.R" "$INVALID_OUT" 'BIOTYPE_ENABLED <- "yes"'
if PATH="${FAKE_BIN}:$PATH" LOCAL_QC_CHILD_CALL_LOG="$CALL_LOG" \
    "$REAL_RSCRIPT" "$RUNNER" --config "${TEST_TMP}/invalid_config.R" >"${TEST_TMP}/invalid.out" 2>&1; then
    echo "Expected invalid module switch config to fail" >&2
    exit 1
fi
assert_contains 'BIOTYPE_ENABLED must be TRUE or FALSE.' "${TEST_TMP}/invalid.out"
[[ ! -s "$CALL_LOG" ]] || {
    echo "Invalid config should not invoke child scripts" >&2
    exit 1
}

printf '%s\n' \
    'OUT_ROOT <- "/test/missing_counts_out"' \
    > "${TEST_TMP}/missing_counts_config.R"
if "$REAL_RSCRIPT" "$RUNNER" --config "${TEST_TMP}/missing_counts_config.R" >"${TEST_TMP}/missing_counts.out" 2>&1; then
    echo "Expected config without COUNTS to fail" >&2
    exit 1
fi
assert_contains 'Config is missing required setting: COUNTS' "${TEST_TMP}/missing_counts.out"

# The rendered report omits disabled sections when given a manifest.
REPORT_RESULTS="${TEST_TMP}/report_results"
REPORT_OUT="${TEST_TMP}/report_out"
mkdir -p "$REPORT_RESULTS/run_logs"
printf '%s\n' \
    $'module\tenabled\treason' \
    $'normalization\tTRUE\tMandatory module.' \
    $'biotype\tFALSE\tDisabled' \
    $'pca_ssgsea\tFALSE\tDisabled' \
    $'sex_inference\tFALSE\tDisabled' \
    $'hpc_metric\tFALSE\tDisabled' \
    $'report\tTRUE\tMandatory module.' \
    > "${REPORT_RESULTS}/run_logs/local_qc_module_manifest.tsv"
"$REAL_RSCRIPT" "$REPORT_SCRIPT" \
    --results "$REPORT_RESULTS" \
    --template "$REPORT_TEMPLATE" \
    --out "$REPORT_OUT" \
    --manifest "${REPORT_RESULTS}/run_logs/local_qc_module_manifest.tsv" >/dev/null
REPORT_HTML="${REPORT_OUT}/local_qc_report.html"
assert_not_contains 'Biotype distribution QC' "$REPORT_HTML"
assert_not_contains 'PCA on top variable genes' "$REPORT_HTML"
assert_not_contains 'Gene-set contamination scores' "$REPORT_HTML"
assert_not_contains 'Sex chromosome QC' "$REPORT_HTML"
assert_not_contains 'HPC QC metric correlation with expression PCs' "$REPORT_HTML"
assert_contains 'Normalization and log-transformation QC' "$REPORT_HTML"

echo "Local module switch tests passed."
