#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/environment.sh"
load_kraken_config "${1:-${SCRIPT_DIR}/config.sh}"

umask 0022

command -v sbatch >/dev/null || { echo 'ERROR: sbatch is unavailable.' >&2; exit 1; }
bash "${SCRIPT_DIR}/00_setup_directories.sh" "$KRAKEN_CONFIG"
bash "${SCRIPT_DIR}/01_check_environment.sh" "$KRAKEN_CONFIG"

submit_stage() {
    local stage="$1" dependency="${2:-}" submission job_id
    local arguments=(--parsable
        "--output=${KRAKEN_LOG_DIR}/${stage}_%j.out"
        "--error=${KRAKEN_LOG_DIR}/${stage}_%j.err")
    if [[ -n "$dependency" ]]; then
        arguments+=("--dependency=afterok:${dependency}" --kill-on-invalid-dep=yes)
    fi
    submission="$(sbatch "${arguments[@]}" "${SCRIPT_DIR}/${stage}.sbatch" "$KRAKEN_CONFIG" "$SCRIPT_DIR")" || return 1
    job_id="${submission%%;*}"
    [[ "$job_id" =~ ^[0-9]+$ ]] || { echo "ERROR: Invalid sbatch response: $submission" >&2; return 1; }
    printf '%s\t%s\t%s\n' "$(utc_now)" "$stage" "$job_id" >> "${KRAKEN_LOG_DIR}/submitted_jobs.tsv" || return 1
    printf '%s\n' "$job_id"
}

download_submission="$(submit_stage 02_download_database)"
download_job_id="${download_submission%%;*}"

build_submission="$(submit_stage 03_build_database "$download_job_id")"
build_job_id="${build_submission%%;*}"

validation_submission="$(submit_stage 04_validate_database "$build_job_id")"
validation_job_id="${validation_submission%%;*}"

printf '\nSubmitted Kraken2 database workflow:\n'
printf '  Download job:    %s\n' "${download_job_id}"
printf '  Build job:       %s (after successful download)\n' "${build_job_id}"
printf '  Validation job:  %s (after successful build)\n' "${validation_job_id}"
printf '\nMonitor with:\n'
printf '  squeue -u %s\n' "$(whoami)"
printf '  tail -f %s/02_download_database_%s.out\n' "${KRAKEN_LOG_DIR}" "${download_job_id}"
printf 'Job log: %s/submitted_jobs.tsv\n' "$KRAKEN_LOG_DIR"
