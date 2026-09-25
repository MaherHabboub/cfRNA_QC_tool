#!/usr/bin/env bash
#SBATCH --job-name=k2_db_build
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=3-00:00:00

set -euo pipefail

# Slurm executes a spooled copy: submit_all passes the original bundle directory.
SCRIPT_DIR="${2:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "${SCRIPT_DIR}/environment.sh"
load_kraken_config "${1:-${SCRIPT_DIR}/config.sh}"

umask 0022

load_kraken_modules

if [[ ! -f "${KRAKEN_DB_DIR}/.taxonomy_download_complete" ]]
then
    printf 'ERROR: taxonomy download is not marked complete.\n' >&2
    exit 1
fi

for ref_library in "${KRAKEN_LIBRARIES[@]}"
do
    completion_marker="${KRAKEN_DB_DIR}/.library_${ref_library}_download_complete"
    if [[ ! -f "${completion_marker}" ]]
    then
        printf 'ERROR: library download is not marked complete: %s\n' "${ref_library}" >&2
        exit 1
    fi
done

if [[ -f "${KRAKEN_DB_DIR}/.database_build_complete" ]] \
    && [[ -s "${KRAKEN_DB_DIR}/hash.k2d" ]] \
    && [[ -s "${KRAKEN_DB_DIR}/opts.k2d" ]] \
    && [[ -s "${KRAKEN_DB_DIR}/taxo.k2d" ]]
then
    printf 'Database is already marked complete; nothing to build.\n'
    exit 0
fi

rm -f "${KRAKEN_DB_DIR}/.database_build_complete" "${KRAKEN_DB_DIR}/.database_validation_complete"
record_stage_time build_started_utc
record_stage_environment build

printf 'Kraken2 database build started: %s\n' "$(utc_now)"
printf 'Host: %s\n' "$(hostname)"
printf 'Database: %s\n' "${KRAKEN_DB_DIR}"
printf 'Threads: %s\n' "${SLURM_CPUS_PER_TASK:-32}"
printf 'Maximum compact-hash size: uncapped (full database)\n'
printf 'k-mer length: %s\n' "${KRAKEN_KMER_LENGTH}"
printf 'Minimizer length: %s\n' "${KRAKEN_MINIMIZER_LENGTH}"
printf 'Minimizer spaces: %s\n' "${KRAKEN_MINIMIZER_SPACES}"

module list 2>&1
free -h
df -h "${KRAKEN_ROOT}"

/usr/bin/time -v kraken2-build \
    --build \
    --threads "${SLURM_CPUS_PER_TASK:-32}" \
    --kmer-len "${KRAKEN_KMER_LENGTH}" \
    --minimizer-len "${KRAKEN_MINIMIZER_LENGTH}" \
    --minimizer-spaces "${KRAKEN_MINIMIZER_SPACES}" \
    --db "${KRAKEN_DB_DIR}"

for required_database_file in hash.k2d opts.k2d taxo.k2d
do
    if [[ ! -s "${KRAKEN_DB_DIR}/${required_database_file}" ]]
    then
        printf 'ERROR: expected database file is missing or empty: %s\n' \
            "${required_database_file}" >&2
        exit 1
    fi
done


printf 'Database build completed: %s\n' "$(utc_now)"
du -sh "${KRAKEN_DB_DIR}"
free -h
df -h "${KRAKEN_ROOT}"

record_stage_time build_completed_utc
touch "${KRAKEN_DB_DIR}/.database_build_complete"
