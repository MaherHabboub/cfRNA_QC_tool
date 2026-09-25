#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/environment.sh"
load_kraken_config "${1:-${SCRIPT_DIR}/config.sh}"

umask 0022

load_kraken_modules

printf 'Loaded modules:\n'
module list 2>&1

required_commands=(kraken2 kraken2-build kraken2-inspect rsync wget dustmasker gzip sha256sum awk xargs readlink)

for required_command in "${required_commands[@]}"
do
    if ! command -v "${required_command}" >/dev/null 2>&1
    then
        printf 'ERROR: required command is unavailable: %s\n' "${required_command}" >&2
        if [[ "${required_command}" == "dustmasker" ]]
        then
            printf 'Run module spider BLAST+ and load an available BLAST+ module before building.\n' >&2
        fi
        exit 1
    fi
    command -v "${required_command}"
done

printf '\nKraken2 versions:\n'
kraken2 --version
kraken2-build --version

build_help="$(kraken2-build --help 2>&1)"
if ! grep -q -- '--skip-maps' <<< "$build_help"
then
    printf 'ERROR: this kraken2-build does not support --skip-maps.\n' >&2
    exit 1
fi

kraken_bin_dir="$(dirname "$(readlink -f "$(command -v kraken2-build)")")"
for helper in scan_fasta_file.pl mask_low_complexity.sh; do
    [[ -x "$kraken_bin_dir/$helper" ]] || {
        printf 'ERROR: Kraken helper is unavailable: %s\n' "$kraken_bin_dir/$helper" >&2
        exit 1
    }
done

printf '\nEnvironment check passed.\n'
printf 'Database path: %s\n' "${KRAKEN_DB_DIR}"
