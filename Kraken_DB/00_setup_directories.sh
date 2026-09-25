#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/environment.sh"
load_kraken_config "${1:-${SCRIPT_DIR}/config.sh}"

umask 0022

mkdir -p "${KRAKEN_ROOT}/database"
mkdir -p "${KRAKEN_LOG_DIR}"
mkdir -p "${KRAKEN_DB_DIR}"

apply_kraken_permissions

printf 'Kraken root:     %s\n' "${KRAKEN_ROOT}"
printf 'Database path:  %s\n' "${KRAKEN_DB_DIR}"
printf 'Log directory:  %s\n' "${KRAKEN_LOG_DIR}"
printf 'Hash size cap:  none (full database)\n'
printf 'Shared group:   %s\n' "${KRAKEN_SHARED_GROUP:-default ownership}"
printf 'Permissions:    readable data; group inheritance when configured\n'

df -h "${KRAKEN_ROOT}"
