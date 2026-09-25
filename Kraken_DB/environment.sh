#!/usr/bin/env bash
# Shared software environment and workflow helpers. Source; do not submit.
KRAKEN_ENV_MODULE="env/software/doduo"
KRAKEN_SOFTWARE_MODULE="Kraken2/2.1.3-gompi-2023a"
# Add a cluster-specific BLAST+ module here if dustmasker is unavailable.
KRAKEN_EXTRA_MODULES=()

load_kraken_modules() {
    module --force purge
    module load "$KRAKEN_ENV_MODULE"
    module load "$KRAKEN_SOFTWARE_MODULE"
    local extra_module
    for extra_module in ${KRAKEN_EXTRA_MODULES[@]+"${KRAKEN_EXTRA_MODULES[@]}"}; do
        module load "$extra_module"
    done
}

load_kraken_config() {
    local requested="$1" setting library
    [[ -r "$requested" ]] || { echo "ERROR: Unreadable config: $requested" >&2; return 1; }
    KRAKEN_CONFIG="$(cd "$(dirname "$requested")" && pwd)/$(basename "$requested")"
    source "$KRAKEN_CONFIG"
    [[ "${KRAKEN_ROOT:-}" == /* && "$KRAKEN_ROOT" != / && "$KRAKEN_ROOT" != /path/to/* ]] || {
        echo 'ERROR: Set KRAKEN_ROOT to an absolute storage directory, not / or the example path.' >&2; return 1;
    }
    [[ "${KRAKEN_DB_NAME:-}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || {
        echo 'ERROR: KRAKEN_DB_NAME must be a single safe directory name.' >&2; return 1;
    }
    for setting in KRAKEN_KMER_LENGTH KRAKEN_MINIMIZER_LENGTH KRAKEN_MINIMIZER_SPACES; do
        [[ "${!setting:-}" =~ ^[0-9]+$ ]] || { echo "ERROR: Invalid $setting" >&2; return 1; }
    done
    (( KRAKEN_KMER_LENGTH > 0 && KRAKEN_MINIMIZER_LENGTH > 0 &&
       KRAKEN_MINIMIZER_LENGTH <= 31 && KRAKEN_MINIMIZER_LENGTH <= KRAKEN_KMER_LENGTH &&
       KRAKEN_MINIMIZER_SPACES < KRAKEN_MINIMIZER_LENGTH )) || {
        echo 'ERROR: Invalid k-mer/minimizer settings.' >&2; return 1;
    }
    [[ ${#KRAKEN_LIBRARIES[@]} -gt 0 ]] || { echo 'ERROR: Select at least one library.' >&2; return 1; }
    for library in "${KRAKEN_LIBRARIES[@]}"; do
        case "$library" in archaea|viral|fungi|human|bacteria) ;;
            *) echo "ERROR: Unsupported library: $library" >&2; return 1;;
        esac
    done
    KRAKEN_SHARED_GROUP="${KRAKEN_SHARED_GROUP:-}"
    if [[ -n "$KRAKEN_SHARED_GROUP" ]] && ! id -nG | tr ' ' '\n' | grep -Fx -- "$KRAKEN_SHARED_GROUP" >/dev/null; then
        echo "ERROR: User is not a member of $KRAKEN_SHARED_GROUP" >&2; return 1
    fi
    KRAKEN_DB_DIR="${KRAKEN_ROOT%/}/database/${KRAKEN_DB_NAME}"
    KRAKEN_LOG_DIR="${KRAKEN_ROOT%/}/logs/${KRAKEN_DB_NAME}"
}

apply_kraken_permissions() {
    local directory
    for directory in "$KRAKEN_DB_DIR" "$KRAKEN_LOG_DIR"; do
        if [[ -n "$KRAKEN_SHARED_GROUP" ]]; then
            chgrp -R "$KRAKEN_SHARED_GROUP" "$directory"
            find "$directory" -type d -exec chmod 2755 {} +
        else
            find "$directory" -type d -exec chmod 0755 {} +
        fi
        find "$directory" -type f -exec chmod 0644 {} +
    done
}

utc_now() { date -u '+%Y-%m-%dT%H:%M:%SZ'; }

record_stage_time() {
    local field="$1" directory="${KRAKEN_DB_DIR}/provenance"
    mkdir -p "$directory"
    if [[ ! -s "$directory/$field" ]]; then
        utc_now > "$directory/$field.tmp.$$"
        mv "$directory/$field.tmp.$$" "$directory/$field"
    fi
}

record_stage_environment() {
    local stage="$1" directory="${KRAKEN_DB_DIR}/provenance"
    mkdir -p "$directory"
    kraken2 --version > "$directory/${stage}_kraken2_version.txt.tmp.$$" 2>&1
    mv "$directory/${stage}_kraken2_version.txt.tmp.$$" "$directory/${stage}_kraken2_version.txt"
    module list > "$directory/${stage}_modules.txt.tmp.$$" 2>&1
    mv "$directory/${stage}_modules.txt.tmp.$$" "$directory/${stage}_modules.txt"
    {
        printf 'field\tvalue\n'
        printf 'included_libraries\t%s\n' "${KRAKEN_LIBRARIES[*]}"
        printf 'kmer_length\t%s\n' "$KRAKEN_KMER_LENGTH"
        printf 'minimizer_length\t%s\n' "$KRAKEN_MINIMIZER_LENGTH"
        printf 'minimizer_spaces\t%s\n' "$KRAKEN_MINIMIZER_SPACES"
        printf 'database_hash_size_limit\tuncapped\n'
        printf 'low_complexity_masking\tenabled\n'
    } > "$directory/${stage}_settings.tsv.tmp.$$"
    mv "$directory/${stage}_settings.tsv.tmp.$$" "$directory/${stage}_settings.tsv"
}
