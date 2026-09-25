#!/usr/bin/env bash
#SBATCH --job-name=k2_db_validate
#SBATCH --cpus-per-task=4
#SBATCH --mem=200G
#SBATCH --time=1-00:00:00

set -euo pipefail

# Slurm executes a spooled copy: submit_all passes the original bundle directory.
SCRIPT_DIR="${2:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "${SCRIPT_DIR}/environment.sh"
load_kraken_config "${1:-${SCRIPT_DIR}/config.sh}"

umask 0022

# Invalidate previous success before any operation that can fail.
rm -f "${KRAKEN_DB_DIR}/.database_validation_complete"
load_kraken_modules
record_stage_time validation_started_utc
record_stage_environment validation

if [[ ! -f "${KRAKEN_DB_DIR}/.database_build_complete" ]]
then
    printf 'ERROR: database build is not marked complete.\n' >&2
    exit 1
fi

for required_database_file in hash.k2d opts.k2d taxo.k2d
do
    if [[ ! -s "${KRAKEN_DB_DIR}/${required_database_file}" ]]
    then
        printf 'ERROR: expected database file is missing or empty: %s\n' \
            "${required_database_file}" >&2
        exit 1
    fi
done

printf 'Database validation started: %s\n' "$(utc_now)"
printf 'Database: %s\n' "${KRAKEN_DB_DIR}"

/usr/bin/time -v kraken2-inspect \
    --db "${KRAKEN_DB_DIR}" \
    > "${KRAKEN_DB_DIR}/inspect.txt"

# Validate the taxa corresponding to the configured libraries.
for library in "${KRAKEN_LIBRARIES[@]}"; do
    case "$library" in
        bacteria) taxid=2;; archaea) taxid=2157;; viral) taxid=10239;;
        fungi) taxid=4751;; human) taxid=9606;;
    esac
    if ! awk -F '\t' -v required_taxid="${taxid}" '$5 == required_taxid {found=1} END {exit !found}' \
        "${KRAKEN_DB_DIR}/inspect.txt"
    then
        printf 'ERROR: expected taxon was not found: %s (taxid %s)\n' \
            "$library" "${taxid}" >&2
        exit 1
    fi
done

awk -F '\t' '$5 == 2 || $5 == 2157 || $5 == 10239 || $5 == 4751 || $5 == 9606' \
    "${KRAKEN_DB_DIR}/inspect.txt" \
    > "${KRAKEN_DB_DIR}/expected_taxa.inspect.tsv"

# Relative paths keep checksums usable from inside the database directory.
(
    cd "$KRAKEN_DB_DIR"
    sha256sum hash.k2d opts.k2d taxo.k2d > checksums.sha256.tmp
    mv checksums.sha256.tmp checksums.sha256
    find library taxonomy -type f \
        \( -name 'assembly_summary.txt' -o -name 'manifest.txt' \
        -o -name 'https_manifest.tsv' -o -name 'names.dmp' -o -name 'nodes.dmp' \) \
        -print0 > provenance/reference_files.list
    while IFS= read -r -d '' reference_metadata_file; do
        sha256sum "$reference_metadata_file"
    done < provenance/reference_files.list > reference_metadata_checksums.sha256.tmp
    mv reference_metadata_checksums.sha256.tmp reference_metadata_checksums.sha256
)

{
    printf 'field\tvalue\n'
    printf 'database_name\t%s\n' "${KRAKEN_DB_NAME}"
    printf 'database_path\t%s\n' "${KRAKEN_DB_DIR}"
    for field in download_started_utc download_completed_utc build_started_utc build_completed_utc validation_started_utc; do
        value="unknown (legacy build)"
        [[ ! -s "$KRAKEN_DB_DIR/provenance/$field" ]] || value="$(cat "$KRAKEN_DB_DIR/provenance/$field")"
        printf '%s\t%s\n' "$field" "$value"
    done
    printf 'kraken2_module\t%s\n' "${KRAKEN_SOFTWARE_MODULE}"
    printf 'included_libraries\t%s\n' "${KRAKEN_LIBRARIES[*]}"
    printf 'database_hash_size_limit\tuncapped\n'
    printf 'kmer_length\t%s\n' "${KRAKEN_KMER_LENGTH}"
    printf 'minimizer_length\t%s\n' "${KRAKEN_MINIMIZER_LENGTH}"
    printf 'minimizer_spaces\t%s\n' "${KRAKEN_MINIMIZER_SPACES}"
    printf 'fast_build\tfalse\n'
    printf 'low_complexity_masking\tenabled\n'
    printf 'database_size_bytes\t%s\n' "$(du -sb "${KRAKEN_DB_DIR}" | awk '{print $1}')"
} > "${KRAKEN_DB_DIR}/database_manifest.tsv.tmp"

module list > "${KRAKEN_DB_DIR}/module_list.txt" 2>&1

printf '\nValidated major taxa:\n'
cat "${KRAKEN_DB_DIR}/expected_taxa.inspect.tsv"

printf '\nDatabase files:\n'
ls -lh \
    "${KRAKEN_DB_DIR}/hash.k2d" \
    "${KRAKEN_DB_DIR}/opts.k2d" \
    "${KRAKEN_DB_DIR}/taxo.k2d"

printf '\nChecksums:\n'
cat "${KRAKEN_DB_DIR}/checksums.sha256"

printf '\nApplying final group ownership and permissions...\n'
apply_kraken_permissions

# Only now publish successful validation. Keep the first successful date on reruns.
record_stage_time validation_completed_utc
printf 'validation_completed_utc\t%s\n' "$(cat "$KRAKEN_DB_DIR/provenance/validation_completed_utc")" >> "$KRAKEN_DB_DIR/database_manifest.tsv.tmp"
printf 'last_validation_utc\t%s\n' "$(utc_now)" >> "$KRAKEN_DB_DIR/database_manifest.tsv.tmp"
printf 'kraken2_version\t%s\n' "$(tr '\n\t' '  ' < "$KRAKEN_DB_DIR/provenance/build_kraken2_version.txt" 2>/dev/null || printf unknown)" >> "$KRAKEN_DB_DIR/database_manifest.tsv.tmp"
mv "$KRAKEN_DB_DIR/database_manifest.tsv.tmp" "$KRAKEN_DB_DIR/database_manifest.tsv"
touch "$KRAKEN_DB_DIR/.database_validation_complete"
printf '\nDatabase validation completed successfully: %s\n' "$(utc_now)"
