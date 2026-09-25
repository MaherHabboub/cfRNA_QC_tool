#!/usr/bin/env bash
#SBATCH --job-name=k2_db_download
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=3-00:00:00

set -euo pipefail

# Slurm executes a spooled copy: submit_all passes the original bundle directory.
SCRIPT_DIR="${2:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "${SCRIPT_DIR}/environment.sh"
load_kraken_config "${1:-${SCRIPT_DIR}/config.sh}"

umask 0022

load_kraken_modules

mkdir -p "${KRAKEN_DB_DIR}"

download_bacteria_https() {
    local bacteria_dir="${KRAKEN_DB_DIR}/library/bacteria"
    local download_dir="${bacteria_dir}/https_downloads"
    local assembly_summary="${bacteria_dir}/assembly_summary.txt"
    local manifest="${bacteria_dir}/https_manifest.tsv"
    local manifest_marker="${bacteria_dir}/.https_manifest_complete"
    local downloads_marker="${bacteria_dir}/.https_downloads_complete"
    local assembled_marker="${bacteria_dir}/.https_library_assembled"
    local prepared_marker="${bacteria_dir}/.https_library_prepared"
    local library_file="${bacteria_dir}/library.fna"
    local library_tmp="${bacteria_dir}/library.fna.tmp"
    local worker_count="${SLURM_CPUS_PER_TASK:-16}"

    # Limit simultaneous NCBI transfers even when more CPUs were requested.
    if (( worker_count > 8 ))
    then
        worker_count=8
    fi

    mkdir -p "${bacteria_dir}" "${download_dir}"

    if [[ ! -f "${manifest_marker}" ]]
    then
        printf '\nDownloading the bacterial assembly summary over HTTPS...\n'
        wget \
            --https-only \
            --timeout=120 \
            --tries=20 \
            --retry-connrefused \
            --waitretry=10 \
            -O "${assembly_summary}.tmp" \
            "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt"

        if [[ ! -s "${assembly_summary}.tmp" ]]
        then
            printf 'ERROR: downloaded bacterial assembly summary is empty.\n' >&2
            exit 1
        fi

        mv "${assembly_summary}.tmp" "${assembly_summary}"

        # Match Kraken2's normal selection: RefSeq entries whose assembly level
        # is Complete Genome or Chromosome and that have a downloadable path.
        {
            printf 'taxid\turl\tfilename\n'
            awk -F '\t' '
                BEGIN { OFS="\t" }
                /^#/ { next }
                ($12 == "Complete Genome" || $12 == "Chromosome") && $20 != "na" {
                    n = split($20, path_parts, "/")
                    base = path_parts[n]
                    url = $20 "/" base "_genomic.fna.gz"
                    sub(/^ftp:/, "https:", url)
                    print $6, url, base "_genomic.fna.gz"
                }
            ' "${assembly_summary}" | LC_ALL=C sort -t $'\t' -k3,3
        } > "${manifest}.tmp"

        if [[ "$(wc -l < "${manifest}.tmp")" -le 1 ]]
        then
            printf 'ERROR: no bacterial assemblies were selected.\n' >&2
            exit 1
        fi

        mv "${manifest}.tmp" "${manifest}"
        touch "${manifest_marker}"
    else
        printf '\nExisting frozen bacterial HTTPS manifest found; reusing it.\n'
    fi

    printf 'Bacterial assemblies in manifest: %s\n' "$(( $(wc -l < "${manifest}") - 1 ))"
    printf 'Parallel HTTPS downloads: %s\n' "${worker_count}"

    # Salvage any complete files left by the previous FTP attempts.
    if [[ -d "${bacteria_dir}/all" ]]
    then
        printf '\nChecking previous FTP files for reusable complete downloads...\n'
        while IFS= read -r -d '' previous_file
        do
            if gzip -t "${previous_file}" 2>/dev/null
            then
                mv -n "${previous_file}" "${download_dir}/"
            fi
        done < <(
            find "${bacteria_dir}/all" \
                -maxdepth 1 \
                -type f \
                -name '*_genomic.fna.gz' \
                -print0
        )
    fi

    if [[ ! -f "${downloads_marker}" ]]
    then
        printf '\nDownloading missing bacterial genomes over restartable HTTPS...\n'

        export KRAKEN_BACTERIA_HTTPS_DIR="${download_dir}"

        download_one_bacterial_genome() {
            local url="$1"
            local filename="$2"
            local destination="${KRAKEN_BACTERIA_HTTPS_DIR}/${filename}"
            local attempt

            if [[ -s "${destination}" ]] && gzip -t "${destination}" 2>/dev/null
            then
                printf 'Already complete: %s\n' "${filename}"
                return 0
            fi

            for attempt in $(seq 1 12)
            do
                printf 'Downloading %s (attempt %s/12)\n' "${filename}" "${attempt}"

                if wget \
                    --continue \
                    --https-only \
                    --connect-timeout=30 \
                    --read-timeout=120 \
                    --timeout=120 \
                    --tries=3 \
                    --retry-connrefused \
                    --retry-on-http-error=429,500,502,503,504 \
                    --waitretry=10 \
                    --no-verbose \
                    --directory-prefix="${KRAKEN_BACTERIA_HTTPS_DIR}" \
                    "${url}"
                then
                    if gzip -t "${destination}" 2>/dev/null
                    then
                        printf 'Completed: %s\n' "${filename}"
                        return 0
                    fi

                    printf 'Downloaded file failed gzip validation; restarting it: %s\n' \
                        "${filename}" >&2
                    rm -f "${destination}"
                fi

                sleep "$(( attempt * 15 ))"
            done

            printf 'ERROR: failed to download a valid copy of %s\n' "${filename}" >&2
            return 1
        }

        export -f download_one_bacterial_genome

        tail -n +2 "${manifest}" \
            | cut -f2,3 \
            | xargs -r -n2 -P "${worker_count}" \
                bash -c 'download_one_bacterial_genome "$1" "$2"' _

        printf '\nValidating all downloaded bacterial genome archives...\n'
        while IFS=$'\t' read -r taxid url filename
        do
            if [[ ! -s "${download_dir}/${filename}" ]] \
                || ! gzip -t "${download_dir}/${filename}" 2>/dev/null
            then
                printf 'ERROR: missing or invalid bacterial genome: %s\n' \
                    "${filename}" >&2
                exit 1
            fi
        done < <(tail -n +2 "${manifest}")

        touch "${downloads_marker}"
    else
        printf '\nAll bacterial HTTPS downloads are already marked complete.\n'
    fi

    if [[ ! -f "${assembled_marker}" ]]
    then
        printf '\nCreating Kraken2 bacterial library.fna...\n'
        rm -f "${library_tmp}"

        while IFS=$'\t' read -r taxid url filename
        do
            gzip -cd "${download_dir}/${filename}" \
                | awk -v taxid="${taxid}" '
                    /^>/ { sub(/^>/, ">kraken:taxid|" taxid "|") }
                    { print }
                ' >> "${library_tmp}"
        done < <(tail -n +2 "${manifest}")

        if [[ ! -s "${library_tmp}" ]]
        then
            printf 'ERROR: assembled bacterial library is empty.\n' >&2
            exit 1
        fi

        mv "${library_tmp}" "${library_file}"
        rm -f "${library_file}.masked"

        touch "${assembled_marker}"
    else
        printf '\nKraken2 bacterial library.fna is already assembled.\n'
    fi

    if [[ ! -s "${library_file}" ]]
    then
        printf 'ERROR: bacterial library assembly marker exists but library.fna is missing.\n' >&2
        exit 1
    fi

    if [[ ! -f "${prepared_marker}" ]]
    then

        kraken_bin_dir="$(dirname "$(readlink -f "$(command -v kraken2-build)")")"

        printf 'Creating bacterial preliminary sequence map...\n'
        "${kraken_bin_dir}/scan_fasta_file.pl" \
            "${library_file}" \
            > "${bacteria_dir}/prelim_map.txt"

        printf 'Masking low-complexity bacterial sequence...\n'
        "${kraken_bin_dir}/mask_low_complexity.sh" "${bacteria_dir}"

        touch "${prepared_marker}"
    else
        printf '\nKraken2 bacterial library is already prepared.\n'
    fi

    if [[ ! -s "${library_file}" ]] \
        || [[ ! -s "${bacteria_dir}/prelim_map.txt" ]] \
        || [[ ! -f "${library_file}.masked" ]]
    then
        printf 'ERROR: bacterial Kraken2 library preparation is incomplete.\n' >&2
        exit 1
    fi

    printf '\nBacterial HTTPS library completed successfully.\n'
    du -sh "${bacteria_dir}"
}

record_stage_time download_started_utc
if [[ ! -s "${KRAKEN_DB_DIR}/provenance/download_completed_utc" ]]; then
    record_stage_environment download
fi

printf 'Kraken2 database download started: %s\n' "$(utc_now)"
printf 'Host: %s\n' "$(hostname)"
printf 'Database: %s\n' "${KRAKEN_DB_DIR}"
printf 'Threads: %s\n' "${SLURM_CPUS_PER_TASK:-16}"
df -h "${KRAKEN_ROOT}"

if [[ ! -f "${KRAKEN_DB_DIR}/.taxonomy_download_complete" ]]
then
    printf '\nDownloading NCBI taxonomy...\n'
    kraken2-build \
        --download-taxonomy \
        --skip-maps \
        --use-ftp \
        --db "${KRAKEN_DB_DIR}"
    touch "${KRAKEN_DB_DIR}/.taxonomy_download_complete"
else
    printf '\nNCBI taxonomy already completed; skipping.\n'
fi

for ref_library in "${KRAKEN_LIBRARIES[@]}"
do
    completion_marker="${KRAKEN_DB_DIR}/.library_${ref_library}_download_complete"

    if [[ -f "${completion_marker}" ]]
    then
        printf '\nLibrary already completed; skipping: %s\n' "${ref_library}"
        continue
    fi

    printf '\nDownloading Kraken2 library: %s\n' "${ref_library}"

    if [[ "${ref_library}" == "bacteria" ]]
    then
        download_bacteria_https
    else
        # These libraries are much smaller and previously completed on Doduo
        # using Kraken2's FTP mode. rsync is blocked from the compute nodes.
        kraken2-build \
            --download-library "${ref_library}" \
            --threads "${SLURM_CPUS_PER_TASK:-16}" \
            --use-ftp \
            --db "${KRAKEN_DB_DIR}"
    fi

    touch "${completion_marker}"
done

printf '\nAll requested libraries downloaded successfully.\n'
printf 'Download completed: %s\n' "$(utc_now)"
du -sh "${KRAKEN_DB_DIR}"
df -h "${KRAKEN_ROOT}"

record_stage_time download_completed_utc
