#!/usr/bin/env bash
#SBATCH --job-name=move_kraken
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=1-00:00:00
#SBATCH --output=/scratch/gent/vo/000/gvo00027/projects/MHB/move_kraken_%j.out
#SBATCH --error=/scratch/gent/vo/000/gvo00027/projects/MHB/move_kraken_%j.err

set -euo pipefail

SOURCE="/data/gent/vo/000/gvo00027/tools/kraken"
DESTINATION="/scratch/gent/vo/000/gvo00027/projects/MHB/kraken"
SHARED_GROUP="gvo00027"

if [[ ! -d "${SOURCE}" ]]
then
    printf 'ERROR: source directory does not exist: %s\n' "${SOURCE}" >&2
    exit 1
fi

mkdir -p "${DESTINATION}"

printf 'Copy started: %s\n' "$(date --iso-8601=seconds)"
printf 'Source: %s\n' "${SOURCE}"
printf 'Destination: %s\n' "${DESTINATION}"

rsync \
    -aH \
    --partial \
    --info=progress2 \
    "${SOURCE}/" \
    "${DESTINATION}/"

printf '\nChecking whether anything still differs...\n'

verification_file="$(mktemp)"

rsync \
    -aHn \
    --delete \
    --itemize-changes \
    "${SOURCE}/" \
    "${DESTINATION}/" \
    > "${verification_file}"

if [[ -s "${verification_file}" ]]
then
    printf 'ERROR: source and destination still differ:\n' >&2
    cat "${verification_file}" >&2
    exit 1
fi

rm -f "${verification_file}"

printf '\nApplying group ownership and permissions...\n'

chgrp -R "${SHARED_GROUP}" "${DESTINATION}"
find "${DESTINATION}" -type d -exec chmod 2755 {} +
find "${DESTINATION}" -type f -exec chmod 0755 {} +

printf '\nSource size:\n'
du -sh "${SOURCE}"

printf 'Destination size:\n'
du -sh "${DESTINATION}"

printf '\nCopy and verification completed: %s\n' \
    "$(date --iso-8601=seconds)"

printf 'The original directory has intentionally not been deleted.\n'