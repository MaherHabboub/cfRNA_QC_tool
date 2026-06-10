#!/bin/bash

set -euo pipefail

CONFIG="${1:-}"

if [[ -z "$CONFIG" ]]; then
    echo "Usage: bash submit_make_bins_only.sh path/to/config.sh"
    exit 1
fi

if [[ ! -f "$CONFIG" ]]; then
    echo "ERROR: Config file not found: $CONFIG"
    exit 1
fi

module swap cluster/doduo || true
module load env/software/doduo || true

CONFIG="$(readlink -f "$CONFIG")"

TOOL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULE_DIR="${TOOL_DIR}/modules"
WRAPPER_DIR="${TOOL_DIR}/wrappers"
LOG_DIR="${TOOL_DIR}/logs"

mkdir -p "$WRAPPER_DIR" "$LOG_DIR"

source "$CONFIG"

: "${OUTDIR:?ERROR: OUTDIR not set in config}"

mkdir -p "${OUTDIR}/logs"

WRAPPER="${WRAPPER_DIR}/02_make_dropoff_bins_only.sh"

cat > "$WRAPPER" <<EOF
#!/bin/bash

#SBATCH -J make_bins_only
#SBATCH -D ${TOOL_DIR}
#SBATCH -t 03:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=FAIL
#SBATCH -o ${LOG_DIR}/02_make_bins_only_%j.out
#SBATCH -e ${LOG_DIR}/02_make_bins_only_%j.err

set -euo pipefail

echo "============================================================"
echo "Running Make_Dropoff_Bins.sh only"
echo "Config: ${CONFIG}"
echo "Started: \$(date)"
echo "Host: \$(hostname)"
echo "============================================================"

bash "${MODULE_DIR}/Make_Dropoff_Bins.sh" "${CONFIG}"

echo "============================================================"
echo "Finished Make_Dropoff_Bins.sh"
echo "Finished: \$(date)"
echo "============================================================"
EOF

chmod +x "$WRAPPER"

echo "Submitting:"
echo "$WRAPPER"

sbatch "$WRAPPER"