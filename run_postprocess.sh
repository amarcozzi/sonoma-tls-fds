#!/bin/bash
#
# Submit the whole postprocessing pipeline on Ceres:
#   1. resolve the Python environment ONCE on the login node (array tasks
#      must not race to build .venv on a shared filesystem)
#   2. one array task per simulation in identifiers.txt
#   3. an aggregation job chained with --dependency=afterany
#
# Usage:
#   ./run_postprocess.sh                  # all of identifiers.txt
#   ./run_postprocess.sh c4_p10 c8_p25    # just these
#   WINDOW=steady ./run_postprocess.sh    # summary table over the steady window
#
# Environment overrides: SIM_ROOT (default simulations), OUT_DIR (default
# postprocess), WINDOW (default roi), CONCURRENT (default 10).

set -euo pipefail

cd "$(dirname "$0")"
REPO="$(pwd)"

SIM_ROOT="${SIM_ROOT:-simulations}"
OUT_DIR="${OUT_DIR:-postprocess}"
WINDOW="${WINDOW:-roi}"
CONCURRENT="${CONCURRENT:-10}"
export SIM_ROOT OUT_DIR WINDOW
export UV_CACHE_DIR="${UV_CACHE_DIR:-${REPO}/.uv-cache}"

module purge
module load uv

# If simulation ids were given, run only those via a scratch identifier list.
LIST=identifiers.txt
if [ "$#" -gt 0 ]; then
    LIST="$(mktemp "${REPO}/.identifiers_subset.XXXXXX")"
    printf '%s\n' "$@" | sort > "${LIST}"
    echo "Restricting to $# simulation(s) from ${LIST}"
fi

N=$(grep -cve '^[[:space:]]*$' "${LIST}")
if [ "${N}" -eq 0 ]; then
    echo "ERROR: ${LIST} lists no simulations" >&2
    exit 1
fi

mkdir -p logs_postprocess "${OUT_DIR}"

echo "Resolving the Python environment (uv sync)..."
uv sync
uv run --no-sync python -c "import fdsreader, numpy, scipy, matplotlib; \
print(f'fdsreader {fdsreader.__version__}, numpy {numpy.__version__}, \
scipy {scipy.__version__}, matplotlib {matplotlib.__version__}')"

echo "Submitting ${N} postprocessing task(s), ${CONCURRENT} concurrent..."
ARRAY_JOB=$(sbatch --parsable \
    --array="0-$((N - 1))%${CONCURRENT}" \
    --export=ALL,SIM_ROOT="${SIM_ROOT}",OUT_DIR="${OUT_DIR}",IDENT_LIST="${LIST}" \
    submit_postprocess.sh)
echo "  array job ${ARRAY_JOB}"

AGG_JOB=$(sbatch --parsable \
    --dependency=afterany:"${ARRAY_JOB}" \
    --export=ALL,OUT_DIR="${OUT_DIR}",WINDOW="${WINDOW}" \
    submit_aggregate.sh)
echo "  aggregate job ${AGG_JOB} (runs after the array, pass or fail)"

cat <<EOF

Watch:    squeue -u \$USER
Logs:     logs_postprocess/
Results:  ${OUT_DIR}/metrics_summary.csv
          ${OUT_DIR}/metrics_long.csv
          ${OUT_DIR}/metrics_dictionary.csv
          ${OUT_DIR}/comparison_*.png

To pull only the tables and figures back to your laptop:
  rsync -av --include='*/' --include='*.csv' --include='*.png' --exclude='*' \\
        ceres:${REPO}/${OUT_DIR}/ ./${OUT_DIR}/
EOF
