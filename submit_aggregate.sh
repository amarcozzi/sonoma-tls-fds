#!/bin/bash

#SBATCH -J FDS_Sonoma_Agg
#SBATCH --nodes=1
#SBATCH -n 1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH -t 0:30:00
#SBATCH -A umontana_fire_modeling
#SBATCH -o logs_postprocess/aggregate.log
#SBATCH -e logs_postprocess/aggregate.log

# Merges every per-simulation metrics.json written by submit_postprocess.sh
# into the cross-simulation tables and figures. Chained after the array job
# with --dependency=afterany, so it still runs (and reports the failures) if
# some simulations error out.

set -uo pipefail

module purge
module load uv

cd "${SLURM_SUBMIT_DIR}"

OUT_DIR="${OUT_DIR:-postprocess}"
WINDOW="${WINDOW:-roi}"

export MPLBACKEND=Agg
export MPLCONFIGDIR="${TMPDIR:-/tmp}/mpl-${SLURM_JOB_ID}"
export UV_CACHE_DIR="${UV_CACHE_DIR:-${SLURM_SUBMIT_DIR}/.uv-cache}"
mkdir -p "${MPLCONFIGDIR}"

uv run --no-sync python scripts/aggregate_metrics.py \
    --out "${OUT_DIR}" --window "${WINDOW}"
STATUS=$?

# the other windows cost nothing and save a re-run when 'roi' turns out to be
# too short a sample for the upper percentiles
for w in roi steady full; do
    if [ "$w" != "${WINDOW}" ]; then
        uv run --no-sync python scripts/aggregate_metrics.py \
            --out "${OUT_DIR}" --window "$w" --no-figures > /dev/null 2>&1
        cp "${OUT_DIR}/metrics_summary.csv" "${OUT_DIR}/metrics_summary_${w}.csv"
    fi
done
# leave the requested window as the canonical metrics_summary.csv
uv run --no-sync python scripts/aggregate_metrics.py \
    --out "${OUT_DIR}" --window "${WINDOW}" --no-figures > /dev/null 2>&1

echo "Aggregation finished with status ${STATUS}"
echo "Tables and figures: ${OUT_DIR}"
ls -la "${OUT_DIR}"/*.csv "${OUT_DIR}"/*.png 2>/dev/null

exit ${STATUS}
