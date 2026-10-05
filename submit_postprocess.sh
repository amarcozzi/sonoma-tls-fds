#!/bin/bash

#SBATCH -J FDS_Sonoma_Post
#SBATCH --nodes=1
#SBATCH -n 1
#SBATCH --cpus-per-task=4         # numpy/scipy only; the work is I/O bound
#SBATCH --mem=32G                 # peak is one assembled HRRPUV slice (~250 MB)
#SBATCH -t 4:00:00
#SBATCH -A umontana_fire_modeling
#SBATCH --array=0-5%10            # overridden by run_postprocess.sh
#SBATCH -o logs_postprocess/%a.log
#SBATCH -e logs_postprocess/%a.log

# One array task per simulation. Reads the raw FDS output in place and writes
# a few MB of CSV/PNG per run into $OUT_DIR, so nothing large leaves Ceres.
# Submit via ./run_postprocess.sh (which sizes the array from identifiers.txt
# and chains the aggregation job).

set -uo pipefail

module purge
module load uv

cd "${SLURM_SUBMIT_DIR}"

SIM_ROOT="${SIM_ROOT:-simulations}"
OUT_DIR="${OUT_DIR:-postprocess}"

# matplotlib and uv both want writable caches; keep them off the home quota
export MPLBACKEND=Agg
export MPLCONFIGDIR="${TMPDIR:-/tmp}/mpl-${SLURM_JOB_ID}"
export UV_CACHE_DIR="${UV_CACHE_DIR:-${SLURM_SUBMIT_DIR}/.uv-cache}"
mkdir -p "${MPLCONFIGDIR}"

mapfile -t SIM_IDS < <(sort "${IDENT_LIST:-identifiers.txt}")
CURRENT_SIM_ID=${SIM_IDS[$SLURM_ARRAY_TASK_ID]}
SIM_DIR="${SIM_ROOT}/${CURRENT_SIM_ID}"

echo "Task ${SLURM_ARRAY_TASK_ID}: postprocessing ${CURRENT_SIM_ID} from ${SIM_DIR}"

if [ ! -d "${SIM_DIR}" ]; then
    echo "ERROR: ${SIM_DIR} does not exist" >&2
    exit 1
fi
if ! compgen -G "${SIM_DIR}/*.smv" > /dev/null; then
    echo "ERROR: no .smv in ${SIM_DIR} — FDS never produced output here" >&2
    exit 1
fi

# --no-sync: the environment was resolved once on the login node by
# run_postprocess.sh, so array tasks must not race to rebuild .venv
uv run --no-sync python scripts/postprocess.py "${SIM_DIR}" --out "${OUT_DIR}"
STATUS=$?

echo "Task ${SLURM_ARRAY_TASK_ID} (${CURRENT_SIM_ID}) finished with status ${STATUS}"
exit ${STATUS}
