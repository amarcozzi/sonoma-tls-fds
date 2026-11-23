#!/bin/bash

#SBATCH -J FDS_Sonoma_Array
#SBATCH --nodes=1
#SBATCH -n 1
#SBATCH --cpus-per-task=80        # 80 meshes -> 80 MPI ranks (pure-MPI build, 1 core/rank)
#SBATCH --hint=nomultithread      # 1 rank per PHYSICAL core; needs >48 cores, so lands on ceres24/25 (128-core, 2.2T)
#SBATCH --mem-per-cpu=4G          # 80 x 4G = 320G (peak use was 69-154 GB on the 64-mesh domain)
#SBATCH -t 5-0
#SBATCH -A umontana_fire_modeling
#SBATCH --array=0-5               # size this to the case list on the command line, e.g. --array=0-20
#SBATCH -o logs/%A_%a.log
#SBATCH -e logs/%A_%a.log

# Usage:
#   sbatch --array=0-$(( $(wc -l < identifiers.txt) - 1 )) submit.sh
#   sbatch --array=0-5 --export=ALL,IDENT_LIST=batch1.txt submit.sh
# IDENT_LIST is read when each task starts, so give every batch its own list
# file instead of rewriting identifiers.txt while tasks are still queued.

FDS_BIN=/project/umontana_fire_modeling/anthony.marcozzi/fds/Build/impi_intel_linux/fds_impi_intel_linux_6.11.2-0-g8e4c2c4a47

# --- Environment Setup ---
module purge
module load intel

# Read the sorted list of simulation IDs into a bash array
mapfile -t SIM_IDS < <(grep -v '^[[:space:]]*$' "${IDENT_LIST:-identifiers.txt}" | sort)

# Get the identifier for the current Slurm array task
CURRENT_SIM_ID=${SIM_IDS[$SLURM_ARRAY_TASK_ID]}
if [ -z "${CURRENT_SIM_ID}" ]; then
    echo "ERROR: no simulation at index ${SLURM_ARRAY_TASK_ID} in ${IDENT_LIST:-identifiers.txt}" >&2
    exit 1
fi

# Define the simulation directory for this specific task
SIM_DIR="simulations/${CURRENT_SIM_ID}"
echo "Starting Slurm Task ID: ${SLURM_ARRAY_TASK_ID}, Simulation: ${CURRENT_SIM_ID}"

# Change into the simulation directory to run the job there.
cd "${SIM_DIR}" || exit 1

# Launch the FDS simulation: one rank per mesh (80), pinned to physical cores.
mpirun -n 80 "${FDS_BIN}" input.fds

echo "Job ${CURRENT_SIM_ID} finished."
