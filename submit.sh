#!/bin/bash

#SBATCH -J FDS_Sonoma_Array
#SBATCH --nodes=1
#SBATCH -n 1
#SBATCH --cpus-per-task=64        # 64 meshes -> 64 MPI ranks (pure-MPI build, 1 core/rank)
#SBATCH --hint=nomultithread      # 1 rank per PHYSICAL core; needs >48 cores, so lands on ceres24/25 (128-core, 2.2T)
#SBATCH --mem-per-cpu=4G          # 64 x 4G = 256G
#SBATCH -t 5-0
#SBATCH -A umontana_fire_modeling
#SBATCH --array=0-3%4
#SBATCH -o logs/%a.log
#SBATCH -e logs/%a.log

# --- Environment Setup ---
module purge
module load intel

# Read the sorted list of simulation IDs from a file into a bash array
mapfile -t SIM_IDS < <(sort identifiers.txt)

# Get the identifier for the current Slurm array task
CURRENT_SIM_ID=${SIM_IDS[$SLURM_ARRAY_TASK_ID]}

# Define the simulation directory for this specific task
SIM_DIR="simulations/${CURRENT_SIM_ID}"
echo "Starting Slurm Task ID: ${SLURM_ARRAY_TASK_ID}, Simulation: ${CURRENT_SIM_ID}"

# Change into the simulation directory to run the job there.
cd ${SIM_DIR}

# Launch the FDS simulation: one rank per mesh (64), pinned to physical cores.
mpirun -n 64 /project/umontana_fire_modeling/anthony.marcozzi/fds/Build/impi_intel_linux/fds_impi_intel_linux input.fds

echo "Job ${CURRENT_SIM_ID} finished."