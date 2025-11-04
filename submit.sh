#!/bin/bash

#SBATCH -J FDS_Sonoma_Array
#SBATCH --nodes=1
#SBATCH -n 1
#SBATCH --cpus-per-task=81
#SBATCH --mem-per-cpu=2G
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

# Launch the FDS simulation.
mpirun -n 81 /90daydata/umontana_fire_modeling/anthony.marcozzi/fds/Build/impi_intel_linux/fds_impi_intel_linux input.fds

echo "Job ${CURRENT_SIM_ID} finished."