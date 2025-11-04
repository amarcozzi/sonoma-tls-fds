#!/bin/bash

#SBATCH -J FDS_Sonoma
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=25
#SBATCH --mem-per-cpu=2G
#SBATCH -t 5-0
#SBATCH -A umontana_fire_modeling
#SBATCH --partition=atlas
#SBATCH --array=0-25%1

module load intel-oneapi-compilers intel-oneapi-mpi

export FI_PROVIDER=verbs
export I_MPI_FABRICS=shm:ofi

# Read the sorted list of simulation IDs from a file into a bash array
mapfile -t SIM_IDS < <(sort identifiers.txt)

# Get the identifier for the current Slurm array task
CURRENT_SIM_ID=${SIM_IDS[$SLURM_ARRAY_TASK_ID]}

# Define the simulation directory
SIM_DIR="simulations/${CURRENT_SIM_ID}"

echo "Slurm Task ID: ${SLURM_ARRAY_TASK_ID}, Simulation: ${CURRENT_SIM_ID}"

# Navigate to the correct simulation directory and execute
cd ${SIM_DIR} && srun ~/fds/Build/fds_impi_intel_linux input.fds > fds_run.log 2>&1

echo "Job ${CURRENT_SIM_ID} finished."