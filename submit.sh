#!/bin/bash

#SBATCH -J FDS_Sonoma_Array
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=25        # 4 nodes * 25 tasks/node = 100 MPI processes
#SBATCH --mem-per-cpu=2G
#SBATCH -t 5-0
#SBATCH -A umontana_fire_modeling
#SBATCH --array=0-25%4

# --- Environment Setup ---
module load intel

# export FI_PROVIDER=verbs
# export I_MPI_FABRICS=shm:ofi

# Read the sorted list of simulation IDs from a file into a bash array
mapfile -t SIM_IDS < <(sort identifiers.txt)

# Get the identifier for the current Slurm array task
CURRENT_SIM_ID=${SIM_IDS[$SLURM_ARRAY_TASK_ID]}

# Define the simulation directory for this specific task
SIM_DIR="simulations/${CURRENT_SIM_ID}"
cd ${SIM_DIR}

# Define a unique log file for this task inside its simulation directory
LOG_FILE="fds_run.log"

echo "Starting Slurm Task ID: ${SLURM_ARRAY_TASK_ID}, Simulation: ${CURRENT_SIM_ID}" > ${LOG_FILE}

# Launch the single, parallel FDS simulation.
# srun automatically uses the resources allocated to this specific array task (--nodes=4, --ntasks=100)
srun /90daydata/umontana_fire_modeling/anthony.marcozzi/fds/Build/impi_intel_linux/fds_impi_intel_linux input.fds >> ${LOG_FILE} 2>&1

echo "Job ${CURRENT_SIM_ID} finished." >> ${LOG_FILE}