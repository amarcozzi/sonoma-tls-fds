#!/bin/bash

#SBATCH -J FDS_Sonoma_Array
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=25        # 4 nodes * 25 tasks/node = 100 MPI processes
#SBATCH --mem-per-cpu=2G
#SBATCH -t 5-0
#SBATCH -A umontana_fire_modeling
#SBATCH --partition=atlas
#SBATCH --array=0-25%4

# --- Environment Setup ---
module load intel-oneapi-compilers intel-oneapi-mkl

source /project/umontana_fire_modeling/anthony.marcozzi/miniforge3/etc/profile.d/conda.sh
conda activate fds

export FI_PROVIDER=verbs
export I_MPI_FABRICS=shm:ofi


# Read the sorted list of simulation IDs from a file into a bash array
mapfile -t SIM_IDS < <(sort identifiers.txt)

# Get the identifier for the current Slurm array task
CURRENT_SIM_ID=${SIM_IDS[$SLURM_ARRAY_TASK_ID]}

# Define the simulation directory for this specific task
SIM_DIR="simulations/${CURRENT_SIM_ID}"

# Define a unique log file for this task inside its simulation directory
LOG_FILE="${SIM_DIR}/fds_run_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"

echo "Starting Slurm Task ID: ${SLURM_ARRAY_TASK_ID}, Simulation: ${CURRENT_SIM_ID}" > ${LOG_FILE}

# Navigate to the correct simulation directory
cd ${SIM_DIR}

# Launch the single, parallel FDS simulation.
# srun automatically uses the resources allocated to this specific array task (--nodes=4, --ntasks=100)
srun ~/fds/Build/fds_impi_intel_linux input.fds >> ${LOG_FILE} 2>&1

echo "Job ${CURRENT_SIM_ID} finished." >> ${LOG_FILE}