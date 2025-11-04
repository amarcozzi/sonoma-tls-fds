#!/bin/bash

#SBATCH -J FDS_Sonoma_Array
#SBATCH --nodes=16                  # Total nodes for the array (4 simulations * 4 nodes each)
#SBATCH --ntasks-per-node=25        # Tasks per node for each simulation
#SBATCH --mem-per-cpu=2G
#SBATCH -t 5-0
#SBATCH -A umontana_fire_modeling
#SBATCH --partition=atlas
#SBATCH --array=0-25%4              # Run 26 jobs (0-25), with a maximum of 4 running at once

# --- Environment Setup ---
module load intel-oneapi-compilers intel-oneapi-mkl

source /project/umontana_fire_modeling/anthony.marcozzi/miniforge3/etc/profile.d/conda.sh
conda activate fds

export FI_PROVIDER=verbs
export I_MPI_FABRICS=shm:ofi

# --- Job Array Logic ---

# Read the sorted list of simulation IDs from a file into a bash array
mapfile -t SIM_IDS < <(sort identifiers.txt)

# Get the identifier for the current Slurm array task
CURRENT_SIM_ID=${SIM_IDS[$SLURM_ARRAY_TASK_ID]}

# Define the simulation directory for this specific task
SIM_DIR="simulations/${CURRENT_SIM_ID}"

# Define a unique log file for this task inside its simulation directory
LOG_FILE="${SIM_DIR}/fds_run_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"

echo "========================================================" > ${LOG_FILE}
echo "Slurm Job Array ID:  $SLURM_ARRAY_JOB_ID" >> ${LOG_FILE}
echo "Slurm Task ID:       $SLURM_ARRAY_TASK_ID" >> ${LOG_FILE}
echo "Running Simulation:  ${CURRENT_SIM_ID}" >> ${LOG_FILE}
echo "Working Directory:   $(pwd)/${SIM_DIR}" >> ${LOG_FILE}
echo "Timestamp:           $(date)" >> ${LOG_FILE}
echo "========================================================" >> ${LOG_FILE}

# --- Execution ---
# Navigate to the correct simulation directory
cd ${SIM_DIR}

# Launch a single, multi-node FDS simulation as a job step.
# This srun command carves out the necessary resources (4 nodes, 100 tasks)
# from the total 16-node allocation provided to the parent job array.
# The '--exclusive' flag is crucial to ensure tasks don't share nodes.
srun --exclusive \
     --nodes=4 \
     --ntasks=100 \
     --ntasks-per-node=25 \
     ~/fds/Build/fds_impi_intel_linux input.fds >> ${LOG_FILE} 2>&1

echo "Job ${CURRENT_SIM_ID} finished." >> ${LOG_FILE}