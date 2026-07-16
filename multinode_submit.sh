#!/bin/bash

#SBATCH -J FDS_Sonoma_Array
#SBATCH --nodes=3                 # spread each sim's 81 meshes across 3 nodes
#SBATCH --ntasks=81              # 1 MPI rank per mesh (9x9 MULT grid = 81 meshes)
#SBATCH --ntasks-per-node=27     # 81 ranks / 3 nodes = 27 ranks per node (even split)
#SBATCH --cpus-per-task=1        # FDS non-OpenMP build: 1 core per rank
#SBATCH --mem-per-cpu=2G         # 27 ranks x 2G = 54G per node
#SBATCH -t 5-0
#SBATCH -A umontana_fire_modeling
#SBATCH --array=0-3%2            # 4 sims, 2 concurrent -> up to 6 nodes at once
#SBATCH -o logs/%a.log
#SBATCH -e logs/%a.log

# --- Environment Setup ---
module purge
module load intel

# Make Intel MPI's Hydra launcher place ranks according to the SLURM
# allocation, so meshes are distributed across all allocated nodes.
export I_MPI_HYDRA_BOOTSTRAP=slurm

# Read the sorted list of simulation IDs from a file into a bash array
mapfile -t SIM_IDS < <(sort identifiers.txt)

# Get the identifier for the current Slurm array task
CURRENT_SIM_ID=${SIM_IDS[$SLURM_ARRAY_TASK_ID]}

# Define the simulation directory for this specific task
SIM_DIR="simulations/${CURRENT_SIM_ID}"
echo "Starting Slurm Task ID: ${SLURM_ARRAY_TASK_ID}, Simulation: ${CURRENT_SIM_ID}"
echo "Allocation: ${SLURM_NTASKS} ranks across ${SLURM_NNODES} nodes (${SLURM_JOB_NODELIST})"

# Change into the simulation directory to run the job there.
cd ${SIM_DIR}

# Launch FDS: one rank per mesh, distributed across all allocated nodes.
mpirun -n ${SLURM_NTASKS} /project/umontana_fire_modeling/anthony.marcozzi/fds/Build/impi_intel_linux/fds_impi_intel_linux input.fds

echo "Job ${CURRENT_SIM_ID} finished."
