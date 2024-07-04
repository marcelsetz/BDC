#!/bin/bash
#SBATCH --job-name=assignment4
#SBATCH --time=0-01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --partition=assemblix
#SBATCH --nodelist=assemblix2019

# Load necessary modules
module load anaconda
module load mpi

# Activate your conda environment
source activate /../.conda

# Run a simple test to check MPI execution
mpiexec -n 4 python -c "from mpi4py import MPI; print('Hello from rank', MPI.COMM_WORLD.Get_rank())"