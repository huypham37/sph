#!/bin/bash
#PBS -N sph_simulation
#PBS -l nodes=1:ppn=32
#PBS -l walltime=01:00:00
# Change the above to match your cluster's job scheduler syntax
# Common alternatives:
# SLURM: #SBATCH --ntasks=1 --cpus-per-task=32 --time=01:00:00
# LSF: #BSUB -n 32 -W 60

# Setup environment
cd ${PBS_O_WORKDIR:-$(pwd)}
source ./setup_env.sh

# Run simulation with performance tracking
cd "/home/coder/project/sph/build"
echo "Starting SPH simulation with ${OMP_NUM_THREADS} threads at $(date)"
time ./sph_simulation

echo "Simulation completed at $(date)"
