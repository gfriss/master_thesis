#!/bin/bash
# SLURM resource specifications
# (use an extra '#' in front of SBATCH to comment-out any unused options)
#SBATCH --job-name=test   # shows up in the output of 'squeue'
#SBATCH --time=1:59:59       # specify the requested wall-time
#SBATCH --partition=astro_devel  # specify the partition to run on
#SBATCH --nodes=1              # number of nodes allocated for this job
#SBATCH --ntasks-per-node=3    # number of MPI ranks per node
##SBATCH --cpus-per-task=1       # number of OpenMP threads per MPI rank
# nodes*ntasks-per-node = number of CPUs you want to use
##SBATCH --mail-type=ALL

module load astro
module load python/anaconda3

rm -rf /lustre/astro/gfriss/test_many   # this is only for testing purposes
mkdir /lustre/astro/gfriss/test_many    # to have a fresh folder each try

srun python test_many_sim.py # running the python script
