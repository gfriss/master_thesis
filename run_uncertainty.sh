#!/bin/bash
# SLURM resource specifications
# (use an extra '#' in front of SBATCH to comment-out any unused options)
#SBATCH --job-name=uncer   # shows up in the output of 'squeue'
#SBATCH --time=11:59:59       # specify the requested wall-time
#SBATCH --partition=astro2_short  # specify the partition to run on
#SBATCH --nodes=1              # number of nodes allocated for this job
#SBATCH --ntasks-per-node=20    # number of MPI ranks per node
##SBATCH --cpus-per-task=1       # number of OpenMP threads per MPI rank
# nodes*ntasks-per-node = number of CPUs you want to use
#SBATCH --mail-type=ALL

srun python uncertainty.py
