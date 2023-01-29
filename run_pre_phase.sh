#!/bin/bash
# SLURM resource specifications
# (use an extra '#' in front of SBATCH to comment-out any unused options)
#SBATCH --job-name=pre_phase   # shows up in the output of 'squeue'
#SBATCH --time=1:59:59       # specify the requested wall-time
#SBATCH --partition=astro2_short  # specify the partition to run on
##SBATCH --partition=astro_devel  # specify the partition to run on
#SBATCH --nodes=1              # number of nodes allocated for this job
#SBATCH --ntasks-per-node=5    # number of MPI ranks per node
##SBATCH --cpus-per-task=1       # number of OpenMP threads per MPI rank
#SBATCH --mail-type=ALL

[ -d /lustre/astro/gfriss/pre_phase ] && rm -rf /lustre/astro/gfriss/pre_phase
mkdir /lustre/astro/gfriss/pre_phase

srun python pre_phase.py # running the python script
