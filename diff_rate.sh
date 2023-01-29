#!/bin/bash
# SLURM resource specifications
# (use an extra '#' in front of SBATCH to comment-out any unused options)
#SBATCH --job-name=diff_rates   # shows up in the output of 'squeue'
#SBATCH --time=3-23:59:59       # specify the requested wall-time
#SBATCH --partition=astro_long  # specify the partition to run on
#SBATCH --nodes=1              # number of nodes allocated for this job
#SBATCH --ntasks-per-node=20    # number of MPI ranks per node
##SBATCH --cpus-per-task=1       # number of OpenMP threads per MPI rank
#SBATCH --mail-type=ALL

module load astro
module load python/anaconda3

##[ -d /lustre/astro/gfriss/diff_rate ] && rm -rf /lustre/astro/gfriss/diff_rate
##mkdir /lustre/astro/gfriss/diff_rate

##srun python diff_rate_sim.py # running the python script

srun python uncertainty.py