#!/bin/bash
# script to run 1 simulation, launched from /lustre/hpc/astro/gfriss
# folder with parameter files is here and given as input
#SBATCH --job-name=pnautilus_test    	# shows up in the output of 'squeue'
#SBATCH --partition=astro_devel   	   	# specify the partition to run on
#SBATCH --time=0:04:59       # specify the requested wall-time
#SBATCH --nodes=1               # number of nodes allocated for this job
##SBATCH --ntasks-per-node=20    # number of MPI ranks per node

module load astro
module load python/anaconda3

cd /lustre/astro/gfriss/test

/lustre/astro/gfriss/pnautilus/pnautilus
/lustre/astro/gfriss/pnautilus/pnautilus_outputs
