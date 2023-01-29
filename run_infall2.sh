#!/bin/bash
#SBATCH --job-name=infall    	# shows up in the output of 'squeue'
#SBATCH --partition=astro_short   	   	# specify the partition to run on
#SBATCH --time=0-23:59:59       # specify the requested wall-time
#SBATCH --nodes=1               # number of nodes allocated for this job
##SBATCH --ntasks-per-node=20    # number of MPI ranks per node

#cd /lustre/astro/gfriss/infall_own
#/lustre/astro/gfriss/pnautilus/scripts/pnautilus-clean.sh
#/lustre/astro/gfriss/pnautilus/pnautilus
#/lustre/astro/gfriss/pnautilus/pnautilus_outputs
#rm abundances.*.out

cd /lustre/astro/gfriss/infall_own_with_pre
/lustre/astro/gfriss/pnautilus/scripts/pnautilus-clean.sh
/lustre/astro/gfriss/pnautilus/pnautilus
/lustre/astro/gfriss/pnautilus/pnautilus_outputs
/lustre/astro/gfriss/pnautilus/pnautilus_rates
rm abundances.*.out
rm rates.*.out

#cd /lustre/astro/gfriss/noburst_infall_own
#/lustre/astro/gfriss/pnautilus/scripts/pnautilus-clean.sh
#/lustre/astro/gfriss/pnautilus/pnautilus
#/lustre/astro/gfriss/pnautilus/pnautilus_outputs
#rm abundances.*.out

cd /lustre/astro/gfriss/noburst_infall_own_with_pre
/lustre/astro/gfriss/pnautilus/scripts/pnautilus-clean.sh
/lustre/astro/gfriss/pnautilus/pnautilus
/lustre/astro/gfriss/pnautilus/pnautilus_outputs
/lustre/astro/gfriss/pnautilus/pnautilus_rates
rm abundances.*.out
rm rates.*.out