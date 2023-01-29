import os
import subprocess # to run bash commands like in a terminal
# this is how to use it:
# subprocess.run(['everything', 'you', 'would', 'type', 'word', 'by', 'word'], check = True, text = True)
import new_rates as nr # this is a python program written by me to change the rates
from mpi4py.MPI import COMM_WORLD as CW # for paralellisation

rank = CW.Get_rank()
size = CW.Get_size()

nsim = 12
sim_per_rank = int(nsim / size) # this is needed to distribure the tasks between tha CPUs

main_folder = '/lustre/astro/gfriss/test_many' # place to store outputs

for i in range(rank*sim_per_rank, (rank+1)*sim_per_rank):   # paralellisation itself, it spreads the task between the CPUs
                                                            # this is the magic, after this just think of it as a normal, sequential loop
    if i < 10:
        new_folder = os.path.join(main_folder, 'sim_0' + str(i)) # my folder naming conventions
    else:
        new_folder = os.path.join(main_folder, 'sim_' + str(i))
    #subprocess.run(['mkdir', new_folder], check = True, text = True)
    subprocess.run(['cp', '-r', '/lustre/astro/gfriss/og_inputs', new_folder], check = True, text = True)
    nr.write(new_folder)
    # input files done
    # run the simulation in the given folder
    subprocess.check_call(['run_sim.sh', new_folder]) # this calls a bash script that runs nautilus in the given folder

