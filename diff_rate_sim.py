import os
import subprocess 
import new_rates as nr
from mpi4py.MPI import COMM_WORLD as CW 

rank = CW.Get_rank()
size = CW.Get_size()

nsim = 2000
sim_per_rank = int(nsim / size) 

main_folder = '/lustre/astro/gfriss/diff_rate'

for i in range(rank*sim_per_rank, (rank+1)*sim_per_rank):
    if i < 10:
        new_folder = os.path.join(main_folder, 'sim_000' + str(i))
    elif i >= 10 and i < 100:
        new_folder = os.path.join(main_folder, 'sim_00' + str(i))
    elif i >= 100 and i < 1000:
        new_folder = os.path.join(main_folder, 'sim_0' + str(i)) 
    else:
        new_folder = os.path.join(main_folder, 'sim_' + str(i))
    subprocess.run(['cp', '-r', '/lustre/astro/gfriss/og_inputs', new_folder], check = True, text = True)
    if i != 0: #let's keep the original sim as the 0th sim
        nr.write(new_folder)
    # input files done
    # run the simulation in the given folder
    subprocess.check_call(['run_sim.sh', new_folder])

