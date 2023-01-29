#%%
import numpy as np
#%%
import os
import subprocess 
import struct_maker as sm
from mpi4py.MPI import COMM_WORLD as CW 

rank = CW.Get_rank()
size = CW.Get_size()

nsim = 15
sim_per_rank = int(nsim / size) 

main_folder = '/lustre/astro/gfriss/pre_phase'
T_p = np.logspace(2, 6, nsim-1)
T_p = np.concatenate((np.array([0.]), T_p))

for i in range(rank*sim_per_rank, (rank+1)*sim_per_rank):
    if i < 10:
        new_folder = os.path.join(main_folder, 'sim_0' + str(i))
    else:
        new_folder = os.path.join(main_folder, 'sim_' + str(i))
    
    subprocess.run(['cp', '-r', '/lustre/astro/gfriss/og_inputs', new_folder], check = True, text = True)
    if i != 0: # keep i=0 as the original
        sm.write_w_pre(new_folder, T_p[i])

    # input files done
    # run the simulation in the given folder
    subprocess.check_call(['run_sim.sh', new_folder])

