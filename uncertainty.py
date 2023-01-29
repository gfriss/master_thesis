#%%
from json.tool import main
from matplotlib import tight_layout
import numpy as np
import matplotlib.pyplot as plt
import os
import plot_reset as pr
from mpi4py.MPI import COMM_WORLD as CW
from scipy.integrate import cumtrapz

rank = CW.Get_rank()
size = CW.Get_size()
N = 2000
sim_per_rank = int(N / size)

def name_for_plot(name):
    new_name = ''
    surface = False
    mantle = False
    for letter in name:
        if letter == 'J':
            surface = True
        elif letter == 'K':
            mantle = True
        elif letter.isnumeric() == True:
            new_name += '$_'+letter+'$'
        elif letter == '+' or letter == '-':
            if letter[-1] == '$':
                new_name = new_name[:-1] + '^' + letter + '$'
            else:
                new_name += '$^'+letter+'$'
        else:
            new_name += letter
    if surface == True:
        new_name += '$_s$'
    elif mantle == True:
        new_name += '$_m$'
    return new_name

def mean_uncertain(main_folder):
    data1, data2 = {}, {}
    species = []
    path = os.path.join(main_folder, os.listdir(main_folder)[0], 'ab')
    for file in os.listdir(path):
        spec = file[:-3]
        species.append(spec)
    interval1 = int(N * 0.683) # for 1 sigma
    interval2 = int(N * 0.954) # for 2 sigma
    
    for s in species:
        data1[s] = one_spec(s, main_folder, interval1)
        data2[s] = one_spec(s, main_folder, interval2)
    return data1, data2

def one_spec(spec, main_folder, interval):
    X_s = []
    for folder in os.listdir(main_folder)[rank*sim_per_rank:(rank+1)*sim_per_rank]:
        sim = os.path.join(main_folder, folder)
        #if os.path.isdir(sim) == True:
        path = os.path.join(sim, 'ab')
        for file in os.listdir(path):
            if file[:-3] == spec:
                with open(os.path.join(path, file)) as f:
                    # reads in the abundance (2nd cloumn) at the last timestep (1st column)
                    X_s.append(float(f.readlines()[-1].split()[-1]))
    X_spec = CW.gather(X_s, root = 0)
    if rank == 0:
        X_spec = np.array(X_spec).flatten()
        logX_spec = np.sort(np.log10(X_spec))
        X_mean = np.mean(X_spec)

        logX_diff = logX_spec[-1] - logX_spec[0]
        j = 0   # this will tell the index-interval+index values for 
                # logXmin and logXmax, respectively
        for i in range(N - interval):
            logX_diff_new = logX_spec[interval+i] - logX_spec[i]
            if logX_diff_new < logX_diff:
                logX_diff = logX_diff_new
                j = i

        logX_min = logX_spec[j]
        logX_max = logX_spec[interval + j]
        delta_logX = 0.5 * (logX_max - logX_min)
        return [X_mean, delta_logX]

def plot_uncertainty(data, fsize = (8,6), figname = None):
    delta_logX = np.array(list(data.values()))[:, 1]
    fig, ax = plt.subplots(figsize = fsize)

    ax.hist(delta_logX, bins = 20, color = 'r', edgecolor = 'k', linewidth = 1.5)
    ax.axvline(np.mean(delta_logX), linestyle = '--', color = 'gray', label = r'$\mu$ = {:.2f}'.format(np.mean(delta_logX)))
    ax.axvline(np.median(delta_logX), linestyle = '-.', color = 'darkorange', label = 'm = {:.2f}'.format(np.median(delta_logX)))
    if '1sigma' in figname:
        ax.set_xlabel(r'$1\sigma_{log(X)}$')
    elif '2sigma' in figname:
        ax.set_xlabel(r'$2\sigma_{log(X)}$')
    ax.set_ylabel('# species')
    ax.legend()
    if figname != None:
        fig.savefig(figname)
        
def trim_statistical(data, fsize = (14,6), figname = None):
    average_logX = np.log10(np.array(list(data.values()))[:, 0])
    delta_logX = np.array(list(data.values()))[:, 1]
    max_logX = average_logX + delta_logX
    n, be, _ = plt.hist(delta_logX, bins = 100, density = True)
    histrange = (0, be[-1])
    bw = be[1] - be[0]
    bc = be[:-1] + bw/2
    I = cumtrapz(n, bc, initial = 0.)
    
    fig, ax = plt.subplots(ncols = 3, figsize = fsize, sharey = True, tight_layout = True)
    ax = ax.flatten()
    # species with rel. abundance < 1e-20 (even 1e-14 according to Jes so set it to 1e-18) are insignificant
    delta_trimmed = delta_logX[max_logX > -18]
    ax[0].hist(delta_trimmed, bins = 20, range = histrange, color = 'r', edgecolor = 'k', linewidth = 1.5)
    ax[0].axvline(np.mean(delta_trimmed), linestyle = '--', color = 'gray', label = r'$\mu$ = {:.2f}'.format(np.mean(delta_trimmed)))
    ax[0].axvline(np.median(delta_trimmed), linestyle = '-.', color = 'darkorange', label = 'm = {:.2f}'.format(np.median(delta_trimmed)))
    if '1sigma' in figname:
        ax[0].set_xlabel(r'$1\sigma_{log(X)}$')
    elif '2sigma' in figname:
        ax[0].set_xlabel(r'$2\sigma_{log(X)}$')
    ax[0].set_ylabel('# species')
    ax[0].legend()
    
    # throw out species whose error values have a p value < 0.01
    trim_01 = np.max(bc[I < 1 - 0.01]) + bw/2
    delta_trimmed_01 = delta_logX[delta_logX < trim_01]
    ax[1].hist(delta_trimmed_01, bins = 20, range = histrange, color = 'r', edgecolor = 'k', linewidth = 1.5)
    ax[1].axvline(np.mean(delta_trimmed_01), linestyle = '--', color = 'gray', label = r'$\mu$ = {:.2f}'.format(np.mean(delta_trimmed_01)))
    ax[1].axvline(np.median(delta_trimmed_01), linestyle = '-.', color = 'darkorange', label = 'm = {:.2f}'.format(np.median(delta_trimmed_01)))
    if '1sigma' in figname:
        ax[1].set_xlabel(r'$1\sigma_{log(X)}$')
    elif '2sigma' in figname:
        ax[1].set_xlabel(r'$2\sigma_{log(X)}$')
    ax[1].legend()
    
    # and p value < 0.05
    trim_05 = np.max(bc[I < 1 - 0.05]) + bw/2
    delta_trimmed_05 = delta_logX[delta_logX < trim_05]
    ax[2].hist(delta_trimmed_05, bins = 20, range = histrange, color = 'r', edgecolor = 'k', linewidth = 1.5)
    ax[2].axvline(np.mean(delta_trimmed_05), linestyle = '--', color = 'gray', label = r'$\mu$ = {:.2f}'.format(np.mean(delta_trimmed_05)))
    ax[2].axvline(np.median(delta_trimmed_05), linestyle = '-.', color = 'darkorange', label = 'm = {:.2f}'.format(np.median(delta_trimmed_05)))
    if '1sigma' in figname:
        ax[2].set_xlabel(r'$1\sigma_{log(X)}$')
    elif '2sigma' in figname:
        ax[2].set_xlabel(r'$2\sigma_{log(X)}$')
    ax[2].legend()
    
    if figname != None:
        fig.savefig(figname, bbox_inches = 'tight')

def print_large(data, main_folder, fname = 'large_error.txt'):
    with open(os.path.join(main_folder, fname), mode = 'w') as f:
        for spec in data.keys():
            if data[spec][1] > 2.:
                f.write(spec.ljust(13) + '\t{:.4e}'.format(data[spec][0]) + '\t{:.4f}\n'.format(data[spec][1]))

def print_small(data, main_folder, fname = 'small_error.txt'):
    with open(os.path.join(main_folder, fname), mode = 'w') as f:
        for spec in data.keys():
            if data[spec][1] <= 2.:
                f.write(spec.ljust(13) + '\t{:.4e}'.format(data[spec][0]) + '\t{:.4f}\n'.format(data[spec][1]))

def plot_convergence(data, species, main_folder, fsize = (8,6), figname = None):
    X_final = {}
    X_sim = {}
    X_sim_all = {}
    X_rel = {}
    if rank == 0:
        fig, ax = plt.subplots(figsize = fsize)
        n = np.arange(1, N+1)
    for spec in species:
        X_sim[spec] = []
        for folder in os.listdir(main_folder)[rank*sim_per_rank:(rank+1)*sim_per_rank]:
            sim = os.path.join(main_folder, folder)
            #if os.path.isdir(sim) == True:
            path = os.path.join(sim, 'ab')
            for file in os.listdir(path):
                if file[:-3] == spec:
                    with open(os.path.join(path, file)) as f:
                        # reads in the abundance (2nd cloumn) at the last timestep (1st column)
                        X_sim[spec].append(float(f.readlines()[-1].split()[-1]))

        X_sim_all[spec] = CW.gather(X_sim[spec], root = 0)
        if rank ==  0:
            X_final[spec] = data[spec][0]
            X_sim_all[spec] = np.array(X_sim_all[spec]).flatten()
            X_rel[spec] = np.cumsum(X_sim_all[spec]) / n / X_final[spec]

            ax.plot(n, X_rel[spec], label = name_for_plot(spec))
    if rank == 0:
        ax.set_xlabel('# simulation')
        ax.set_ylabel(r'<$X_n$>/<$X_f$>')
        ax.legend(loc = 'lower right', ncol = 2)
        fig.savefig(figname)

def collect_species(main_folder, species):
    data = {}    
    for s in species:
        data[s] = collect_one(s, main_folder)
    return data

def collect_one(spec, main_folder):
    interval1 = int(N * 0.683) # for 1 sigma
    interval2 = int(N * 0.954) # for 2 sigma

    # time is same for all, so just get it from one
    path_to_ab = os.path.join(main_folder, os.listdir(main_folder)[0], 'ab')
    path_to_time = os.path.join(path_to_ab, os.listdir(path_to_ab)[0])
    t = np.loadtxt(path_to_time, comments = '!')[:, 0]
    X_s = {}
    for T in t:
        X_s[T] = []
    
    for folder in sorted(os.listdir(main_folder))[rank*sim_per_rank:(rank+1)*sim_per_rank]:
        sim = os.path.join(main_folder, folder)
        path = os.path.join(sim, 'ab')
        for file in os.listdir(path):
            if file[:-3] == spec:
                X = np.loadtxt(os.path.join(path, file), comments = '!')[:, 1]
                for i in range(len(t)):
                    X_s[t[i]].append(X[i])
                break  # so that it doesn't go through unnecessary files
    X_spec = {}
    for T in t:
        X_spec[T] = CW.gather(X_s[T], root = 0)

    if rank == 0:    
        X_mean = []
        X_min1, X_min2 = [], []
        X_max1, X_max2 = [], []
        for key in X_spec.keys():
            X_time = np.array(X_spec[key]).flatten()
            X_mean.append(np.mean(X_time))
            logX_spec = np.sort(np.log10(X_time))
            logX_diff1 = logX_spec[-1] - logX_spec[0]
            logX_diff2 = np.copy(logX_diff1)
            j1 = 0   # this will tell the index-interval+index values for 
            j2 = 0   # logXmin and logXmax, respectively
            for i in range(N - interval1):
                logX_diff_new = logX_spec[interval1+i] - logX_spec[i]
                if logX_diff_new < logX_diff1:
                    logX_diff1 = logX_diff_new
                    j1 = i
            for i in range(N - interval2):
                logX_diff_new = logX_spec[interval2+i] - logX_spec[i]
                if logX_diff_new < logX_diff2:
                    logX_diff2 = logX_diff_new
                    j2 = i

            X_min1.append(10**logX_spec[j1])
            X_max1.append(10**logX_spec[interval1 + j1])
            X_min2.append(10**logX_spec[j2])
            X_max2.append(10**logX_spec[interval2 + j2])
        output = {}
        output['time'] = t
        output['mean'] = np.array(X_mean)
        output['min1'] = np.array(X_min1)
        output['max1'] = np.array(X_max1)
        output['min2'] = np.array(X_min2)
        output['max2'] = np.array(X_max2)
        return output

def plot_species(data, og_path = 'sim_n5_t50', fsize = (14,6), figname = None):
    tim = data[list(data.keys())[0]]['time']
    fig, ax = plt.subplots(ncols = 2, sharey = True, figsize = fsize, tight_layout = True)
    ax = ax.flatten()
    for spec in data.keys():
        path_spec = os.path.join(og_path, 'ab', spec+'.ab')
        X_og = np.loadtxt(path_spec, comments = '!')[:, 1]
        p = ax[0].semilogy(tim, data[spec]['mean'])
        c = p[0].get_color()
        ax[0].semilogy(tim, data[spec]['min1'], linestyle = '--', color = c, alpha = 0.8)
        ax[0].semilogy(tim, data[spec]['max1'], linestyle = '--', color = c, alpha = 0.8)
        ax[0].semilogy(tim, data[spec]['min2'], linestyle = ':', color = c, alpha = 0.6)
        ax[0].semilogy(tim, data[spec]['max2'], linestyle = ':', color = c, alpha = 0.6)
        ax[1].semilogy(tim, X_og, color = c, label = name_for_plot(spec))
    ax[0].set_xlabel('t [yr]')
    ax[1].set_xlabel('t [yr]')
    ax[0].set_ylabel(r'n(X)/n(H$_2$)')
    ax[1].legend(loc = 10, bbox_to_anchor = (1.2, 0.5))
    if figname != None:
        fig.savefig(figname, bbox_inches = 'tight')


#%%
pr.reset_plt(20, 20)
sims = '/lustre/astro/gfriss/diff_rate'
out = '/lustre/astro/gfriss/diff_rate_out'
#out = '/lustre/hpc/astro/gfriss'
#%%

spec1 = ['JCO', 'JHCO', 'JH2CO', 'JCH3O', 'JCH2OH', 'JCH3OH']
spec2 = ['JHCOOH', 'JHNCO', 'JNH2CHO', 'JHCOOCH3']
spec3 = ['CO', 'N2', 'HCO+', 'H3+', 'N2H+']
d1 = collect_species(sims, spec1)
d2 = collect_species(sims, spec2)
d3 = collect_species(sims, spec3)

D1, D2 = mean_uncertain(sims)
if rank == 0:
    '''plot_uncertainty(D1, figname = os.path.join(out, 'uncertainties_1sigma.pdf'))
    print_large(D1, main_folder = out, fname = 'large_error_1sigma.txt')
    trim_statistical(D1, figname = os.path.join(out, 'trimmed_uncertainties_1sigma.pdf'))
    plot_uncertainty(D2, figname = os.path.join(out, 'uncertainties_2sigma.pdf'))
    print_large(D2, main_folder = out, fname = 'large_error_2sigma.txt')
    trim_statistical(D2, figname = os.path.join(out, 'trimmed_uncertainties_2sigma.pdf'))
    print_small(D1, main_folder = out, fname = 'small_error_1sigma.txt')
    print_small(D2, main_folder = out, fname = 'small_error_2sigma.txt')
    '''
    plot_species(d1, figname = os.path.join(out, 'methanol_chain_w_uncer.pdf'))
    plot_species(d2, figname = os.path.join(out, 'COMs_w_uncer.pdf'))
    plot_species(d3, figname = os.path.join(out, 'main5_w_uncer.pdf'))
plot_convergence(D1, ['JCO', 'JHCO', 'JH2CO', 'JCH3O', 'JCH2OH', 'JCH3OH'], sims, figname = os.path.join(out, 'convergence_methanol_chain_1sigma.pdf'))
plot_convergence(D1, ['JHCOOH', 'JHNCO', 'JNH2CHO', 'JHCOOCH3'], sims, figname = os.path.join(out, 'convergence_other_1sigma.pdf'))
plot_convergence(D1, ['CO', 'N2', 'HCO+', 'H3+', 'N2H+'], sims, figname = os.path.join(out, 'convergence_5main_1sigma.pdf'))
