#%%
from json.tool import main
from matplotlib import tight_layout
import numpy as np
import matplotlib.pyplot as plt
import os
import plot_reset as pr
from scipy.integrate import cumtrapz

N = 2000

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
    
    for folder in sorted(os.listdir(main_folder)):
        sim = os.path.join(main_folder, folder)
        path = os.path.join(sim, 'ab')
        for file in os.listdir(path):
            if file[:-3] == spec:
                X = np.loadtxt(os.path.join(path, file), comments = '!')[:, 1]
                for i in range(len(t)):
                    X_s[t[i]].append(X[i])
                break  # so that it doesn't go through unnecessary files
    
    X_mean = []
    X_min1, X_min2 = [], []
    X_max1, X_max2 = [], []
    for key in X_s.keys():
        X_time = np.array(X_s[key]).flatten()
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

def plot_species(data, og_path, fsize = (14,6), figname = None):
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
        ax[1].semilogy(tim, X_og, color = c, label = spec)
    ax[0].set_xlabel('t [yr]')
    ax[1].set_xlabel('t [yr]')
    ax[0].set_ylabel(r'n(X)/n(H$_2$)')
    ax[1].legend(loc = 10, bbox_to_anchor = (1.2, 0.5))
    if figname != None:
        fig.savefig(figname, bbox_inches = 'tight')


#%%
pr.reset_plt(20, 20)
sims = '/lustre/astro/gfriss/diff_rate'
#out = '/lustre/astro/gfriss/diff_rate_out'
out = '/lustre/hpc/astro/gfriss'
og_sim = 'sim_n5_t50'
#%%

spec1 = ['JCO', 'JHCO', 'JH2CO', 'JCH3O', 'JCH2OH', 'JCH3OH']
spec2 = ['JHCOOH', 'JHNCO', 'JNH2CHO', 'JHCOOCH3']
spec3 = ['CO', 'N2', 'HCO+', 'H3+', 'N2H+']
d1 = collect_species(sims, spec1)
d2 = collect_species(sims, spec2)
d3 = collect_species(sims, spec3)

plot_species(d1, og_sim, figname = os.path.join(out, 'methanol_chain_w_uncer.pdf'))
plot_species(d2, og_sim, figname = os.path.join(out, 'COMs_w_uncer.pdf'))
plot_species(d3, og_sim, figname = os.path.join(out, 'main5_w_uncer.pdf'))
