#%%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import os
import plot_reset as pr

def get_gas(path, species):
    path_gas = os.path.join(path, 'ab', species+'.ab')
    X = np.loadtxt(path_gas, comments = '!')[:, 1]
    return X

def add_ices(path, species, time = False):
    path_j = os.path.join(path, 'ab', 'J'+species+'.ab')
    path_k = os.path.join(path, 'ab', 'K'+species+'.ab')
    jspec = np.loadtxt(path_j, comments = '!')[:, 1]
    kspec = np.loadtxt(path_k, comments = '!')[:, 1]
    if time == True:
        t = np.loadtxt(path_j, comments = '!')[:, 0]
        return t, jspec + kspec
    else:
        return jspec + kspec

def ice_wrt_water_param(folder, species, co = False, observations = None, fsize = (8,6), figname = None):
    ''' Compare abundances of water (or methanol) and species for every 
        simulation for certain parameter with observations.'''
    '''if typ == 'dens': # from many_sim_plot 
        n = np.logspace(4, 7, 15)
        lab = []
        for rho in n:
            lab.append('{:.3e}'.format(rho) + r' cm$^{-3}$')
    elif typ == 'Tb':
        Tb = np.linspace(10, 50, 15)
        lab = []
        for tb in Tb:
            lab.append('{:.1f} K'.format(tb))
    elif typ == 'tq':
        tauq = np.logspace(2.69, 4, 15)
        lab = []
        for tq in tauq:
            lab.append('{:.3e} yr'.format(tq))
    elif typ == 'pre':
        ...
    # in case this is needed, put typ back as argument
    '''
    fig, ax = plt.subplots(figsize = fsize, tight_layout = True)
    cmap = plt.get_cmap('Reds_r')
    N = len(os.listdir(folder))
    i = 0
    for sim in sorted(os.listdir(folder), reverse = True):
        path_ab = os.path.join(folder, sim)
        if os.path.isdir(path_ab) == True:
            if co == True:
                t, X_water = add_ices(path_ab, 'CO', time = True)
            else:
                t, X_water = add_ices(path_ab, 'H2O', time = True)
            X_spec = add_ices(path_ab, species)
            if i == 0 and observations != None:
                for key in observations.keys(): # label all the or only in caption?
                    if key == 'range':
                        y1 = observations[key][0]*np.ones(np.shape(t))
                        y2 = observations[key][1]*np.ones(np.shape(t))
                        ax.fill_between(t, y1, y2, facecolor = 'darkgray')
                    elif key == 'mean':
                        ax.axhline(observations[key], xmin = 0.045, xmax = 0.955, color = 'k', linestyle = ':')
                    elif key == 'comet':
                        y1 = observations[key][0]*np.ones(np.shape(t))
                        y2 = observations[key][1]*np.ones(np.shape(t))
                        ax.fill_between(t, y1, y2, facecolor = 'saddlebrown', alpha = 0.8)
            c = cmap(i/N)
            ax.semilogy(t, X_spec/X_water, color = c)#, label = lab[i])
            i += 1
    
    ax.set_ylabel(r'n(X)/n(H$_2$O)')
    ax.set_xlabel('t [yr]')
    #ax.legend(loc = 10, bbox_to_anchor = (1.2, 0.5))
    if figname != None:
        fig.savefig(figname)
    plt.close(fig)

def gas_vs_ice_param(folder, species, observations = None, fsize = (8,6), figname = None):
    ''' Compare abundances of gas and ice of given species for 
        every simulation for certain parameter with observations.'''
    fig, ax = plt.subplots(figsize = fsize, tight_layout = True)
    cmap = plt.get_cmap('Reds_r')
    N = len(os.listdir(folder))
    i = 0
    for sim in sorted(os.listdir(folder), reverse = True):
        path_ab = os.path.join(folder, sim)
        if os.path.isdir(path_ab) == True: 
            t, X_ice = add_ices(path_ab, species, time = True)
            X_gas = get_gas(path_ab, species)
            if i == 0 and observations != None:
                for key in observations.keys(): # label all the or only in caption?
                    if key == 'range':
                        y1 = observations[key][0]*np.ones(np.shape(t))
                        y2 = observations[key][1]*np.ones(np.shape(t))
                        ax.fill_between(t, y1, y2, facecolor = 'darkgray')
                    elif key == 'mean':
                        ax.axhline(observations[key], xmin = 0.045, xmax = 0.955, color = 'k', linestyle = ':')
                    elif key == 'comet':
                        y1 = observations[key][0]*np.ones(np.shape(t))
                        y2 = observations[key][1]*np.ones(np.shape(t))
                        ax.fill_between(t, y1, y2, facecolor = 'saddlebrown', alpha = 0.8)
            c = cmap(i/N)
            ax.semilogy(t, X_gas/X_ice, color = c)#, label = lab[i])
            i += 1
    
    ax.set_ylabel(r'n(X$_{gas}$)/n(X$_{ice}$)')
    ax.set_xlabel('t [yr]')
    #ax.legend(loc = 10, bbox_to_anchor = (1.2, 0.5))
    if figname != None:
        fig.savefig(figname)
    plt.close(fig)
#%%
pr.reset_plt(20, 20)
# %%
obs_methanol_water = {'range':[0.05, 0.32], 'mean':np.mean([0.28, 0.32, 0.06, 0.11, 0.14, 0.17, 0.09, 0.06, 0.08]), 'comet':[0.002, 0.07]}
obs_co_water = {'range':[0.04, 0.73], 'mean':np.mean([0.65, 0.07, 0.21, 0.25, 0.052, 0.11, 0.73]), 'comet':[0.004, 0.3]}
obs_co2_water = {'range':[0.12, 0.39], 'mean':np.mean([0.2, 0.19, 0.28, 0.26]), 'comet':[0.1,0.24]}
obs_nh3_water = {'range':[0.03,0.1], 'mean':0.06, 'comet':[0.002,0.014]}
obs_ch4_water = {'range':[0.01,0.11], 'mean':0.045, 'comet':[0.004,0.016]}
obs_h2co_water = {'range':[0.02,0.07], 'mean':0.05, 'comet':[0.0011,0.01]}
obs_methanol_gas_ice = {'range':[1.4e-4,3.7e-3], 'mean':np.mean([1.4e-4,3.7e-3])}
obs_co_gas_ice = {'range':[1.,6.], 'mean':3.5}
obs_methanol_co = {'range':[0.12,3.17], 'mean':np.mean([0.32, 0.42, 0.12, 2.06, 0.23, 0.89, 1.73, 0.2, 3.17, 0.24, 0.27])}

n5_t50 = 'sim_n5_t50'

p_pre = '/lustre/astro/gfriss/pre_phase'
ice_wrt_water_param(p_pre, 'CH3OH', observations = obs_methanol_water, figname = 'ice_plots/methanol_water_diff_pre.pdf')
ice_wrt_water_param(p_pre, 'CO', observations = obs_co_water, figname = 'ice_plots/CO_water_diff_pre.pdf')
ice_wrt_water_param(p_pre, 'CO2', observations = obs_co2_water, figname = 'ice_plots/CO2_water_diff_pre.pdf')
ice_wrt_water_param(p_pre, 'NH3', observations = obs_nh3_water, figname = 'ice_plots/NH3_water_diff_pre.pdf')
ice_wrt_water_param(p_pre, 'CH4', observations = obs_ch4_water, figname = 'ice_plots/CH4_water_diff_pre.pdf')
ice_wrt_water_param(p_pre, 'H2CO', observations = obs_h2co_water, figname = 'ice_plots/H2CO_water_diff_pre.pdf')
ice_wrt_water_param(p_pre, 'CH3OH', co = True, observations = obs_methanol_co, figname = 'ice_plots/methanol_co_diff_pre.pdf')
gas_vs_ice_param(p_pre, 'CH3OH', observations = obs_methanol_gas_ice, figname = 'ice_plots/methanol_gas_ice_diff_pre.pdf')
gas_vs_ice_param(p_pre, 'CO', observations = obs_co_gas_ice, figname = 'ice_plots/CO_gas_ice_diff_pre.pdf')
