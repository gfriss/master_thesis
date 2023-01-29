#%%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import os
import plot_reset as pr

def name_for_plot(name):
    new_name = ''
    surface = False
    for letter in name:
        if letter == 'J':
            surface = True
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
    return new_name

main_folder = '/lustre/astro/gfriss'
def comp_spec(fold, species, ice = True):
    ''' Finding the given species also in their ice form (if needed) since the
        temperature are low for the single point model.
        Returns the changed conditions and the corresponding abundaces.'''
    X = {}
    if ice == True:
        X_ice = {}
        jspecies, kspecies = [], []
        for spec in species:
            jspecies.append('J'+spec)
            kspecies.append('K'+spec)
            X['J'+spec] = []
            X_ice['J'+spec] = []
            X_ice['K'+spec] = []
    for spec in species:
        X[spec] = []
    for folder in sorted(os.listdir(os.path.join(main_folder,fold))):
        sim = os.path.join(main_folder,fold,folder)
        if os.path.isdir(sim) == True:
            path = os.path.join(main_folder,fold,folder,'ab')
            for file in os.listdir(path):
                if file[:-3] in species:
                    with open(path + '/' + file) as f:
                        X[file[:-3]].append(float(f.readlines()[-1].split()[-1]))
                if ice == True and file[:-3] in jspecies:
                    with open(path + '/' + file) as f:
                        X_ice[file[:-3]].append(float(f.readlines()[-1].split()[-1]))
                elif ice == True and file[:-3] in kspecies:
                    with open(path + '/' + file) as f:
                        X_ice[file[:-3]].append(float(f.readlines()[-1].split()[-1]))

    if ice == True:
        for spec in species:
            X['J'+spec] = list(np.array(X_ice['J'+spec]) + np.array(X_ice['K'+spec]))
    # parameter values either from folder name or the input file or the jupyter notebook (many_sim)
    # let's go with the jupyter notebook version so that this can get the species
    # -> the values will be according to the parameter value order in the notebook

    return X

def plot(param, X, xlab, fsize = (8,6), xscale = 'log', yscale = 'log', figname = None):
    fig, ax = plt.subplots(figsize = fsize)
    colour = {}
    for spec in X.keys():
        if spec[0] == 'J':
            p = ax.plot(param, X[spec], label = name_for_plot(spec[1:]))
            colour[spec[1:]] = p[0].get_color()
        else:
            ax.plot(param, X[spec], color = colour[spec], linestyle = ':')
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    ax.set_xlabel(xlab)
    ax.set_ylabel(r'n(X)/n(H$_2$)')
    ax.legend(loc = 10, bbox_to_anchor = (1.2, 0.5))
    if figname != None:
        fig.savefig(figname, bbox_inches = 'tight')

pr.reset_plt(20, 20)

spec1 = ['CO', 'HCO', 'H2CO', 'CH3O', 'CH2OH', 'CH3OH']
spec2 = ['HCOOCH3', 'HCOOH', 'NH2CHO', 'HNCO']
nsim= 15
T_p = np.logspace(2, 6, nsim-1)
T_p = np.concatenate((np.array([0.]), T_p))
#%%
X = comp_spec('pre_phase', spec1)
plot(T_p, X, r'$\tau_{pre}$ [yr]', figname = os.path.join(main_folder, 'pre_phase', 'pre1.pdf'))

X = comp_spec('pre_phase', spec2)
plot(T_p, X, r'$\tau_{pre}$ [yr]', figname = os.path.join(main_folder, 'pre_phase', 'pre2.pdf'))

X = comp_spec('pre_12K', spec1)
plot(T_p, X, r'$\tau_{pre}$ [yr]', figname = os.path.join(main_folder, 'pre_12K', 'pre1_12K.pdf'))

X = comp_spec('pre_12K', spec2)
plot(T_p, X, r'$\tau_{pre}$ [yr]', figname = os.path.join(main_folder, 'pre_12K', 'pre2_12K.pdf'))
#%%
X = comp_spec('pre_15K', spec1)
plot(T_p, X, r'$\tau_{pre}$ [yr]', figname = os.path.join(main_folder, 'pre_15K', 'pre1_15K.pdf'))

X = comp_spec('pre_15K', spec2)
plot(T_p, X, r'$\tau_{pre}$ [yr]', figname = os.path.join(main_folder, 'pre_15K', 'pre2_15K.pdf'))