#%%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import os

from requests import get
# structure of rates.out:
# r1 r2 r3 p1 p2 p3 p4 p5 rate_t1 rate_t2.... id noidea
# if e.g. r3 is not existant then it says XXX
# aim: find given species and the corresponding production/destruction rates

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

def filt_xxx(side):
    ind = []
    for i in range(len(side)):
        if side[i] == 'XXX':
            ind.append(i)
    return np.delete(side, ind)

def filt_swap(reacts, prods):
    if len(reacts) == 1 and len(prods) == 1:
        if reacts[0][1:] == prods[0][1:]:
            return True # meaning we have a swapping reaction
    else:
        return False

def prod_or_dest(line, species, l, no_swap = True):
    ''' Decides whether the given reaction (line in the file) is a 
        production or destruction reaction for the given species.
        For the ices, it ignores (if wanted) the swapping reacions 
        to avoid bias like the built in version in pnautilus.'''
    values = line.split()[:-1] # we don't need the last element (itype)
    reactants = filt_xxx(values[:3])
    products = filt_xxx(values[3:8])
    if len(values[8:]) != l + 1: # 5-digit ID is together with last rate value (+1 for the ID)
        ID = values[-1][-5:]
        values[-1] = values[-1][:-5]
        values.append(ID)
    ID = values[-1]
    production = True
    swap = filt_swap(reactants, products)
    if species not in products:
        production = False
    if no_swap == True and swap == True:
        production = 'swap'
        rates = np.array([])
        reaction = ''
    else:
        rates = np.array(values[8: -1], dtype = float)
        reaction = ID + ': ' + ' + '.join(reactants) + ' -> ' + ' + '.join(products)
    return production, reaction, rates
    
def find_spec(folder, species, no_swap = True):
    p, rp = [], []
    d, rd = [], []
    file = os.path.join(folder,'rates.out')
    with open(file) as f:
        n = 0
        for row in f:
            if n == 0:
                t = np.array(row.split()[2:], dtype = float)
                n += 1
            if species in row.split() and '*' not in row:                  # too low rates are labeled with *
                prod, react, rate = prod_or_dest(row, species, l = len(t), no_swap = no_swap) # instead of 0, skip them
                if prod == True: # production reaction
                    p.append(rate)
                    rp.append(react)
                elif prod == 'swap': 
                    continue
                elif prod == False: # destruction reaction
                    d.append(rate)
                    rd.append(react)
    return t, np.array(p), rp, np.array(d), rd

def plot_rates(t, rate, reaction, folder, typ, spec, \
               lim = 0.05, fsize = (14, 6), yscale = 'log', xscale = 'linear', \
               yrange0 = (1e-5, 1.1), yrange1 = (1e-20, 1e-8), figname = None):
    ''' Plots the relative and original rates for the processes that are
        over the set limit (lim). The labelling is done with ID and the
        corresponding reactions are written to [prod/dest]_[spec].dat'''
    file = os.path.join(folder,typ + '_' + spec + '.dat')
    fig, ax = plt.subplots(ncols = 2, figsize = fsize, tight_layout = True)
    ax = ax.flatten()
    summa = np.sum(rate, 0)
    rel = rate / summa
    ind = []
    num = 0
    with open(file, 'w') as f:
        for i in range(len(rel)):
            if np.any(rel[i] > lim):
                ind.append(i)
                f.write(reaction[i] + '\n')
                lab = reaction[i].split()[0][:-1] # label is the ID, without the ':'
                if num % 2 == 1:
                    ax[0].plot(t, rel[i], '--', label = lab)
                    ax[1].plot(t, rate[i], '--')
                else:
                    ax[0].plot(t, rel[i], label = lab)
                    ax[1].plot(t, rate[i])
                num += 1
    for i in range(len(ax)):
#         ax[i].set_yscale(yscale)
#         ax[i].set_xscale(xscale)
        ax[i].set_xlabel('t [yr]')
    ax[0].set_yscale(yscale)
    ax[0].set_xscale(xscale)
    ax[0].set_ylim(yrange0)
    ax[1].set_yscale('log')
    ax[1].set_xscale('linear')
    ax[1].set_ylim(yrange1)
    ax[0].set_ylabel('Relative ' + typ + ' rate of ' + name_for_plot(spec))
    ax[1].set_ylabel(r'Actual rate [cm$^{-3}$ s$^{-1}$]')
    fig.legend(loc = 'center', bbox_to_anchor = (0.525, 0), ncol = 5)
    print(np.array(reaction)[ind])
    if figname != None:
        fig.savefig(figname, bbox_inches = 'tight')
    
def plot_sum(simulation, spec, no_swap = True, figname = None, \
               fsize = (8, 6), yscale = 'log', xscale = 'linear', version = 1):
    ''' Plots the overall production and destruction rate for 
    	both the burst and no burst cases.'''
    	
    if version == 1:
        t, pb, _, db, _ = find_spec('remake_visser/sim_' + simulation + '/', spec, no_swap = no_swap)
        _, pnb, _, dnb, _ = find_spec('sim_noburst/' + simulation + '/', spec, no_swap = no_swap)
    elif version == 3:
        t, pb, _, db, _ = find_spec('v3_sim/' + simulation + '/', spec, no_swap = no_swap)
        _, pnb, _, dnb, _ = find_spec('v3_sim/nb_' + simulation + '/', spec, no_swap = no_swap)
    
    fig, ax = plt.subplots(ncols = 1, figsize = fsize, tight_layout = True)
    
    summa_pb = np.sum(pb, 0)
    summa_db = np.sum(db, 0)
    summa_pnb = np.sum(pnb, 0)
    summa_dnb = np.sum(dnb, 0)
    ax.plot(t, summa_pb, 'r')
    ax.plot(t, summa_pnb, 'r--')    
    ax.plot(t, summa_db, 'k')
    ax.plot(t, summa_dnb, 'k--')
    ax.set_xlabel('t [yr]')
    ax.set_yscale(yscale)
    ax.set_xscale(xscale)
    ax.set_ylabel('Rate')
    line_b = mlines.Line2D([], [], linestyle = '-', color = 'gray', label = 'Burst')
    line_nb = mlines.Line2D([], [], linestyle = '--', color = 'gray', label = 'No burst')
    line_prod = mlines.Line2D([], [], linestyle = '-', color = 'r', label = 'Production')
    line_dest = mlines.Line2D([], [], linestyle = '-', color = 'k', label = 'Destruction')
    ax.legend(loc = 'lower right', handles = [line_b, line_nb, line_prod, line_dest])
    if figname != None:
        fig.savefig(figname, bbox_inches = 'tight')
        
def get_species(species, folder, all_ice = False):
    '''abundances in dict so that the ices can be compared'''
    X = {}
    j = 0
    for spec in species:
        filename = os.path.join(folder,'ab',spec + '.ab')
        with open (filename) as f:
            x = []
            for row in f:
                if row[0] != '!':
                    x.append(float(row.split()[1]))
        if spec[1:] == species[j - 1][1:] and spec != species[0]: # in case we look at ices in 3 phase
                                                                  # also, not circulating with the first species
            X[species[j-1]] += np.array(x)
        else:
            X[spec] = np.array(x)
        j += 1
    return X
        
# def plot_sum_all(simulation, species, no_swap = True, figname = None, \
#                fsize = (14, 10), yscale = 'log', xscale = 'linear', version = 1):
#     ''' Plots the overall production and destruction rate for 
#         both the burst and no burst cases for all species.'''
    
#     t, summa_pb, summa_db = {}, {}, {}
#     summa_pnb, summa_dnb = {}, {}
#     for spec in species:
#         if version == 1:
#             t[spec], pb, _, db, _ = find_spec('remake_visser/sim_' + simulation + '/', spec, no_swap = no_swap)
#             _, pnb, _, dnb, _ = find_spec('sim_noburst/' + simulation + '/', spec, no_swap = no_swap)
#         elif version == 3:
#             t[spec], pb, _, db, _ = find_spec('v3_sim/' + simulation + '/', spec, no_swap = no_swap)
#             _, pnb, _, dnb, _ = find_spec('v3_sim/nb_' + simulation + '/', spec, no_swap = no_swap)

#         summa_pb[spec] = np.sum(pb, 0)
#         summa_db[spec] = np.sum(db, 0)
#         summa_pnb[spec] = np.sum(pnb, 0)
#         summa_dnb[spec] = np.sum(dnb, 0)
    
#     fig, ax = plt.subplots(ncols = 3, nrows = 2, figsize = fsize, tight_layout = True, sharey = True)
#     ax = ax.flatten()
    
#     for i in range(len(ax)):
#         ax[i].plot(t[species[i]], summa_pb[species[i]], 'r')
#         ax[i].plot(t[species[i]], summa_pnb[species[i]], 'r--')
#         ax[i].plot(t[species[i]], summa_db[species[i]], 'k')
#         ax[i].plot(t[species[i]], summa_dnb[species[i]], 'k--')
#         ax[i].set_yscale(yscale)
#         ax[i].set_xscale(xscale)
#         ax[i].set_ylim(1e-19, )
#         if i > 2:
#             ax[i].set_xlabel('t [yr]')
#         if i == 0 or i == 3:
#             ax[i].set_ylabel('Rate')
#     line_b = mlines.Line2D([], [], linestyle = '-', color = 'gray', label = 'Burst')
#     line_nb = mlines.Line2D([], [], linestyle = '--', color = 'gray', label = 'No burst')
#     line_prod = mlines.Line2D([], [], linestyle = '-', color = 'r', label = 'Production')
#     line_dest = mlines.Line2D([], [], linestyle = '-', color = 'k', label = 'Destruction')
#     ax[0].legend(loc = 'lower right', handles = [line_b, line_nb, line_prod, line_dest])
#     if figname != None:
#         fig.savefig(figname, bbox_inches = 'tight')

def plot_sum_all(species, no_swap = True, figname = None, overall = True, \
               fsize = (14, 36), yscale = 'log', xscale = 'linear', version = 1):
    ''' Plots the overall production and destruction rate for 
        both the burst and no burst cases for all species and 2 simulations.'''
    sim1 = 'n5_t5'
    sim2 = 'n5_t50'
    summa_pb1, summa_db1 = {}, {}
    summa_pnb1, summa_dnb1 = {}, {}
    summa_pb2, summa_db2 = {}, {}
    summa_pnb2, summa_dnb2 = {}, {}
    if version == 1:
        X_b1 = get_species(species, 'remake_visser/sim_' + sim1 + '/')
        X_nb1 = get_species(species, 'sim_noburst/' + sim1 + '/')
        X_b2 = get_species(species, 'remake_visser/sim_' + sim2 + '/')
        X_nb2 = get_species(species, 'sim_noburst/' + sim2 + '/')
    elif version == 3:
        X_b1 = get_species(species, 'v3_sim/' + sim1 + '/')
        X_nb1 = get_species(species, 'v3_sim/nb_' + sim1 + '/')
        X_b2 = get_species(species, 'v3_sim/' + sim2 + '/')
        X_nb2 = get_species(species, 'v3_sim/nb_' + sim2 + '/')
        
    for spec in species:
        j = 0
        if version == 1:
            if j == 0:
                t1, pb1, _, db1, _ = find_spec('remake_visser/sim_' + sim1 + '/', spec, no_swap = no_swap)
                t2, pb2, _, db2, _ = find_spec('remake_visser/sim_' + sim2 + '/', spec, no_swap = no_swap)
            else:
                _, pb1, _, db1, _ = find_spec('remake_visser/sim_' + sim1 + '/', spec, no_swap = no_swap)
                _, pb2, _, db2, _ = find_spec('remake_visser/sim_' + sim2 + '/', spec, no_swap = no_swap)
            _, pnb1, _, dnb1, _ = find_spec('sim_noburst/' + sim1 + '/', spec, no_swap = no_swap)
            _, pnb2, _, dnb2, _ = find_spec('sim_noburst/' + sim2 + '/', spec, no_swap = no_swap)
        elif version == 3:
            if j == 0:
                t1, pb1, _, db1, _ = find_spec('v3_sim/' + sim1 + '/', spec, no_swap = no_swap)
                t2, pb2, _, db2, _ = find_spec('v3_sim/' + sim2 + '/', spec, no_swap = no_swap)
            else:
                _, pb1, _, db1, _ = find_spec('v3_sim/' + sim1 + '/', spec, no_swap = no_swap)
                _, pb2, _, db2, _ = find_spec('v3_sim/' + sim2 + '/', spec, no_swap = no_swap)
            _, pnb1, _, dnb1, _ = find_spec('v3_sim/nb_' + sim1 + '/', spec, no_swap = no_swap)
            _, pnb2, _, dnb2, _ = find_spec('v3_sim/nb_' + sim2 + '/', spec, no_swap = no_swap)
        j += 1
        summa_pb1[spec] = np.sum(pb1, 0)
        summa_db1[spec] = np.sum(db1, 0)
        summa_pnb1[spec] = np.sum(pnb1, 0)
        summa_dnb1[spec] = np.sum(dnb1, 0)
        summa_pb2[spec] = np.sum(pb2, 0)
        summa_db2[spec] = np.sum(db2, 0)
        summa_pnb2[spec] = np.sum(pnb2, 0)
        summa_dnb2[spec] = np.sum(dnb2, 0)
    
    fig, ax = plt.subplots(ncols = 4, nrows = 6, figsize = fsize, tight_layout = True)
    ax = ax.flatten()
    
    p1_ind = [0, 4, 8, 12, 16, 20]
    x1_ind = [1, 5, 9, 13, 17, 21]
    p2_ind = [2, 6, 10, 14, 18, 22]
    x2_ind = [3, 7, 11, 15, 19, 23]
    for i in p1_ind:
        ind = i//4
        if overall == True:
            ax[i].plot(t1, summa_pb1[species[ind]] - summa_db1[species[ind]], 'r')
            ax[i].plot(t1, summa_pnb1[species[ind]] - summa_dnb1[species[ind]], 'r--')
        else:
            ax[i].plot(t1, summa_pb1[species[ind]], 'r')
            ax[i].plot(t1, summa_pnb1[species[ind]], 'r--')
            ax[i].plot(t1, summa_db1[species[ind]], 'k')
            ax[i].plot(t1, summa_dnb1[species[ind]], 'k--')
        ax[i].set_ylabel(r'Rate [cm$^{-3}$ s$^{-1}$]')
        ax[i].set_yscale(yscale)
        ax[i].set_xscale(xscale)
    for i in x1_ind:
        ind = i//4
        ax[i].plot(t1, X_b1[species[ind]], 'g', label = name_for_plot(species[ind]))
        ax[i].plot(t1, X_nb1[species[ind]], 'g--')
        ax[i].set_yscale(yscale)
        ax[i].set_xscale(xscale)
        ax[i].legend(loc = 'lower right')
        ax[i].set_ylabel(r'n(X)/n(H$_2$)')
    for i in p2_ind:
        ind = i//4
        if overall == True:
            ax[i].plot(t1, summa_pb2[species[ind]] - summa_db2[species[ind]], 'r')
            ax[i].plot(t1, summa_pnb2[species[ind]] - summa_dnb2[species[ind]], 'r--')
        else:
            ax[i].plot(t2, summa_pb2[species[ind]], 'r')
            ax[i].plot(t2, summa_pnb2[species[ind]], 'r--')
            ax[i].plot(t2, summa_db2[species[ind]], 'k')
            ax[i].plot(t2, summa_dnb2[species[ind]], 'k--')
        ax[i].set_yscale(yscale)
        ax[i].set_xscale(xscale)
        ax[i].set_ylabel(r'Rate [cm$^{-3}$ s$^{-1}$]')
    for i in x2_ind:
        ind = i//4
        ax[i].plot(t2, X_b2[species[ind]], 'g')
        ax[i].plot(t2, X_nb2[species[ind]], 'g--')
        ax[i].set_yscale(yscale)
        ax[i].set_xscale(xscale)
        ax[i].set_ylabel(r'n(X)/n(H$_2$)')
    
    line_b = mlines.Line2D([], [], linestyle = '-', color = 'gray', label = 'Burst')
    line_nb = mlines.Line2D([], [], linestyle = '--', color = 'gray', label = 'No burst')
    if overall == True:
        ax[0].legend(loc = 'lower right', handles = [line_b, line_nb])
    else:    
        line_prod = mlines.Line2D([], [], linestyle = '-', color = 'r', label = 'Production')
        line_dest = mlines.Line2D([], [], linestyle = '-', color = 'k', label = 'Destruction')
        ax[0].legend(loc = 'lower right', handles = [line_b, line_nb, line_prod, line_dest])
    ax[20].set_xlabel('t [yr]')
    ax[21].set_xlabel('t [yr]')
    ax[22].set_xlabel('t [yr]')
    ax[23].set_xlabel('t [yr]')
    if figname != None:
        fig.savefig(figname, bbox_inches = 'tight')



#%%

T, R, Rp, D, Rd = find_spec('/lustre/astro/gfriss/infall_own_with_pre', 'CH3CHO')
plot_rates(T, R, Rp, '/lustre/astro/gfriss/infall_own_with_pre', 'prod', 'CH3CHO', figname='/lustre/astro/gfriss/infall_own_with_pre/prod_CH3CHO.pdf')
plot_rates(T, D, Rd, '/lustre/astro/gfriss/infall_own_with_pre', 'dest', 'CH3CHO', figname='/lustre/astro/gfriss/infall_own_with_pre/dest_CH3CHO.pdf')

T, R, Rp, D, Rd = find_spec('/lustre/astro/gfriss/noburst_infall_own_with_pre', 'CH3CHO')
plot_rates(T, R, Rp, '/lustre/astro/gfriss/noburst_infall_own_with_pre', 'prod', 'CH3CHO', figname='/lustre/astro/gfriss/noburst_infall_own_with_pre/prod_CH3CHO_nb.pdf')
plot_rates(T, D, Rd, '/lustre/astro/gfriss/noburst_infall_own_with_pre', 'dest', 'CH3CHO', figname='/lustre/astro/gfriss/noburst_infall_own_with_pre/dest_CH3CHO_nb.pdf')