import numpy as np
import os

ksi = 5e-17 # cr rate is 1/s
Av = 15
T = 10 # calculate uncertainties for T=10K
kb = 1.38e-16 # Boltzmann in cgs
a = 1e-5 # grain radius in cm
amu = 1.66e-24 # atomic mass in g
cr_duration = 1e-5 # peak duration after CR hits dust
cr_T = 70 # peak temperature afte CR hits dust
Fe_ion = 3e-14
n_s = 1.5e15 # site density on grains in 1/cm2
N_s = n_s * 4 * np.pi * a**2 # number of sites on 1 grain
nh2 = 1e5 # h2 density in 1/cm3, overall gas density basically
chem_barr_thick = 1e-8 # chemical barrier thickness in cm
hbar = 1.055e-27 # hbar in cgs


def rate(a, b, c, typ):
    if typ == '0':
        return a * (T/300)**b # from ode_solver.f90
    elif typ == '1':
        return a * ksi
    elif typ == '2':
        return a * np.exp(-c*Av)
    elif typ == '3':
        return a * (T/300)**b * np.exp(-c/T)
    elif typ == '4':
        return a * b * (0.62 + 0.4767*c*np.sqrt(300/T))
    elif typ == '5':
        return a * b * (1. + 0.0967*c*np.sqrt(300/T) + (c**2/10.526)*300/T)
    elif typ == '10':
        return a # according to ode_solver.f90, these two aren't relevant for us since surface reactions are enabled
    elif typ == '11':
        return a
    
def new_a(k, b, c, typ):
    if typ == '0':
        new = k / (T/300)**b
    elif typ == '1':
        new = k / ksi
    elif typ == '2':
        new = k / np.exp(-c*Av)
    elif typ == '3':
        new = k / (T/300)**b * np.exp(-c/T)
    elif typ == '4':
        new = k / b * (0.62 + 0.4767*c*np.sqrt(300/T))
    elif typ == '5':
        new = k / b * (1. + 0.0967*c*np.sqrt(300/T) + (c**2/10.526)*300/T)
    elif typ == '10':
        new = k
    elif typ == '11':
        new = k

    if new == float('inf') or new == float('-inf'): # just to make sure
        return 1e-13
    else:
        return max([new, 1e-13]) # 1e-13 basically turns it off, smallest than every in the original

def new_params(line):
    values = line.split()
    F0_typ = values[-8] # logn or NA ()
    F0 = 2. # say it's 2 in default (in case NA)
    if F0_typ == 'logn' and float(values[-10]) != 0:
        F0 = float(values[-10])
    typ = values[-4]
    alpha = float(values[-13])
    beta = float(values[-12])
    gamma = float(values[-11])
    k0 = rate(alpha, beta, gamma, typ)
    N = np.random.normal(0, 1)
    k = np.exp(np.log(k0) + np.log(F0)*N)
    return '{:.3e}'.format(new_a(k, beta, gamma, typ))

def self_check(line, a):
    b = line.split()[-12]
    c = line.split()[-11]
    if a == b:
        print('Beta is tha same as alpha.')
    elif a == c:
        print('Gamma is tha same as alpha.')

def write(new_path):
    original_path = '/lustre/astro/gfriss/og_inputs'    
    original_gas = os.path.join(original_path, 'gas_reactions.in')
    original_grain = os.path.join(original_path, 'grain_reactions.in')
    with open(original_gas, mode = 'r') as og:
        rows = og.readlines()
        with open(os.path.join(new_path,'gas_reactions.in'), mode = 'w') as f:
            for row in rows:
                if row[0] == '!':
                    f.write(row)
                elif row == rows[-1]:
                    new_alpha = new_params(row)
                    old_alpha = row.split()[-13]
                    self_check(row, old_alpha)
                    new_row = row.replace(old_alpha, new_alpha).rstrip()
                    f.write(new_row)
                else:
                    new_alpha = new_params(row)
                    old_alpha = row.split()[-13]
                    self_check(row, old_alpha)
                    new_row = row.replace(old_alpha, new_alpha)
                    f.write(new_row)
    with open(original_grain, mode = 'r') as og:
        rows = og.readlines()
        with open(os.path.join(new_path, 'grain_reactions.in'), mode = 'w') as f:
            for row in rows:
                if row[0] == '!':
                    f.write(row)
                elif row == rows[-1]:
                    new_alpha = new_params(row)
                    old_alpha = row.split()[-13]
                    self_check(row, old_alpha)
                    new_row = row.replace(old_alpha, new_alpha).rstrip()
                    f.write(new_row)
                else:
                    new_alpha = new_params(row)
                    old_alpha = row.split()[-13]
                    self_check(row, old_alpha)
                    new_row = row.replace(old_alpha, new_alpha)
                    f.write(new_row)


