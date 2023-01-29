import numpy as np
import os
import math
from sympy import FU

ksi = 5e-17 # cr rate is 1/s
Av = 15
T = 10 # calculate uncertainties for T=10K
kb = 1.38e-16 # Boltzmann in cgs
a_dust = 1e-5 # grain radius in cm
amu = 1.66e-24 # atomic mass in g
cr_duration = 1e-5 # peak duration after CR hits dust
cr_T = 70 # peak temperature afte CR hits dust
Fe_ion = 3e-14
n_s = 1.5e15 # site density on grains in 1/cm2
N_s = n_s * 4 * np.pi * a_dust**2 # number of sites on 1 grain
nh2 = 1e5 # h2 density in 1/cm3, overall gas density basically
chem_barr_thick = 1e-8 # chemical barrier thickness in cm
hbar = 1.055e-27 # hbar in cgs
FUV = 1.
UVCR = ksi/1.3e-17
gotdn = 100.
eb_file = 'pnautius_own/example_simulation/surface_parameters.in'
br_f = 'branching_ratios.txt'

def get_stick(species):
    if species[-1] == '+' or species[-1] == '-':
        return 0
    else:
        return 1

def get_masses(species):
    mass = []
    for spec in species:
        ms = 0
        for i in range(spec):
            if spec[i] == 'H':
                if i+1 < len(spec) and spec[i+1] == 'e':
                    ms += 4
                else:
                    ms += 1
            elif spec[i] == 'C':
                if i+1 < len(spec) and spec[i+1] == 'l':
                    ms += 31
                else:
                    ms += 12
            elif spec[i] == 'N':
                if i+1 < len(spec) and spec[i+1] == 'a':
                    ms += 23
                else:
                    ms += 14
            elif spec[i] == 'O':
                ms += 16
            elif spec[i] == 'S':
                if i+1 < len(spec) and spec[i+1] == 'i':
                    ms += 28
                else:
                    ms += 32
            elif spec[i] == 'F':
                if i+1 < len(spec) and spec[i+1] == 'e':
                    ms += 56
                else:
                    ms += 19
            elif spec[i] == 'M':
                ms += 24
            elif spec[i] == 'P':
                ms += 31
        mass.append(ms)
    return mass

def get_Eb(spec):
    with open() as f:
        for row in f:
            if row.split()[0] == spec:
                val = row.split()
                eb = float(val[3])
                ed = float(val[2])
                if eb == 0 and spec[0] == 'J':
                    return ed / 0.4
                elif eb == 0 and spec[0] == 'K':
                    return ed / 0.8
                else:
                    return eb

def get_ED(spec):
    with open() as f:
        for row in f:
            if row.split()[0] == spec:
                return float(row.split()[2])

def vib_freq(spec, E_b, mass):
    return np.sqrt(2*kb/(np.pi**2*amu) * n_s * E_b/mass)

def calc_barr(species, Ea):
    m = get_masses(species) # get both reactant's masses
    redmas = m[0]*m[1]/np.sum(m)
    surf_react_proba = 2 * (chem_barr_thick/hbar) * np.sqrt(2*amu*redmas*kb*Ea)
    if Ea < surf_react_proba:
        return np.exp(-Ea/T)
    else:
        return np.exp(-surf_react_proba)

def calc_diff(species):
    m = get_masses(species)
    Eb1 = get_Eb(species[0])
    Eb2 = get_Eb(species[1])
    diff_bar1 = get_ED(species[0])
    diff_bar2 = get_ED(species[1])
    nu1 = vib_freq(species[0], Eb1, m[0])
    nu2 = vib_freq(species[1], Eb2, m[1])
    return (nu1*np.exp(-diff_bar1/T) + nu2*np.exp(-diff_bar2/T)) / N_s

def get_r_from_br(line):
    items = line.split()
    r_br = []
    for item in items:
        if item == '->':
            break
        else:
            r_br.append(item)
    return r_br

def comp_r(r1, r2):
    rcheck = False
    for item in r1:
        rcheck = item in r2
        if rcheck == False:
            break
    return rcheck

def rate(a, b, c, typ, reactants):
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
    elif typ == '14':
        with open(br_f) as B:
            for row in B:
                rea = get_r_from_br(row)
                if comp_r(reactants, rea) == True:
                    BR = float(row.split()[-2])
                    Ea = float(row.split()[-1])
        BARR = calc_barr(reactants, Ea)
        DIFF = calc_diff(reactants)
        return a * BR * BARR * DIFF * gotdn / nh2
    elif typ == '15':
        m = get_masses(reactants)
        Eb = get_Eb(reactants[0])
        nu = vib_freq(reactants[0], Eb, m[0])
        return a * 1 * nu * np.exp(-Eb/T)
        # branching ratio is assumet to be one (no info, should not change much in the lognoral dist)
    elif typ == '16':
        m = get_masses(reactants)
        Eb = get_Eb(reactants[0])
        nu = vib_freq(reactants[0], Eb, m[0])
        return a * 1 * (ksi/1e3-17) * nu * Fe_ion * cr_duration * np.exp(-Eb/cr_T)
    elif typ == '17' or typ == '18':
        return a * ksi * nh2 * 2
    elif typ == '19' or typ == '20':
        return a * np.exp(-c*Av) * FUV
    elif typ == '66':
        return a * FUV * 1e8 * np.exp(-2*Av) / n_s
    elif typ == '67':
        return a * 1e4 * UVCR / n_s
    elif typ == '99': # this is only if species is neutral
        m = get_masses(reactants)
        S = get_stick(reactants[0])
        if S == 1:
            pre = np.pi * a_dust**2 * np.sqrt(8*kb/(np.pi*amu)) / np.sqrt(m)
            return a * 1 * pre * np.sqrt(T) * nh2 / gotdn
            # assume branching ratio is one (did not find for these reactions)
        else:
            return a

def new_a(k, b, c, typ, reactants):
    if typ == '0':
        new = k / (T/300)**b
    elif typ == '1':
        new = k / ksi
    elif typ == '2':
        new = k / np.exp(-c*Av)
    elif typ == '3':
        new = k / ((T/300)**b * np.exp(-c/T))
    elif typ == '4':
        new = k / (b * (0.62 + 0.4767*c*np.sqrt(300/T)))
    elif typ == '5':
        new = k / (b * (1. + 0.0967*c*np.sqrt(300/T) + (c**2/10.526)*300/T))
    elif typ == '10':
        new = k
    elif typ == '11':
        new = k
    elif typ == '14':
        with open(br_f) as B:
            for row in B:
                rea = get_r_from_br(row)
                if comp_r(reactants, rea) == True:
                    BR = float(row.split()[-2])
                    Ea = float(row.split()[-1])
        BARR = calc_barr(reactants, Ea)
        DIFF = calc_diff(reactants)
        new = k / (BR * BARR * DIFF * gotdn / nh2)
    elif typ == '15':
        m = get_masses(reactants)
        Eb = get_Eb(reactants[0])
        nu = vib_freq(reactants[0], Eb, m[0])
        new = k / (1 * nu * np.exp(-Eb/T))
        # branching ratio is assumed to be 1
    elif typ == '16':
        m = get_masses(reactants)
        Eb = get_Eb(reactants[0])
        nu = vib_freq(reactants[0], Eb, m[0])
        new = k / (1 * (ksi/1e3-17) * nu * Fe_ion * cr_duration * np.exp(-Eb/cr_T))
    elif typ == '17' or typ == '18':
        new = k / (ksi * nh2 * 2)
    elif typ == '19' or typ == '20':
        new = k / (np.exp(-c*Av) * FUV)
    elif typ == '66':
        new = k / (FUV * 1e8 * np.exp(-2*Av) / n_s)
    elif typ == '67':
        new = k / (1e4 * UVCR / n_s)
    elif typ == '99': # this is only if species is neutral
        m = get_masses(reactants)
        S = get_stick(reactants[0])
        if S == 1:
            pre = np.pi * a_dust**2 * np.sqrt(8*kb/(np.pi*amu)) / np.sqrt(m)
            new = k / (1 * pre * np.sqrt(T) * nh2 / gotdn)
            # one is branching ration, did not findfor these reactions, so assume
        else:
            new = k

    if new == float('inf') or new == float('-inf') or math.isnan(new) == True: # just to make sure
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
    rrrr = line[:23].split()
    k0 = rate(alpha, beta, gamma, typ, rrrr)
    N = np.random.normal(0, 1)
    k = np.exp(np.log(k0) + np.log(F0)*N)
    return '{:.3e}'.format(new_a(k, beta, gamma, typ, rrrr))            

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


