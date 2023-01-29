#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import os
#%%
a = 7.7e-15
n0 = 3e3
Av0 = 2
# ----------Phase I-----------
def t_from_n(n, year = True):
    t_s = (3/a) * ( np.sqrt((n/n0)**(1/3)-1)/(n**(1/3)) + np.arctan(np.sqrt((n/n0)**(1/3)-1))/(n0**(1/3)) )
    if year == True:
        return t_s / (365.25*24*3600)
    else:
        return t_s

def Av_from_n(n):
    return Av0 * (n/n0)**(2/3)

def temp_dust(av):
    Temp_D = []
    for val in av:
        if val < 10:
            Temp_D.append(18.67 - 1.637*val + 0.07515*val**2 - 0.001492*val**3 + 0.316)
        else:
            Temp_D.append( ( (43*val**(-0.56)-77*val**(-1.28))**(5.6) + (7.9*val**(-0.089))**(5.6)*(1.8-0.089*val**(0.5)+7.9e-5*val**(1.5)) + (6.2-0.0031*val) )**(1/5.6) )
    return np.array(Temp_D) 

def add_zerotime(n, time, av, temp_d):
    n_new = np.concatenate((np.array([n0]), n))
    time_new = np.concatenate((np.array([0.]), time))
    temp_d_new = np.concatenate((np.array([temp_d[0]]), temp_d))
    av_new = np.concatenate((np.array([Av0]), av))
    return n_new, time_new, av_new, temp_d_new
#%%
#------------Phase II--------------
nf = 1e7
tc = 1.15e6 # contraction time in years for M
M = 5 # mass in solar masses
Mdot = M / tc
#t2 = np.logspace(0, np.log10(tc), 512)
#Msol = t2 * Mdot
Rsol = 2. # reasoable value? should it vary?
G = 6.67e-8 # grav const in cgs

def time_rate_ep(t0, tf, tq = 1e4, Int = 42):
    t_last = t0 + 100
    time = np.linspace(t0, t_last, 10)
    mrate = np.ones(np.shape(time)) * Mdot * Int
    stop = False
    while stop == False:
        t_last_new = np.min([t_last + tq, tf])
        if t_last_new == tf:
            stop = True
        t_new = np.linspace(t_last + 10, t_last_new, 100) # no time is dupicated this way
        time = np.concatenate((time, t_new))
        mrate = np.concatenate((mrate, np.ones(np.shape(t_new))*Mdot))
        t_last = t_last_new
        if t_last != tf:
            t_last_new = t_last + 110
            t_new = np.linspace(t_last + 10, t_last_new, 10) # still a 100 years long
            time = np.concatenate((time, t_new))
            mrate = np.concatenate((mrate, np.ones(np.shape(t_new))*Mdot*Int))
            t_last = t_last_new
    return time, mrate

def lum(mrate, time):
    Msol = time * Mdot # approcimation instead of integration, should be fine
    M_cgs = Msol * 1.99e33
    Mdot_cgs = mrate * 1.99e33 / (365.25*24*3600)
    R_cgs = Rsol * 6.96e10
    L_cgs = G * M_cgs * Mdot_cgs / R_cgs
    return L_cgs / 3.9e33 # returns in solar luminosity

def Temp_ep(temp_d0, lum_ep, T_lim = 400, tf = 200):
    const = (tf-temp_d0) / (tc**2)
    tim = np.concatenate((np.array([0.]),np.logspace(0, 7.7, 1000))) # by lazy trial and error
    T_int = temp_d0 + const*(tim**2)
    L_int = lum(Mdot, tim)
    T_L = interp1d(L_int, T_int, kind = 'quadratic')
    temp_ep = T_L(lum_ep)
    temp_ep[temp_ep > T_lim] = T_lim
    return temp_ep
#%%
#--------Write to file-----------
def write(folder, time, logAv, logn, logTg, logTd):
    with open(os.path.join(folder,'structure_evolution.dat'), mode = 'w') as f:
        f.write('! time\tlog(Av)\tlog(n)\tlog(Tg)\tlog(Td)\n')
        f.write('! (yr)\tlog(mag)\tlog(cm-3)\tlog(K)\tlog(K)\n')
        #f.write('{:.3e} {:.3e} {:.3e} {:.3e} {:.3e}\n'.format(0, logAv, logn, logTg[0], logTd[0]))
        for i in range(len(time)):
            f.write('{:.7e} {:.3e} {:.3e} {:.3e} {:.3e}\n'.format(time[i], logAv[i], logn[i], logTg[i], logTd[i]))

def write_ab(folder_1, folder_2):
    ab = os.path.join(folder_1, 'ab')
    f_ab = os.path.join(folder_2, 'abundances.in')
    with open(f_ab, mode = 'w') as f:
        for file in os.listdir(ab):
            spec = file[:-3]
            x = np.loadtxt(os.path.join(ab, file), comments = '!')[-1,-1]
            if x > 1e-40:
                x = str(x)
                x = x.replace('E', 'D')
                x = x.replace('e', 'D')
                row = spec.ljust(12) + '= ' + x + '\n'
                f.write(row)
#%%

#N = np.linspace(3e3, 1e7, 1000)
N_int = np.logspace(np.log10(n0), 7, 512)
t_int = np.concatenate((np.array([0.]), t_from_n(N_int)))
n_t = interp1d(t_int, np.concatenate((np.array([n0]), N_int)))
t_t_1 = np.linspace(0, t_int[-1]*5/6, 100)
t_t_2 = t_int[-1] + t_int[-1]*5/6 - np.logspace(np.log10(t_int[-1]), np.log10(t_int[-1]*5/6), 412)
N1_1 = n_t(t_t_1)
N1_2 = n_t(t_t_2)
t1 = np.concatenate((t_t_1, t_t_2))
N1 = np.concatenate((N1_1, N1_2))
Av1 = Av_from_n(N1)
Td1 = temp_dust(Av1)
N1, t1, Av1, Td1 = add_zerotime(N1, t1, Av1, Td1)
Tg1 = np.ones(np.shape(N1)) * 10 # gas is held at 10K in phase 1
#%%
t2, Mdot_ep = time_rate_ep(0, tc)
#L = lum(Mdot, t2)
L_ep = lum(Mdot_ep, t2)
Td2 = Temp_ep(Td1[-1], L_ep)
Tg2 = np.copy(Td2)
Tg2[Tg2 < 10] = 10. # they are only coupled onces Td reaches 10 K
N2 = np.ones(np.shape(t2)) * N1[-1]
Av2 = np.ones(np.shape(t2)) * Av1[-1]
'''t2 += t1[-1]
t = np.concatenate((t1, t2))
Av = np.concatenate((Av1, Av2))
N = np.concatenate((N1, N2))
Tg = np.concatenate((Tg1, Tg2))
Td = np.concatenate((Td1, Td2))
'''
#%%
write('/lustre/astro/gfriss/infall_1', t1, np.log10(Av1), np.log10(N1), np.log10(Tg1), np.log10(Td1))
write('/lustre/astro/gfriss/infall_2', t2, np.log10(Av2), np.log10(N2), np.log10(Tg2), np.log10(Td2))
#%%
write_ab('/lustre/astro/gfriss/infall_1', '/lustre/astro/gfriss/infall_2')
# %%
