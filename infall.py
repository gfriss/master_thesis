#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.integrate import cumtrapz
import os
import plot_reset as pr
#%%
au_to_cm = 1.49597871e13
yr_to_s = 365.24*24*3600
M = 1. * 1.99e33 # envelope mass in g
r0 = 1e4 * au_to_cm # starting radius in cm (10000 AU)
rf = 10. * au_to_cm # final radius in cm (10 AU)
G = 6.674e-8 # grav const in cgs
m_h = 1.67e-24 # hydrogen mass in g

def make_time_T0(t0, tf, tq = 1e4, Tq = 10., Tb = 30.):
    ''' Prepares the time and temperature porfiles. The latter relies on Visser: 
        a factor of 240 change in L changes T0 from 10 K to 30 K'''
    t_last = t0 + 100
    time = np.linspace(t0, t_last, 10)
    temp_0 = np.ones(np.shape(time)) * Tb
    stop = False
    while stop == False:
        t_last_new = np.min([t_last + tq, tf])
        if t_last_new == tf:
            stop = True
        t_new = np.linspace(t_last + 10, t_last_new, 1000) # no time is dupicated this way
        time = np.concatenate((time, t_new))
        temp_0_new = np.ones(np.shape(t_new)) * Tq
        temp_0 = np.concatenate((temp_0, temp_0_new))
        t_last = t_last_new
        if t_last != tf:
            t_last_new = t_last + 110
            t_new = np.linspace(t_last + 10, t_last_new, 10) # still a 100 years long
            time = np.concatenate((time, t_new))
            temp_0_new = np.ones(np.shape(t_new)) * Tb
            temp_0 = np.concatenate((temp_0, temp_0_new))
            t_last = t_last_new
    return time, temp_0

def calc_n0(R0, Rf, p, M_env):
    ''' using the results of manual integration of the mass that gives M_env
        uses cgs!!!'''
    prefactor = 4 * np.pi * m_h / (R0**(-p) * (3-p))
    upper = R0**(3-p)
    lower = Rf**(3-p)
    return M_env / (prefactor * (upper - lower))

def rad_t(time, R0 = r0, Rf = rf, p = 1.5, M_env = M):
    ''' Gives the radius as a functionof time. Derived by inverting the isothermal-sphere collapse time,
        e.g. free-fall time with changing density. CGS again!!'''
    N0 = calc_n0(R0, Rf, p, M_env)
    prefactor = (32*G*m_h*N0 / (3*np.pi))**(1/p)
    return R0 - R0 * prefactor * (time**(2/p))

def make_n(R, R0 = r0, Rf = rf, p = 1.5, M_env = M):
    n0 = calc_n0(R0, Rf, p, M_env)
    return n0 * (R/R0)**(-p)

def temp(R, R0, T0):
    return T0 * (R/R0)**(-0.4)

def vis_ext(N, R, R0 = r0, Av0 = 10):
    N_col = -cumtrapz(N, x = R) # negative to adjust for the decreasing distance
    vis = Av0 + N_col / 2.21e21
    return np.concatenate((np.array([Av0]), vis))
    #return Av0 + (N*(R0-R) / 2.2e21)

#%%
#--------Write to file-----------
def write(folder, time, logAv, logn, logTg):
    with open(os.path.join(folder,'structure_evolution.dat'), mode = 'w') as f:
        f.write('! time\tlog(Av)\tlog(n)\tlog(Tg)\n')
        f.write('! (yr)\tlog(mag)\tlog(cm-3)\tlog(K)\n')
        #f.write('{:.3e} {:.3e} {:.3e} {:.3e} {:.3e}\n'.format(0, logAv, logn, logTg[0], logTd[0]))
        for i in range(len(time)):
            f.write('{:.7e} {:.3e} {:.3e} {:.3e}\n'.format(time[i], logAv[i], logn[i], logTg[i]))
#%%

t, T0 = make_time_T0(0, 2.4978e5) # enf time by trial to go down until ~10 AU
r = rad_t(t*yr_to_s)# / au_to_cm # radius in au calculated using time in s
n = make_n(r)
T = temp(r, r0, T0)
Av = vis_ext(n, r)
r = r/au_to_cm

pr.reset_plt(20,20)
fig, ax = plt.subplots(ncols = 2, nrows = 2, sharex = True, figsize = (14,12), tight_layout = True)
ax = ax.flatten()
ax[0].semilogy(t+1e6, r, 'r')
ax[0].set_ylabel('r[AU]')
ax[1].semilogy(t+1e6, n, 'r')
ax[1].set_ylabel(r'n[cm$^{-3}$]')
ax[2].plot(t+1e6, T, 'r')
ax[2].set_xlabel('t[yr]')
ax[2].set_ylabel('T[K]')
ax[3].plot(t+1e6, Av, 'r')
ax[3].set_xlabel('t[yr]')
ax[3].set_ylabel(r'A$_V$')
fig.savefig('infall_struct.pdf')

write('/lustre/astro/gfriss/infall_own', t, np.log10(Av), np.log10(n), np.log10(T))
# %%
t_pre_final = 1e6
t_pre = np.concatenate((np.array([0.]), np.geomspace(1, t_pre_final, 1000)))
n_pre = np.ones(np.shape(t_pre)) * calc_n0(r0, rf, 1.5, M)
Av_pre = np.ones(np.shape(t_pre)) * 10
T_pre = np.ones(np.shape(t_pre)) * 10

t = np.concatenate((t_pre[:-1], t_pre[-1]+t)) # continuity
n = np.concatenate((n_pre[:-1], n))
Av = np.concatenate((Av_pre[:-1], Av))
T = np.concatenate((T_pre[:-1], T))

write('/lustre/astro/gfriss/infall_own_with_pre', t, np.log10(Av), np.log10(n), np.log10(T))

#%%
#------------no burst---------------
t, T0 = make_time_T0(0, 2.4978e5, Tb = 10.) # enf time by trial to go down until ~10 AU
r = rad_t(t*yr_to_s)# / au_to_cm # radius in au calculated using time in s
n = make_n(r)
T = temp(r, r0, T0)
Av = vis_ext(n, r)
r = r/au_to_cm

write('/lustre/astro/gfriss/noburst_infall_own', t, np.log10(Av), np.log10(n), np.log10(T))
# %%
t_pre_final = 1e6
t_pre = np.concatenate((np.array([0.]), np.geomspace(1, t_pre_final, 1000)))
n_pre = np.ones(np.shape(t_pre)) * calc_n0(r0, rf, 1.5, M)
Av_pre = np.ones(np.shape(t_pre)) * 10
T_pre = np.ones(np.shape(t_pre)) * 10

t = np.concatenate((t_pre[:-1], t_pre[-1]+t)) # continuity
n = np.concatenate((n_pre[:-1], n))
Av = np.concatenate((Av_pre[:-1], Av))
T = np.concatenate((T_pre[:-1], T))

write('/lustre/astro/gfriss/noburst_infall_own_with_pre', t, np.log10(Av), np.log10(n), np.log10(T))

# %%
