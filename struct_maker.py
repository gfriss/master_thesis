import numpy as np
import os

def make_time_temp(t0, tf, tq, nq, T0, Tb, v = 1):
    if v == 1:
        t_last = t0 + 100
        time = np.linspace(t0, t_last, 10)
        logtemp = np.log10(Tb) * np.ones(np.shape(time))
    else:
        t_last = t0
        time = np.array([])
        logtemp = np.array([])
    for i in range(nq):
        t_last_new = np.min([t_last + tq, tf])
        t_new = np.linspace(t_last + 10, t_last_new, 100) # no time is dupicated this way
        logtemp_new = np.log10(T0) * np.ones(np.shape(t_new))
        time = np.concatenate((time, t_new))
        logtemp = np.concatenate((logtemp, logtemp_new))
        t_last = t_last_new
        if t_last != tf:
            t_last_new = t_last + 110
            t_new = np.linspace(t_last + 10, t_last_new, 10) # still a 100 years long
            logtemp_new = np.log10(Tb) * np.ones(np.shape(t_new))
            time = np.concatenate((time, t_new))
            t_last = t_last_new
            logtemp = np.concatenate((logtemp, logtemp_new))
        else:
            break
    return time, logtemp
        
def write(folder, logn, tq, logAv, Tb, t0 = 1e0, nq = 5, T0 = 10, v = 1, tf = None):
    if tf == None:
        t_f = tq*nq# + 5*100
        n_q = nq
    else:
        t_f = tf
        n_q = int(t_f // tq) + 1
    t, logT = make_time_temp(t0, t_f, tq, n_q, T0, Tb, v = v)
    with open(folder+'structure_evolution.dat', mode = 'w') as f:
        f.write('! time\tlog(Av)\tlog(n)\tlog(T)\n')
        f.write('! (yr)\tlog(mag)\tlog(cm-3)\tlog(K)\n')
        f.write('{:.3e} {:.3e} {:.3e} {:.3e}\n'.format(0, logAv, logn, logT[0]))
        for i in range(len(t)):
            f.write('{:.3e} {:.3e} {:.3e} {:.3e}\n'.format(t[i], logAv, logn, logT[i]))

def make_w_pre(t_pre, tq, nq, T0, Tb, Tpre, n_pre = 64):
    t_last = t_pre
    time = np.linspace(0, t_pre, n_pre)
    time = np.logspace(-1, np.log10(t_pre), n_pre)
    logtemp = np.ones(np.shape(time)) * np.log10(Tpre)
    tf = t_pre + (tq + 100)*nq
    for i in range(nq):
        if t_last != tf:
            t_last_new = t_last + 110
            t_new = np.linspace(t_last + 10, t_last_new, 10) # still a 100 years long
            logtemp_new = np.log10(Tb) * np.ones(np.shape(t_new))
            time = np.concatenate((time, t_new))
            t_last = t_last_new
            logtemp = np.concatenate((logtemp, logtemp_new))
        else:
            break
        t_last_new = np.min([t_last + tq, tf])
        t_new = np.linspace(t_last + 10, t_last_new, 100) # no time is dupicated this way
        logtemp_new = np.log10(T0) * np.ones(np.shape(t_new))
        time = np.concatenate((time, t_new))
        logtemp = np.concatenate((logtemp, logtemp_new))
        t_last = t_last_new
        
    return time, logtemp

def write_w_pre(folder, t_pre, Tpre = 10, n_pre = 64, tq = 1e4, nq = 5, T0 = 10, Tb = 30, logn_pre = 4, logn = 5, logAv = np.log10(15)):
    t, logT = make_w_pre(t_pre, tq, nq, T0, Tb, Tpre)
    with open(os.path.join(folder, 'structure_evolution.dat'), mode = 'w') as f:
        f.write('! time\tlog(Av)\tlog(n)\tlog(T)\n')
        f.write('! (yr)\tlog(mag)\tlog(cm-3)\tlog(K)\n')
        f.write('{:.3e} {:.3e} {:.3e} {:.3e}\n'.format(0, logAv, logn, logT[0]))
        for i in range(len(t)):
            if i < n_pre: # prestellar phase has lower density
                f.write('{:.3e} {:.3e} {:.3e} {:.3e}\n'.format(t[i], logAv, logn_pre, logT[i]))
            else:
                f.write('{:.3e} {:.3e} {:.3e} {:.3e}\n'.format(t[i], logAv, logn, logT[i]))
