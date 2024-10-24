import numpy as np 

from scipy import interpolate
from scipy.optimize import brentq


import matplotlib
matplotlib.use('TKAgg')

import matplotlib.pyplot as plt

import functools
from scipy.integrate import quad
import streamlit as st

import itertools

#https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]#, idx

def _sLor(f, H, b, c):
    return H / (1 + (f / b) ** c)

def _sinc(x):
    return np.sinc(x/np.pi)

def bgModel(nu, theta, nuNyq, n_comps, n_gauss, individual=True):
    """
    Background model value at a given frequency 'nu'
    """
    white_noise = theta["Pn"].values
    sc = _sinc(np.pi * nu / (2 * nuNyq)) ** 2
    
    comps = []
    comp_names = []
    model = np.zeros(len(nu))
    model_no_osc = np.zeros(len(nu))
    if "H1" in theta.columns:
        height_var = "H"
    elif "A1" in theta.columns:
        height_var = "A"
    for i in range(n_comps):
        H = theta[height_var+f"{i+1}"].values
        b = theta[f"b{i+1}"].values
        if f'c{i+1}' in theta.columns:
            c = theta[f"c{i+1}"].values
        else:
            c = 4
        comp = _sLor(nu, H, b, c)
        model += comp
        model_no_osc += comp
        comps.append(comp)  
        comp_names.append(f'gran_comp_{i+1}')  


    for i in range(n_gauss):
        Henv = theta[f"Pg{i+1}"].values if i > 0 else theta["Pg"].values
        numax = theta[f"numax{i+1}"].values if i > 0 else theta["numax"].values
        sigmaEnv = theta[f"sigmaEnv{i+1}"].values if i > 0 else theta["sigmaEnv"].values
        gauss = Henv * np.exp(-((nu - numax) ** 2) / (2 * sigmaEnv ** 2))
        model += gauss
        comps.append(gauss)
        comp_names.append(f'gauss_comp_{i+1}')  

    model *= sc
    model_no_osc *= sc
    model += white_noise
    model_no_osc += white_noise
    comps.append(white_noise*np.ones(len(nu)))
    comp_names.append("white_noise")

    if individual == True:
        return comps, comp_names, model, model_no_osc
    return model, model_no_osc

def lorentzian(f, amp, lw, freq):
    height = 2*amp**2/(np.pi*2*lw)
    x = (2/lw) * (f - freq)
    return height / (1 + x**2)

def sinc_sq(f, amp, lw, freq):
    height = amp**2 / (np.pi * (f[1]-f[0]))
    return height*np.sinc((f-freq)/(f[1]-f[0]))**2

def model(f, row):
    if np.isfinite(row['linewidth']):
        return lorentzian(f, row['amplitude'], row['linewidth'], row['frequency'])
    else:
        #print(row['linewidth'])
        return sinc_sq(f, row['amplitude'], row['linewidth'], row['frequency'])

def construct_peaksmodel(pds, peaks):
    fit_model = np.ones(len(pds))
    for idx, i in peaks.iterrows():
        #plt.plot(pds['frequency'], model(pds['frequency'], i))
        #plt.show()
        #fit_model02 += model(pds['frequency'], i)
        fit_model += model(pds['frequency'].values, i)

    return fit_model

def construct_MLEmodel(pds, peaks):
    fit_model02 = np.zeros(len(pds))
    fit_model1 = np.zeros(len(pds))
    for idx, i in peaks.loc[(peaks['l'] == 0) | (peaks['l'] == 2), ].iterrows():
        #plt.plot(pds['frequency'], model(pds['frequency'], i))
        #plt.show()
        #fit_model02 += model(pds['frequency'], i)
        fit_model02 += model(pds['frequency'].values, i)
        
    for idx, i in peaks.loc[(peaks['l'] != 0) & (peaks['l'] != 2), ].iterrows():
        #plt.plot(pds['frequency'], model(pds['frequency'], i))
        #plt.show()
        fit_model1 += model(pds['frequency'].values, i)
    return fit_model02, fit_model1

def echelle(freq, power, dnu, fmin=0., fmax=None, offset=0.0):
    # This is a slightly modified version of Dan Hey's fantastic echelle code
    # https://github.com/danhey/echelle/blob/master/echelle/echelle.py
    """Calculates the echelle diagram. Use this function if you want to do
    some more custom plotting.
    
    Parameters
    ----------
    freq : array-like
        Frequency values
    power : array-like
        Power values for every frequency
    dnu : float
        Value of deltanu
    fmin : float, optional
        Minimum frequency to calculate the echelle at, by default 0.
    fmax : float, optional
        Maximum frequency to calculate the echelle at. If none is supplied, 
        will default to the maximum frequency passed in `freq`, by default None
    offset : float, optional
        An offset to apply to the echelle diagram, by default 0.0
    
    Returns
    -------
    array-like
        The x, y, and z values of the echelle diagram.
    """
    if fmax is None:
        fmax = freq[-1]

    fmin = fmin - offset
    fmax = fmax - offset
    freq = freq - offset

    if fmin <= 0.0:
        fmin = 0.0
    else:
        fmin = fmin - (fmin % dnu)

    # trim data
    index = (freq>=fmin) & (freq<=fmax)
    trimx = freq[index]

    samplinginterval = np.median(trimx[1:-1] - trimx[0:-2])# * 0.1
    xp = np.arange(fmin,fmax+dnu,samplinginterval)
    yp = np.interp(xp, freq, power)

    n_stack = int(np.ceil((fmax-fmin)/dnu))
    n_element = int(np.ceil(dnu/samplinginterval))

    morerow = 2
    arr = np.arange(1,n_stack) * dnu
    arr2 = np.array([arr,arr])
    yn = np.reshape(arr2,len(arr)*2,order="F")
    yn = np.insert(yn,0,0.0)
    yn = np.append(yn,n_stack*dnu) + fmin + offset

    xn = np.arange(1,n_element+1)/n_element * dnu
    #print(yn)
    z = np.zeros([n_stack*morerow,n_element])
    for i in range(n_stack):
        for j in range(i*morerow,(i+1)*morerow):
            z[j,:] = yp[n_element*(i):n_element*(i+1)]
    return xn, yn, z

def l0_from_UP(N, eps_p, alpha_, n_max, DeltaNu):
    # Theoretical radial mode frequencies from Universal Pattern
    return ((N + eps_p + (alpha_/2) * (N - n_max)**2) * DeltaNu)

def l1_nominal_p_freqs(freqs_zero, deltanu, d1=None):
    if d1 is not None:
        return freqs_zero + deltanu/2 + d1
    else:
        #d1 = 0.0553 - 0.036*np.log(deltanu)
        d1 =  -0.056 + -0.002*np.log10(deltanu)
        return freqs_zero + deltanu/2 + d1

def all_mixed_l1_freqs(delta_nu, nu_zero, nu_p, DPi1, eps_g, 
                       coupling, return_order=True, calc_zeta=True):

    l1_freqs = []
    l1_g_freqs = []
    zeta = []
    order = []
    N_g = []

    for i in range(len(nu_p)):

        if nu_p[i] > nu_zero[-1]:
            radial = np.array([nu_zero[-1], nu_zero[-1] + delta_nu])
        else:
            radial = np.array([nu_zero[i], nu_zero[i+1]])
            

        tmp, tmp_g, tmp_ng = find_mixed_l1_freqs(delta_nu, radial, nu_p[i], 
                                                 DPi1, eps_g, coupling)
        if calc_zeta == True:
            tmp_zeta = zeta_Mosser(tmp, nu_p[i], delta_nu, DPi1, coupling, eps_g)
            zeta.append(tmp_zeta)
        order.append([i]*len(tmp))
        l1_freqs.append(tmp)
        l1_g_freqs.append(tmp_g)
        N_g.append(tmp_ng)

    if calc_zeta and return_order:
        return np.array(list(itertools.chain(*l1_freqs))), \
               np.array(list(itertools.chain(*zeta))), \
               np.array(list(itertools.chain(*order))), \
               np.array(list(itertools.chain(*l1_g_freqs))), \
               np.array(list(itertools.chain(*N_g)))
    elif calc_zeta:
        return np.array(list(itertools.chain(*l1_freqs))), \
               np.array(list(itertools.chain(*zeta)))
    else:
        return np.array(list(itertools.chain(*l1_freqs)))


def find_mixed_l1_freqs(delta_nu, nu_zero, nu_p, DPi1, eps_g, coupling):
    """
    Find all mixed modes in a given radial order
    """

    nmin = np.floor(1 / (DPi1*1e-6 * nu_zero[1]) - (1/2) - eps_g)
    nmax = np.floor(1 / (DPi1*1e-6 * nu_zero[0]) - (1/2) - eps_g)

    N_modes = (delta_nu * 1e-6) / (DPi1 * (nu_p*1e-6)**2)
    
    # Changed to this to be exactly the same as Mosser 2018

    # For some reason have to do this to ensure that find all solutions!
    N = np.arange(nmin, nmax + 2, 1)
  
    frequencies, g_mode_freqs, N_g = find_mixed_l1_freq(delta_nu, nu_zero, nu_p, 
                                                        DPi1, eps_g, coupling, N)

    return np.sort(frequencies[np.isfinite(frequencies)]), np.sort(g_mode_freqs[np.isfinite(g_mode_freqs)]), np.sort(N_g[np.isfinite(N_g)])


def find_mixed_l1_freq(delta_nu, pzero, pone, DPi1, eps_g, coupling, N):
    """
    Find individual mixed mode
    """
    def opt_funcM(nu, nu_g):
        theta_p = (np.pi / (pzero[1]-pzero[0])) * (nu - pone)
        theta_g = np.pi/DPi1 * 1e6 * (1/nu - 1/nu_g) + np.pi/2
        y = np.tan(theta_p) - coupling * np.tan(theta_g)
        return y

    nu_g = 1 / (DPi1*1e-6 * (N + 1/2 + eps_g))
    # Go +/- 1/2 * DPi1 away from g-mode period
    lower_bound = 1 / (DPi1*1e-6 * (N     + 1/2 + 1/2 + eps_g)) + np.finfo(float).eps * 1e4
    upper_bound = 1 / (DPi1*1e-6 * (N - 1 + 1/2 + 1/2 + eps_g)) - np.finfo(float).eps * 1e4

    f = np.linspace(pzero[0], pzero[1], 10000)#[1:-1]
    dnu = np.diff(pzero)
    solns = []
    for i in range(len(nu_g)):
        if (upper_bound[i] > pzero[1]):
            upper_bound[i] = pzero[1]
        elif (lower_bound[i] < pzero[0]):
            lower_bound[i] = pzero[0]
        if (upper_bound[i] < lower_bound[i]) or (lower_bound[i] > upper_bound[i]):
            pass
        else:
            ff = np.linspace(lower_bound[i], upper_bound[i], 1000)
            y = opt_funcM(ff, nu_g[i])
    
            idx = np.where(np.diff(np.sign(y)) > 0)[0]
            if len(idx) == 0:
                soln = np.array([])
            else:
                soln = ff[idx]

            solns = np.append(solns, soln)
            theta_p = (np.pi / (pzero[1]-pzero[0])) * (solns - pone)
            # Approximate pure g-mode frequencies and radial orders
            g_period = 1/(solns*1e-6) - DPi1/np.pi * np.arctan2(np.tan(theta_p), coupling) 
            n_g = np.floor(g_period / DPi1 - eps_g - 1/2)
        
    return solns, 1e6/g_period, n_g
    
"""
def find_mixed_l1_freqs(DeltaNu, nu_p, DPi1, eps_g, coupling):
    #l1_freqs = []
    #for i in range(len(freq_zero)):
    #    l1_freqs.append(find_mixed_l1_freq(freq, delta_nu, freq_zero[i], 
    #                                                      period_spacing, 
    #                                                      epsilon_g, coupling))
    #return np.array(l1_freqs)

    nmin = np.ceil(1 / (DPi1*1e-6 * (nu_p + (DeltaNu/2))) - (1/2) - eps_g)
    nmax = np.ceil(1 / (DPi1*1e-6 * (nu_p - (DeltaNu/2))) - (1/2) - eps_g)
    #nmin -= 10
    nmax += 2
    #st.write(nmin, nmax)
    frequencies = []
    for i in np.arange(nmin, nmax, 1):
        tmp = find_mixed_l1_freq(DeltaNu, nu_p, DPi1, eps_g, coupling, i)
        frequencies = np.append(frequencies, tmp)
#    print(frequencies)
    return np.sort(frequencies[np.isfinite(frequencies)])

def find_mixed_l1_freqs_alt(DeltaNu, nu_p, DPi1, eps_g, coupling):
    # Calculate the mixed mode frequencies ...
    nu = np.arange(nu_p - DeltaNu, nu_p + 1 * DeltaNu, 1e-5)
    nu *= 1e-6
    lhs = np.pi * (nu - (nu_p * 1e-6)) / (DeltaNu*1e-6)
    rhs = np.arctan(coupling * np.tan(np.pi/(DPi1 * nu) - eps_g))
    mixed1 = np.zeros(100)
    counter = 0
    for i in np.arange(0, nu.size-1):
        if lhs[i] - rhs[i] < 0 and lhs[i+1] - rhs[i+1] > 0:
            mixed1[counter] = nu[i]
            counter += 1
    mixed1 = mixed1[:counter]
    # add in the rotational splitting ...
    mixed1 *= 1e6
    return mixed1

def all_mixed_l1_freqs(DeltaNu, nu_p, DPi1, eps_g, coupling):
    #l1_freqs = []
    #for i in range(len(freq_zero)):
    #    l1_freqs.append(find_mixed_l1_freq(freq, delta_nu, freq_zero[i], 
    #                                                      period_spacing, 
    #                                                      epsilon_g, coupling))
    #return np.array(l1_freqs)

    l1_freqs = []
    for i in range(len(nu_p)):
        tmp = find_mixed_l1_freqs(DeltaNu, nu_p[i], 
                                  DPi1, eps_g, coupling)
        l1_freqs.append(tmp)

    return np.array(list(itertools.chain(*l1_freqs)))
"""

def zeta_Mosser(freq, nu_p, DeltaNu, DPi1, coupling, eps_g):
    # Mosser et al. 2018: A&A, AA/2018/32777
    N = DeltaNu / (DPi1*1e-6 * freq**2)

    theta_p = (np.pi / DeltaNu) * (freq - nu_p)
    b = 1 + (coupling / N) / ((coupling**2) * np.cos(theta_p)**2 + np.sin(theta_p)**2)
    res = 1/b

    return res

def zeta_Deheuvels(freq, nu_p, DeltaNu, DPi1, coupling, eps_g):
    # Deheuvels et al. (2015) <http://dx.doi.org/10.1051/0004-6361/201526449>
    a1 = np.cos(np.pi * ((1 / (freq * DPi1*1e-6)) - eps_g))**2
    a2 = np.cos(np.pi * ((freq - nu_p) / DeltaNu))**2
    a3 = (freq**2 * DPi1*1e-6) / (coupling * DeltaNu)
    b = 1 + a1 * a3 / a2
    return 1/b

def calc_zeta(freq_zero, nu_p, DeltaNu, DPi1, coupling, eps_g):
    zeta = []
    #l1_freqs = all_mixed_l1_freqs(DeltaNu, nu_p,
    #                              DPi1, eps_g, coupling)
    l1_freqs = []
    for i in range(len(nu_p)):
        if nu_p[i] != nu_p[-1]:
            radial = [freq_zero[i], freq_zero[i+1]]
        else:
            radial = [freq_zero[-1], freq_zero[-1]+DeltaNu]
        tmp_l1_freqs, _, _ = find_mixed_l1_freqs(DeltaNu, radial, nu_p[i], DPi1, eps_g, coupling)
        zeta = np.append(zeta, zeta_Mosser(tmp_l1_freqs, nu_p[i], DeltaNu, DPi1,
                                coupling, eps_g))
        l1_freqs = np.append(l1_freqs, tmp_l1_freqs)
    return l1_freqs, zeta

#@st.cache
def zeta_interp(freq, freq_zero, nu_p, DeltaNu, 
                DPi1, coupling, eps_g,
                numDPi1=100, DPi1_range=[0.99, 1.01],  return_func=True):
    # Interpolate zeta function
    l1_freqs = []
    zeta = []
    DPi1_vals = np.linspace(DPi1_range[0]*DPi1, DPi1_range[1]*DPi1, numDPi1)

    for i in range(len(DPi1_vals)):
        #print(DPi1_vals[i])
        tmp_l1_freqs, tmp_zeta = calc_zeta(freq_zero, nu_p, DeltaNu,
                                           DPi1_vals[i], coupling, eps_g)
        l1_freqs = np.append(l1_freqs, tmp_l1_freqs)
        zeta = np.append(zeta, tmp_zeta)

    l1_freqs = l1_freqs.ravel()
    zeta = zeta.ravel()

    idx = np.argsort(l1_freqs)
    l1_freqs = l1_freqs[idx]
    zeta = zeta[idx]

    zeta_fun = interpolate.interp1d(l1_freqs, zeta, bounds_error=False)
    #interp_zeta = np.interp(freq, l1_freqs, zeta)
    interp_zeta = zeta_fun(freq)
    if return_func:
        return l1_freqs, interp_zeta, zeta_fun
    return l1_freqs, interp_zeta

def stretched_pds(pds, freq_zero, DeltaNu, 
                  DPi1, coupling, eps_g, 
                  numDPi1=100, DPi1_range=[0.99, 1.01], oversample=1):
    # Compute frequency bin-width
    bw = pds.frequency[1]-pds.frequency[0]
    # Compute interpolated zeta across entire frequency range
    nom_l1_freqs = l1_nominal_p_freqs(freq_zero, DeltaNu)
    l1_freqs, zz, zeta_fun = zeta_interp(pds.frequency.values, freq_zero, 
                                    nom_l1_freqs, DeltaNu,
                                    DPi1, coupling, eps_g, numDPi1, DPi1_range)
    # Compute dtau
    if oversample > 1:
        new_freq = np.arange(pds.frequency.min(), pds.frequency.max(), bw/oversample)
        #dtau = 1 / (zz*(pds.frequency.values)**2) * 1e6
        dtau = 1 / (zeta_fun(new_freq)*(new_freq)**2) * 1e6
    else:
        new_freq = pds.frequency.values
        dtau = 1 / (zeta_fun(new_freq)*(new_freq)**2) * 1e6
   
    dtau[np.isnan(dtau)] = 0
    # Compute tau
    tau = -np.cumsum(dtau)*(bw/oversample)
    tau -= np.min(tau)

    #plt.plot(pds.frequency, dtau)
    #plt.plot(pds.frequency, tau)
    #plt.show()

    # Place tau into seconds
    #tau *= 1e6
    # Compute tau values of l1 frequencies to shift tau
    #l1_tau = np.interp(l1_freqs, pds.frequency.values, tau)
    l1_tau = np.interp(l1_freqs, new_freq, tau)
    l1_x = ((l1_tau + DPi1/2) % DPi1) - DPi1/2
    #l1_x = l1_tau % DPi1
    tau_shift = np.median(l1_x)
    #st.write(tau_shift)
    #l1_tau = l1_tau - tau_shift + DPi1
    # Compute l1 zeta
    #l1_zeta = np.interp(l1_freqs, pds.frequency.values, zz)

    tau = tau - tau_shift + DPi1
    return new_freq, tau, zeta_fun#, pds.power

def peaks_stretched_period(frequency, pds_frequency, tau):
    assert len(tau) == len(pds_frequency)
    return np.interp(frequency, pds_frequency, tau)

def l1_rot_from_zeta(pds, nu_0, nu_m, drot, zeta_fun):
    
    # Upper and lower limits for integration
    # Minimum value of nu_0 and nu_m or maximum value
    lower_limit = nu_0 if nu_0 < nu_m else nu_m
    upper_limit = nu_0 if nu_0 > nu_m else nu_m

    # Integrate zeta over that range
    res = quad(zeta_fun, lower_limit, upper_limit)
    int_zeta = res[0] / (nu_m - nu_0)
    return nu_0 + drot * int_zeta

def l1_rot_from_zeta_iter(pds, nu_0, nu_m, drot, zeta_fun, tol, max_iters=50, curr_iter=1):
    # Compute rotational splitting iteratively
    # If no rotational splitting then return nu_m, or nu_0?
    #st.write(curr_iter, max_iters)
    #st.write(drot)
    if drot == 0:
        return nu_m
    if curr_iter >= max_iters:
        st.write("Maximum number of iterations reached without convergence")
        return nu_m
    
    nu_m_new = l1_rot_from_zeta(pds, nu_0, nu_m, drot, zeta_fun)

    #st.write(abs(nu_m_new - nu_m))
    if abs(nu_m_new - nu_m) < tol:
        return nu_m_new
    else:
        return l1_rot_from_zeta_iter(pds, nu_0, nu_m_new, drot, zeta_fun,
                                     tol, max_iters, curr_iter+1)

def l1_theoretical_rot_M(pds, l1_m0_freqs, drot, zeta_fun, max_iters=50, tol=1e-4):

    l_mp1_freqs = []
    l_mn1_freqs = []
    for i in range(len(l1_m0_freqs)):
        tmp_p1 = l1_rot_from_zeta_iter(pds, l1_m0_freqs[i], l1_m0_freqs[i]+drot,
                              drot, zeta_fun, tol, max_iters)
        tmp_n1 = l1_rot_from_zeta_iter(pds, l1_m0_freqs[i], l1_m0_freqs[i]-drot,
                              drot, zeta_fun, tol, max_iters)
        l_mp1_freqs = np.append(l_mp1_freqs, tmp_p1)
        l_mn1_freqs = np.append(l_mn1_freqs, tmp_n1)

        #st.write(l1_m0_freqs[i], tmp_p1, tmp_p1-l1_m0_freqs[i])
        #st.write(l1_m0_freqs[i], tmp_n1, l1_m0_freqs[i]-tmp_n1)
        #sys.exit()
    return l_mp1_freqs, l_mn1_freqs

#https://stackoverflow.com/questions/29166353/how-do-you-add-error-bars-to-bokeh-plots-in-python
def errorbar(fig, x, y, xerr=None, yerr=None, color='red', 
             point_kwargs={}, error_kwargs={}):

  fig.circle(x, y, color=color, **point_kwargs)

  if np.all(xerr != None):
      x_err_x = []
      x_err_y = []
      for px, py, err in zip(x, y, xerr):
          x_err_x.append((px - err, px + err))
          x_err_y.append((py, py))
      fig.multi_line(x_err_x, x_err_y, color=color, **error_kwargs)

  if np.all(yerr != None):
      y_err_x = []
      y_err_y = []
      for px, py, err in zip(x, y, yerr):
          y_err_x.append((px, px))
          y_err_y.append((py - err, py + err))
      fig.multi_line(y_err_x, y_err_y, color=color, **error_kwargs)

def visualise_stretched_echelle(psd_bgr, peaks, summary, session):

    #st.write(psd_bgr)
    
    from sloscillations import frequencies, mixed_modes_utils
    from taco.rotation import rotation_utils

    orig_DPi1 = summary['DeltaPi1'].item()
    orig_numax = summary['numax'].item()
    orig_eps_p = summary['eps_p'].item()
    orig_alpha = summary['alpha'].item()
    #'', summary.columns
    orig_DeltaNu = summary['DeltaNu'].item()
    orig_coupling = summary['coupling'].item()
    orig_eps_g = summary['eps_g'].item()

    orig_splitting = 0.0

    # Set up sliders
    reset = st.sidebar.button('Reset parameters')
    if reset:
        session.run_id += 1
    save = st.sidebar.button('Save parameters')
    if save:
        st.write('Doesn\'t do anything at the moment, but will do soon!')
    
    DeltaNu = st.sidebar.slider('Δν (μHz)', min_value=orig_DeltaNu*0.9, max_value=orig_DeltaNu*1.1, value=orig_DeltaNu, key=session.run_id)
    DPi1 = st.sidebar.slider('ΔΠ₁ (s)', min_value=orig_DPi1*0.9, max_value=orig_DPi1*1.1, value=orig_DPi1, key=session.run_id)
    coupling = st.sidebar.slider('q', min_value=0.0, max_value=0.8, value=orig_coupling, key=session.run_id)
    eps_g = st.sidebar.slider('ε', min_value=-1.0, max_value=1.0, value=float(orig_eps_g), key=session.run_id)
    splitting = st.sidebar.slider('δνₛ (μHz)', min_value=0., max_value=500.0/1e3, value=orig_splitting, key=session.run_id)

    peaks['x'] = ((peaks['frequency'] % DeltaNu - summary['eps_p'].values) / summary['DeltaNu'].values) % 1
    #st.write(peaks.loc[peaks['l'] == 0, ])
    #st.write([(np.minimum(np.min(peaks.loc[peaks['l'] == 0, 'x']), np.min(peaks.loc[peaks['l'] == 2, 'x'])) - 0.05) % 1,
    #           (np.maximum(np.max(peaks.loc[peaks['l'] == 0, 'x']), np.max(peaks.loc[peaks['l'] == 2, 'x'])) + 0.05) % 1]
    #)

    l1_peaks = rotation_utils.prepare_l1_peaks(peaks, summary)

    # Create frequencies instance
    freqs = frequencies.Frequencies(frequency=psd_bgr.frequency.values,
                            numax=orig_numax, 
                            delta_nu=DeltaNu, 
                            epsilon_p=orig_eps_p,
                            alpha=orig_alpha)
    
    params = {'calc_l0': True, # Compute radial mode properties
            'calc_l2': True, # Compute l=2 mode properties
            'calc_l3': False, # Don't need to calculate l=3 theoretical freqs
            'calc_nom_l1': True, # Compute nominal l=1 p-mode properties
            'calc_mixed': True, # Don't compute mixed modes (as not needed)
            'calc_rot': True, # Don't compute rotation
            'DPi1': DPi1,
            'coupling': coupling,
            'eps_g': eps_g,
            'split_core': splitting,
            'split_env': 0.0,
            'l': 1, # Mixed modes are dipole mixed modes
            }

    # Make computation - in our case this is for the computation of zeta
    freqs(params)
    freqs.generate_tau_values()

    new_peaks_tau = mixed_modes_utils.peaks_stretched_period(l1_peaks.frequency.values, 
                                                                psd_bgr.frequency.values, 
                                                                freqs.tau)
    new_peaks_zeta = mixed_modes_utils.peaks_stretched_period(l1_peaks.frequency.values, 
                                                                psd_bgr.frequency.values, 
                                                                freqs.zeta)
    #st.write(freqs.l1_mixed_freqs_p1, freqs.l1_mixed_freqs_n1)
    # Plot stretched echelle
    model_freqs = np.c_[freqs.l1_mixed_freqs, freqs.l1_mixed_freqs_p1, freqs.l1_mixed_freqs_n1]
    model_tau = np.c_[freqs.l1_mixed_tau, freqs.l1_mixed_tau_p1, freqs.l1_mixed_tau_n1]
    real_tau = new_peaks_tau
    real_freqs = l1_peaks.frequency.values
    heights = l1_peaks.amplitude*10
    y_real = (real_tau - DPi1*(1/2 + (-eps_g)))  % DPi1 - DPi1/2
    # The shift is already accounted for in the calculation of the theoretical frequencies, so don't add it in here
    y_theo = (model_tau[:,0] - DPi1/2) % DPi1  - DPi1/2
    y_theo_p1 = (model_tau[:,1] - DPi1/2) % DPi1  - DPi1/2
    y_theo_n1 = (model_tau[:,2] - DPi1/2) % DPi1  - DPi1/2

    
    # plt.scatter(y_real, real_freqs, s=heights)
    # plt.scatter(y_theo, model_freqs[:,0], marker='x')
    # plt.scatter(y_theo_p1, model_freqs[:,1], marker='x')
    # plt.scatter(y_theo_n1, model_freqs[:,2], marker='x')

    p = figure(
        title="",
        x_axis_label="Stretched Period (s)",
        y_axis_label="Frequency (μHz)",
        match_aspect=True,
        x_range=(-DPi1/2, DPi1/2),
        y_range=(psd_bgr.frequency.min(), 
                    psd_bgr.frequency.max()),
        plot_height=400,
        plot_width=600)  


    p.triangle(y_real, real_freqs, size=10*l1_peaks.amplitude, 
                color='red', alpha=0.5)

    p.line(y_theo,
                model_freqs[:,0],
                #angle=np.pi, size=10, 
                #alpha=0.5, 
                color="blue", legend_label=r'l=1, m=0')

    if splitting > 0.0:

        p.line(y_theo_n1,
                    model_freqs[:,2],
                    #angle=-np.pi/2, size=10, alpha=0.5, 
                    color="cyan", 
                    legend_label=r'l=1, m=-1')
        p.line(y_theo_p1,
                    model_freqs[:,1], 
                    #angle=np.pi/2, size=10, alpha=0.5, 
                    color="purple", legend_label=r'l=1, m=+1')

    p.legend.click_policy = 'hide'
    st.bokeh_chart(p)