import sys

sys.path.insert(0, 'libs/sloscillations')

import site

site.addsitedir('../../libs/sloscillations')

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
from astropy.timeseries import LombScargle
from joblib import Parallel, delayed
from loguru import logger
from scipy.stats import chi2
from sloscillations import frequencies, mixed_modes_utils
from tqdm import tqdm


#' It is assumed that DeltaNu is in μHz
def DeltaPi1_from_DeltaNu_RGB(DeltaNu):
    # Compute Period spacing (in s) from deltanu
    return 60 + 1.7*DeltaNu

def Lor_model(pds, peak):
    return peak.height / (1 + ((pds.frequency.values - peak.frequency)/peak.linewidth)**2)

def sinc2_model(pds, peak):
    deltanu = np.mean(np.diff(pds.frequency.values))
    return peak.height * np.sinc((pds.frequency.values - peak.frequency)/deltanu)**2

def fit_model(pds, peaks):
    model = np.ones_like(pds.frequency.values)

    for i in range(len(peaks)):
        if np.isfinite(peaks.linewidth.iloc[i]):
            model += Lor_model(pds, peaks.iloc[i,])
        else:
            model += sinc2_model(pds, peaks.iloc[i, ])
    return model

def DeltaNu_from_numax(numax):
    return 0.276 * numax ** 0.751

def coupling_from_dnu_RGB(dnu):
    """
    Based on RGB coupling values in Mosser et al.(2017)
    using a linear fit
    """
    return 0.104 + 0.0024*dnu

def density_mixed_modes(DeltaNu, numax):
    dpi1 = DeltaPi1_from_DeltaNu_RGB(DeltaNu) * 1e-6
    return DeltaNu / (dpi1 * (numax ** 2))


@logger.catch
def DPi1_from_stretched_PDS(DPi1, q, freqs, pds, return_max=False, plot=False, search_range=None):
    """
    Determine the period spacing from the power spectrum of the stretched power
    spectrum.

    Parameters
    -------------
    DPi1 : float
        Current DPi1 guess

    q : float
        Current coupling guess

    freqs : sloscillations object
        Sloscillations class object used to compute the stretched period for
        a given DPi1, q combination.

    pds : dataframe
        Dataframe containing the power spectrum (frequency and pds columns must
        be present).

    return_max : boolean (optional)
        Whether or not to also return the value of the maximum of the power
        spectrum of the stretched power spectrum.

    plot : boolean (optional)
        Whether or not to plot the power spectrum of the stretched power spectrum.

    search_range : list (optional)
        Search range in period spacing in which to search for possible peak. If
        given then should be a list (or array) containing the lower and upper
        bounds.

    Returns - dependent on set keywords
    -------
    float
        Candidate period spacing value
    float
        Power value corresponding to candidate period spacing value (if
        return_max is True)
    boolean
        Whether candidate peak is significant

    """
    params = {'calc_l0': True, # Compute radial mode properties
            'calc_l2': True, # Compute l=2 mode properties
            'calc_l3': False, # Don't need to calculate l=3 theoretical freqs
            'calc_nom_l1': True, # Compute nominal l=1 p-mode properties
            'calc_mixed': False, # Don't compute mixed modes (as not needed)
            'calc_rot': False, # Don't compute rotation
            'DPi1': DPi1,
            'coupling': q,
            'eps_g': 0.0, # Epsilon_g isn't needed for computation of tau due to chosen formulation of zeta
            'l': 1, # Mixed modes are dipole mixed modes
            }
    # Make computation - in our case this is for the computation of zeta
    freqs(params)


    # Compute tau from the zeta value just computed
    new_frequency, tau, zeta = mixed_modes_utils.stretched_pds(pds.frequency.values,
                                                               freqs.zeta)
    if tau.max() > 0:
        # If the search range isn't given then default to frequency range
        # corresponding to period range of 20-400s
        #if search_range is None:
        f = np.arange(1/(400.), 1/(20.), 0.1/tau.max())
        #else:
        #    f = np.arange(1/search_range[1], 1/search_range[0], 0.1/tau.max())

        # Set up Lomb-Scargle periodogram calculation
        ls = LombScargle(tau, pds.power.values)
        PSD_LS = ls.power(f)

        # Compute background noise level for significance level
        # This means we can normalise the periodogram to SNR which makes significance
        # level calculation much easier

        noise = np.nanmedian(PSD_LS) / (1 - 1/9)**3

        # Cut down period and power arrays to search range if given
        if search_range is None:
            cut_f = f
            cut_PSD_LS = PSD_LS
        else:
            cut_PSD_LS = PSD_LS[(1/f > search_range[0]) & (1/f < search_range[1])]
            cut_f = f[(1/f > search_range[0]) & (1/f < search_range[1])]


        # Report significance level of highest peak in search area
        sig_prob = chi2.cdf(np.nanmax(cut_PSD_LS/noise), df=2)

        # Highest peak
        highest_peak_period = (1/cut_f)[cut_PSD_LS == cut_PSD_LS.max()].item()
        highest_peak = np.nanmax(cut_PSD_LS/noise)

        #logger.debug(f"Highest peak found at {highest_peak_period} seconds, with power {highest_peak} and significance {sig_prob}")

        if plot:
            plt.plot(1/cut_f, cut_PSD_LS/noise)
            plt.axhline(5.99, linestyle='--', color='r', label=r'FAP 95.0%')
            plt.xlabel(r'$\Delta\Pi_{1}$ (s)', fontsize=18)
            plt.ylabel(r'Power (arb. units)', fontsize=18)
            plt.legend(loc='best')
            plt.show()

        # Return candidate period spacing value as period spacing corresponding to
        # highest peak in chosen search range.
        if return_max:
            # 9.21 is required level for chi-squared 2 d.o.f with FAP of 99%
            # 7.378 is required level for chi-squared 2 d.o.f with FAP of 97.5%
            # 5.99 is required level for chi-squared 2 d.o.f with FAP of 95.0%
            return highest_peak_period, highest_peak, sig_prob #, highest_peak > 5.99 #7.378 #9.210
        else:
            return highest_peak_period, sig_prob #highest_peak > 5.99 #7.378 #9.210
    else:
        highest_peak_period = 0
        highest_peak = 0
        sig_prob = False
        if return_max:
            # 9.21 is required level for chi-squared 2 d.o.f with FAP of 99%
            # 7.378 is required level for chi-squared 2 d.o.f with FAP of 97.5%
            # 5.99 is required level for chi-squared 2 d.o.f with FAP of 95.0%
            return highest_peak_period, highest_peak, sig_prob #, highest_peak > 5.99 #7.378 #9.210
        else:
            return highest_peak_period, sig_prob #highest_peak > 5.99 #7.378 #9.210

def DPi1_iteration(maxiters, conv_thresh, sig_thresh, DPi10, q0, freqs, pds, return_max, plot, search_range):
    """
                    Iteration over period spacing from stretched PDS

    Input parameters:
    -------------------
        maxiter:      number of maximum iterations
        conv_thresh:  Minimum DPi1 difference with previous iteration to claim convergence of iterations.
        sig_thresh:   Minimum DPi1 significance value
        DPi10:        Initial DPi1 guess
        q0:           Initial q value
        freqs:        Artificial frequencies for creation of stretched power spectrum using values determined from TACO for this star
        pds:          PDS of potential dipole modes (e.g. PDS after removing l = 0,2,3 modes.)
        return_max:   Return value of highest peak of the LS of strecthed PDS
        plot:         Boolean
        search_range: Range of DPi1 to search for RGB or RC stars
    """

    converged = False
    curr_DPi1 = DPi10

    print(f"Starting DPi1: {curr_DPi1}")

    estimates_DPi1, value_DPi1, significance_DPi1 = [], [], []
    looped = False

    for i in range(maxiters):
        if curr_DPi1 > 0:
            old_DPi1 = curr_DPi1
            new_DPi1, val, sig = DPi1_from_stretched_PDS(curr_DPi1, q0,
                                                        freqs,
                                                        pds,
                                                        return_max = return_max,
                                                        plot = plot,
                                                        search_range = search_range)

            curr_DPi1 = new_DPi1
            curr_val = val
            curr_sig = sig

            estimates_DPi1.append(curr_DPi1)
            value_DPi1.append(curr_val)
            significance_DPi1.append(curr_sig)

            if (abs(curr_DPi1 - old_DPi1) < conv_thresh) and (max(np.array(significance_DPi1)) >= sig_thresh) and (i == maxiters - 1):

                converged = True

                # Retrive Dpi1 with highest significance.
                highest_sig = np.array(significance_DPi1) == max(np.array(significance_DPi1))

                curr_DPi1 = np.unique(np.array(estimates_DPi1)[highest_sig])[0]
                curr_val = np.unique(np.array(value_DPi1)[highest_sig])[0]
                curr_sig = max(np.array(significance_DPi1))

            else:
                if looped == False:
                    if (i == maxiters - 1):
                        maxiters *= 10
                        looped = True
                    else:
                        looped = False
                else:
                    return False, np.nan, np.nan, np.nan

    if converged:
        return converged, curr_DPi1, curr_val, curr_sig

    else:
        return False, np.nan, np.nan, np.nan


def peak_bag_period_spacing(pds, peaks, data,
        maxiters = 10, conv_thresh = 1e-3, sig_thresh = 0.9,
        niters = 5, dpi_only = False, ncores = 1, plot = False):

    """
    Asymptotic period spacing (ΔΠ1) from the periodogram of the 'stretched' PDS.

    Parameters:
        pds(pandas.DataFrame): Periodogram
            Columns:
                Name: frequency, dtype: float[micro-Hertz]
                Name: power, dtype: float
        data(pandas.DataFrame): Summary data
            Columns:
                Name: numax, dtype: float
                Name: sigmaEnv, dtype: float
        maxiters(int): Maximum number of iterations to use in the self-consistent calculation of ΔΠ1
        niters(int): Number of iterations to repeat the calculation of ΔΠ1, q and ε_g
        dpi_only(bool): Only infer the period spacing and don't calculate tau or q
        ncores(int): Number of cores to use for a parallel calculation
        plot(bool): Plot diagnostic plots

    Returns:
        pds(pandas.DataFrame): Periodogram
        peaks(pandas.DataFrame): Identified peaks
        data(pandas.DataFrame): Summary data
    """

    #################################
    #    Keep peaks with AIC > 10
    #################################
    AIC_cut = 2.0
    peaks = peaks.loc[(peaks.AIC.values > AIC_cut) | (peaks.AIC1.values > AIC_cut), ].sort_values(by = 'frequency')


    # Adding in FLAG for DeltaPi1 computation
    # 0 - computation succeeded
    # 1 - DeltaNu < 3
    # 2 - No significant peak detected in RC PS of stretched PS and delta nu too low to obtain period spacing for RGB star
    # 3 - No significant peak detected in either PS of stretched power spectrum
    # 4 - Algorithm didn't converge after specified number of iterations
    # 5 - Algorithm couldn't compute set of artificial frequencies with sloscillations
    # 6 - Not enought l = 1 detected candidates (might be suppresed dipole modes)
    # 7 - Star has more mixed-modes than predicted (might have rotational splittings)

    flag = 0

    if (data.DeltaNu.values < 3.0):
        data['DeltaPi1'] = np.nan
        data['coupling'] = np.nan
        data['eps_g'] = np.nan
        data['DeltaPi1_val'] = np.nan
        data['DeltaPi1_Flag'] = 1
        data['DeltaPi1_sig'] = np.nan
        print('Delta nu too low to obtain period spacing for any star')
        flag = 1
        return(pds, peaks, flag, data)

    # Filter peaks file to be within +/-3 sigmaEnv of numax
    # NOTE: This is already done from PeaksFind on!!!

    peaks = peaks.loc[abs(peaks.frequency.values - data.numax.values) < 3 * data.sigmaEnv.values, ]
    pds = pds.loc[abs(pds.frequency.values - data.numax.values) < 3 * data.sigmaEnv.values, ]

    # Split the peaks in the l=0,2,3 peaks (which have been already identified)
    # and the rest, which should hopefully be unidentified l=3
    l023_peaks = peaks.loc[(peaks.l == 0) | (peaks.l == 2) | (peaks.l == 3), ]
    #l023_peaks = peaks.loc[(peaks.l == 0) | (peaks.l == 2), ]
    l0_peaks = peaks.loc[(peaks.l==0), ]
    l1_peaks = peaks.loc[(peaks.l == 1) | (np.isfinite(peaks.l) == False)]


    # Don't proceed if there are less than 3 potential dipole modes per radial order.
    # To avoid computing period spacing for stars with suppressed dipoles
    nl1_flag = False
    l0_freqs = l0_peaks.frequency.values
    l1_freqs = l1_peaks.frequency.values
    for i in range(0, len(l0_freqs) - 1):
        l1_n = l1_peaks[(l0_freqs[i] <= l1_freqs) & (l1_freqs <= l0_freqs[i + 1])]
        if len(l1_n) < 2:
            nl1_flag = True

    if nl1_flag:
        data['DeltaPi1'] = np.nan
        data['coupling'] = np.nan
        data['eps_g'] = np.nan
        data['DeltaPi1_val'] = np.nan
        data['DeltaPi1_Flag'] = 6
        data['DeltaPi1_sig'] = np.nan
        print('Not enought dipole modes detected to obtain period spacing of this star')
        flag = 6
        return(pds, peaks, flag, data)



    # Initial ΔΠ1 values
    RC_init = 300
    RGB_init = DeltaPi1_from_DeltaNu_RGB(data.DeltaNu.values[0]
                                         if np.isfinite(data.DeltaNu.values[0])
                                         else DeltaNu_from_numax(data.numax.values[0]))

    # Compute ranges of period spacing to search for RC and RGB stars
    RC_DPi_range = np.array([100, 400])

    if (data.DeltaNu.values >= 30.0): # Use a more constrained range before SGs (Dnu ~ 30)
        RGB_DPi_range = np.array([45, 400])
    else:
        RGB_DPi_range = np.array([45, 150])
        #RGB_DPi_range = np.array([RGB_init*0.6, RGB_init*1.4])

    #RGB_DPi_range = np.array([RGB_init*0.6, RGB_init*1.4]) # Too restrictive
    #RGB_DPi_range = np.array([RGB_init*0.8, RGB_init*1.2])


    # Initial coupling values
    q_RC = 0.3
    q_RGB = coupling_from_dnu_RGB(data.DeltaNu.values[0]
                                  if np.isfinite(data.DeltaNu.values[0])
                                  else DeltaNu_from_numax(data.numax.values[0]))

    # We DON'T estimate epsilon_g here and instead do that in the next module
    # when searching for rotation

    # Remove l=0,2 frequencies
    pds_l023_removed = pds.assign(power = pds.power / fit_model(pds, l023_peaks))


    bw = np.mean(np.diff(pds.frequency.values))
    ratio = np.sum(fit_model(pds, l1_peaks)-1) / np.sum(fit_model(pds, l0_peaks)-1)
    data['visibility_ratio'] = ratio
    print(f"Visibility ratio (l=1/l=0): {ratio}")

    # Create artificial frequencies for creation of stretched power spectrum using values determined from TACO for this star
    freqs = frequencies.Frequencies(frequency=pds_l023_removed.frequency.values,
                                    numax=data.numax.values,
                                    delta_nu=data.DeltaNu.values if np.isfinite(data.DeltaNu.values) else None,
                                    epsilon_p=data.eps_p.values if np.isfinite(data.eps_p.values) else None,
                                    alpha=data.alpha.values if np.isfinite(data.alpha.values) else None)


    ## Estimate ΔΠ1 self-consistently with some initial guesses using both a high and low initial guess
    _, RC_test_maximum, RC_sig = DPi1_from_stretched_PDS(RC_init, q_RC,
                                                         freqs,
                                                         pds_l023_removed,
                                                         return_max = True,
                                                         plot = plot,
                                                         search_range = RC_DPi_range)


    if (RC_sig == False) and (data.DeltaNu.values < 4.0): # Shouldn't it be dnu < 3.0?
        data['DeltaPi1'] = np.nan
        data['coupling'] = np.nan
        data['eps_g'] = np.nan
        data['DeltaPi1_val'] = np.nan
        data['DeltaPi1_Flag'] = 2
        data['DeltaPi1_sig'] = np.nan
        print('No significant peak detected in power spectrum and delta nu too low to obtain period spacing for RGB star')
        flag = 2
        return(pds, peaks, flag, data)

    _, RGB_test_maximum, RGB_sig = DPi1_from_stretched_PDS(RGB_init, q_RGB,
                                                           freqs,
                                                           pds_l023_removed,
                                                           return_max = True,
                                                           plot = plot,
                                                           search_range = RGB_DPi_range)

    #############################################
    #    New threshold to find RGB or RC
    # Based on central Dnu and eps_p relation
    # Also corrects very high eps_p values
    #############################################
    DNu = data.Central_DeltaNu.values[0]
    eps_p = data.Central_eps_p.values[0]

    if (eps_p >= 1.4) and (DNu < 10.0):
        eps_p = eps_p - 1.0

    eps_p_threshold = (0.592 * np.log10(DNu)) + 0.6  # Computed using current TACO results from ideal_sample v10.05 (2024)


    if (RC_sig == False) and (RGB_sig == False):
        data['DeltaPi1'] = np.nan
        data['coupling'] = np.nan
        data['eps_g'] = np.nan
        data['DeltaPi1_val'] = np.nan
        data['DeltaPi1_sig'] = np.nan
        data['DeltaPi1_Flag'] = 3
        print('No significant peak detected in power spectrum of stretched power spectrum')
        flag = 3
        return(pds, peaks, flag, data)

    # If significant peak only found in one of two test spectra
    elif ((RC_sig < sig_thresh) and (RGB_sig >= sig_thresh)) or ((data['DeltaNu'].values > 10.0)): # Take delta_nu > 10uHz as definitely RGB to avoid secondary clump stars from v10.05 TACO results.
        DPi1_init = RGB_init
        q_init = q_RGB
        search_range = RGB_DPi_range


    elif (RC_sig >= sig_thresh) and (RGB_sig < sig_thresh):
        DPi1_init = RC_init
        q_init = q_RC
        search_range = RC_DPi_range

    # Otherwise proceed as normally would
    else:
        #######################################
        #  New constrains using eps_p and dnu
        #######################################
        if (RC_test_maximum > RGB_test_maximum) and (eps_p < eps_p_threshold) and (DNu < 10.0):
            DPi1_init = RC_init
            q_init = q_RC
            search_range = RC_DPi_range

        else:
            DPi1_init = RGB_init
            q_init = q_RGB
            search_range = RGB_DPi_range


    ###################################
    #         Estimating DPi1
    ###################################

    curr_DPi1 = DPi1_init
    converged, curr_DPi1, curr_val, curr_sig = DPi1_iteration(maxiters=maxiters, conv_thresh=conv_thresh, sig_thresh=sig_thresh,
                                                              DPi10=curr_DPi1, q0=q_init, freqs=freqs, pds=pds_l023_removed,
                                                              return_max = True, plot = plot, search_range = search_range)


    #####################################################
    #  Iterate again if DPi1 > 360s or didn't converge
    #####################################################
    if converged == False or curr_DPi1 >= 360.0:
        prev_DPi1 = curr_DPi1
        prev_sig = curr_sig

        if curr_DPi1 > 100:
            print('Testing RGB instead, because it didnt converge')
            DPi1_init = RGB_init
            q_init = q_RGB
            search_range = RGB_DPi_range

        else:
            if data.DeltaNu.values < 9.0:
                print('Testing RC instead, because it didnt converge')
                DPi1_init = RC_init
                q_init = q_RC
                search_range = RC_DPi_range

            else:
                print('Testing RGB again')
                DPi1_init = prev_DPi1
                q_init = q_RGB
                search_range = RGB_DPi_range


        converged, curr_DPi1, curr_val, curr_sig = DPi1_iteration(maxiters=maxiters, conv_thresh=conv_thresh, sig_thresh=sig_thresh,
                                                              DPi10=DPi1_init, q0=q_init, freqs=freqs, pds=pds_l023_removed,
                                                              return_max = True, plot = plot, search_range = search_range)


        if prev_DPi1 >= 360.0 and converged == True:
            if (prev_sig < curr_sig) & (curr_DPi1 <= 360):
                print('Found a better estimation.')

            else:
                converged = False


    N = density_mixed_modes(data.DeltaNu.values, data.numax.values)
    if len(l1_freqs) / (l0_peaks.n.values.max() - l0_peaks.n.values.min()) > 1.5 * N[0]:
        print('Found more peaks than expected. ', curr_DPi1)


    if not converged:
        data['DeltaPi1'] = np.nan
        data['coupling'] = np.nan
        data['eps_g'] = np.nan
        data['DeltaPi1_val'] = np.nan
        data['DeltaPi1_Flag'] = 4
        data['DeltaPi1_sig'] = np.nan
        print('Algorithm did not converge, no period spacing found')
        flag = 4
        return(pds, peaks, flag, data)

    print(f"DPi1: {curr_DPi1}")

    if dpi_only:
        data['DeltaPi1'] = curr_DPi1
        data['coupling'] = np.nan
        data['eps_g'] = np.nan
        data['DeltaPi1_val'] = curr_val
        data['DeltaPi1_Flag'] = 0
        data['DeltaPi1_sig'] = curr_sig
        print('Find DPi1 only flag set')
        return(pds, peaks, flag, data)


    # Find optimal coupling value
    if curr_DPi1 < 100.0:
        q_range = np.linspace(0.01, 0.3, 200)
        search_range = np.array([curr_DPi1 - 0.05 * curr_DPi1,
                                 curr_DPi1 + 0.05 * curr_DPi1])
        #RGB_DPi_range #np.array([20, 120])
    else:
        q_range = np.linspace(0.01, 0.7, 1000)
        #search_range = RC_DPi_range #np.array([100, 400])
        search_range = np.array([curr_DPi1 - 0.05 * curr_DPi1,
                                 curr_DPi1 + 0.05 * curr_DPi1])
    with Parallel(n_jobs=1) as parallel:
    #aprun = ParallelExecutor(use_bar='tqdm', n_jobs=4)
        results = parallel(delayed(DPi1_from_stretched_PDS)(curr_DPi1, q_range[j], freqs, pds_l023_removed, return_max=True, search_range=search_range, plot=False) for j in tqdm(range(len(q_range)), total=len(q_range)))
    results = np.array(results)
    # Multiply by significance to act as a mask, zeros all non-significant peaks - doesn't work now as return significance levels explicitly
    #results[:,1] *= results[:,2]



    # plt.plot(q_range, results[:,0])
    # plt.plot(q_range, results[:,2]*100)
    # plt.show()

    if np.all(results[:,1] < 0.9):
        print("No significant peaks found (significance level > 0.9), cannot accurately measure the coupling!")
        q_best = np.nan
    else:
        q_best = q_range[results[:,1]==np.max(results[:,1])]

    dpi_best = results[results[:,1] == np.max(results[:,1]), 0]
    sig_best = results[results[:,1] == np.max(results[:,1]), 2]

    print(f"Best parameters q={q_best}, DPi={dpi_best} at significance level {sig_best}")

    data['DeltaPi1'] = dpi_best
    data['coupling'] = q_best
    data['eps_g'] = 0
    data['DeltaPi1_sig'] = sig_best


    # Update aritifical frequencies to compute zeta and tau for all frequencies
    params = {'calc_l0': True,
            'calc_l2': True,
            'calc_l3': False, # Don't need to calculate l=3 theoretical freqs
            'calc_nom_l1': True,
            'calc_mixed': True,
            'calc_rot': False,
            'DPi1': dpi_best,
            'coupling': q_best,
            'eps_g': 0.0,
            'l': 1,
            }

    try:
        freqs(params)

    except (ValueError, RuntimeError, TypeError, NameError):
        print("Problems to compute delta Pi")
        flag = 5
        return (pds, peaks, flag, data)

    # Compute stretched_pds
    _, tau, zeta = mixed_modes_utils.stretched_pds(pds.frequency.values,
                                                          freqs.zeta)

    # PDS tau and zeta values
    pds['tau'] = tau
    pds['zeta'] = zeta

    # Values for peaks
    peaks['tau'] = mixed_modes_utils.peaks_stretched_period(peaks.frequency.values,
                                                            pds.frequency.values,
                                                            tau)

    # This function can also be used with zeta as it just builds a simple interpolation function
    peaks['zeta'] = mixed_modes_utils.peaks_stretched_period(peaks.frequency.values,
                                                             pds.frequency.values,
                                                             zeta)

    pds = pds.reset_index(drop = True)
    return(pds, peaks, flag, data)
