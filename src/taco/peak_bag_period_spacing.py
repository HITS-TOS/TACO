import site
import sys
site.addsitedir('../../libs/sloscillations')

import argparse
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy

from astropy.timeseries import LombScargle
from pathlib import Path
from sloscillations import frequencies, mixed_modes_utils
from tqdm import tqdm

#from peakBagRotation import validate_arguments

from scipy.stats import chi2

from joblib import Parallel, delayed

from loguru import logger


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
    noise = np.median(PSD_LS) / (1 - 1/9)**3

    # Cut down period and power arrays to search range if given
    if search_range is None:
        cut_f = f
        cut_PSD_LS = PSD_LS
    else:
        cut_PSD_LS = PSD_LS[(1/f > search_range[0]) & (1/f < search_range[1])]
        cut_f = f[(1/f > search_range[0]) & (1/f < search_range[1])]


    # Report significance level of highest peak in search area
    sig_prob = chi2.cdf(np.max(cut_PSD_LS/noise), df=2)

    # Highest peak
    highest_peak_period = (1/cut_f)[cut_PSD_LS == cut_PSD_LS.max()].item()
    highest_peak = np.max(cut_PSD_LS/noise)

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

def validate_arguments(argv: argparse.Namespace) -> None:
    """
    Raises an `IOError` if file given in argv doesn't exist.

    Parameters
    ----------
    argv : NameSpace
        argv namespace of inputs to file.
    Returns
    -------
    None
    """
    # Loop over argparser arguments
    for arg in vars(argv):
        # Only look if attribute is a string, i.e. ignore ncores as not a file
        # name
        if type(getattr(argv, arg)) == str:
            # If file doesn't exists then raise IOError
            if not Path(getattr(argv, arg)).is_file():
                raise IOError(f"{arg} file {getattr(argv, arg)} does not exist!")

if __name__=="__main__":

    # Get the filenames and options
    parser = argparse.ArgumentParser(description="Asymptotic period spacing (ΔΠ1) from the periodogram of the 'stretched' PDS.")
    parser.add_argument('--peaks', dest='peaks', default = "peaksMLE.csv",
                        help = "File name of the identified peaks")
    parser.add_argument('--summary', dest='summary', default='summary.csv',
                        help = "File name of the csv file with summary values. It must contain numax and sigmaEnv")
    parser.add_argument("--pds", dest='pds', default = "pds_bgr.csv",
                        help = "File name of the csv file with the background-subtracted PDS It must contain the columns 'frequency' and 'power")
    parser.add_argument("--full_pds", dest='full_pds', default = "pds.csv",
                        help = "File name of the csv file with the original PDS It must contain the columns 'frequency' and 'power")
    parser.add_argument("--maxiters", dest='maxiters', default = 10, type=int,
                        help = "Maximum number of iterations to use in the self-consistent calculation of ΔΠ1")
    parser.add_argument("--niters", dest='niters', default = 5, type=int,
                        help = "Number of iterations to repeat the calculation of ΔΠ1, q and ε_g.")
    parser.add_argument("--dpi_only", dest='dpi_only', default=False, type=bool,
                        help = "Only infer the period spacing and don't calculate tau or q.")
    parser.add_argument("--ncores", dest='ncores', default = 1, type=int,
                        help = "Number of cores to use for a parallel calculation")
    parser.add_argument("--plot", dest="plot", default=False, type=bool,
                        help = "Plot diagnostic plots.")
    argv = parser.parse_args()

    # Check files exist by validating arguments in argparser
    validate_arguments(argv)

    # Read in summary file
    summary = pd.read_csv(argv.summary)

    # Adding in FLAG for DeltaPi1 computation
    # 0 - computation succeeded
    # 1 - DeltaNu < 3
    # 2 - No significant peak detected in RC PS of stretched PS and delta nu too low to obtain period spacing for RGB star
    # 3 - No significant peak detected in either PS of stretched power spectrum
    # 4 - Algorithm didn't converge after specified number of iterations

    if (summary.DeltaNu.values < 3.0):
        summary['DeltaPi1'] = np.nan
        summary['coupling'] = np.nan
        summary['eps_g'] = np.nan
        summary['DeltaPi1_val'] = np.nan
        summary['DeltaPi1_Flag'] = 1
        summary['DeltaPi1_sig'] = np.nan
        summary.to_csv(argv.summary, index=False)
        sys.exit('Delta nu too low to obtain period spacing for any star, exiting ...')

    # Read in and filter peaks file to be within +/-3 sigmaEnv of numax
    peaks = pd.read_csv(argv.peaks)
    peaks = peaks.loc[abs(peaks.frequency.values - summary.numax.values) < 3*summary.sigmaEnv.values, ]

    # Read in pds
    pds = pd.read_csv(argv.pds)
    pds = pds.loc[abs(pds.frequency.values - summary.numax.values) < 3*summary.sigmaEnv.values, ]


    #plt.plot(pds.frequency, pds.power)
    #plt.show()

    # Read in original pds incase we need it
    pds_full = pd.read_csv(argv.full_pds)

    # Split the peaks in the l=0,2,3 peaks (which have been already identified)
    # and the rest, which should hopefully be unidentified l=3
    l023_peaks = peaks.loc[(peaks.l == 0) | (peaks.l == 2) | (peaks.l == 3), ]
    l0_peaks = peaks.loc[(peaks.l==0), ]
    l1_peaks = peaks.loc[(peaks.l == 1) | (np.isfinite(peaks.l) == False)]

    #if len(l1_peaks) < 3:
    #    sys.exit('Not enough peaks to determine period spacing')

    #################################

    # Initial ΔΠ1 values
    RC_init = 300
    RGB_init = DeltaPi1_from_DeltaNu_RGB(summary.DeltaNu.values
                                         if np.isfinite(summary.DeltaNu.values)
                                         else DeltaNu_from_numax(summary.numax.values))

    # Compute ranges of period spacing to search for RC and RGB stars
    RC_DPi_range = np.array([100, 400])
    RGB_DPi_range = np.array([20, 120])
    RGB_DPi_range = np.array([RGB_init*0.6, RGB_init*1.4])

    #RGB_DPi_range = np.array([RGB_init*0.8, RGB_init*1.2])

    # If, for some reason, the RGB initial value is below our range then change
    # it slightly so it is
    if RGB_init < np.min(RGB_DPi_range):
        RGB_init = np.min(RGB_DPi_range) + 1
    elif RGB_init > np.max(RGB_DPi_range):
        RGB_init = np.max(RGB_DPi_range) - 1

    # Initial coupling values
    q_RC = 0.3
    q_RGB = 0.2

    # How much leeway we leave while estimating ΔΠ1 self-consistently. A larger range produces faster convergence but reduced stability
    search_mult = np.array([0.9, 1.1])

    # We DON'T estimate epsilon_g here and instead do that in the next module
    # when searching for rotation

    # Remove l=0,2 frequencies
    pds_l023_removed = pds.assign(power = pds.power / fit_model(pds, l023_peaks))

    #plt.plot(pds_l023_removed.frequency, pds_l023_removed.power)
    #plt.show()

    bw = np.mean(np.diff(pds.frequency.values))
    ratio = np.sum(fit_model(pds, l1_peaks)-1) / np.sum(fit_model(pds, l0_peaks)-1)
    summary['visibility_ratio'] = ratio
    print(f"Visibility ratio (l=1/l=0): {ratio}")

    # Create artificial frequencies for creation of stretched power spectrum using values determined from TACO for this star
    freqs = frequencies.Frequencies(frequency=pds_l023_removed.frequency.values,
                                    numax=summary.numax.values,
                                    delta_nu=summary.DeltaNu.values if np.isfinite(summary.DeltaNu.values) else None,
                                    epsilon_p=summary.eps_p.values if np.isfinite(summary.eps_p.values) else None,
                                    alpha=summary.alpha.values if np.isfinite(summary.alpha.values) else None)
    ## Estimate ΔΠ1 self-consistently with some initial guesses using both a high and low initial guess

    RC_test_dpi, RC_test_maximum, RC_sig = DPi1_from_stretched_PDS(RC_init, q_RC,
                                                            freqs,
                                                            pds_l023_removed,
                                                            return_max=True,
                                                            plot=argv.plot,
                                                            search_range = RC_DPi_range)


    if (RC_sig == False) and (summary.DeltaNu.values < 4.0):
        summary['DeltaPi1'] = np.nan
        summary['coupling'] = np.nan
        summary['eps_g'] = np.nan
        summary['DeltaPi1_val'] = np.nan
        summary['DeltaPi1_Flag'] = 2
        summary['DeltaPi1_sig'] = np.nan
        summary.to_csv(argv.summary, index=False)
        sys.exit('No significant peak detected in power spectrum and delta nu too low to obtain period spacing for RGB star, exiting ...')

    RGB_test_dpi, RGB_test_maximum, RGB_sig = DPi1_from_stretched_PDS(RGB_init, q_RGB,
                                                                freqs,
                                                                pds_l023_removed,
                                                                return_max=True,
                                                                plot=argv.plot,
                                                                search_range = RGB_DPi_range)
    # If no significant peak detected in PSxPS then exit
    #print("SIGS: ", RC_sig, RGB_sig)

    if (RC_sig == False) and (RGB_sig == False):
        summary['DeltaPi1'] = np.nan
        summary['coupling'] = np.nan
        summary['eps_g'] = np.nan
        summary['DeltaPi1_val'] = np.nan
        summary['DeltaPi1_sig'] = np.nan
        summary['DeltaPi1_Flag'] = 3
        summary.to_csv(argv.summary, index=False)
        sys.exit('No significant peak detected in power spectrum of stretched power spectrum, exiting ...')

    # If significant peak only found in one of two test spectra
    elif (RC_sig < 0.9) and (RGB_sig >= 0.9) or (summary['DeltaNu'].values > 13.0): # Take delta_nu > 13uHz as definitely RGB to avoid secondary clump stars
        DPi1_init = RGB_init
        q_init = q_RGB
        search_range = RGB_DPi_range
    elif (RC_sig >= 0.9) and (RGB_sig < 0.9):
        DPi1_init = RC_init
        q_init = q_RC
        search_range = RC_DPi_range
    # Otherwise proceed as normally would
    else:
        if RC_test_maximum > RGB_test_maximum:
            DPi1_init = RC_init
            q_init = q_RC
            search_range = RC_DPi_range
        else:
            DPi1_init = RGB_init
            q_init = q_RGB
            search_range = RGB_DPi_range

    curr_DPi1 = DPi1_init

    print(f"Starting DPi1: {curr_DPi1}")

    converged = False
    for i in range(argv.maxiters):

        old_DPi1 = curr_DPi1
        new_DPi1, val, sig = DPi1_from_stretched_PDS(curr_DPi1, q_init,
                                                     freqs,
                                                     pds_l023_removed,
                                                     return_max=True,
                                                     plot=argv.plot,
                                                     search_range=search_range)
        if sig > 0.9:
            curr_DPi1 = new_DPi1
            curr_val = val
            curr_sig = sig
            print(f"Current DPi1: {curr_DPi1} at significance level {curr_sig}.")

            if abs(curr_DPi1 - old_DPi1) < 1e-3:
                converged = True
                break

    if converged == False:
        summary['DeltaPi1'] = np.nan
        summary['coupling'] = np.nan
        summary['eps_g'] = np.nan
        summary['DeltaPi1_val'] = np.nan
        summary['DeltaPi1_Flag'] = 4
        summary['DeltaPi1_sig'] = np.nan
        summary.to_csv(argv.summary, index=False)

        sys.exit('Algorithm did not converge, no period spacing found, exiting ...')

    #if (curr_DPi1 < 50) and (abs(curr_DPi1 - DPi1_init) > 15.0):
    #    summary['DeltaPi1'] = curr_DPi1
    #    summary['coupling'] = np.nan
    #    summary['eps_g'] = np.nan
    #    summary['DeltaPi1_val'] = curr_val
    #    summary['DeltaPi1_Flag'] = 5
    #    summary.to_csv(argv.summary, index=False)
    #    sys.exit('DPi1 value is too low, star needs to be checked, exiting ...')



    print(f"DPi1: {curr_DPi1}")

    if argv.dpi_only:
        summary['DeltaPi1'] = curr_DPi1
        summary['coupling'] = np.nan
        summary['eps_g'] = np.nan
        summary['DeltaPi1_val'] = curr_val
        summary['DeltaPi1_Flag'] = 0
        summary['DeltaPi1_sig'] = curr_sig
        summary.to_csv(argv.summary, index=False)
        sys.exit('Find DPi1 only flag set, exiting ...')


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

    summary['DeltaPi1'] = dpi_best
    summary['coupling'] = q_best
    summary['eps_g'] = 0
    summary['DeltaPi1_sig'] = sig_best

    if np.all(results[:,1] == 0):
        summary.to_csv(argv.summary, index=False)

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
    freqs(params)

    # Compute stretched_pds
    new_freq, tau, zeta = mixed_modes_utils.stretched_pds(pds.frequency.values,
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

    # Write out results
    summary.to_csv(argv.summary, index=False)
    pds.to_csv(argv.pds, index=False)
    peaks.to_csv(argv.peaks, index=False)
