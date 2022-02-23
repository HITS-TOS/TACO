import site
import sys
site.addsitedir('../../')

import argparse
import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd 
import scipy

from lib.rotation import refine_rotation_parameters, rotation_utils

from pathlib import Path
from sklearn.neighbors import KDTree 
from sloscillations import frequencies
from tqdm import tqdm

from scipy.spatial import distance_matrix
from scipy.optimize import linear_sum_assignment

from typing import List, Optional

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
    parser = argparse.ArgumentParser(description="Identification of the azimuthal order of l=1 mixed modes.")
    parser.add_argument('--peaks', dest='peaks', default = "peaksMLE.csv",
                        help = "File name of the identified peaks")
    parser.add_argument('--summary', dest='summary', default='summary.csv',
                        help = "File name of the csv file with summary values. It must contain numax and sigmaEnv")
    parser.add_argument("--pds", dest='pds', default = "pds_bgr.csv",
                        help = "File name of the csv file with the background-subtracted PDS It must contain the columns 'frequency' and 'power")
    parser.add_argument("--ncores", dest='ncores', default = 1, type=int,
                        help = "Number of cores to use for a parallel calculation")
    parser.add_argument("--cut", dest='cut', default=0.7, type=float, 
                        help = "Cut in zeta for rotational splitting finder.")
    argv = parser.parse_args()
    
    # Check files exist by validating arguments in argparser
    validate_arguments(argv)

    # Read in summary file
    summary = pd.read_csv(argv.summary)
    
    if (np.isfinite(summary['DeltaPi1'].values) == False) or (np.isfinite(summary['coupling'].values) == False):
        sys.exit('No valid period spacing or coupling value found, cannot compute rotational splitting, exiting ....')

    # Read in power spectrum
    pds = pd.read_csv(argv.pds)
    # Only keep pds around oscillations
    pds = pds.loc[abs(pds['frequency'].values - summary['numax'].values) < 3 * summary['sigmaEnv'].values, ]

    # Read in and filter peaks file to be within +/-3 sigmaEnv of numax
    peaks = pd.read_csv(argv.peaks)
    # Compute "reduced frequency" for modes
    peaks['x'] = ((peaks['frequency'] % summary['DeltaNu'].values - summary['eps_p'].values) / summary['DeltaNu'].values) % 1
    peaks.loc[~np.isfinite(peaks['l']), 'l'] = 1

    data_n_max = np.floor((summary.numax.values / summary.DeltaNu.values) - summary.eps_p.values).item()
    # -1 in lower range here in case have mixed modes in radial order below lowest detected radial mode
    radial_order_range = [peaks.n.min() - data_n_max - 1, peaks.n.max() - data_n_max ]

    l1_peaks = rotation_utils.prepare_l1_peaks(peaks, summary)

    # For aliasing, use this as range in ΔΠ1
    alias_spacing = rotation_utils.compute_alias_spacing(summary.numax.item(), summary.DeltaPi1.values.item())

    dpi_range = np.linspace(summary.DeltaPi1.values - 2.0*alias_spacing,
                            summary.DeltaPi1.values + 2.0*alias_spacing, 50)

    # Set up rotational splitting range (in uHz)
    if summary.DeltaNu.values > 5:
        split_range = np.logspace(np.log(0.05), np.log(1.0), 100, base=np.exp(1))
    else:
        split_range = np.logspace(np.log(0.01), np.log(0.5), 100, base=np.exp(1))

    # Add no rotational splitting to start of array
    split_range = np.insert(split_range, 0, 0)

    # Estimate rotational splitting
    refined_dpi, refined_split, refined_eps_g, n_ridges = refine_rotation_parameters.find_rotational_splitting(pds.frequency.values, summary.numax.values, l1_peaks, 
                                                                        dpi_range, split_range, summary.DeltaNu.values, summary.eps_p.values, 
                                                                        summary.alpha.values, radial_order_range, summary.coupling.values, 
                                                                        real_heights=None, plot=False, return_distances=False)

    print(f"Best δν_rot {refined_split}")
    print(f"Best ε_g {refined_eps_g}")
    print(f"Best ΔΠ_1 {refined_dpi}")      

    if refined_dpi < 100:
        coupling_range = np.linspace(0.01, 0.5, 100)
    else:
        coupling_range = np.linspace(0.1, 0.7, 100)

    # Refined coupling
    refined_coupling, (dists, _) = refine_rotation_parameters.refine_coupling(pds.frequency.values, summary.numax.values, l1_peaks, 
                                                    coupling_range, refined_dpi, refined_eps_g, refined_split, 
                                                summary.DeltaNu.values, summary.eps_p.values, 
                                                    summary.alpha.values, radial_order_range, summary.coupling.values, 
                                                    real_heights=None, plot=False, return_distances=True)
    print(f"Best q {refined_coupling}")

    # Refining d01
    d01_range = np.linspace(-0.05*summary.DeltaNu.values, 0.05*summary.DeltaNu.values, 100)
    new_d01, (dists, _) = refine_rotation_parameters. refined_nominal_pmodes(pds.frequency.values, summary.numax.values, l1_peaks, 
                                                    d01_range, refined_dpi, refined_coupling, refined_eps_g, refined_split, 
                                                summary.DeltaNu.values, summary.eps_p.values, 
                                                    summary.alpha.values, radial_order_range, 
                                                    real_heights=None, plot=False, return_distances=True)
    print(f"Best d01 {new_d01}")

    # Plot result
    rotation_utils.plot_results(pds, summary, l1_peaks, radial_order_range, refined_dpi, refined_coupling, refined_eps_g, refined_split, new_d01)

    summary['DeltaPi1'] = refined_dpi
    summary['coupling'] = refined_coupling
    summary['eps_g'] = refined_eps_g
    summary['rot_split'] = refined_split
    #summary['d01'] = np.nan
    #summary['radial_order_range'] = np.nan
    #summary.at[0, 'd01'] = list(new_d01)
    #summary.at[0, 'radial_order_range'] = list(radial_order_range)
    summary.to_csv(argv.summary, index=False)
