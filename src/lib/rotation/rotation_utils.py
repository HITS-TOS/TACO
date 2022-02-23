import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from tqdm import tqdm

from sklearn.neighbors import KDTree 
from sloscillations import frequencies, mixed_modes_utils

import itertools
from . import refine_rotation_parameters

from typing import List, Optional, Tuple

def prepare_l1_peaks(peaks: pd.DataFrame, summary: pd.DataFrame,
                     AIC_cut: Optional[float] = 0.0, height_cut: Optional[float] = 0.0) -> pd.DataFrame:
    """
    Extract the mixed modes from the peaks dataframe.
    
    Parameters
    ----------
    peaks: pd.DataFrame
        Dataframe containing the detected peaks and parameters.
        
    summary: pd.DataFrame
        Dataframe containing the global stellar information.
    
    AIC_cut: Optional[float] = 0.0
        Cut to make in the Akaike Information Criterion if desired.
        
    height_cut: Optional[float] = 0.0
        Cut to make in the mode height if desired.
        
    Outputs
    -------
    pd.DataFrame
        Dataframe containing the mixed mode peaks and associated mode parameters.
    """
    
    # Don't want to include any modes near l=0 or 2s, this is why this and the step in the next cell is performed.
    x_range = [(np.minimum(np.min(peaks.loc[peaks['l'] == 0, 'x']), np.min(peaks.loc[peaks['l'] == 2, 'x'])) - 0.05) % 1,
               (np.maximum(np.max(peaks.loc[peaks['l'] == 0, 'x']), np.max(peaks.loc[peaks['l'] == 2, 'x'])) + 0.05) % 1]
    
    l1_peaks = peaks.loc[(peaks.l == 1) | ~np.isfinite(peaks.l), ]
    l1_peaks['x'] = ((l1_peaks['frequency'] % summary['DeltaNu'].values - summary['eps_p'].values) / summary['DeltaNu'].values) % 1
    #l1_peaks = l1_peaks.loc[(l1_peaks['x'] > 0.2) & (l1_peaks['x'] < 0.85), ]

    #if x_range[0] < x_range[1]:
    #    l1_peaks = l1_peaks.loc[(l1_peaks['x']  x_range[1]) & (l1_peaks['x'] > x_range[0]), ]
    #else:
    #    l1_peaks = l1_peaks.loc[(l1_peaks['x'] > x_range[1]) | (l1_peaks['x'] < x_range[0]), ]

    l1_peaks = l1_peaks.loc[(l1_peaks['height'] > height_cut), ]
    l1_peaks = l1_peaks.loc[(l1_peaks['AIC'] > AIC_cut), ]

    return l1_peaks

def compute_AIC(distances: np.ndarray) -> np.ndarray:
    """
    Compute the Akaike Information Criterion.
    
    Parameters
    ----------
    distances: np.ndarray
        Distances matrix from rotational splitting computation.
        
    Returns
    -------
    np.ndarray
        Akaike Information Criterion for each model (i.e. 1, 2 or 3 ridges).
    """
    return 2*np.array([1, 2, 3]) + 2*np.min(distances, axis=0)

def how_many_ridges(AIC: np.ndarray, return_idx: bool = False) -> float:
    """
    Chose how many ridges there are given the input AIC.
    
    Parameters
    ----------
    AIC: np.ndarray
        AIC values for each number of ridges (1, 2 or 3).
    
    Returns
    -------
    float
        Number of ridges from AIC.
    """
    if return_idx:
        return np.argmin(AIC, axis=0) + 1, np.argmin(AIC, axis=0)
    else:
        return np.argmin(AIC, axis=0) + 1

def compute_alias_spacing(numax: float, DPi1: float) -> float:
    """
    Compute the difference in seconds of an alias caused by the wrong radial order in seconds
    """
    return (numax * 1e-6) * (DPi1)**2

def plot_stretched_echelle(model_freqs: np.ndarray, model_tau: np.ndarray, 
                           real_freqs: np.ndarray, real_tau: np.ndarray, 
                           DPi1: float, shift: float, heights: Optional[np.ndarray]=None):
    """
    Utility function to plot the stretched echelle for a given set of real and theoretical frequencies
    
    Parameters
    ----------
    model_freqs: np.ndarray
        Theoretical mixed mode frequencies.
        
    model_tau: np.ndarray
        Theoretical mixed mode tau values.
        
    real_freqs: np.ndarray
        Detected peak frequencies.
        
    real_tau: np.ndarray
        Tau values of detected peaks.
        
    DPi1: float
        Period spacing by whcih to compute the stretched echelle.
        
    shift: float
        The computed shift to ensure that the stretched echelle is properly aligned.
        
    heights: Optional[np.ndarray] = None
        Heights of detected peaks. Optional and if given are used to set the size of the 
        markers of the real frequencies in the stretched echelle.
    """
    
    y_real = (real_tau - DPi1*(1/2 + shift))  % DPi1 - DPi1/2
    # The shift is already accounted for in the calculation of the theoretical frequencies, so don't add it in here
    y_theo = (model_tau[:,0] - DPi1/2) % DPi1  - DPi1/2
    y_theo_p1 = (model_tau[:,1] - DPi1/2) % DPi1  - DPi1/2
    y_theo_n1 = (model_tau[:,2] - DPi1/2) % DPi1  - DPi1/2

    if heights is None:
        plt.scatter(y_real, real_freqs)
    else:
        plt.scatter(y_real, real_freqs, s=heights)
    plt.scatter(y_theo, model_freqs[:,0], marker='x')
    plt.scatter(y_theo_p1, model_freqs[:,1], marker='x')
    plt.scatter(y_theo_n1, model_freqs[:,2], marker='x')
    
def plot_results(pds, summary, l1_peaks, radial_order_range, dpi, coupling, eps_g, split, d01=None, use_heights=True):
    
    freqs = frequencies.Frequencies(frequency=pds.frequency.values,
                                numax=summary.numax.values, 
                                delta_nu=summary.DeltaNu.values if np.isfinite(summary.DeltaNu.values) else None, 
                                epsilon_p=summary.eps_p.values if np.isfinite(summary.eps_p.values) else None,
                                alpha=summary.alpha.values if np.isfinite(summary.alpha.values) else None,
                                radial_order_range=radial_order_range)
    params = {'calc_l0': True, # Compute radial mode properties
            'calc_l2': True, # Compute l=2 mode properties
            'calc_l3': False, # Don't need to calculate l=3 theoretical freqs
            'calc_nom_l1': True, # Compute nominal l=1 p-mode properties
            'calc_mixed': True, # Don't compute mixed modes (as not needed)
            'calc_rot': True, # Don't compute rotation
            'DPi1': dpi,
            'coupling': coupling,
            'eps_g': eps_g,
            'split_core': split,
            'split_env': 0.0,
            'l': 1, # Mixed modes are dipole mixed modes
            }
    if d01 is not None:
        params["d01"] = d01

    # Make computation - in our case this is for the computation of zeta
    freqs(params)

    #print(freqs.DPi1)
    freqs.generate_tau_values()

    new_peaks_tau = mixed_modes_utils.peaks_stretched_period(l1_peaks.frequency.values, 
                                                                pds.frequency.values, 
                                                                freqs.tau)
    new_peaks_zeta = mixed_modes_utils.peaks_stretched_period(l1_peaks.frequency.values, 
                                                                pds.frequency.values, 
                                                                freqs.zeta)
    new_peaks_tau -= freqs.shift*freqs.DPi1
    plt.figure(figsize=(12,8))
    plot_stretched_echelle(np.c_[freqs.l1_mixed_freqs, freqs.l1_mixed_freqs_p1, freqs.l1_mixed_freqs_n1], 
                           np.c_[freqs.l1_mixed_tau, freqs.l1_mixed_tau_p1, freqs.l1_mixed_tau_n1], 
                           l1_peaks.frequency.values, new_peaks_tau, freqs.DPi1, shift=0.0, heights=l1_peaks.height*5)
    # Don't include shift in plot of echelle as shift already included
    for i in freqs.l1_nom_freqs:
        plt.axhline(i, color='k', linestyle='--', alpha=0.5)
    plt.xlabel(r'$\tau\;\mathrm{mod}\;\Delta\Pi_{1}$ (s)', fontsize=18);
    plt.ylabel(r'Frequency ($\mu$Hz)', fontsize=18);
    
    return freqs, new_peaks_tau
