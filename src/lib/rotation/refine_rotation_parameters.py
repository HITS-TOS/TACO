import matplotlib.pyplot as plt
import numpy as np
import sklearn

from sklearn.neighbors import KDTree 
from sloscillations import frequencies, mixed_modes_utils

from tqdm import tqdm
from typing import List, Optional, Tuple, Union

from . import inner_computations
from . import rotation_utils

from copy import deepcopy

#from numba import njit

def find_rotational_splitting(pds_frequency, numax, l1_peaks, dpi_range, split_range, DeltaNu, eps_p, alpha, 
                                radial_order_range, coupling, real_heights=None, plot=False, return_distances=False):
    """
    DOCSTRING NEEDED!
    """
    full_distances, full_splitting, full_shift, full_matches = inner_computations.compute_rotational_splitting(pds_frequency, numax, l1_peaks.frequency.values, dpi_range, split_range, 
                                                                                              DeltaNu, eps_p, alpha, radial_order_range, 
                                                                                              coupling, real_heights=real_heights)#l1_peaks.height.values)

    AIC = rotation_utils.compute_AIC(full_distances)
    n_ridges, ridge_idx = rotation_utils.how_many_ridges(AIC, return_idx=True)
    print(f"AIC values for each ridge (1, 2, 3): {AIC}")
    print(f"Found {n_ridges} ridge(s)")

    # Extract dpi, splitting, eps_g
    refined_dpi = dpi_range[np.argmin(full_distances, axis=0)][ridge_idx]
    refined_splitting = full_splitting[np.argmin(full_distances, axis=0), [0,1,2]][ridge_idx]
    refined_eps_g = -1*full_shift[np.argmin(full_distances, axis=0), [0,1,2]][ridge_idx]

    # If also want distances etc.
    if return_distances:
        return refined_dpi, refined_splitting, refined_eps_g, n_ridges, (full_distances, full_splitting, full_shift, full_matches)
    return refined_dpi, refined_splitting, refined_eps_g, n_ridges


def refine_coupling(pds_frequency, numax, l1_peaks, coupling_range, dpi, eps_g, split, DeltaNu, eps_p, alpha,
                    radial_order_range, coupling, 
                    real_heights=None, plot=False, return_distances=False):
    
    full_distances_c = np.zeros_like(coupling_range)
    full_shifts_c = np.zeros_like(coupling_range)
    full_matches_c = np.zeros_like(coupling_range)

    for i in tqdm(range(len(coupling_range)), total=len(coupling_range)):

            
        freqs = frequencies.Frequencies(frequency=pds_frequency,
                                    numax=numax, 
                                    delta_nu=DeltaNu, 
                                    epsilon_p=eps_p,
                                    alpha=alpha,
                                    radial_order_range=radial_order_range)
        
        params = {'calc_l0': True, # Compute radial mode properties
            'calc_l2': True, # Compute l=2 mode properties
            'calc_l3': False, # Don't need to calculate l=3 theoretical freqs
            'calc_nom_l1': True, # Compute nominal l=1 p-mode properties
            'calc_mixed': True, # Don't compute mixed modes (as not needed)
            'calc_rot': True, # Don't compute rotation
            'DPi1': dpi,
            'coupling': coupling_range[i],
            'eps_g': eps_g,
            'split_core': split,
            'l': 1, # Mixed modes are dipole mixed modes
            }

        # Make computation - in our case this is for the computation of zeta
        freqs(params)
        freqs.generate_tau_values()

        # Compute tau from the zeta value just computed
        new_peaks_tau = mixed_modes_utils.peaks_stretched_period(l1_peaks.frequency.values, 
                                                                    pds_frequency, 
                                                                    freqs.tau)

        # The shift is already accounted for in the calculation of the theoretical frequencies, so don't add it in here
        y_theo = (freqs.l1_mixed_tau - freqs.DPi1/2) % freqs.DPi1  - freqs.DPi1/2
        y_theo_p1 = (freqs.l1_mixed_tau_p1 - freqs.DPi1/2) % freqs.DPi1  - freqs.DPi1/2
        y_theo_n1 = (freqs.l1_mixed_tau_n1 - freqs.DPi1/2) % freqs.DPi1  - freqs.DPi1/2

        X = np.c_[np.r_[freqs.l1_mixed_tau, 
                        freqs.l1_mixed_tau_p1, 
                        freqs.l1_mixed_tau_n1], 
                  np.r_[freqs.l1_mixed_freqs, 
                        freqs.l1_mixed_freqs_p1, 
                        freqs.l1_mixed_freqs_n1]]
        tree = KDTree(X)

        shift, distances, matches = inner_computations.compute_shift(X, new_peaks_tau, l1_peaks.frequency.values, freqs.DPi1, freqs.shift, tree)
        full_distances_c[i] = np.min(distances)
        full_shifts_c[i] = shift[np.argmin(distances)]
        full_matches_c[i] = matches[np.argmin(distances)]
        
    if return_distances:
        return coupling_range[np.argmin(full_distances_c)] , (full_distances_c, full_shifts_c, full_matches_c)
    return coupling_range[np.argmin(full_distances_c)] 

def refined_nominal_pmodes(pds_frequency, numax, l1_peaks, d01_range, dpi, coupling, eps_g, split, 
                           DeltaNu, eps_p, alpha, radial_order_range, 
                           agg_fn = np.mean,
                           real_heights=None, plot=False, return_distances=False):
    """
    
    """

    # Baseline d01
    freqs = frequencies.Frequencies(frequency=pds_frequency,
                        numax=numax, 
                        delta_nu=DeltaNu, 
                        epsilon_p=eps_p,
                        alpha=alpha,
                        radial_order_range=radial_order_range)

    params = {'calc_l0': True, # Compute radial mode properties
            'calc_l2': True, # Compute l=2 mode properties
            'calc_l3': False, # Don't need to calculate l=3 theoretical freqs
            'calc_nom_l1': True, # Compute nominal l=1 p-mode properties
            'calc_mixed': False, # Don't compute mixed modes (as not needed)
            'calc_rot': False, # Don't compute rotation
            'DPi1': dpi,
            'coupling': coupling,
            'eps_g': eps_g,
            'split_core': split,
            'l': 1, # Mixed modes are dipole mixed modes
            }
    freqs(params)
    full_distances_d = np.zeros([len(d01_range), len(freqs.n)])
    full_shifts_d = np.zeros([len(d01_range), len(freqs.n)])
    full_matches_d = np.zeros([len(d01_range), len(freqs.n)])

  
    original_d01 = deepcopy(freqs.d01)

    for k in range(len(freqs.n)):

        for i in tqdm(range(len(d01_range)), total=len(d01_range)):

            freqs = frequencies.Frequencies(frequency=pds_frequency,
                                    numax=numax, 
                                    delta_nu=DeltaNu, 
                                    epsilon_p=eps_p,
                                    alpha=alpha,
                                    radial_order_range=radial_order_range)

            new_d01 = deepcopy(original_d01)
            new_d01[k] = d01_range[i]

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
                'l': 1, # Mixed modes are dipole mixed modes
                'd01': new_d01
                }

            # Make computation - in our case this is for the computation of zeta
            freqs(params)
            freqs.generate_tau_values()
            
            # Compute tau from the zeta value just computed
            new_peaks_tau = mixed_modes_utils.peaks_stretched_period(l1_peaks.frequency.values, 
                                                                        pds_frequency, 
                                                                        freqs.tau)

            # The shift is already accounted for in the calculation of the theoretical frequencies, so don't add it in here
            y_theo = (freqs.l1_mixed_tau - freqs.DPi1/2) % freqs.DPi1  - freqs.DPi1/2
            y_theo_p1 = (freqs.l1_mixed_tau_p1 - freqs.DPi1/2) % freqs.DPi1  - freqs.DPi1/2
            y_theo_n1 = (freqs.l1_mixed_tau_n1 - freqs.DPi1/2) % freqs.DPi1  - freqs.DPi1/2

            X = np.c_[np.r_[freqs.l1_mixed_tau, 
                            freqs.l1_mixed_tau_p1, 
                            freqs.l1_mixed_tau_n1], 
                      np.r_[freqs.l1_mixed_freqs, 
                            freqs.l1_mixed_freqs_p1, 
                            freqs.l1_mixed_freqs_n1]]
            tree = KDTree(X)

            shift, distances, matches = inner_computations.compute_shift(X, new_peaks_tau, l1_peaks.frequency.values, freqs.DPi1, freqs.shift, tree, agg_fn=agg_fn)#, l1_peaks.height.values)
            full_distances_d[i,k] = np.min(distances)
            full_shifts_d[i,k] = shift[np.argmin(distances)]  
            full_matches_d[i,k] = matches[np.argmin(distances)]
            
    if plot:
        plt.plot(d01_range, full_distances_d)
        plt.ylabel("Distance", fontsize=18)
        plt.xlabel(r'$\delta\nu_{0,1}$ ($\mu$Hz)', fontsize=18)
        plt.show()
        
    if return_distances:
        return d01_range[np.argmin(full_distances_d, axis=0)].squeeze(), (full_distances_d, full_shifts_d, full_matches_d)
    return d01_range[np.argmin(full_distances_d, axis=0)].squeeze()