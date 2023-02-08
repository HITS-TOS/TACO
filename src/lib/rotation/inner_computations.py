import numpy as np
import sklearn

from sklearn.neighbors import KDTree 
from libs.sloscillations.sloscillations import frequencies, mixed_modes_utils

from tqdm import tqdm
from typing import List, Optional, Tuple, Union

#from numba import njit

def compute_shift(X_theo: np.ndarray, real_tau: np.ndarray, real_freqs: np.ndarray, 
                  dpi1: float, alpha_tau: float, 
                  tree: sklearn.neighbors._kd_tree.KDTree, 
                  agg_fn=np.median,
                  heights: Optional[np.ndarray] = None) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute the 'shift' between a set of theoretical mixed mode frequencies and a set of 
    detected peaks.
    
    Parameters
    -----------
    X_theo: np.ndarray
        Array with dimension (N,2) where N is the number of theoretical mixed modes. The array
        contains the tau values of the theoretical mixed modes in the first columns and the 
        mixed mode frequencies in the second.
    
    real_tau: np.ndarray
        Array containing the predicted tau values of the detected frequencies.
        
    real_freqs: np.ndarray
        Array containing the detected frequencies.
        
    dpi1: float
        Period spacing used for the mixed mode computation.
        
    alpha_tau: float
        Tau shift computed for the mixed modes.
    
    tree: sklearn.neighbors._kd_tree.KDTree
        KDTree used to determine the distances between the theoretical mixed modes
        and the detected frequencies.
    
    heights: np.ndarray (Optional)
        Heights of the detected frequencies. If given these are used to computed a height-weighted
        average of the distances.
        
    Outputs
    -------
    shift: np.ndarray
        Array of the shift values used in the computation.
        
    distances: np.ndarray
        Array of minimum distance between each detected frequency and the theoretical mixed mode. 
        This corresponds to the closest theoretical frequency to each detected frequency.
    
    """
    shift = np.linspace(-0.5, 0.5, 100)
    distances = np.zeros(100)
    matches = np.zeros(100)
    for i in range(len(shift)):
        y_real = (real_tau - dpi1*(alpha_tau + shift[i]))
        dist, idx = tree.query(np.c_[y_real, real_freqs], k=1)
        matches[i] = len(np.unique(idx)) / len(real_freqs)
        # Can we make this a weighted average!
        if heights is None:
            distances[i] = agg_fn(dist)
        else:
            weights = heights / np.sum(heights)
            distances[i] = np.average(dist.ravel(), weights=weights)
    return shift, distances, matches

def inner_loop(model_freqs: np.ndarray, model_tau: np.ndarray, 
               real_freqs: np.ndarray, real_tau: np.ndarray, 
               freqs: np.ndarray, agg_fn = np.median, heights: Optional[np.ndarray] = None, 
               n_comps: Optional[int] = 3) -> Tuple[float, float]:
    """
    The inner computational loop in the determination of the rotational splitting.
    
    Parameters
    ----------
    model_freqs: np.array
        Array containing the theoretical mixed mode frequencies.
        
    model_tau: np.array
        Array containing the theoretical mixed mode tau values. 
    
    real_tau: np.array
        Array containing the predicted tau values of the detected frequencies.
        
    real_freqs: np.array
        Array containing the detected frequencies.   
        
    freqs: sloscillations.frequency object
        Contains all the information about the theoretical mixed modes.
        
    heights: np.array (Optional)
        Heights of the detected frequencies. If given these are used to computed a height-weighted
        average of the distances.
        
    Outputs
    -------
    float:
        distance between best fitting theoretical frequencies and observed frequencies.
    float:
        shift corresponding to smallest distance.
    """

    
    metric = 'euclidean'
    if n_comps == 1:
        X = np.c_[model_tau[:,0], 
                  model_freqs[:,0]]
        tree = KDTree(X, metric=metric)

        shift, distances, matches = compute_shift(X, real_tau, real_freqs, freqs.DPi1, freqs.shift, tree, agg_fn, heights)

        return np.min(distances), shift[np.argmin(distances)], matches[np.argmin(distances)]
    
    elif n_comps == 2:
        
        # 90 degree case
        X = np.c_[np.r_[ 
                        model_tau[:,1] - freqs.shift * freqs.DPi1, 
                        model_tau[:,2] - freqs.shift * freqs.DPi1], 
                  np.r_[ 
                        model_freqs[:,1], 
                        model_freqs[:,2]]]
        tree = KDTree(X, metric=metric)

        shift, distances, matches = compute_shift(X, real_tau, real_freqs, freqs.DPi1, freqs.shift, tree, agg_fn, heights)

        return np.min(distances), shift[np.argmin(distances)], matches[np.argmin(distances)]
       
    elif n_comps == 3:

        # 45 degree case
        X = np.c_[np.r_[model_tau[:,0], 
                        model_tau[:,1] - freqs.shift * freqs.DPi1, 
                        model_tau[:,2] - freqs.shift * freqs.DPi1], 
                  np.r_[model_freqs[:,0], 
                        model_freqs[:,1], 
                        model_freqs[:,2]]]
        tree = KDTree(X, metric=metric)

        shift, distances, matches = compute_shift(X, real_tau, real_freqs, freqs.DPi1, freqs.shift, tree, agg_fn, heights)
   
        
        return np.min(distances), shift[np.argmin(distances)], matches[np.argmin(distances)]

    else:
        sys.exit("Incorrect number of components given!")
        
def compute_rotational_splitting(frequency: np.ndarray, numax: float, real_frequencies: np.ndarray, 
                                 dpi_range: np.ndarray, split_range: np.ndarray, 
                                 delta_nu: Optional[float] = None, 
                                 epsilon_p: Optional[float] = None, 
                                 alpha: Optional[float] = None, 
                                 radial_order_range: Optional[Union[np.ndarray, List]] = None, 
                                 coupling: Optional[float] = None, 
                                 real_heights: Optional[np.ndarray] = None) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Find the rotational splitting that best fits the observed data.
    
    Parameters
    ----------
    frequency: np.ndarray
        Frequency array of power spectrum.
        
    numax: float
        Numax value of star in question.
            
    real_frequencies: np.array
        Array containing the detected mixed mode frequencies.   
    
    dpi_range: np.ndarray
        Period spacing values to use in finding rotational splitting.
        This allows some flexibility around original calculated value.
        
    split_range: np.ndarray
        Rotational splitting values to use.
        
    delta_nu: Optional[float] = None
        Delta nu value. If not given then computed from numax-dnu scaling relation.
        
    epsilon_p: Optional[float] = None
        Epsilon_p value. If not given then computed from scaling relation.
        
    alpha: Optional[float] = None
        Alpha value (curvature). If not given then computed from scaling relation.
        
    radial_order_range: Optional[Union[np.ndarray, List]] = None
        Range of radial orders to compute theoretical radial modes over.
        
    coupling: Optional[float] = None
        Mixed mode coupling value. If not given then set to 0.15 (for RGB). TODO: Should be evolutionary state dependent!
            
    real_heights: Optional[np.array] = None
        Heights of the detected frequencies. If given these are used to computed a height-weighted
        average of the distances.
        
    Outputs
    -------
    full_distance np.ndarray:
        Closest distances of theoretical mixed modes to real frequencies. Size is relate to dpi_range.
        
    full_splitting np.ndarray:
        Best rotational splitting for a given dpi_range value.
        
    full_shift: np.ndarray:
        Best shift for a given dpi_range value.
    """
    
    # With the full (slow) looping it would look like this
    full_distances = np.zeros([len(dpi_range), 3])
    full_splitting = np.zeros_like(full_distances)
    full_shift = np.zeros_like(full_distances)
    full_matches = np.zeros_like(full_distances)

    for i in tqdm(range(len(dpi_range)), total=len(dpi_range)):

        distances = np.zeros([len(split_range), 3])
        shifts = np.zeros([len(split_range), 3])
        matches = np.zeros_like(shifts)

        freqs = frequencies.Frequencies(frequency=frequency,
                                    numax=numax, 
                                    delta_nu = delta_nu, 
                                    epsilon_p = epsilon_p,
                                    alpha = alpha,
                                    radial_order_range = radial_order_range)

        # Precompute mixed mode frequencies
        params = {'calc_l0': True, # Compute radial mode properties
                'calc_l2': True, # Compute l=2 mode properties
                'calc_l3': False, # Don't need to calculate l=3 theoretical freqs
                'calc_nom_l1': True, # Compute nominal l=1 p-mode properties
                'calc_mixed': True, # Don't compute mixed modes (as not needed)
                'calc_rot': False, # Don't compute rotation
                'DPi1': dpi_range[i],
                'coupling': coupling,
                'eps_g': 0.0, # Epsilon_g isn't needed for computation of tau due to chosen formulation of zeta
                'split_core': 0.0,
                'split_env': 0.0,
                'l': 1, # Mixed modes are dipole mixed modes
                }

        # Make computation - in our case this is for the computation of zeta
        freqs(params)

        # Compute m=0 tau values
        freqs.generate_tau_values()#params['DPi1'], params['coupling'], params['eps_g'])

        # Compute tau from the zeta value just computed
        new_peaks_tau = mixed_modes_utils.peaks_stretched_period(real_frequencies, 
                                                                 freqs.frequency, 
                                                                 freqs.tau)
        
        # Do m=0 component outside of loop
        distances[:, 0], shifts[:, 0], matches[:,0] = inner_loop(freqs.l1_mixed_freqs[:,None], 
                                                                   freqs.l1_mixed_tau[:,None], 
                                                                   real_frequencies,
                                                                   new_peaks_tau, freqs, 
                                                                   heights=real_heights, n_comps=1)

        # If splitting = 0.0 is included in given range
        if split_range[0] == 0.0:
            # Cases of no rotational splitting are identical to one ridge case, irrespective of number of ridges
            # All single ridge values are the same (as splitting isn't changing), therefore just assign first value.
            distances[0, 1], shifts[0, 1], matches[0, 1] = distances[0, 0].copy(), shifts[0, 0].copy(), matches[0,0].copy()
            distances[0, 2], shifts[0, 2], matches[0, 2] = distances[0, 0].copy(), shifts[0, 0].copy(), matches[0,0].copy()

            # Don't run for case of no splitting
            for j in range(1, len(split_range)):        
                freqs_p1 = freqs.l1_mixed_freqs + freqs.l1_zeta * split_range[j]
                freqs_n1 = freqs.l1_mixed_freqs - freqs.l1_zeta * split_range[j]

                tau_p1 = mixed_modes_utils.peaks_stretched_period(freqs_p1, freqs.frequency, freqs.tau)
                tau_n1 = mixed_modes_utils.peaks_stretched_period(freqs_n1, freqs.frequency, freqs.tau)

                model_freqs = np.c_[freqs.l1_mixed_freqs, freqs_p1, freqs_n1]
                model_tau = np.c_[freqs.l1_mixed_tau, tau_p1, tau_n1]



                distances[j, 1], shifts[j, 1], matches[j, 1] = inner_loop(model_freqs, model_tau, 
                                                        real_frequencies, 
                                                        new_peaks_tau, freqs, 
                                                        heights=real_heights, n_comps=2)   
                distances[j, 2], shifts[j, 2], matches[j, 2] = inner_loop(model_freqs, model_tau, 
                                                        real_frequencies, 
                                                        new_peaks_tau, freqs, 
                                                        heights=real_heights, n_comps=3)   
        else: # Case when splitting = 0uHz not included
            for j in range(len(split_range)):        
                freqs_p1 = freqs.l1_mixed_freqs + freqs.l1_zeta * split_range[j]
                freqs_n1 = freqs.l1_mixed_freqs - freqs.l1_zeta * split_range[j]

                tau_p1 = mixed_modes_utils.peaks_stretched_period(freqs_p1, freqs.frequency, freqs.tau)
                tau_n1 = mixed_modes_utils.peaks_stretched_period(freqs_n1, freqs.frequency, freqs.tau)

                model_freqs = np.c_[freqs.l1_mixed_freqs, freqs_p1, freqs_n1]
                model_tau = np.c_[freqs.l1_mixed_tau, tau_p1, tau_n1]



                distances[j, 1], shifts[j, 1], matches[j, 1] = inner_loop(model_freqs, model_tau, 
                                                        real_frequencies, 
                                                        new_peaks_tau, freqs, 
                                                        heights=real_heights, n_comps=2)   
                distances[j, 2], shifts[j, 2], matches[j, 2] = inner_loop(model_freqs, model_tau, 
                                                        real_frequencies, 
                                                        new_peaks_tau, freqs, 
                                                        heights=real_heights, n_comps=3)  


        full_distances[i, :] = np.min(distances, axis=0)
        full_splitting[i,:] = split_range[np.argmin(distances, axis=0)]
        full_shift[i,:] = shifts[np.argmin(distances, axis=0), [0,1,2]]
        full_matches[i,:] = matches[np.argmin(distances, axis=0), [0,1,2]]
    
    return full_distances, full_splitting, full_shift, full_matches
