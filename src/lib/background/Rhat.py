# Potential scale reduction factor by Gelman and Rubin 1992

import numpy as np

def Rhat(chains):
    """
    Calculate Rhat (Gelman and Rubin 1992, Statistical Science, 7, 457) for
    a given set of chains and posterior. This should always be greater than 1.
    It is near 1 when the chains have converged.
    """
    n = chains.shape[1]  # number of iterations
    m = chains.shape[0]  # number of chains (walkers)
    global_mean = np.mean(chains)
    chains_mean = np.mean(chains, 1)
    Bn = np.sum((chains_mean - global_mean) ** 2) / (m - 1)
    W = np.sum([np.sum([(chains[j, t] - chains_mean[j]) ** 2
                        for t in range(n)])
                for j in range(m)]) / (m * (n - 1))
    Rhat = ((m + 1) * (n - 1) / (m * n)) + ((m + 1) * Bn / (m * W)) - ((n - 1) / (m * n))
    return Rhat
