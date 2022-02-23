## Background fit plot class for Kepler targets in the Long-Cadence mode

import argparse
import corner
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy

from KeplerLCBgFit import *
import PDSBgFit

def _sLor(nu, A, b, c):
    return A / (1 + (nu / b) ** c)

def _gauss(nu, A, b, c):
    return A * np.exp(-((nu - b)**2) / (2 * c**2))

def _sinc(x):
    return np.sinc(x/np.pi)

def bin_pds(pds, bins):
    assert(type(bins) is int)
    if (bins == -1):
        pds = self.pds
    else:
        assert(bins > 0)
        bin_pds  = scipy.stats.binned_statistic(
            x         = pds["frequency"],
            values    = pds["power"],
            statistic = 'mean',
            bins      = bins)
        bin_edges = bin_pds[1]
        bin_width = (bin_edges[1] - bin_edges[0])
        pds = {"frequency": bin_edges[1:] - bin_width/2,
               "power": bin_pds[0]}
    return pds

def _compute_model(nu, params):

    sc = _sinc(np.pi * nu / (2 * nu.max())) ** 2
    white = np.median(params['Pn'])
    comp1 = _sLor(nu, np.median(params['A1']),
                      np.median(params['b1']),
                      4)
    comp2 = _sLor(nu, np.median(params['A2']),
                      np.median(params['b2']),
                      4)
    comp3 = _sLor(nu, np.median(params['A3']),
                      np.median(params['b3']),
                      4)
    gaussian = _gauss(nu, np.median(params['Pg']),
                          np.median(params['numax']),
                          np.median(params['sigmaEnv']))

    bg = sc * comp1
    bg = bg + sc * comp2
    bg = bg + sc * comp3
    bg = bg + sc * gaussian
    bg = bg + white
    return bg, (comp1, comp2, comp3, gaussian, white)

def corner_plot(posteriors):

    param_names = list(posteriors)
    print(posteriors.head())
    corner.corner(posteriors.values, 
                  labels=param_names,
                  quantiles=[0.16, 0.5, 0.84],
                  show_titles=False)
                  
    plt.tight_layout()
    plt.show()

def plot_model(pds, posteriors):

    model, indiv_comps = _compute_model(pds['frequency'], posteriors)
    binned_pds = bin_pds(pds, 300)
    plt.plot(pds['frequency'], pds['power'], 'k')
    plt.plot(binned_pds['frequency'], binned_pds['power'], 'b')
    plt.plot(pds['frequency'], model, color='g', lw=2)
    plt.plot(pds['frequency'], indiv_comps[0] * _sinc(np.pi * pds['frequency'] / (2 * pds['frequency'].max())) ** 2, color='C0', linestyle='--')
    plt.plot(pds['frequency'], indiv_comps[1] * _sinc(np.pi * pds['frequency'] / (2 * pds['frequency'].max())) ** 2, color='C0', linestyle='--')
    plt.plot(pds['frequency'], indiv_comps[2] * _sinc(np.pi * pds['frequency'] / (2 * pds['frequency'].max())) ** 2, color='C0', linestyle='--')
    plt.plot(pds['frequency'], indiv_comps[3] * _sinc(np.pi * pds['frequency'] / (2 * pds['frequency'].max())) ** 2, color='r', linestyle='--')
    plt.plot(pds['frequency'], indiv_comps[4]*np.ones_like(pds['frequency']), 
                               color='C1', linestyle='--')
    plt.xscale('log')
    plt.xlim(0.1, 283)
    plt.ylim(pds['power'].min() * 0.1, pds['power'].max() * 10)
    plt.yscale('log')
    plt.xlabel(r'Frequency ($\mu$Hz)', fontsize=18)
    plt.ylabel(r'Power (ppm$^{2}\mu$Hz$^{-1}$)', fontsize=18)
    plt.show()

#############
# Entry point
#############
def main(argv):

    pds = pd.read_csv(argv.pds)
    posteriors = pd.read_csv(argv.posterior)
    #print(posteriors.head())
    plot_model(pds, posteriors)
    #corner_plot(posteriors)

if __name__=="__main__":

    parser = argparse.ArgumentParser(description="Make an MCMC background fit.")
    parser.add_argument('--pds', dest='pds', default='pds.csv',
                        help="A csv file with the periodogram on which to make the fit")
    parser.add_argument('--posterior', dest='posterior', default='pds_fit_posterior.csv',
                        help="Destination on which to write the MCMC posterior (collapsed chains)")
    argv = parser.parse_args()
    main(argv)
