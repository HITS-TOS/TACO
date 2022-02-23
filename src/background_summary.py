## Extract some useful quantities from the MCMC performed by background_fit.py.
## Estimate the background function of the PDS and save a background-removed PDS


import argparse
import emcee
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pathlib import Path

def validate_arguments(argv):
    """Raises an `IOError` if file given in argv doesn't exist.
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
        # pds_bgr file and oversampled version don't exist yet so don't check 
        if (arg != 'pds_bgr') and (arg != 'ofac_bgr') and (arg != 'cov') and (arg != 'quantiles'):
            if type(getattr(argv, arg)) == str:
                # If file doesn't exists then raise IOError
                if not Path(getattr(argv, arg)).is_file():
                    raise IOError(f"{arg} file {getattr(argv, arg)} does not exist!")

def Harvey(x, a, b, c=4):
    return a/(1 + (x/b)**c)

def Gaussian(x, a, b, c):
    return a * np.exp(-(x-b)**2/(2*c**2))

def eta_sq(x, nyq):
    return np.sinc(x/(2*nyq))**2

if __name__=="__main__":
    
    # Get the filenames and options
    parser = argparse.ArgumentParser(description="Take the chains of an MCMC background fit, extract soome statistics and produce a background-removed PDS. The background parameters are estimated as the median of the posterior distribution for each parameter. The summary file is updated with these median values")
    parser.add_argument("--pds", dest="pds", default = "pds.csv",
                        help = "File name of the PDS in csv format with headers 'frequency' and 'power'",
                        type = str)
    parser.add_argument("--ofac_pds", dest="ofac_pds", default = "ofac_pds.csv",
                    help = "File name of the oversampled PDS in csv format with headers 'frequency' and 'power'",
                    type = str)
    parser.add_argument("--posterior", dest="posterior", default = "pds_fit_posterior.h5",
                    help = "A csv file cointaining the collapsed background fit chains",
                    type = str)
    parser.add_argument("--summary", dest="summary", default = "summary.csv",
                    help = "A csv file on which to append the background fit parameter estimates",
                    type = str)
    parser.add_argument("--pds_bgr", dest="pds_bgr", default = "pds_bgr.csv",
                    help = "File name on which to save the background-corrected PDS",
                    type = str)
    parser.add_argument("--ofac_bgr", dest="ofac_bgr", default = "ofac_pds_bgr.csv",
                    help = "File name on which to save the background-corrected oversampled PDS",
                    type = str)
    parser.add_argument("--bkg_names", dest="bkg_names", default = "background_parameter_names.txt",
                    help = "File name where the background parameters are kept.",
                    type = str)
    parser.add_argument("--cov", dest="cov", default = "pds_fit_covariance.csv",
                    help = "File on which to save the covariance matrix of the MCMC chains",
                    type = str)
    parser.add_argument("--quantiles", dest="quantiles", default = "pds_fit_quantiles.csv",
                    help = "File on which to save the 16-50-84 quantiles of the MCMC chains",
                    type = str)

    argv = parser.parse_args()

    # Check files exist by validating arguments in argparser
    validate_arguments(argv)

    # Read in all the files
    pds = pd.read_csv(argv.pds)
    ofac_pds = pd.read_csv(argv.ofac_pds)
    summary = pd.read_csv(argv.summary)
    bkg_names = np.loadtxt(argv.bkg_names, dtype=str)

    # Burnin kept as default value of 1000
    # NEEDS TO BE MADE A VARIABLE INCASE DEFAULT VALUE ISN'T SUITABLE
    burnin = 1000

    # Read in chains
    reader = emcee.backends.HDFBackend(argv.posterior, read_only=True)

    # Flattened chains and log-probability
    flatchain = reader.get_chain(discard=burnin, flat=True)
    lnprob = reader.get_log_prob(discard=burnin, flat=True)
    posteriors = np.c_[flatchain, lnprob]

    # Summary of posterior distribution
    df = pd.DataFrame(data=posteriors, columns=bkg_names)

    bkg_summary = df.describe(percentiles=[0.16, 0.50, 0.84]).T
    bkg_summary = bkg_summary.drop(['count', 'min', 'max'], axis=1)
    bkg_summary.reset_index(level=0, inplace=True)
    bkg_summary.rename(columns={'index': 'parameter', 'std': 'sd',
                                '16%': 'Q16', '50%': 'Q50', '84%': 'Q84'},
                       inplace=True)

    # Correlation matrix of background parameters
    """
    ## Correlations of background parameters
    fit.corr <- cor(as.matrix(fit.posterior))
    fit.corr.order <- hclust(dist(fit.corr), method = "complete")$order
    fit.corr <- fit.corr[fit.corr.order, fit.corr.order]
    """

    # Remove the background
    nuNyq = summary['nuNyq'].values
    bg_parameters = bkg_summary['Q50'].values

    bg = eta_sq(pds['frequency'], nuNyq) * (Harvey(pds['frequency'], bg_parameters[1], bg_parameters[2], 4) + 
                                            Harvey(pds['frequency'], bg_parameters[3], bg_parameters[4], 4) + 
                                            Harvey(pds['frequency'], bg_parameters[5], bg_parameters[6], 4)) + bg_parameters[0]
    pds_bgr = pds.assign(power = pds['power'] / bg)

    ofac_bg = eta_sq(ofac_pds['frequency'], nuNyq) * (Harvey(ofac_pds['frequency'], bg_parameters[1], bg_parameters[2], 4) + 
                                                      Harvey(ofac_pds['frequency'], bg_parameters[3], bg_parameters[4], 4) + 
                                                      Harvey(ofac_pds['frequency'], bg_parameters[5], bg_parameters[6], 4)) + bg_parameters[0]
    ofac_bgr = ofac_pds.assign(power = ofac_pds['power'] / ofac_bg)

    full_bg = eta_sq(pds['frequency'], nuNyq) * (Harvey(pds['frequency'], bg_parameters[1], bg_parameters[2], 4) + 
                                                 Harvey(pds['frequency'], bg_parameters[3], bg_parameters[4], 4) + 
                                                 Harvey(pds['frequency'], bg_parameters[5], bg_parameters[6], 4) +
                                                 Gaussian(pds['frequency'], 
                                                          bg_parameters[7], 
                                                          bg_parameters[8],
                                                          bg_parameters[9])) + bg_parameters[0] 
    idx_closest_numax = np.where(abs(pds['frequency'].values - bg_parameters[8]) == np.min(np.abs(pds['frequency'].values - bg_parameters[8])))

    # Update summary file
    summary['Hmax'] = full_bg.values[idx_closest_numax]
    summary['Bmax'] = bg.values[idx_closest_numax]
    summary['HBR'] = full_bg.values[idx_closest_numax] / bg.values[idx_closest_numax]

    for idx, row in bkg_summary.iterrows():
        summary[row['parameter']] = row['Q50']

    # Write out files
    summary.to_csv(argv.summary, index=False)
    bkg_summary.to_csv(argv.quantiles, index=False)
    pds_bgr.to_csv(argv.pds_bgr, index=False)
    ofac_bgr.to_csv(argv.ofac_bgr, index=False)