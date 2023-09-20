# Background fitting using emcee's MCMC procedure
import ast
import os
import argparse
import emcee
import numpy as np
import pandas as pd
import lib.background.mESS as mESS
import lib.background.KeplerLCBgFit

import matplotlib.pyplot as plt

from lib.background.KeplerLCBgFit import *


## ===========
## Entry point
## ===========

def main(argv):
    kw = {"nwalkers" : argv.nwalkers,
          "nwarmup"  : argv.nwarmup,
          "minsteps" : argv.minsteps,
          "maxsteps" : argv.maxsteps,
          "bins"     : argv.bins,    # If this is None or -1, we don't bin the pds
          "save_posteriors": argv.save_posteriors,
          "backend_filename"  : argv.posterior}

    # Set random number generator seed for reproducible results
    if argv.seed >= 0:
        np.random.seed(argv.seed)

    # Fetch background model
    bkg_model = getattr(lib.background.KeplerLCBgFit, argv.bkg_model)

    pds = pd.read_csv(argv.pds)
    numax0 = pd.read_csv(argv.summary)["numax0"][0]
    nuNyq = pd.read_csv(argv.summary)["nuNyq"][0]
    bgFit = bkg_model(pds, numax0, nuNyq, logfile=argv.logfile)
    minESS = mESS.minESS(bgFit.ndim, alpha=0.05, eps=0.1)
    i=0
    done_p = False
    flag = 0
    print("flag")
    print(flag)
    
    try:
        print("Starting initial MCMC with binned PDS. Number of bins: %s" % kw["bins"])
        bgFit.MCMC(bgFit.bg_params, **kw)  # MCMC with binned PDS
        print("Finished initial MCMC with binned PDS")
    except (ValueError, RuntimeError, TypeError, NameError):
        print("Background fit did not converge")
        flag = 1
        return

    print(flag)
    chain_i = np.argmax(bgFit.MCMC_sampler.lnprobability[:,-1])
    theta0 = bgFit.MCMC_sampler.chain[chain_i,-1,:]

    #print(bgFit.theta_to_dict(theta0))
    #sys.exit()
    while (i < 5 and not done_p):
        if (bgFit.MCMCp["mixed_p"]
            and bgFit.MCMCp["converged_p"]
            and bgFit.MCMCp["mESS"] >= minESS):
            done_p = True
            print("MCMC with binned PDS done.")
        else:
            chain_i = np.argmax(bgFit.MCMC_sampler.lnprobability[:,-1])
            theta0 = bgFit.MCMC_sampler.chain[chain_i,-1,:]
            #Pn, A1, b1, A2, b2, A3, b3, Pg, numax, sigmaEnv = theta0
            iguess = bgFit.theta_to_dict(theta0)
            #iguess = {"Pn":Pn, "A1":A1, "b1":b1, "A2":A2, "b2":b2, "A3":A3, "b3":b3,
            #          "Pg":Pg, "numax":numax, "sigmaEnv":sigmaEnv}
            print("Attempting again (%s/5) the MCMC with binned PDS. Number of bins: %s" %
                  (i, kw["bins"]))
            bgFit.MCMC(iguess, **kw)
            i = i+1
    #bgFit.write_chains(
    #    filename = argv.chains, nwarmup = bgFit.MCMCp["nwarmup"])
    #bgFit.write_posterior(
    #    filename = argv.posterior, nwarmup = bgFit.MCMCp["nwarmup"])
    #bgFit.write_chains(
    #    filename = argv.chains + "_full")
    if not done_p:
        print("Giving up")
        flag_file = os.path.join(os.path.dirname(argv.pds), "BACKGROUND_FIT_FLAG")
        open(flag_file, "a").close()
        flag = 1
        #summary = pd.read_csv(argv.summary)
        #ofac_pds = pd.read_csv(argv.ofac_pds)
    else:
        # Create background_subtracted spectra
        #summary = pd.read_csv(argv.summary)
        # Read in pds and oversampled pds
        #ofac_pds = pd.read_csv(argv.ofac_pds)

        if argv.save_posteriors:
            # Read in chains
            reader = emcee.backends.HDFBackend(argv.posterior, read_only=True)

            #print(reader.get_autocorr_time())
            # Flattened chains and log-probability
            flatchain = reader.get_chain(discard=argv.nwarmup, flat=True)
            lnprob = reader.get_log_prob(discard=argv.nwarmup, flat=True)
        else:
            # Flattened chains and log-probability         
            #tau = bgFit.MCMC_sampler.get_autocorr_time()
            #print(f"Autocorrelation time: {tau}")
            #print(argv.nwarmup)
            flatchain = bgFit.MCMC_sampler.get_chain(discard=argv.nwarmup, flat=True)
            lnprob = bgFit.MCMC_sampler.get_log_prob(discard=argv.nwarmup, flat=True)

        posteriors = np.c_[flatchain, lnprob]

        # Summary of posterior distribution
        names = bgFit.bg_param_names()
        names.append("lnprob")
        df = pd.DataFrame(data=posteriors, columns=names)

        bkg_summary = df.describe(percentiles=[0.16, 0.50, 0.84]).T
        bkg_summary = bkg_summary.drop(['count', 'min', 'max'], axis=1)

        bkg_summary.reset_index(level=0, inplace=True)
        bkg_summary.rename(columns={'index': 'parameter', 'std': 'sd',
                                    '16%': 'Q16', '50%': 'Q50', '84%': 'Q84'},
                           inplace=True)

        # Remove the background
        nuNyq = summary['nuNyq'].values
        bg_parameters = bkg_summary['Q50']

        numax_index = bgFit.bg_param_names().index("numax")

        # Go to -1 because of logprob added at end
        model = bgFit.bgModel(theta=bg_parameters[:-1], nu=pds.frequency, no_osc=True)
        full_model = bgFit.bgModel(bg_parameters[:-1], pds.frequency)
        ofac_model = bgFit.bgModel(bg_parameters[:-1], ofac_pds.frequency, no_osc=True)

        pds_bgr = pds.assign(power = pds['power'] / model)
        ofac_pds_bgr = ofac_pds.assign(power = ofac_pds['power'] / ofac_model)

        idx_closest_numax = np.where(abs(pds['frequency'].values - bg_parameters[numax_index]) == np.min(np.abs(pds['frequency'].values - bg_parameters[numax_index])))[0]

        # Update summary file
        summary['Hmax'] = full_model[idx_closest_numax].values
        summary['Bmax'] = model[idx_closest_numax].values
        summary['HBR'] = full_model[idx_closest_numax].values / model[idx_closest_numax].values

        for idx, row in bkg_summary.iterrows():
            summary[row['parameter']] = row['Q50']

        # Write out files
        summary.to_csv(argv.summary, index=False)
        bkg_summary.to_csv(argv.quantiles, index=False)
        pds_bgr.to_csv(argv.pds_bgr, index=False)
        ofac_pds_bgr.to_csv(argv.ofac_bgr, index=False)
        


## ======================
## Command-line arguments
## ======================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Make an MCMC background fit.")
    parser.add_argument('--pds', dest='pds', default='pds.csv',
                        help="A csv file with the periodogram on which to make the fit")
    parser.add_argument('--ofac_pds', dest='ofac_pds', default='ofac_pds.csv',
                        help="A csv file with the oversampled periodogram, used when creating background-subtracted spectrum")
    parser.add_argument('--summary', dest='summary', default='summary.csv',
                        help="A csv file containing 'numax0' (initial numax estimation) and 'nuNyq' (Nyquist frequency).")
    parser.add_argument('--save_posteriors', dest='save_posteriors', action='store_true',
                        help="Whether or not MCMC posteriors (full chains) should be saved or not.")
    parser.add_argument('--posterior', dest='posterior', default='pds_fit_posterior.h5',
                        help="Destination on which to write the MCMC posteriors (full chains)")
    parser.add_argument('--logfile', dest='logfile', default='pds_fit.log',
                        help="Location on which to write the MCMC log")
    # MCMC parameters
    parser.add_argument('--nwalkers', dest='nwalkers', default=50, type=int,
                        help="Number of walkers (chains) for the MCMC fit")
    parser.add_argument('--nwarmup', dest='nwarmup', default=1000, type=int,
                        help="Number of steps for the MCMC warmup")
    parser.add_argument('--minsteps', dest='minsteps', default=2000, type=int,
                        help="Minimum number of steps for the MCMC estimation")
    parser.add_argument('--maxsteps', dest='maxsteps', default=5000, type=int,
                        help="Maximum number of steps for the whole MCMC run")
    parser.add_argument('--bins', dest='bins', default=-1, type=int,
                        help="The PDS will be binned by this number of bins for the MCMC. Setting it to -1 will not bin it.")
    parser.add_argument('--bkg_model', dest='bkg_model', default="KeplerBg3Comp", type=str,
                        help="The background model to fit to the data.")
    parser.add_argument("--pds_bgr", dest="pds_bgr", default = "pds_bgr.csv",
                    help = "File name on which to save the background-corrected PDS",
                    type = str)
    parser.add_argument("--ofac_bgr", dest="ofac_bgr", default = "ofac_pds_bgr.csv",
                    help = "File name on which to save the background-corrected oversampled PDS",
                    type = str)
    parser.add_argument("--out_backg_model", dest="out_backg_model", default = "out_backg_model.pkl",
                    help = "Save the background model out to file.",
                    type = str)
    parser.add_argument("--quantiles", dest="quantiles", default = "pds_fit_quantiles.csv",
                    help = "File on which to save the 16-50-84 quantiles of the MCMC chains",
                    type = str)
    parser.add_argument("--seed", dest="seed", default = -1, type=int,
                    help = "Random number generator seed for reproducible results (Default: -1).")
    argv = parser.parse_args()
    main(argv)
else:
    class argv:
        pds       = "pds.csv"
        summary   = "summary.csv"
        posterior = 'pds_fit_posterior.h5'
        logfile   = 'pds_fit.log'
        nwalkers  = 50
        nwarmup   = 1000
        minsteps  = 2000
        maxsteps  = 5000
        bins      = -1
        seed      = -1
