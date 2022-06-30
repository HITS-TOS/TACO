import os

import emcee
import lib.background.KeplerLCBgFit
import lib.background.mESS as mESS
import numpy as np
import pandas as pd


def background_fit(pds):
    """
    Background fitting using emcee's MCMC procedure

    Parameters:
        pds(pandas.DataFrame):Periodogram
            Columns:
                Name: frequency, dtype: float[micro-Hertz]
                Name: power, dtype: float
    """

    kw = {"nwalkers" : argv.nwalkers,
          "nwarmup"  : argv.nwarmup,
          "minsteps" : argv.minsteps,
          "maxsteps" : argv.maxsteps,
          "bins"     : argv.bins,    # If this is None or -1, we don't bin the pds
          "save_posteriors": argv.save_posteriors,
          "backend_filename"  : argv.posterior}

    # Fetch background model
    bkg_model = getattr(lib.background.KeplerLCBgFit, argv.bkg_model)

    numax0 = pd.read_csv(argv.summary)["numax0"][0]
    nuNyq = pd.read_csv(argv.summary)["nuNyq"][0]
    bgFit = bkg_model(pds, numax0, nuNyq, logfile=argv.logfile)
    minESS = mESS.minESS(bgFit.ndim, alpha=0.05, eps=0.1)
    i=0
    done_p = False

    print("Starting initial MCMC with binned PDS. Number of bins: %s" % kw["bins"])
    bgFit.MCMC(bgFit.bg_params, **kw)  # MCMC with binned PDS
    print("Finished initial MCMC with binned PDS")

    chain_i = np.argmax(bgFit.MCMC_sampler.lnprobability[:,-1])
    theta0 = bgFit.MCMC_sampler.chain[chain_i,-1,:]

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

    if not done_p:
        print("Giving up")
        flag_file = os.path.join(os.path.dirname(argv.pds), "BACKGROUND_FIT_FLAG")
        open(flag_file, "a").close()
        return
    else:
        # Create background_subtracted spectra
        summary = pd.read_csv(argv.summary)
        # Read in pds and oversampled pds
        ofac_pds = pd.read_csv(argv.ofac_pds)

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

        Hmax = full_model[idx_closest_numax].values
        Bmax = model[idx_closest_numax].values
        HBR = Hmax / Bmax

        for idx, row in bkg_summary.iterrows():
            summary[row['parameter']] = row['Q50']

        summary.to_csv(argv.summary, index=False)
        bkg_summary.to_csv(argv.quantiles, index=False)
        pds_bgr.to_csv(argv.pds_bgr, index=False)
        ofac_pds_bgr.to_csv(argv.ofac_bgr, index=False)        

    return Hmax, Bmax, HBR
