import os

import emcee
import lib.background.KeplerLCBgFit
import lib.background.mESS as mESS
import numpy as np
import pandas as pd

class Settings(object):
    """ Turns kwargs dictionary into a settings object """
    def _setup_attrs(self):
        self.bins = 300
        self.bkg_model = 'KeplerBg3Comp'
        self.logfile = 'pds_fit.log'
        self.maxsteps = 5000
        self.minsteps = 2000
        self.nwalkers = 50
        self.nwarmup = 1000
        self.output_backg_model = 'out_backg_model.pkl'
        self.output_ofac_bgr = 'ofac_pds_bgr.csv'
        self.output_pds_pgr = 'pds_bgr.csv'
        self.output_quantiles = 'pds_fit_quantiles.csv'
        self.posterior = 'pds_fit_posterior.h5'
        self.save_posteriors = False

    def __init__(self, **kwargs):
        self._setup_attrs()
        assert kwargs.keys() <= self.__dict__.keys()
        self.__dict__.update(kwargs)
    
    def get_mcmc_settings(self):
        """ MCMC subset of settings """
        kwargs = {k: self.__dict__[k] for k in (
            'bins',
            'posterior',
            'save_posteriors',
            'maxsteps',
            'minsteps',
            'nwalkers',
            'nwarmup'
        )}
        # Rename keyword
        kwargs['backend_filename'] = kwargs.pop('posterior')
        return kwargs


def background_fit(pds, numax0, nuquist, **kwargs):
    """
    Background fitting using emcee's MCMC procedure

    Parameters:
        pds(pandas.DataFrame):Periodogram
            Columns:
                Name: frequency, dtype: float[micro-Hertz]
                Name: power, dtype: float
        numax0(float):
        nuquist(float):Nyquist frequency
        kwargs(dict):

    Returns:
        Hmax(float):
        Bmax(float):
        HBR(float):
    """

    settings = Settings(**kwargs)
    print(settings.bkg_model)

    # Fetch background model
    bkg_model = getattr(lib.background.KeplerLCBgFit, settings.bkg_model)

    bg_fit = bkg_model(pds, numax0, nuquist, logfile = settings.logfile)
    minESS = mESS.minESS(bg_fit.ndim, alpha=0.05, eps=0.1)
    i = 0
    done_p = False

    print("Starting initial MCMC with binned PDS. Number of bins:", settings.bins)
    bg_fit.MCMC(bg_fit.bg_params, **settings.get_mcmc_settings())  # MCMC with binned PDS
    print("Finished initial MCMC with binned PDS")

    chain_i = np.argmax(bg_fit.MCMC_sampler.lnprobability[:,-1])
    theta0 = bg_fit.MCMC_sampler.chain[chain_i,-1,:]

    # while (i < 5 and not done_p):
    #     if (bg_fit.MCMCp["mixed_p"]
    #         and bg_fit.MCMCp["converged_p"]
    #         and bg_fit.MCMCp["mESS"] >= minESS):
    #         done_p = True
    #         print("MCMC with binned PDS done.")
    #     else:
    #         chain_i = np.argmax(bg_fit.MCMC_sampler.lnprobability[:,-1])
    #         theta0 = bg_fit.MCMC_sampler.chain[chain_i,-1,:]
    #         #Pn, A1, b1, A2, b2, A3, b3, Pg, numax, sigmaEnv = theta0
    #         iguess = bg_fit.theta_to_dict(theta0)
    #         #iguess = {"Pn":Pn, "A1":A1, "b1":b1, "A2":A2, "b2":b2, "A3":A3, "b3":b3,
    #         #          "Pg":Pg, "numax":numax, "sigmaEnv":sigmaEnv}
    #         print("Attempting again (%s/5) the MCMC with binned PDS. Number of bins: %s" %
    #               (i, kw["bins"]))
    #         bg_fit.MCMC(iguess, **kw)
    #         i = i+1

    # if not done_p:
    #     print("Giving up")
    #     flag_file = os.path.join(os.path.dirname(kwargs.pds), "BACKGROUND_FIT_FLAG")
    #     open(flag_file, "a").close()
    #     return
    # else:
    #     # Create background_subtracted spectra
    #     summary = pd.read_csv(kwargs.summary)
    #     # Read in pds and oversampled pds
    #     ofac_pds = pd.read_csv(kwargs.ofac_pds)

    #     if kwargs.save_posteriors:
    #         # Read in chains
    #         reader = emcee.backends.HDFBackend(kwargs.posterior, read_only=True)

    #         #print(reader.get_autocorr_time())
    #         # Flattened chains and log-probability
    #         flatchain = reader.get_chain(discard=kwargs.nwarmup, flat=True)
    #         lnprob = reader.get_log_prob(discard=kwargs.nwarmup, flat=True)
    #     else:
    #         # Flattened chains and log-probability         
    #         #tau = bg_fit.MCMC_sampler.get_autocorr_time()
    #         #print(f"Autocorrelation time: {tau}")
    #         #print(kwargs.nwarmup)
    #         flatchain = bg_fit.MCMC_sampler.get_chain(discard=kwargs.nwarmup, flat=True)
    #         lnprob = bg_fit.MCMC_sampler.get_log_prob(discard=kwargs.nwarmup, flat=True)

    #     posteriors = np.c_[flatchain, lnprob]

    #     # Summary of posterior distribution
    #     names = bg_fit.bg_param_names()
    #     names.append("lnprob")
    #     df = pd.DataFrame(data=posteriors, columns=names)

    #     bkg_summary = df.describe(percentiles=[0.16, 0.50, 0.84]).T
    #     bkg_summary = bkg_summary.drop(['count', 'min', 'max'], axis=1)

    #     bkg_summary.reset_index(level=0, inplace=True)
    #     bkg_summary.rename(columns={'index': 'parameter', 'std': 'sd',
    #                                 '16%': 'Q16', '50%': 'Q50', '84%': 'Q84'},
    #                        inplace=True)

    #     # Remove the background
    #     nuquist = summary['nuNyq'].values
    #     bg_parameters = bkg_summary['Q50']

    #     numax_index = bg_fit.bg_param_names().index("numax")

    #     # Go to -1 because of logprob added at end
    #     model = bg_fit.bgModel(theta=bg_parameters[:-1], nu=pds.frequency, no_osc=True)
    #     full_model = bg_fit.bgModel(bg_parameters[:-1], pds.frequency)
    #     ofac_model = bg_fit.bgModel(bg_parameters[:-1], ofac_pds.frequency, no_osc=True)

    #     pds_bgr = pds.assign(power = pds['power'] / model)
    #     ofac_pds_bgr = ofac_pds.assign(power = ofac_pds['power'] / ofac_model)

    #     idx_closest_numax = np.where(abs(pds['frequency'].values - bg_parameters[numax_index]) == np.min(np.abs(pds['frequency'].values - bg_parameters[numax_index])))[0]

    #     Hmax = full_model[idx_closest_numax].values
    #     Bmax = model[idx_closest_numax].values
    #     HBR = Hmax / Bmax

    #     for idx, row in bkg_summary.iterrows():
    #         summary[row['parameter']] = row['Q50']

    #     bkg_summary.to_csv(kwargs.quantiles, index=False)
    #     pds_bgr.to_csv(kwargs.pds_bgr, index=False)
    #     ofac_pds_bgr.to_csv(kwargs.ofac_bgr, index=False)        

    # return Hmax, Bmax, HBR
    return 1.0, 1.0, 1.0