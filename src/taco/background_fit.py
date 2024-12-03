from pathlib import Path

import math
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
        self.output_ofac_pds_bgr = 'ofac_pds_bgr.csv'
        self.output_pds_bgr = 'pds_bgr.csv'
        self.output_quantiles = 'pds_fit_quantiles.csv'
        self.posterior = 'pds_fit_posterior.h5'
        self.save_posteriors = False
        self.seed = -1

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


def background_fit(pds, ofac_pds, data, output = '', output_directory = '', **kwargs):
    """
    Background fitting using emcee's MCMC procedure

    Parameters:
        pds(pandas.DataFrame):Periodogram
            Columns:
                Name: frequency, dtype: float[micro-Hertz]
                Name: power, dtype: float
        ofac_pds(pandas.DataFrame):Oversampled periodogram
            Columns:
                Name: frequency, dtype: float[micro-Hertz]
                Name: power, dtype: float
        data(pandas.DataFrame):Summary data
            Columns:
                initial_numax(float):Initial numax estimation
                nuNyq(float):Nyquist frequency
        kwargs(dict):

    Returns:
        pds_bgr(pandas.DataFrame): Periodogram, background corrected
        ofac_pds_bgr(pandas.DataFrame): Oversampled periodogram, background corrected
        data(pandas.DataFrame): Summary data
            Columns:
                Hmax(float): Power of background model at numax [ppm^2/micro-Hertz]
                Bmax(float): Power of background model (without gaussian envelope) at numax [ppm^2/micro-Hertz]
                HBR(float): Hmax over Bmax ratio
                Pn (float): White noise [ppm^2/micro-Hertz]
                Pn_sd (float): White noise standard deviation [ppm^2/micro-Hertz]
                Ai (float): Height of background component i [ppm^2/micro-Hertz]
                Ai_sd (float): Ai standard deviation [ppm^2/micro-Hertz]
                bi (float): Characteristic frequency of background component i [micro-Hertz]
                bi_sd (float): bi standard deviation [micro-Hertz]
                Pg (float): Gaussian envelope height [ppm^2/micro-Hertz]
                Pg_sd (float): Pg standard deviation [ppm^2/micro-Hertz]
                numax (float): Numax fit [micro-Hertz]
                numax_sd (float): numax standard deviation [micro-Hertz]
                sigmaEnv (float): [micro-Hertz]
                sigmaEnv_sd (float): [micro-Hertz]
        flag_bgr(int): Background fit flag
            Values:
                0: Background fit converged
                1: Background fit didn't converge
                2: Amplitude of power excess is too low
                3: Background overfitting

    """

    settings = Settings(**kwargs)
    if (max(pds['frequency']) > 10 * settings.bins):
        settings.bins = math.floor(0.1 * (max(pds['frequency'])))
    print(settings.bkg_model)

    # Set random number generator seed for reproducible results
    if settings.seed >= 0:
        np.random.seed(settings.seed)

    # Fetch background model
    bkg_model = getattr(lib.background.KeplerLCBgFit, settings.bkg_model)

    j = 1
    done_q = False

    while (j < 3 and not done_q):
        if (data['initial_numax_sd'][0] < data['initial_numax'][0]/10.0):
            data['initial_numax_sd'][0] = data['initial_numax'][0]/5.0

        if (data['initial_numax_sd'][0] > data['initial_numax'][0]/4.0):
            data['initial_numax_sd'][0] = data['initial_numax'][0]/4.0

        if (j > 1):
            data['initial_numax'][0] = np.nanmean([data['numax_var'][0],data['numax_Morlet'][0],data['numax_CWTMexHat'][0]])
            data['initial_numax_sd'][0] = max([data['numax_var'][0],data['numax_Morlet'][0],data['numax_CWTMexHat'][0]]) - np.nanmin([data['numax_var'][0],data['numax_Morlet'][0],data['numax_CWTMexHat'][0]])

        bg_fit = bkg_model(pds, data['initial_numax'].iloc[0], data['initial_numax_sd'].iloc[0], data['nuNyq'].iloc[0], logfile = Path(output_directory,settings.logfile))
        minESS = mESS.minESS(bg_fit.ndim, alpha=0.05, eps=0.1)
        i = 0
        done_p = False
        flag = 0

        try:
            print("Starting initial MCMC with binned PDS. Number of bins:", settings.bins)
            bg_fit.MCMC(bg_fit.bg_params, output_directory, **settings.get_mcmc_settings())  # MCMC with binned PDS
            print("Finished initial MCMC with binned PDS")
        except (ValueError, RuntimeError, TypeError, NameError):
            print("Background fit did not converge")
            flag = 1
            return None, None, data, flag

        chain_i = np.argmax(bg_fit.MCMC_sampler.get_log_prob()[-1,:]) #Check if needed
        theta0 = bg_fit.MCMC_sampler.get_chain()[-1,chain_i,:]

        while (i < 3 and not done_p):
            if (bg_fit.MCMCp["mixed_p"]
                and bg_fit.MCMCp["converged_p"]
                and bg_fit.MCMCp["mESS"] >= minESS):
                done_p = True
                print("MCMC with binned PDS done.")
            else:
                chain_i = np.argmax(bg_fit.MCMC_sampler.get_log_prob()[-1,:])
                theta0 = bg_fit.MCMC_sampler.get_chain()[-1,chain_i,:]
                #Pn, A1, b1, A2, b2, A3, b3, Pg, numax, sigmaEnv = theta0
                iguess = bg_fit.theta_to_dict(theta0)
                #iguess = {"Pn":Pn, "A1":A1, "b1":b1, "A2":A2, "b2":b2, "A3":A3, "b3":b3,
                #          "Pg":Pg, "numax":numax, "sigmaEnv":sigmaEnv}
                try:
                    print("Attempting again (%s/5) the MCMC with binned PDS. Number of bins: %s" %
                    (i, settings.bins))
                    print("Starting MCMC with binned PDS. Number of bins:", settings.bins)
                    bg_fit.MCMC(iguess, output_directory, **settings.get_mcmc_settings())
                    print("Finished MCMC with binned PDS")
                    i += 1

                except (ValueError, RuntimeError, TypeError, NameError):
                    print("Background fit did not converge")
                    flag = 1
                    return None, None, data, flag



        if not done_p:
            j+=1
        else:
            # Create background_subtracted spectra

            if settings.save_posteriors:
                # Read in chains
                # reader = emcee.backends.HDFBackend(settings.posterior, read_only=True)
                reader = emcee.backends.HDFBackend(Path(output_directory, settings.posterior), read_only=True)
                #print(reader.get_autocorr_time())
                # Flattened chains and log-probability
                flatchain = reader.get_chain(discard=settings.nwarmup, flat=True)
                lnprob = reader.get_log_prob(discard=settings.nwarmup, flat=True)
            else:
                # Flattened chains and log-probability
                #tau = bg_fit.MCMC_sampler.get_autocorr_time()
                #print(f"Autocorrelation time: {tau}")
                flatchain = bg_fit.MCMC_sampler.get_chain(discard=settings.nwarmup, flat=True)
                lnprob = bg_fit.MCMC_sampler.get_log_prob(discard=settings.nwarmup, flat=True)

            posteriors = np.c_[flatchain, lnprob]

            # Summary of posterior distribution
            names = bg_fit.bg_param_names()
            names.append("lnprob")
            df = pd.DataFrame(data=posteriors, columns=names)

            bkg_summary = df.describe(percentiles=[0.16, 0.50, 0.84]).T
            bkg_summary = bkg_summary.drop(['count', 'min', 'max'], axis=1)

            bkg_summary.reset_index(level=0, inplace=True)
            bkg_summary.rename(columns={'index': 'parameter', 'std': 'sd',
                                    '16%': 'Q16', '50%': 'Q50', '84%': 'Q84'},
                           inplace=True)

            # Remove the background
            bg_parameters = bkg_summary['Q50']

            numax_index = bg_fit.bg_param_names().index("numax")

            # Go to -1 because of logprob added at end
            model = bg_fit.bgModel(theta=bg_parameters[:-1], nu=pds.frequency, no_osc=True)
            full_model = bg_fit.bgModel(bg_parameters[:-1], pds.frequency)
            ofac_model = bg_fit.bgModel(bg_parameters[:-1], ofac_pds.frequency, no_osc=True)

            pds_bgr = pds.assign(power = pds['power'] / model)
            ofac_pds_bgr = ofac_pds.assign(power = ofac_pds['power'] / ofac_model)

            pds_bgr_full = pds.assign(power = pds['power'] / full_model)

            if np.mean(pds_bgr_full.power.values) < 0.95:
                print('Background fit needs to be checked!')
                flag = 3
                return None, None, None, flag


            if bg_parameters[numax_index] > data['nuNyq'].iloc[0]:
                flag = 4
                print('Numax larger than Nyquist frequency')


            idx_closest_numax = np.where(abs(pds['frequency'].values - bg_parameters[numax_index]) == np.min(np.abs(pds['frequency'].values - bg_parameters[numax_index])))[0]

            data['Hmax'] = full_model[idx_closest_numax].values
            data['Bmax'] = model[idx_closest_numax].values
            data['HBR'] = data['Hmax'] / data['Bmax']

            if (data['HBR'][0] > 1.05):
                done_q = True

                for idx, row in bkg_summary.iterrows():
                    data[row['parameter']] = row['Q50']
                    data[row['parameter']+'_sd'] = row['sd']

                bkg_summary.to_csv(Path(output_directory,settings.output_quantiles), index = False)
            else:
                j+=1



    if (done_p and not done_q):
        print("too low amplitude of power excess")
        flag = 2
        for idx, row in bkg_summary.iterrows():
            data[row['parameter']] = row['Q50']
            data[row['parameter']+'_sd'] = row['sd']

        bkg_summary.to_csv(Path(output_directory,settings.output_quantiles), index = False)

    if (not done_p and not done_q):
        print("Giving up")
        flag = 1
        return None, None, None, flag


    return pds_bgr, ofac_pds_bgr, data, flag
