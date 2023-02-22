## Base PDS background fit class

import logging
import csv
import numpy as np
import scipy
import scipy.optimize as op
import emcee
from pathlib import Path
from . import mESS
from sklearn import linear_model
from .Rhat import Rhat

import matplotlib.pyplot as plt


class classproperty(object):
    def __init__(self, func):
        self.func = func
    def __get__(self, owner_self, owner_cls):
        return self.func(owner_cls)


class PDSBgFit(object):
    """
    The following methods must be overriden by any subclass:
        - logPrio
        - bg_param_names
    """

    pds=None
    messages=[]
    _maxPDS = None
    __logfile_default = 'pds_fit.log'
    MCMC_sampler        = None
    MCMCp = {
        # MCMC control parameters
        "step_size"      : 500,         # Advance the sampler this number of steps at a time
        "Rhat_max"       : 1.1,         # Maximum tolerance for Rhat
        "converge_test"  : 1.05,        # Maximum slope of the lnprob chains
        "alpha" : 0.05, "eps" : 0.1,    # If the MCMC run were to be repeated, its estimates will fall
                                        # within 'eps' of the true value for (1-'alpha') repetitions.
        # MCMC diagnostic flags
        "converged_p"    : False,        # Passed the convergence_test
        "mixed_p"        : False,        # Pass the Rhat test
        "nsamples_p"     : False,        # We have enough mESS
        # MCMC results
        "nwarmup"        : 0,            # Warm-up iterations
        "niters"         : 0,            # Total iterations
        "slope"          : np.inf,
        "mESS"           : 0,            # Multivariate Effective Sampling Size
        "Rhat"           : np.inf,       # Potential scale reduction factor
        "chains"         : None          # The actual chains at the end
    }


    @property
    def bg_params(self):
        return self.__bg_params


    @bg_params.setter
    def bg_params(self, params):
        if isinstance(params, dict):
            self.__bg_params = params
        elif isinstance(params, list):
            self.__bg_params = self.theta_to_dict(params)
        else:
            raise ValueError("Argument must be a dict or a list")


    @classmethod
    def bg_param_names(cls):
        """
        A list of strings representing the parameter names.
        """
        raise NotImplementedError


    @property
    def ndim(self):
        return len(self.dict_to_theta(self.bg_params))


    def __init__(self, pds, logfile=None):
        self.pds = pds
        if logfile == None:
            self.logfile = self.__logfile_default
        else:
            self.logfile = logfile
        self._maxPDS = max(pds["power"])
        self.bg_params = {}


    @classmethod
    def theta_to_dict(cls, theta):
        return dict(zip(cls.bg_param_names(), theta))


    @classmethod
    def dict_to_theta(cls, bgdic):
        return [bgdic[i] for i in cls.bg_param_names()]


    def _log_message(self, message):
        if not (self.logfile == None):
            logging.basicConfig(filename=self.logfile,
                                filemode='w',
                                level=logging.INFO,
                                format="%(asctime)s %(levelname)s %(message)s")
            logging.info(message)
        self.messages.append(message)


    # TODO  Throw a warning if the selected "bins" is a bad choice... for some definition of "bad"
    def __bin_pds(self, bins):
        assert(type(bins) is int)
        if (bins == -1):
            pds = self.pds
        else:
            assert(bins > 0)
            bin_pds  = scipy.stats.binned_statistic(
                x         = self.pds["frequency"],
                values    = self.pds["power"],
                statistic = 'mean',
                bins      = bins)
            bin_edges = bin_pds[1]
            bin_width = (bin_edges[1] - bin_edges[0])
            pds = {"frequency": bin_edges[1:] - bin_width/2,
                   "power": bin_pds[0]}
        return pds


    def __update_bg_params(self, theta):
        self.bg_theta  = theta
        self.bg_params = self.theta_to_dict(theta)


    def logProb(self, theta, pds, elements_per_bin = 1):
        lprio = self.logPrio(theta)
        if not np.isfinite(lprio):
            return -np.inf
        return lprio + self.logLikelihood(theta, pds, elements_per_bin)


    def logLikelihood(self, theta, pds, elements_per_bin = 1):
        freqs = pds["frequency"]
        power = pds["power"]
        # Faster likelihood calculation
        M = self.bgModel(theta, freqs)
        res = -1.0 * np.sum(np.log(M) + (power/M))
        if not np.isfinite(res):
            return -np.inf
        return elements_per_bin * res


    def logPrio(self, theta):
        raise NotImplementedError


    def bgModel(self, theta, freqi):
        """
        Background model for the 'pds' using the parameters in 'theta' at frequency 'freqi'.
        """
        raise NotImplementedError


    def MLE(self, iguess, bins=-1):
        """
        MLE optimization using 'iguess' as initial values, If 'bins' is provided
        the pds is binned with that many bins before doing the MLE optimization.
        Sets the result to self.bg_params
        """ #TODO  check if the optimization succeded before setting self.bg_params
        elements_per_bin = round(len(self.pds["frequency"]) / bins)
        pds = self.__bin_pds(bins)
        self._log_message("# Starting relaxation to local optima.")
        nll = lambda *args: -self.logProb(*args)
        MLE_result = op.minimize(nll, self.dict_to_theta(iguess), args=(pds, elements_per_bin))
        self.bg_params = MLE_result["x"].tolist()
        self._log_message("# Finished relaxation to local optima.")


    def __init_walkers_pos(self, iguess, nwalkers):
        theta0 = self.dict_to_theta(iguess)
        #pos = [theta0 * (1 + 1e-4 * np.random.randn(self.ndim)) for i in range(nwalkers)]
        pos = [theta0 + 1e-4 * np.random.randn(self.ndim) for i in range(nwalkers)]
        return pos


    def __get_chains(self, sampler, nwarmup=0, include_logprob=True):
        #chains = np.dstack((sampler.chain, sampler.lnprobability))
        if include_logprob:
            chains = np.dstack((sampler.get_chain(discard=nwarmup), 
                                sampler.get_log_prob(discard=nwarmup)))
        else:
            chains = sampler.get_chain(discard=nwarmup)
        return chains


    def __get_posterior(self, sampler, nwarmup=0):
        samples = self.__get_chains(sampler)
        return samples[:, nwarmup:, :].reshape((-1, self.ndim + 1))


    def __test_convergence(self, chains, par_idx):
        """
        Make a linear fit (iterations, value) and compare with the null model (median)
        using the ratio of the mean squared distances.
        Returns: float = mean_sqrt_null_model / mean_sqrt_linear_fit
        """
        flat_chain = {'iteration': [], 'value': []}
        for w in range(chains.shape[0]):
            for i in range(chains.shape[1]):
                flat_chain["iteration"].append(i)
                flat_chain["value"].append(chains[w, i, par_idx])
        fit_X = np.array(flat_chain['iteration'])
        fit_X = fit_X.reshape((-1, 1))
        fit_y = np.array(flat_chain['value'])
        regr = linear_model.LinearRegression()
        #try:
        regr.fit(fit_X, fit_y)
        null_model = np.median(flat_chain['value'])
        sqrt_null = np.sum([(fit_y[i] - null_model) ** 2 for i in range(fit_y.shape[0])])
        sqrt_fit = np.sum([(fit_y[i] - regr.predict(fit_X[i].reshape(-1,1))) ** 2 for i in range(fit_y.shape[0])])
        return sqrt_null / sqrt_fit



    def MCMC(self, iguess, output_directory, nwalkers, backend_filename, save_posteriors=False, nwarmup=1000, minsteps=1000, maxsteps=10000, bins=-1, alpha=0.05, eps=0.1):
        """
        Make an MCMC parameter estimation.
        Parameters:
            iguess      Initial guesses as a dict. The walkers will be initialized in a
                        small region around them
            nwalkers    Number of walkers to use. Must be more than 2*self.ndim
            nwarmup     Warm-up iterations. These will be discarded.
            minsteps    Minimum number of MCMC steps after warm-uü. Too low values might
                        cause issues with the mESS estimation.
            maxsteps    Maximum number of MCMC steps (on each chain) before giving up
            bins        If provided, the pds will be binned with this number of bins. Set this
                        to -1 to use the full pds without binning
            alpha, eps  If the MCMC run were to be repeated, its estimates will fall within
                        'eps' of the true value for (1-'alpha') number of repetitions.
        """
        ## TODO  Error if walkers <= 2*self.ndim
        # Update parameters
        self.MCMCp["alpha"]       = alpha
        self.MCMCp["eps"]         = eps
        self.MCMCp["converged_p"] = False
        self.MCMCp["mixed_p"]     = False
        self.MCMCp["nsamples_p"]  = False
        self.MCMCp["slope"]       = np.inf
        self.MCMCp["Rhat"]        = np.inf
        self.MCMCp["mESS"]        = 0
        self.MCMCp["nwarmup"]     = 0

        # TODO: Add in a way to use a frequency cut if needed

        #self.pds[["frequency", "power"]] = self.pds[["frequency", "power"]][self.pds["frequency"] >= 2]
        #self.pds.dropna(inplace=True, axis=0)

        if bins == -1:
            elements_per_bin = 1
        else:
            elements_per_bin = round(len(self.pds["frequency"]) / bins)

        pds = self.__bin_pds(bins)

        if save_posteriors:
            # Set up the backend
            backend = emcee.backends.HDFBackend(Path(output_directory,backend_filename))
            # Clear the backend in case the file already exists
            backend.reset(nwalkers, self.ndim)
        else:
            backend = None
        
        # Save out background parameter names - this is hacky because
        # emcee backend doesn't allow writing parameter names into hdf5 file
        # as far as I know
        #np.savetxt("background_parameter_names.txt", self.bg_param_names() + ["lnprob"], fmt="%s")
        
        self.MCMC_sampler = emcee.EnsembleSampler(nwalkers, self.ndim, self.logProb, backend=backend,
                                                      kwargs={"pds" : pds, "elements_per_bin" : elements_per_bin})
        iters = 0
        pos = self.__init_walkers_pos(iguess, nwalkers)

        if bins == -1:
            self._log_message("# Sampling with the full PDS (not using binning.")
        else:
            self._log_message("# Sampling with the PDS binned with %s bins." % bins)
        self._log_message("# Sampler started (for warmup).")
        self._log_message("# Sampling %s iterations as warmup." % nwarmup)
        self.__log_MCMC(header_only=True)
        # Warmup
        #while iters < nwarmup:
        #    pos = self.MCMC_sampler.run_mcmc(pos, self.MCMCp["step_size"])[0]
        #    self.MCMCp["chains"] = self.__get_chains(self.MCMC_sampler)
        #    self.__log_MCMC()
        #    iters = self.MCMC_sampler.chain.shape[1]
        #    self.MCMCp["niters"]  = iters

        for sample in self.MCMC_sampler.sample(pos, iterations=nwarmup, progress=True):
            # Check for convergence every "step_size" steps
            if self.MCMC_sampler.iteration % self.MCMCp["step_size"]:
                continue

            self.MCMCp["chains"] = self.__get_chains(self.MCMC_sampler)
            self.__log_MCMC()
            iters = self.MCMC_sampler.chain.shape[1]
            self.MCMCp["niters"]  = iters

        #sys.exit()
        self.MCMCp["nwarmup"] = self.MCMCp["niters"]
        # Minimum iterations after warmup
        self._log_message("# Sampling for a minimum of %s iterations." % minsteps)

        for sample in self.MCMC_sampler.sample(sample, iterations=minsteps, progress=True):
            # Check for convergence every "step_size" steps
            if self.MCMC_sampler.iteration % self.MCMCp["step_size"]:
                continue   
        
            iters = self.MCMC_sampler.chain.shape[1]
            self.MCMCp["niters"]  = iters
            self.MCMCp["chains"] = self.__get_chains(self.MCMC_sampler).transpose(1, 0, 2)
            self.__log_MCMC()




        self._log_message("# Sampling until reaching convergence, sufficient mixing and sample size.")
        minESS = mESS.minESS(self.ndim, alpha, eps)

        for sample in self.MCMC_sampler.sample(sample, iterations=maxsteps, progress=True):
            # Check for convergence every "step_size" steps
            if self.MCMC_sampler.iteration % self.MCMCp["step_size"]:
                continue   
            
            #fig, axes = plt.subplots(self.ndim, figsize=(10, 7), sharex=True)
            #samples = self.MCMC_sampler.get_chain()
            #labels = ["m", "b", "log(f)"]
            #for i in range(self.ndim):
            #    ax = axes[i]
            #    ax.plot(samples[:, :, i], "k", alpha=0.3)
            #    ax.set_xlim(0, len(samples))
                #ax.set_ylabel(labels[i])
                #ax.yaxis.set_label_coords(-0.1, 0.5)

            #axes[-1].set_xlabel("step number");
            #plt.show()

            iters = self.MCMC_sampler.chain.shape[1]
            self.MCMCp["niters"]  = iters
            self.MCMCp["chains"] = self.__get_chains(self.MCMC_sampler,
                                                     nwarmup=self.MCMCp["nwarmup"])

            # Test for convergence in the chains by looking at a flat log-probability   
            self.MCMCp["slope"] = self.__test_convergence(self.MCMCp["chains"], self.ndim)
            if self.MCMCp["slope"] < self.MCMCp["converge_test"]:
                self.MCMCp["converged_p"] = True
            else:
                self.MCMCp["converged_p"] = False
            # Test for sufficient chain mixing
            self.MCMCp["Rhat"] = Rhat(self.MCMCp["chains"][:, :, self.ndim])
            if self.MCMCp["Rhat"] < self.MCMCp["Rhat_max"]:
                self.MCMCp["mixed_p"] = True
            else:
                self.MCMCp["mixed_p"] = False
            # Test for sufficient multivariate effective sample size when we have enough samples
            if not (self.MCMCp["converged_p"] and self.MCMCp["mixed_p"]):
                self.MCMCp["nwarmup"] += self.MCMCp["step_size"]
                self.MCMCp["mESS"] = 0
            else:
                # Transpose chains for use with multiESS as written for emcee2 not emcee3
                self.MCMCp["mESS"] = np.sum([mESS.multiESS(self.MCMCp["chains"].transpose(1, 0, 2)[w, :, 0:self.ndim])
                                             for w in range(self.MCMCp["chains"].transpose(1, 0, 2).shape[0])])
            self.MCMCp["nsamples_p"] = (minESS <= self.MCMCp["mESS"])
            #print(np.shape(self.MCMCp["chains"]))
            #print(self.MCMCp["converged_p"], self.MCMCp["mixed_p"], self.MCMCp["nsamples_p"])

            #print(minESS, self.MCMCp["mESS"])
            self.__log_MCMC()

            if self.MCMCp["converged_p"] and self.MCMCp["mixed_p"] and self.MCMCp["nsamples_p"]:
                break

        if (self.MCMCp["converged_p"] and self.MCMCp["mixed_p"] and self.MCMCp["nsamples_p"]):
            self._log_message("# Updating background parameter estimations.")
            self.bg_params = self.theta_to_dict(np.median(self.MCMC_sampler.flatchain, axis=0))
        else:
            self._log_message("# NOT updating the background parameter estimations.")
        return



        # Continue until convergence, sufficient mixing and enough sample size.
        #self._log_message("# Sampling until reaching convergence, sufficient mixing and sample size.")
        #minESS = mESS.minESS(self.ndim, alpha, eps)
        #while (((not (self.MCMCp["converged_p"] and self.MCMCp["mixed_p"]))
        #       or not self.MCMCp["nsamples_p"])
        #       and (iters < maxsteps)):
        #    # Advance the MCMC sampler
        #    pos = self.MCMC_sampler.run_mcmc(pos, self.MCMCp["step_size"])[0]
        #    iters = self.MCMC_sampler.chain.shape[1]
        #    self.MCMCp["niters"]  = iters
        #    self.MCMCp["chains"] = self.__get_chains(self.MCMC_sampler,
        #                                             nwarmup=self.MCMCp["nwarmup"])

            # Test for convergence in the chains by looking at a flat log-probability
        #    self.MCMCp["slope"] = self.__test_convergence(self.MCMCp["chains"], self.ndim)
        #    if self.MCMCp["slope"] < self.MCMCp["converge_test"]:
        #        self.MCMCp["converged_p"] = True
        #    else:
        #        self.MCMCp["converged_p"] = False
        #    # Test for sufficient chain mixing
        #    self.MCMCp["Rhat"] = Rhat(self.MCMCp["chains"][:, :, self.ndim])
        #    if self.MCMCp["Rhat"] < self.MCMCp["Rhat_max"]:
        #        self.MCMCp["mixed_p"] = True
        #    else:
        #        self.MCMCp["mixed_p"] = False
        #    # Test for sufficient multivariate effective sample size when we have enough samples
        #    if not (self.MCMCp["converged_p"] and self.MCMCp["mixed_p"]):
        #        self.MCMCp["nwarmup"] += self.MCMCp["step_size"]
        #        self.MCMCp["mESS"] = 0
        #    else:
        #        self.MCMCp["mESS"] = np.sum([mESS.multiESS(self.MCMCp["chains"][w, :, 0:self.ndim])
        #                                     for w in range(self.MCMCp["chains"].shape[0])])
        #    self.MCMCp["nsamples_p"] = (minESS <= self.MCMCp["mESS"])
        #    self.__log_MCMC()
        #if (self.MCMCp["converged_p"] and self.MCMCp["mixed_p"] and self.MCMCp["nsamples_p"]):
        #    self._log_message("# Updating background parameter estimations.")
        #    self.bg_params = self.theta_to_dict(np.median(self.MCMC_sampler.flatchain, axis=0))
        #else:
        #    self._log_message("# NOT updating the background parameter estimations.")
        #return


    def __log_MCMC(self, header_only = False):
        if header_only:
            self._log_message(" niter,niterC, slope,nchains,  Rhat,  mESS")
            return
        self._log_message("{:06d},{:06d},{:6.5f},{:07d},{:6.5f},{:06d}".format(
            self.MCMCp["niters"],                         # iterations
            self.MCMCp["niters"] - self.MCMCp["nwarmup"], # iterations after warmup
            self.MCMCp["slope"],                          # convergence test
            self.MCMCp["chains"].shape[0],                # number of chains
            self.MCMCp["Rhat"],                           # Rhat
            int(self.MCMCp["mESS"])                       # mESS
        ))


    def write_posterior(self, filename, nwarmup=0):
        posterior = self.__get_posterior(self.MCMC_sampler, nwarmup)
        with open(filename, "w") as outcsv:
            csvwriter = csv.writer(outcsv)
            csvwriter.writerow(self.bg_param_names() + ["lnprob"])
            csvwriter.writerows(posterior)
        return


    def write_chains(self, filename, nwarmup=0):
        chains = self.__get_chains(self.MCMC_sampler, nwarmup)
        nwalkers = chains.shape[0]
        with open(filename, "w") as outcsv:
            csvwriter = csv.writer(outcsv)
            header = ["chain.%s.%s" % (walker, par)
                      for walker in range(nwalkers)
                      for par in (self.bg_param_names() + ["lnprob"])]
            csvwriter.writerow(header)
            for n in range(chains.shape[1]):
                csvwriter.writerow([chains[walker, n, par]
                                    for walker in range(nwalkers)
                                    for par in range(self.ndim + 1)])
        return
