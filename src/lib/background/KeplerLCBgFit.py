## Background fit class for Kepler targets in the Long-Cadence mode

import numpy as np
from .PDSBgFit import PDSBgFit, classproperty


def _sLor(nu, A, b, c):
    return A / (1 + (nu / b) ** c)


def _sinc(x):
    return np.sinc(x/np.pi)


class _KeplerLCBgFit(PDSBgFit):
    """
    Abstraction for the whole MCMC background fit routine.
    This is intended to work (for now) for red-giants observed by Kepler
    in the long cadence (LC) mode.
    """

    @classproperty
    def par_rels(cls):
        """
        Empirical exponential relations between background parameters and numax.
        parameter = k*numax^s. First column is 'k', second is 's'
        """
        return {"Pn"       : [1.000000e-00,  0.000000], # Inferred from data, left alone
                "H1"       : [1.000000e-00,  0.000000], # Inferred from data left alone
                "b1"       : [0.5787,  0], # For filter of length 40 days
                "H2"       : [3382, -0.609], # Kallinger (2014)
                "b2"       : [0.317,  0.970], # Kallinger (2014)
                "H3"       : [3382, -0.609], # Kallinger (2014)
                "b3"       : [0.948,  0.992], # Kallinger (2014)
                "Pg"       : [2.03e7, -2.38], # Mosser (2012)
                "numax"    : [1.000000e-00,  1.000000],
                "sigmaEnv" : [0.28,  0.88]} # Mosser (2012) converted from denv to sigma


    def __init__(self, pds, initial_numax, initial_numax_sd, nuNyq, logfile=None):
        """
        We require the PDS and an initial guess at numax: 'initial_numax'
        """
        PDSBgFit.__init__(self, pds, logfile)
        self.initial_numax_sd = initial_numax_sd
        self.initial_numax = initial_numax
        self.nuNyq = nuNyq
        self.bg_params = self.guesses_from_numax(initial_numax)
        self._log_message("# Initial background values set by guesses using numax=%s." % initial_numax)


    @classmethod
    def guess_from_numax(cls, param, numax):
        return cls.par_rels[param][0] * (numax ** cls.par_rels[param][1])


    def guesses_from_numax(self, numax):
        raise NotImplementedError


class KeplerBg3Comp(_KeplerLCBgFit):

    @classmethod
    def bg_param_names(cls):
        return ["Pn", "H1", "b1", "H2", "b2", "H3", "b3", "Pg", "numax", "sigmaEnv"]


    def bgModel(self, theta, nu, no_osc=False):
        """
        Background model value at a given frequency 'nu'
        """
        Pn, H1, b1, H2, b2, H3, b3, Pg, numax, sigmaEnv = theta
        sc = _sinc(np.pi * nu / (2 * self.nuNyq)) ** 2
        bg = Pn
        bg = bg + sc * _sLor(nu, H1, b1, 4)
        bg = bg + sc * _sLor(nu, H2, b2, 4)
        bg = bg + sc * _sLor(nu, H3, b3, 4)
        if not no_osc:
            bg = bg + sc * Pg * np.exp(-((nu - numax) ** 2) / (2 * sigmaEnv ** 2))
        return bg


    def logPrio(self, theta):
        Pn, H1, b1, H2, b2, H3, b3, Pg, numax, sigmaEnv = theta

        # Prior 1
        if not ((Pn > 0) and (Pn < self._maxPDS)):
            #print(1)
            return -np.inf
        # Prior 2
        if not (H1 < 1.2*self._maxPDS):
        #    #print(2)
            return -np.inf
        # Prior 3
        if not (H1 > 0):
            #print("3")
            return -np.inf
        # Prior 4
        if not (H2 > H3):
            #print(4)
            return -np.inf
        # Prior 5
        if not (H3 > 0):
            #print(5)
            return -np.inf
        # Prior 6
        if not (b1 < b2 < b3):
            #print(6)
            return -np.inf
        # Prior 7
        if not ((b1 > 0) and (((numax * 0.9) < b3 < (numax * 1.1)) or ((numax * 0.9) < b2 < (numax * 1.1)))):
            #print(7)
            return -np.inf
        # Prior 8
        if not ((self.initial_numax - 1.1*self.initial_numax_sd) < numax < (self.initial_numax + 1.1*self.initial_numax_sd)):
            #print(8)
            return -np.inf
        # Prior 9
        if not (numax < 1.1*self.nuNyq):
            #print(9)
            return -np.inf
        # Prior 10
        if not ((Pg > 0) and (Pg < 1.2*self._maxPDS)):
            #print(10)
            return -np.inf
        #Pg0 = self.guess_from_numax("Pg", numax)
        #if not ((Pg0 * 0.1) < Pg < (Pg0 * 6)):
        #    #print(8)
        #    return -np.inf
        # Prior 11
        sigmaEnv0 = self.guess_from_numax("sigmaEnv", numax)
        if not ((sigmaEnv0 * 0.3) < sigmaEnv < (sigmaEnv0 * 3.0)):
            #print(11)
            return -np.inf
        return 0.0


    def guesses_from_numax(self, numax):
        """
        Guess the background parameters using only numax. The relations were
        empirically obtained beforehand for Kepler red-giants.
        """
        bg = self.theta_to_dict([self.guess_from_numax(par, numax)
                                 for par in self.bg_param_names()])
        ## Check they are consistent with the priors we define above
        # Convert to correct first guesses
        # First guesses on H2 and H3 from Kallinger are granulation 
        # amplitude, need to convert to height
        bg["H2"] = (2*np.sqrt(2))/np.pi * (bg["H2"]**2/bg["b2"])
        bg["H3"] = (2*np.sqrt(2))/np.pi * (bg["H3"]**2/bg["b3"])

        # Infer H1 from data
        bg["H1"] = 1.2 * bg["H2"]#0.9 * np.max(self.pds["power"])

        # Infer Pn from data - compute using median and corrective factor
        bg["Pn"] = 0.5 * np.mean(self.pds["power"][-10:])
        # Prior 1
        if bg["Pn"] <= 0: bg["Pn"] = 1e-4
        if bg["Pn"] >= self._maxPDS: bg["Pn"] = self._maxPDS - 1e-4
        # Prior 2
        if bg["H1"] >= self._maxPDS: bg["H1"] = self._maxPDS - 1e-4
        # Prior 3
        if bg["H1"] <= bg["H2"]: bg["H2"] = bg["H1"] - 1e-4
        if bg["H2"] <= bg["H3"]: bg["H3"] = bg["H2"] - 1e-4
        # Prior 4
        if bg["H3"] <= 0: bg["H3"] = 1e-4
        # Prior 5
        if bg["b1"] >= bg["b2"]: bg["b2"] = bg["b1"] + 1e-4
        if bg["b2"] >= bg["b3"]: bg["b3"] = bg["b2"] + 1e-4
        # Prior 6
        if bg["b1"] <= 0: bg["b1"] = np.min(self.pds["frequency"])
        if bg["b3"] >= 1.1*self.nuNyq: bg["b3"] = 1.1*self.nuNyq - 1e-4
        # Checking for priors 7-9 here would be redundant
        return bg


class KeplerBg2Comp(_KeplerLCBgFit):

    @classproperty
    def par_rels(cls):
        """
        Empirical exponential relations between background parameters and numax.
        parameter = k*numax^s. First column is 'k', second is 's'
        """
        return {"Pn"       : [1.000000e-00,  0.000000], # Inferred from data, left alone
                "H1"       : [3382, -0.609], # Kallinger (2014)
                "b1"       : [0.317,  0.970], # Kallinger (2014)
                "H2"       : [3382, -0.609], # Kallinger (2014)
                "b2"       : [0.948,  0.992], # Kallinger (2014)
                "Pg"       : [2.03e7, -2.38], # Mosser (2012)
                "numax"    : [1.000000e-00,  1.000000],
                "sigmaEnv" : [0.28,  0.88]} # Mosser (2012) converted from denv to sigma

    @classmethod
    def guess_from_numax(cls, param, numax):
        return cls.par_rels[param][0] * (numax ** cls.par_rels[param][1])

    @classmethod
    def bg_param_names(cls):
        return ["Pn", "H1", "b1", "H2", "b2", "Pg", "numax", "sigmaEnv"]


    def bgModel(self, theta, nu, no_osc=False):
        """
        Background model value at a given frequency 'nu'
        """
        Pn, H1, b1, H2, b2, Pg, numax, sigmaEnv = theta
        sc = _sinc(np.pi * nu / (2 * self.nuNyq)) ** 2
        bg = Pn
        bg = bg + sc * _sLor(nu, H1, b1, 4)
        bg = bg + sc * _sLor(nu, H2, b2, 4)
        if not no_osc:
            bg = bg + sc * Pg * np.exp(-((nu - numax) ** 2) / (2 * sigmaEnv ** 2))
        return bg


    def logPrio(self, theta):
        Pn, H1, b1, H2, b2, Pg, numax, sigmaEnv = theta
        # Prior 1
        if not ((Pn > 0) and (Pn < self._maxPDS)):
            return -np.inf
        # Prior 2
        if not (H1 < self._maxPDS):
            return -np.inf
        # Prior 3
        if not (H1 > H2):
            return -np.inf
        # Prior 4
        if not (H2 > 0):
            return -np.inf
        # Prior 5
        if not (b1 < b2):
            return -np.inf
        # Prior 6
        if not (b1 > 0 and b2 < 1.1 * self.nuNyq):
            return -np.inf
        # Prior 7
        if not ((self.initial_numax * 0.7) < numax < (self.initial_numax * 1.3)):
            return -np.inf
        # Prior 8
        Pg0 = self.guess_from_numax("Pg", self.initial_numax)
        if not ((Pg0 * 0.1) < Pg < (Pg0 * 6)):
            return -np.inf
        # Prior 9
        sigmaEnv0 = self.guess_from_numax("sigmaEnv", self.initial_numax)
        if not ((sigmaEnv0 * 0.3) < sigmaEnv < (sigmaEnv0 * 3.0)):
            return -np.inf
        return 0.0


    def guesses_from_numax(self, numax):
        """
        Guess the background parameters using only numax. The relations were
        empirically obtained beforehand for Kepler red-giants.
        """
        bg = self.theta_to_dict([self.guess_from_numax(par, numax)
                                 for par in self.bg_param_names()])
        ## Check they are consistent with the priors we define above
        # Prior 1
        if bg["Pn"] <= 0: bg["Pn"] = 1e-4
        if bg["Pn"] >= self._maxPDS: bg["Pn"] = self._maxPDS - 1e-4
        # Prior 2
        if bg["H1"] >= self._maxPDS: bg["H1"] = self._maxPDS - 1e-4
        # Prior 3
        if bg["H1"] <= bg["H2"]: bg["H2"] = bg["H1"] - 1e-4
        # Prior 4
        if bg["H2"] <= bg["Pn"]: bg["Pn"] = bg["H2"] - 1e-4
        # Prior 5
        if bg["b1"] >= bg["b2"]: bg["b2"] = bg["b1"] + 1e-4
        # Prior 6
        if bg["b1"] <= 0: bg["b1"] = 1e-4
        if bg["b2"] >= 2*self.nuNyq: bg["b2"] = 2*self.nuNyq - 1e-4
        # Checking for priors 7-9 here would be redundant
        return bg


class KeplerBg3CompExpVar(_KeplerLCBgFit):

    @classmethod
    def bg_param_names(cls):
        return ["Pn", "H1", "b1", "H2", "b2", "c2", "H3", "b3", "c3", "Pg", "numax", "sigmaEnv"]

    def par_rels(self):
        """
        Empirical exponential relations between background parameters and numax.
        parameter = k*numax^s. First column is 'k', second is 's'
        """
        return {"Pn"       : [1.000000e-00,  0.000000], # Inferred from data, left alone
                "H1"       : [1.000000e-00,  0.000000], # Inferred from data left alone
                "b1"       : [0.5787,  0], # For filter of length 40 days
                "H2"       : [3382, -0.609], # Kallinger (2014)
                "b2"       : [0.317,  0.970], # Kallinger (2014)
                "c2"       : [4, 0.0],
                "H3"       : [3382, -0.609], # Kallinger (2014)
                "b3"       : [0.948,  0.992], # Kallinger (2014)
                "c3"       : [4, 0.0],
                "Pg"       : [2.03e7, -2.38], # Mosser (2012)
                "numax"    : [1.000000e-00,  1.000000],
                "sigmaEnv" : [0.28,  0.88]} # Mosser (2012) converted from denv to sigma

    def guess_from_numax(self, param, numax):
        return self.par_rels()[param][0] * (numax ** self.par_rels()[param][1])

    def bgModel(self, theta, nu, no_osc=False):
        """
        Background model value at a given frequency 'nu'
        """
        Pn, H1, b1, H2, b2, c2, H3, b3, c3, Pg, numax, sigmaEnv = theta
        sc = _sinc(np.pi * nu / (2 * self.nuNyq)) ** 2
        bg = Pn
        bg = bg + sc * _sLor(nu, H1, b1, 4)
        bg = bg + sc * _sLor(nu, H2, b2, c2)
        bg = bg + sc * _sLor(nu, H3, b3, c3)
        if not no_osc:
            bg = bg + sc * Pg * np.exp(-((nu - numax) ** 2) / (2 * sigmaEnv ** 2))
        return bg


    def logPrio(self, theta):
        Pn, H1, b1, H2, b2, c2, H3, b3, c3, Pg, numax, sigmaEnv = theta
        # Prior 1
        if not ((c2 >= 2) and (c2 <= 6)):
            return -np.inf
        if not ((c3 >= 2) and (c3 <= 6)):
            return -np.inf
        if not ((Pn > 0) and (Pn < self._maxPDS)):
            #print(1)
            return -np.inf
        # Prior 2
        if not (H1 < 1.2*self._maxPDS):
        #    #print(2)
            return -np.inf
        # Prior 3
        if not (H1 > 0):
            #print("2a")
            return -np.inf
        if not (H2 > H3):
            #print(3)
            return -np.inf
        # Prior 4
        if not (H3 > 0):
            #print(4)
            return -np.inf
        # Prior 5
        if not (b1 < b2 < b3):
            #print(5)
            return -np.inf
        # Prior 6
        if not ((b1 > 0) and (b3 < 1.1*self.nuNyq)):
            #print(6)
            return -np.inf
        # Prior 7
        if not ((self.initial_numax * 0.7) < numax < (self.initial_numax * 1.3)):
            #print(7)
            return -np.inf
        # Prior 8
        Pg0 = self.guess_from_numax("Pg", self.initial_numax)
        if not ((Pg0 * 0.1) < Pg < (Pg0 * 6)):
            #print(8)
            return -np.inf
        # Prior 9
        sigmaEnv0 = self.guess_from_numax("sigmaEnv", self.initial_numax)
        if not ((sigmaEnv0 * 0.3) < sigmaEnv < (sigmaEnv0 * 3.0)):
            #print(9)
            return -np.inf
        return 0.0


    def guesses_from_numax(self, numax):
        """
        Guess the background parameters using only numax. The relations were
        empirically obtained beforehand for Kepler red-giants.
        """
        bg = self.theta_to_dict([self.guess_from_numax(par, numax)
                                 for par in self.bg_param_names()])
        ## Check they are consistent with the priors we define above
        # Convert to correct first guesses
        # First guesses on H2 and H3 from Kallinger are granulation 
        # amplitude, need to convert to height
        bg["H2"] = (2*np.sqrt(2))/np.pi * (bg["H2"]**2/bg["b2"])
        bg["H3"] = (2*np.sqrt(2))/np.pi * (bg["H3"]**2/bg["b3"])

        # Infer H1 from data
        bg["H1"] = 1.2 * bg["H2"]#0.9 * np.max(self.pds["power"])

        # Infer Pn from data - compute using median and corrective factor
        bg["Pn"] = 0.5 * np.mean(self.pds["power"][-10:])
        # Prior 1
        if bg["Pn"] <= 0: bg["Pn"] = 1e-4
        if bg["Pn"] >= self._maxPDS: bg["Pn"] = self._maxPDS - 1e-4
        # Prior 2
        if bg["H1"] >= self._maxPDS: bg["H1"] = self._maxPDS - 1e-4
        # Prior 3
        if bg["H1"] <= bg["H2"]: bg["H2"] = bg["H1"] - 1e-4
        if bg["H2"] <= bg["H3"]: bg["H3"] = bg["H2"] - 1e-4
        # Prior 4
        if bg["H3"] <= 0: bg["H3"] = 1e-4
        # Prior 5
        if bg["b1"] >= bg["b2"]: bg["b2"] = bg["b1"] + 1e-4
        if bg["b2"] >= bg["b3"]: bg["b3"] = bg["b2"] + 1e-4
        # Prior 6
        if bg["b1"] <= 0: bg["b1"] = np.min(self.pds["frequency"])
        if bg["b3"] >= 1.1*self.nuNyq: bg["b3"] = 1.1*self.nuNyq - 1e-4
        # Checking for priors 7-9 here would be redundant
        return bg
