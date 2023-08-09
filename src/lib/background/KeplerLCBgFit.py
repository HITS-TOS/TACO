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
                "A1"       : [1.000000e-00,  0.000000], # Inferred from data left alone
                "b1"       : [0.5787,  0], # For filter of length 40 days
                "A2"       : [3382, -0.609], # Kallinger (2014)
                "b2"       : [0.317,  0.970], # Kallinger (2014)
                "A3"       : [3382, -0.609], # Kallinger (2014)
                "b3"       : [0.948,  0.992], # Kallinger (2014)
                "Pg"       : [2.03e7, -2.38], # Mosser (2012)
                "numax"    : [1.000000e-00,  1.000000],
                "sigmaEnv" : [0.28,  0.88]} # Mosser (2012) converted from denv to sigma


    def __init__(self, pds, numax0, numax0_sd, nuNyq, logfile=None):
        """
        We require the PDS and an initial guess at numax: 'numax0'
        """
        PDSBgFit.__init__(self, pds, logfile)
        self.numax0_sd = numax0_sd
        self.numax0 = numax0
        self.nuNyq = nuNyq
        self.bg_params = self.guesses_from_numax(numax0)
        self._log_message("# Initial background values set by guesses using numax=%s." % numax0)


    @classmethod
    def guess_from_numax(cls, param, numax):
        return cls.par_rels[param][0] * (numax ** cls.par_rels[param][1])


    def guesses_from_numax(self, numax):
        raise NotImplementedError



class KeplerBg3Comp(_KeplerLCBgFit):

    @classmethod
    def bg_param_names(cls):
        return ["Pn", "A1", "b1", "A2", "b2", "A3", "b3", "Pg", "numax", "sigmaEnv"]


    def bgModel(self, theta, nu, no_osc=False):
        """
        Background model value at a given frequency 'nu'
        """
        Pn, A1, b1, A2, b2, A3, b3, Pg, numax, sigmaEnv = theta
        sc = _sinc(np.pi * nu / (2 * self.nuNyq)) ** 2
        bg = Pn
        bg = bg + sc * _sLor(nu, A1, b1, 4)
        bg = bg + sc * _sLor(nu, A2, b2, 4)
        bg = bg + sc * _sLor(nu, A3, b3, 4)
        if not no_osc:
            bg = bg + sc * Pg * np.exp(-((nu - numax) ** 2) / (2 * sigmaEnv ** 2))
        return bg


    def logPrio(self, theta):
        Pn, A1, b1, A2, b2, A3, b3, Pg, numax, sigmaEnv = theta
        # Prior 1
        
        if not ((Pn > 0) and (Pn < self._maxPDS)):
            #print(1)
            return -np.inf
        # Prior 2
        if not (A1 < 1.2*self._maxPDS):
        #    #print(2)
            return -np.inf
        # Prior 3
        if not (A1 > 0):
            #print("2a")
            return -np.inf
        if not (A2 > A3):
            #print(3)
            return -np.inf
        # Prior 4
        if not (A3 > 0):
            #print(4)
            return -np.inf
        # Prior 5
        if not (b1 < b2 < b3):
            #print(5)
            return -np.inf
        # Prior 6
        if not ((b1 > 0) and (((numax * 0.9) < b3 < (numax * 1.1)) or ((numax * 0.9) < b2 < (numax * 1.1)))):
            #print(6)
            return -np.inf
        # Prior 7
        if not ((self.numax0 - 1.1*self.numax0_sd) < numax < (self.numax0 + 1.1*self.numax0_sd)):
            #print(7)
            return -np.inf
        # Prior 8
        if not ((Pg > 0) and (Pg < 1.2*self._maxPDS)):
            #print(8)
            return -np.inf
        #Pg0 = self.guess_from_numax("Pg", numax)
        #if not ((Pg0 * 0.1) < Pg < (Pg0 * 6)):
        #    #print(8)
        #    return -np.inf
        # Prior 9
        sigmaEnv0 = self.guess_from_numax("sigmaEnv", numax)
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
        # First guesses on A2 and A3 from Kallinger are granulation 
        # amplitude, need to convert to height
        bg["A2"] = (2*np.sqrt(2))/np.pi * (bg["A2"]**2/bg["b2"])
        bg["A3"] = (2*np.sqrt(2))/np.pi * (bg["A3"]**2/bg["b3"])

        # Infer A1 from data
        bg["A1"] = 1.2 * bg["A2"]#0.9 * np.max(self.pds["power"])

        # Infer Pn from data - compute using median and corrective factor
        bg["Pn"] = 0.5 * np.mean(self.pds["power"][-10:])
        # Prior 1
        if bg["Pn"] <= 0: bg["Pn"] = 1e-4
        if bg["Pn"] >= self._maxPDS: bg["Pn"] = self._maxPDS - 1e-4
        # Prior 2
        if bg["A1"] >= self._maxPDS: bg["A1"] = self._maxPDS - 1e-4
        # Prior 3
        if bg["A1"] <= bg["A2"]: bg["A2"] = bg["A1"] - 1e-4
        if bg["A2"] <= bg["A3"]: bg["A3"] = bg["A2"] - 1e-4
        # Prior 4
        if bg["A3"] <= 0: bg["A3"] = 1e-4
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
                "A1"       : [3382, -0.609], # Kallinger (2014)
                "b1"       : [0.317,  0.970], # Kallinger (2014)
                "A2"       : [3382, -0.609], # Kallinger (2014)
                "b2"       : [0.948,  0.992], # Kallinger (2014)
                "Pg"       : [2.03e7, -2.38], # Mosser (2012)
                "numax"    : [1.000000e-00,  1.000000],
                "sigmaEnv" : [0.28,  0.88]} # Mosser (2012) converted from denv to sigma

    @classmethod
    def guess_from_numax(cls, param, numax):
        return cls.par_rels[param][0] * (numax ** cls.par_rels[param][1])

    @classmethod
    def bg_param_names(cls):
        return ["Pn", "A1", "b1", "A2", "b2", "Pg", "numax", "sigmaEnv"]


    def bgModel(self, theta, nu, no_osc=False):
        """
        Background model value at a given frequency 'nu'
        """
        Pn, A1, b1, A2, b2, Pg, numax, sigmaEnv = theta
        sc = _sinc(np.pi * nu / (2 * self.nuNyq)) ** 2
        bg = Pn
        bg = bg + sc * _sLor(nu, A1, b1, 4)
        bg = bg + sc * _sLor(nu, A2, b2, 4)
        if not no_osc:
            bg = bg + sc * Pg * np.exp(-((nu - numax) ** 2) / (2 * sigmaEnv ** 2))
        return bg


    def logPrio(self, theta):
        Pn, A1, b1, A2, b2, Pg, numax, sigmaEnv = theta
        # Prior 1
        if not ((Pn > 0) and (Pn < self._maxPDS)):
            return -np.inf
        # Prior 2
        if not (A1 < self._maxPDS):
            return -np.inf
        # Prior 3
        if not (A1 > A2):
            return -np.inf
        # Prior 4
        if not (A2 > 0):
            return -np.inf
        # Prior 5
        if not (b1 < b2):
            return -np.inf
        # Prior 6
        if not (b1 > 0 and b2 < 1.1 * self.nuNyq):
            return -np.inf
        # Prior 7
        if not ((self.numax0 * 0.7) < numax < (self.numax0 * 1.3)):
            return -np.inf
        # Prior 8
        Pg0 = self.guess_from_numax("Pg", self.numax0)
        if not ((Pg0 * 0.1) < Pg < (Pg0 * 6)):
            return -np.inf
        # Prior 9
        sigmaEnv0 = self.guess_from_numax("sigmaEnv", self.numax0)
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
        if bg["A1"] >= self._maxPDS: bg["A1"] = self._maxPDS - 1e-4
        # Prior 3
        if bg["A1"] <= bg["A2"]: bg["A2"] = bg["A1"] - 1e-4
        # Prior 4
        if bg["A2"] <= bg["Pn"]: bg["Pn"] = bg["A2"] - 1e-4
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
        return ["Pn", "A1", "b1", "A2", "b2", "c2", "A3", "b3", "c3", "Pg", "numax", "sigmaEnv"]

    def par_rels(self):
        """
        Empirical exponential relations between background parameters and numax.
        parameter = k*numax^s. First column is 'k', second is 's'
        """
        return {"Pn"       : [1.000000e-00,  0.000000], # Inferred from data, left alone
                "A1"       : [1.000000e-00,  0.000000], # Inferred from data left alone
                "b1"       : [0.5787,  0], # For filter of length 40 days
                "A2"       : [3382, -0.609], # Kallinger (2014)
                "b2"       : [0.317,  0.970], # Kallinger (2014)
                "c2"       : [4, 0.0],
                "A3"       : [3382, -0.609], # Kallinger (2014)
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
        Pn, A1, b1, A2, b2, c2, A3, b3, c3, Pg, numax, sigmaEnv = theta
        sc = _sinc(np.pi * nu / (2 * self.nuNyq)) ** 2
        bg = Pn
        bg = bg + sc * _sLor(nu, A1, b1, 4)
        bg = bg + sc * _sLor(nu, A2, b2, c2)
        bg = bg + sc * _sLor(nu, A3, b3, c3)
        if not no_osc:
            bg = bg + sc * Pg * np.exp(-((nu - numax) ** 2) / (2 * sigmaEnv ** 2))
        return bg


    def logPrio(self, theta):
        Pn, A1, b1, A2, b2, c2, A3, b3, c3, Pg, numax, sigmaEnv = theta
        # Prior 1
        if not ((c2 >= 2) and (c2 <= 6)):
            return -np.inf
        if not ((c3 >= 2) and (c3 <= 6)):
            return -np.inf
        if not ((Pn > 0) and (Pn < self._maxPDS)):
            #print(1)
            return -np.inf
        # Prior 2
        if not (A1 < 1.2*self._maxPDS):
        #    #print(2)
            return -np.inf
        # Prior 3
        if not (A1 > 0):
            #print("2a")
            return -np.inf
        if not (A2 > A3):
            #print(3)
            return -np.inf
        # Prior 4
        if not (A3 > 0):
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
        if not ((self.numax0 * 0.7) < numax < (self.numax0 * 1.3)):
            #print(7)
            return -np.inf
        # Prior 8
        Pg0 = self.guess_from_numax("Pg", self.numax0)
        if not ((Pg0 * 0.1) < Pg < (Pg0 * 6)):
            #print(8)
            return -np.inf
        # Prior 9
        sigmaEnv0 = self.guess_from_numax("sigmaEnv", self.numax0)
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
        # First guesses on A2 and A3 from Kallinger are granulation 
        # amplitude, need to convert to height
        bg["A2"] = (2*np.sqrt(2))/np.pi * (bg["A2"]**2/bg["b2"])
        bg["A3"] = (2*np.sqrt(2))/np.pi * (bg["A3"]**2/bg["b3"])

        # Infer A1 from data
        bg["A1"] = 1.2 * bg["A2"]#0.9 * np.max(self.pds["power"])

        # Infer Pn from data - compute using median and corrective factor
        bg["Pn"] = 0.5 * np.mean(self.pds["power"][-10:])
        # Prior 1
        if bg["Pn"] <= 0: bg["Pn"] = 1e-4
        if bg["Pn"] >= self._maxPDS: bg["Pn"] = self._maxPDS - 1e-4
        # Prior 2
        if bg["A1"] >= self._maxPDS: bg["A1"] = self._maxPDS - 1e-4
        # Prior 3
        if bg["A1"] <= bg["A2"]: bg["A2"] = bg["A1"] - 1e-4
        if bg["A2"] <= bg["A3"]: bg["A3"] = bg["A2"] - 1e-4
        # Prior 4
        if bg["A3"] <= 0: bg["A3"] = 1e-4
        # Prior 5
        if bg["b1"] >= bg["b2"]: bg["b2"] = bg["b1"] + 1e-4
        if bg["b2"] >= bg["b3"]: bg["b3"] = bg["b2"] + 1e-4
        # Prior 6
        if bg["b1"] <= 0: bg["b1"] = np.min(self.pds["frequency"])
        if bg["b3"] >= 1.1*self.nuNyq: bg["b3"] = 1.1*self.nuNyq - 1e-4
        # Checking for priors 7-9 here would be redundant
        return bg
