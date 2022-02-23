## Background fit class for PLATO targets in the Long-Cadence mode

import numpy as np
from PDSBgFit import PDSBgFit, classproperty


def _sLor(nu, A, b, c):
    return A / (1 + (nu / b) ** c)


def _sinc(x):
    return np.sinc(x/np.pi)


class _PLATOLCBgFit(PDSBgFit):
    """
    Abstraction for the whole MCMC background fit routine.
    This is intended to work (for now) for red-giants observed by PLATO.
    """

    nuNyq = 4166.61          # TODO  Check this value

    @classproperty
    def par_rels(cls):
        """
        Empirical exponential relations between background parameters and numax.
        parameter = k*numax^s. First column is 'k', second is 's'
        """
        return {"Pn"       : [1.157062e+03, -0.934038],
                "A1"       : [1.546766e+07, -1.787033],
                "b1"       : [1.342063e+00,  0.392020],
                "A2"       : [9.231790e+06, -2.128955],
                "b2"       : [1.019166e+00,  0.873758],
                "A3"       : [8.893120e+05, -2.091983],
                "b3"       : [6.945823e-01,  1.274016],
                "Pg"       : [1.963711e+07, -2.260536],
                "numax"    : [1.000000e-00,  1.000000],
                "sigmaEnv" : [3.081358e-01,  0.805355]}


    def __init__(self, pds, numax0, logfile=None):
        """
        We require the PDS and an initial guess at numax: 'numax0'
        """
        PDSBgFit.__init__(self, pds, logfile)
        self.numax0 = numax0
        self.bg_params = self.guesses_from_numax(numax0)
        self._log_message("# Initial background values set by guesses using numax=%s." % numax0)


    @classmethod
    def guess_from_numax(cls, param, numax):
        return cls.par_rels[param][0] * (numax ** cls.par_rels[param][1])


    def guesses_from_numax(self, numax):
        raise NotImplementedError



class PLATOBg3Comp(_PLATOLCBgFit):

    @classmethod
    def bg_param_names(cls):
        return ["Pn", "A1", "b1", "A2", "b2", "A3", "b3", "Pg", "numax", "sigmaEnv"]


    def bgModel(self, theta, nu):
        """
        Background model value at a given frequency 'nu'
        """
        Pn, A1, b1, A2, b2, A3, b3, Pg, numax, sigmaEnv = theta
        sc = _sinc(np.pi * nu / (2 * self.nuNyq)) ** 2
        bg = Pn
        bg = bg + sc * _sLor(nu, A1, b1, 4)
        bg = bg + sc * _sLor(nu, A2, b2, 4)
        bg = bg + sc * _sLor(nu, A3, b3, 4)
        bg = bg + sc * Pg * np.exp(-((nu - numax) ** 2) / (2 * sigmaEnv ** 2))
        return bg


    def logPrio(self, theta):
        Pn, A1, b1, A2, b2, A3, b3, Pg, numax, sigmaEnv = theta
        # Prior 1
        if not ((Pn > 0) and (Pn < self._maxPDS)):
            print(1)
            return -np.inf
        # Prior 2
        if not (A1 < self._maxPDS):
            print(2)
            return -np.inf
        # Prior 3
        if not (A1 > A2 > A3):
            print(3)
            return -np.inf
        # Prior 4
        if not (A3 > Pn):
            print(4)
            return -np.inf
        # Prior 5
        if not (b1 < b2 < b3):
            print(5)
            return -np.inf
        # Prior 6
        if not ((b1 > 0) and (b3 < 2*self.nuNyq)):
            print(6)
            return -np.inf
        # Prior 7
        if not ((self.numax0 * 0.7) < numax < (self.numax0 * 1.3)):
            print(7)
            return -np.inf
        # Prior 8
        Pg0 = self.guess_from_numax("Pg", self.numax0)
        if not ((Pg0 * 0.1) < Pg < (Pg0 * 6)):
            print(8)
            return -np.inf
        # Prior 9
        sigmaEnv0 = self.guess_from_numax("sigmaEnv", self.numax0)
        if not ((sigmaEnv0 * 0.3) < sigmaEnv < (sigmaEnv0 * 3.0)):
            print(9)
            return -np.inf
        return 0.0


    def guesses_from_numax(self, numax):
        """
        Guess the background parameters using only numax. The relations were
        empirically obtained beforehand for PLATO red-giants.
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
        if bg["A2"] <= bg["A3"]: bg["A3"] = bg["A2"] - 1e-4
        # Prior 4
        if bg["A3"] <= bg["Pn"]: bg["Pn"] = bg["A3"] - 1e-4
        # Prior 5
        if bg["b1"] >= bg["b2"]: bg["b2"] = bg["b1"] + 1e-4
        if bg["b2"] >= bg["b3"]: bg["b3"] = bg["b2"] + 1e-4
        # Prior 6
        if bg["b1"] <= 0: bg["b1"] = 1e-4
        if bg["b3"] >= 2*self.nuNyq: bg["b3"] = 2*self.nuNyq - 1e-4
        # Checking for priors 7-9 here would be redundant
        return bg



class PLATOBg2Comp(_PLATOLCBgFit):

    @classmethod
    def bg_param_names(cls):
        return ["Pn", "A1", "b1", "A2", "b2", "Pg", "numax", "sigmaEnv"]


    def bgModel(self, theta, nu):
        """
        Background model value at a given frequency 'nu'
        """
        Pn, A1, b1, A2, b2, Pg, numax, sigmaEnv = theta
        sc = _sinc(np.pi * nu / (2 * self.nuNyq)) ** 2
        bg = Pn
        bg = bg + sc * _sLor(nu, A1, b1, 4)
        bg = bg + sc * _sLor(nu, A2, b2, 4)
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
        if not (A2 > Pn):
            return -np.inf
        # Prior 5
        if not (b1 < b2):
            return -np.inf
        # Prior 6
        if not (b1 > 0 and b2 < 2 * self.nuNyq):
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
        empirically obtained beforehand for PLATO red-giants.
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
