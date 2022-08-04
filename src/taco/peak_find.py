import os
from pathlib import Path

import rpy2.robjects as ro
import rpy2.robjects.conversion as cv
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import STAP

def _none2null(none_obj):
    return ro.r("NULL")

def peak_find(pds, oversampled_pds, peaks = None, mixedpeaks = None, snr = 1.2, prob = 0.0001,
              maxlwd = None, removel02 = False, minAIC = 2, navg = 1):
    """
    Find the relevant solar-like oscillations in a background-removed PDS
    We use a tree-map of the local maxima found by a (mexican-hat) wavelet
    transform in the PDS and also look at the residuals to identify "unresolved"
    oscillations.

    Parameters:
        pds(pandas.DataFrame):Periodogram
            Columns:
                Name: frequency, dtype: float[micro-Hertz]
                Name: power, dtype: float
        oversampled_pds(pandas.DataFrame):Oversampled periodogram
            Columns:
                Name: frequency, dtype: float[micro-Hertz]
                Name: power, dtype: float
        peaks(pandas.DataFrame): Identified peaks. It must contain the l=0,2 modes already identified
            Columns:
                Name: frequency, dtype: float[micro-Hertz]
                Name: linewidth, dtype: float
                Name: height, dtype: float
                Name: snr, dtype: float
                Name: AIC, dtype: float
                Name: amplitude, dtype: float
        mixedpeaks(pandas.DataFrame): Mixed mode peaks from peak finding
            Columns:
                Name: frequency, dtype: float[micro-Hertz]
        snr(float): Minimum signal-to-noise ratio (on CWT space) for resolved peaks
        prob(float): Minimum (frequentist) probability threshold for unresolved peaks
        maxlwd(float): Maximum search linewidth for resolved peaks in the CWT search
        removel02(bool): Whether or not the l02 peaks should be divided out before running the CWT search
        minAIC(int): Minimum AIC value for a peak to have in order to be considered significant
        navg(int): Number of power spectra averaged to create current power spectrum
    """

    with open(Path(Path(__file__).parent, 'peak_find.R'), 'r') as f:
        owd = os.getcwd()
        os.chdir(Path(__file__).parents[2])
        peak_find = STAP(f.read(), "peak_find_r")
        os.chdir(owd)

        none_converter = cv.Converter("None converter")
        none_converter.py2rpy.register(type(None), _none2null)

        with localconverter(ro.default_converter + pandas2ri.converter + none_converter):
            r_pds = ro.conversion.py2rpy(pds)
            r_oversampled_pds = ro.conversion.py2rpy(oversampled_pds)
            return peak_find.peak_find_r(r_pds, r_oversampled_pds,
                 peaks = peaks,
                 mixedpeaks = mixedpeaks,
                 snr = snr,
                 prob = prob,
                 maxlwd = maxlwd,
                 removel02 = removel02,
                 minAIC = minAIC,
                 navg = navg)
