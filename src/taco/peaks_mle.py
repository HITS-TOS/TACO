import os
from pathlib import Path

import rpy2.robjects as ro
import rpy2.robjects.conversion as cv
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import STAP

def _none2null(none_obj):
    return ro.r("NULL")

def peaks_mle(pds, peaks, data, mixed_peaks = None, maxlwd = None,
              removel02 = False, minAIC = 2, navg = 1, finalfit = False):
    """
    MLE optimization for the peaks found by peakFind.R (or another way)
    It takes a csv file with columns named frequency, height and linewidth
    and makes an MLE optimization on a PDS. If linewidth is NA, the modes are
    assumed to be sinc^2 functions, otherwise they are assumed as Lorentzians.

    Parameters:
        pds(pandas.DataFrame):Periodogram
            Columns:
                Name: frequency, dtype: float[micro-Hertz]
                Name: power, dtype: float
        peaks(pandas.DataFrame): Identified peaks. It must contain the l=0,2 modes already identified
            Columns:
                Name: frequency, dtype: float[micro-Hertz]
                Name: linewidth, dtype: float
                Name: amplitude, dtype: float
        data(pandas.DataFrame):Summary data
            Columns:
                Name: numax, dtype: float
        mixed_peaks(pandas.DataFrame): Mixed mode peaks from peak finding
            Columns:
                Name: frequency, dtype: float[micro-Hertz]
        maxlwd(float): Maximum search linewidth for resolved peaks in the CWT search
        removel02(bool): Whether or not the l02 peaks should be divided out before running the CWT search
        minAIC(int): Minimum AIC value for a peak to have in order to be considered significant
        navg(int): Number of power spectra averaged to create current power spectrum
        finalfit(bool): Whether or not this is the final MLE optimisation
    """

    with open(Path(Path(__file__).parent, 'peaks_mle.R'), 'r') as f:
        owd = os.getcwd()
        os.chdir(Path(__file__).parents[2])
        peaks_mle = STAP(f.read(), "peaks_mle_r")
        os.chdir(owd)

        none_converter = cv.Converter("None converter")
        none_converter.py2rpy.register(type(None), _none2null)

        with localconverter(ro.default_converter + pandas2ri.converter + none_converter):
            r_pds = ro.conversion.py2rpy(pds)
            r_peaks = ro.conversion.py2rpy(peaks)
            r_data = ro.conversion.py2rpy(data)
            result = peaks_mle.peaks_mle_r(r_pds, r_peaks, r_data,
                 mixed_peaks = mixed_peaks,
                 maxlwd = maxlwd,
                 removel02 = removel02,
                 minAIC = minAIC,
                 navg = navg,
                 finalfit = finalfit)

            peaks = ro.conversion.rpy2py(result[0])
            flag = ro.conversion.rpy2py(result[1])
            data = ro.conversion.rpy2py(result[2])
            return peaks, int(flag[0]), data
