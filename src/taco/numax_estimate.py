import os
from pathlib import Path

import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import STAP


def numax_estimate(pds, data, filter_width = 0.2):
    """
    Make a rough numax estimation using the periodogram, the variance of the time-series
    and the Nyquist frequency.

    Parameters:
        pds(pandas.DataFrame): Periodogram
            Columns:
                Name: frequency, dtype: float[micro-Hertz]
                Name: power, dtype: float
        data(pandas.DataFrame): Summary data
            Columns:
                var(float): Variance of the time-series
                nuNyq(float): Nyquist frequency
        filterwidth(float): The width of the log-median filter used to remove the background
                            for the wavelet numax estimation
    """

    with open(Path(Path(__file__).parent, 'numax_estimate.R'), 'r') as f:
        owd = os.getcwd()
        os.chdir(Path(__file__).parents[2])
        numax_estimate = STAP(f.read(), "numax_estimate_r")
        os.chdir(owd)

        with localconverter(ro.default_converter + pandas2ri.converter):
            r_pds = ro.conversion.py2rpy(pds)
            r_data = ro.conversion.py2rpy(data)
            result = numax_estimate.numax_estimate_r(r_pds, r_data, filter_width)
            
            flag = ro.conversion.rpy2py(result[1])
            data = ro.conversion.rpy2py(result[0])
            return data, flag
