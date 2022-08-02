import os
from pathlib import Path

import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import STAP


def peak_find(pds, oversampled_pds):
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
    """

    with open(Path(Path(__file__).parent, 'peak_find.R'), 'r') as f:
        owd = os.getcwd()
        os.chdir(Path(__file__).parents[2])
        numax_estimate = STAP(f.read(), "peak_find_r")
        os.chdir(owd)

        with localconverter(ro.default_converter + pandas2ri.converter):
            r_pds = ro.conversion.py2rpy(pds)
            r_oversampled_pds = ro.conversion.py2rpy(oversampled_pds)
            return peak_find.peak_find_r(r_pds, r_oversampled_pds)
