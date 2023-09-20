import os
from pathlib import Path

import rpy2.robjects as ro
import rpy2.robjects.conversion as cv
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import STAP

def _none2null(none_obj):
    return ro.r("NULL")

def peak_bag_mode_id02(pds, peaks, data):
    """
    Get Δν and label the l=0,2 peaks

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
    """

    with open(Path(Path(__file__).parent, 'peak_bag_mode_id02.R'), 'r') as f:
        owd = os.getcwd()
        os.chdir(Path(__file__).parents[2])
        peak_bag_mode_id02 = STAP(f.read(), "peak_bag_mode_id02_r")
        os.chdir(owd)

        none_converter = cv.Converter("None converter")
        none_converter.py2rpy.register(type(None), _none2null)

        with localconverter(ro.default_converter + pandas2ri.converter + none_converter):
            r_pds = ro.conversion.py2rpy(pds)
            r_peaks = ro.conversion.py2rpy(peaks)
            r_data = ro.conversion.py2rpy(data)
            result = peak_bag_mode_id02.peak_bag_mode_id02_r(r_pds, r_peaks, r_data)

            peaks = ro.conversion.rpy2py(result[0])
            flag = ro.conversion.rpy2py(result[1])
            data = ro.conversion.rpy2py(result[2])
            return peaks, int(flag[0]), data
