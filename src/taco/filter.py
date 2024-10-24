from pathlib import Path

import numpy as np 
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import STAP


def filter(ts, data = pd.DataFrame.from_dict({"KIC": [float("nan")], "raw_data": [float("nan")]}),
    width = 40, remove_gaps = -1, output = '', output_directory = ''):
    """
    Filtering the lightcurves with a triangular smooth
    I use two rectangular smooths with half the provided width.
    Additionally it interpolates the single-point missing values
    It also calculates the mean and variance of the light curve
    and saves them in the summary file.

    Parameters:
        ts(pandas.DataFrame): Time series with units of days.
            Columns:
                Name: time, dtype: float[days]
                Name: flux, dtype: float
        data(pandas.DataFrame): Summary data
        width(int): Width of the triangular filter
        remove_gaps(int): Remove gaps greater than this value (in days).
                          If set to -1, do not remove gaps
    """
    
    with open(Path(Path(__file__).parent, 'filter.R'), 'r') as f:
        filter = STAP(f.read(), "filter_r")
            
            
        with localconverter(ro.default_converter + pandas2ri.converter):
        #with (ro.default_converter + pandas2ri.converter).context():
            r_ts = ro.conversion.py2rpy(ts)
            r_data = ro.conversion.py2rpy(data)
            
            result = filter.filter_r(r_ts, r_data, width, remove_gaps)
            ts_filtered = ro.conversion.rpy2py(result['filtered'])
            data = ro.conversion.rpy2py(result['data'])
            
            if output:
                ts_filtered.to_csv(Path(output_directory, output), index = False)

            return ts_filtered, data
