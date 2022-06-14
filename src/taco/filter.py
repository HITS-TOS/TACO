import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri

from rpy2.robjects.conversion import localconverter

def filter(ts, width=1, remove_gaps=-1):
    """
    Filtering the lightcurves with a triangular smooth
    I use two rectangular smooths with half the provided width.
    Additionally it interpolates the single-point missing values
    It also calculates the mean and variance of the light curve
    and saves them in the summary file.

    Parameters:
        ts(pandas.DataFrame):Time series with units of days.
            Columns:
                Name: time, dtype: datetime64[days]
                Name: flux, dtype: int64
        width(int):Width of the triangular filter
        remove_gaps(int):Remove gaps greater than this value (in days). If set to -1, do not remove gaps
    """

    base = importr('base')

    with localconverter(ro.default_converter + pandas2ri.converter):
        r_ts = ro.conversion.py2rpy(ts)
