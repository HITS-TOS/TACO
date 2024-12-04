# __init__.py

from .filter import filter
from .pds import calc_pds
from .numax_estimate import numax_estimate
from .background_fit import background_fit
from .peak_find import peak_find
from .peaks_mle import peaks_mle
from .peak_bag_mode_id02 import peak_bag_mode_id02
from .peak_bag_mode_id3 import peak_bag_mode_id3
from .peak_bag_period_spacing import peak_bag_period_spacing
from .cv_method import cv_method

__doc__ = """
TACO - Tools for Automated Characterisation of Oscillations
===========================================================
TACO is an auto-mated modular code that takes timeseries data and analyses these to present the characteristics of
each detected stellar oscillation, i.e. frequency, amplitude, linewidth and mode identification. These
parameters provide essential information needed to determine the internal structure of these stars
through asteroseismology - the study of the internal structure of stars through their global oscillation
modes.
"""

__all__ = [
    "filter",
    "pds",
    "cv_method",
    "numax_estimate",
    "background_fit",
    "peak_find",
    "peaks_mle",
    "peak_bag_mode_id02",
    "peak_bag_period_spacing",
]
