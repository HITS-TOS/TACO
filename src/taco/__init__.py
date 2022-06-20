# __init__.py
from .filter import filter
from .pds import calc_pds

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
]
