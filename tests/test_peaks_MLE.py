import pandas as pd
import pytest
from taco import peaks_MLE

def test_peaks_MLE():
    pds = pd.DataFrame({"frequency": [1, 2],
                        "power": [1, 0]})
    peaks = pd.DataFrame({"frequency": [1],
                          "amplitude": [1],
                          "linewidth": [0.1]})
    data = pd.DataFrame({"numax": [1],
                         "sigmaEnv": [1],
                         "DeltaNu": [1],
                         "eps_p": [1],
                         "gamma0": [1]})
    peaks = peaks_MLE(pds, peaks, data)
    assert "frequency" in peaks

def test_peaks_MLE_001296068():
    """Unit test real data with ID 001296068"""

    pds = pd.read_csv("tests/data/test_peaks_MLE_pds.csv")
    peaks = pd.read_csv("tests/data/test_peaks_MLE_peaks.csv")
    data = pd.read_csv("tests/data/test_peaks_MLE_summary.csv")
    peaks = peaks_MLE(pds, peaks, data)
    assert "frequency" in peaks
