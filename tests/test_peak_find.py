import pandas as pd
import pytest
from taco import peak_find

def test_peak_find():
    pds = pd.DataFrame({"frequency": [1, 2], "power": [1, 0]})
    oversampled_pds = pd.DataFrame({"frequency": [0.5, 1, 1.5, 2], "power": [0, 1, 0, 0]})
    data = pd.DataFrame({"numax": [1], "sigmaEnv": [1]})
    peaks = peak_find(pds, oversampled_pds, data)
    assert "frequency" in peaks
    assert "AIC" in peaks

def test_peak_find_with_peaks():
    pds = pd.DataFrame({"frequency": [1, 2], "power": [1, 0]})
    oversampled_pds = pd.DataFrame({"frequency": [0.5, 1, 1.5, 2], "power": [0, 1, 0, 0]})
    data = pd.DataFrame({"numax": [1], "sigmaEnv": [1]})
    peaks = pd.DataFrame({"frequency": [1], "AIC": [1]})
    mixed_peaks = pd.DataFrame({"frequency": [1]})
    peaks = peak_find(pds, oversampled_pds, data, peaks, mixed_peaks)
    assert "frequency" in peaks

def test_peak_find_001296068():
    """Unit test real data with ID 001296068"""

    pds = pd.read_csv("tests/data/test_peak_find_pds.csv")
    ofac_pds = pd.read_csv("tests/data/test_peak_find_ofac_pds.csv")
    data = pd.read_csv("tests/data/test_peak_find_summary.csv")
    peaks = peak_find(pds, ofac_pds, data, snr = 1.1, prob = 0.0001, minAIC = 2)
    assert "frequency" in peaks
