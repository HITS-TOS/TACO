import pandas as pd
import pytest
from taco import peaks_mle

def test_peaks_mle():
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
    peaks, data = peaks_mle(pds, peaks, data)
    assert "frequency" in peaks

def test_peaks_mle_001296068():
    """Unit test real data with ID 001296068"""

    pds = pd.read_csv("tests/data/test_peaks_mle_pds_bgr.csv")
    peaks = pd.read_csv("tests/data/test_peaks_mle_peaks.csv")
    data = pd.read_csv("tests/data/test_peaks_mle_summary.csv")
    peaks, data = peaks_mle(pds, peaks, data)
    assert data["npeaks"][0] == 31
    assert peaks["frequency"][0] == pytest.approx(34.03745602727537, 1e-4)
