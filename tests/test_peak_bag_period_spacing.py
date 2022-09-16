import pandas as pd
from taco import peak_bag_period_spacing

def test_peak_bag_period_spacing():
    pds = pd.DataFrame({"frequency": [1, 2],
                        "power": [1, 0]})
    peaks = pd.DataFrame({"frequency": [0.5, 1.0],
                          "amplitude": [0.5, 1.0],
                          "linewidth": [0.1, 0.1],
                          "AIC": [1, 1]})
    data = pd.DataFrame({"numax": [1],
                         "sigmaEnv": [1],
                         "DeltaNu": [1],
                         "eps_p": [1],
                         "gamma0": [1]})
    pds, peaks, data = peak_bag_period_spacing(pds, peaks, data)
    assert "DeltaPi1" in data

def test_peak_bag_period_spacing_001296068():
    """Unit test real data with ID 001296068"""

    pds = pd.read_csv("tests/data/test_peak_bag_period_spacing/pds_bgr.csv")
    peaks = pd.read_csv("tests/data/test_peak_bag_period_spacing/peaksMLE.csv")
    data = pd.read_csv("tests/data/test_peak_bag_period_spacing/summary.csv")

    pds, peaks, data = peak_bag_period_spacing(pds, peaks, data)

    pds_reference = pd.read_csv("tests/data/test_peak_bag_period_spacing/result/pds_bgr.csv")
    peaks_reference = pd.read_csv("tests/data/test_peak_bag_period_spacing/result/peaksMLE.csv")
    data_reference = pd.read_csv("tests/data/test_peak_bag_period_spacing/result/summary.csv")

    pd.testing.assert_frame_equal(pds, pds_reference)
    pd.testing.assert_frame_equal(peaks, peaks_reference)
    pd.testing.assert_frame_equal(data, data_reference)
