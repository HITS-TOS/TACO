import pandas as pd
from taco import peak_bag_mode_id02

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
    peaks, data = peak_bag_mode_id02(pds, peaks, data)
    assert "frequency" in peaks
