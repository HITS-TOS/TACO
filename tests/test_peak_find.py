import pandas as pd
import pytest
from taco import peak_find

def test_peak_find():
    pds = pd.DataFrame({"frequency": [1, 2], "power": [1, 0]})
    oversampled_pds = pd.DataFrame({"frequency": [0.5, 1, 1.5, 2], "power": [0, 1, 0, 0]})
    peaks = peak_find(pds, oversampled_pds)
    assert "frequency" in peaks

@pytest.mark.skip(reason="not implemented yet")
def test_peak_find_with_peaks():
    pds = pd.DataFrame({"frequency": [1, 2], "power": [1, 0]})
    oversampled_pds = pd.DataFrame({"frequency": [0.5, 1, 1.5, 2], "power": [0, 1, 0, 0]})
    peaks = pd.DataFrame({"frequency": [1]})
    mixed_peaks = pd.DataFrame({"frequency": [1]})
    peaks = peak_find(pds, oversampled_pds, peaks, mixed_peaks)
    assert "frequency" in peaks
