import pandas as pd
import pytest
from taco.peak_find import *

def test_peak_find():
    pds = pd.DataFrame({"frequency": [1, 2, 3, 4], "power": [1, 2, 3, 4]})
    oversampled_pds = pd.DataFrame({"frequency": [1, 2, 3, 4], "power": [1, 2, 3, 4]})
    peak_find(pds, oversampled_pds)
    assert 1 == 1
