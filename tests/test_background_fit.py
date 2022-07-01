import pandas as pd
import pytest
from taco import background_fit


def test_background_fit():
    pds = pd.DataFrame({"frequency": [1, 2, 3, 4], "power": [1, 2, 3, 4]})
    background_fit(pds, 1.0, 1.0)
    assert True
