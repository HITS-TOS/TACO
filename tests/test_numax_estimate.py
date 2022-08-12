import math
from io import StringIO

import pandas as pd
import pytest
from taco import numax_estimate


def test_numax():
    pds = pd.DataFrame({"frequency": [1, 2, 3, 4], "power": [1, 2, 3, 4]})
    data = pd.DataFrame({"nuNyq": [1], "var": [1]})
    data = numax_estimate(pds, data)
    assert data['numax0'][0] == 1.0
