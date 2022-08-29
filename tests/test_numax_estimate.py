import math
from io import StringIO

import pandas as pd
import pytest
from taco import numax_estimate


def test_numax():
    pds = pd.DataFrame({"frequency": [1, 2, 3, 4], "power": [1, 2, 3, 4]})
    data = pd.DataFrame({"nuNyq": [1], "var": [1]})
    data = numax_estimate(pds, data)
    assert data['numax0_flag'][0] == False
    assert data['numax_var'][0] == 539379.0772482151
    assert data['numax_CWTMexHat'][0] == 539379.0772482151
