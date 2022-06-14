import math

import pandas as pd
import pytest
from taco import filter


def test_filter():
    ts = pd.DataFrame({"time": [1, 2, 3, 4], "flux": [1, -1, math.inf, -1]})
    ts_filtered = filter.filter(ts, width=1)
    assert len(ts_filtered) == 3
