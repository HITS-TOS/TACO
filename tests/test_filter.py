import random

import pandas as pd
import pytest
from taco import filter


def test_filter():
    ts = pd.DataFrame({"time": [1, 2, 3, 4], "flux": [1, -1, 1, -1]})
    filter.filter(ts)
