import random

import pandas as pd
import pytest
from taco import pds


def test_value_error():
    ts = pd.DataFrame({"time": [1, 2, 3, 4], "flux": [1, -1, 1, -1]})
    with pytest.raises(ValueError):
        pds.calc_pds(ts)


testdata = [
    (
        pd.DataFrame(
            {"time": range(0, 50), "flux": random.sample(range(10, 3000), 50)}
        ),
        pytest.approx(5.6689, 0.0001),
    ),
    (
        pd.DataFrame({"time": range(0, 10), "flux": [1,-1,1,-1,1,-1,1,-1,1,-1]}),
        pytest.approx(5.1440, 0.0001),
    )
]


@pytest.mark.parametrize("ts,expected", testdata)
def test_pds(ts, expected):
    _, nyquist = pds.calc_pds(ts)
    assert nyquist == expected
