import random

import pandas as pd
import pytest
from taco import pds


def test_zero_division():
    with pytest.raises(ZeroDivisionError):
        1 / 0

def test_pds():
    # ts = pd.DataFrame({"time": [1.0, 2.0, 3.0, 4.0], "flux": [15.0, -3.0, 2.0, 8.0]})
    # ts = pd.DataFrame({"time": [54953.5392, 54953.5596, 54953.5801, 54953.6005],
    #                    "flux": [156.9757, -68.3652, 333.7183, 556.5434]})
    ts = pd.DataFrame({"time": range(0, 50),
                       "flux": random.sample(range(10, 3000), 50)})

    print(ts)
    _, nyquist = pds.calc_pds(ts)
    assert nyquist == 10.0
