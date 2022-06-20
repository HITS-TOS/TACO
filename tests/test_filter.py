import math
from io import StringIO

import pandas as pd
import pytest
from taco import filter


def test_filter():
    ts = pd.DataFrame({"0": [1, 2, 3, 4], "1": [1, -1, math.inf, -1]})
    ts_filtered,_ = filter(ts)
    assert len(ts_filtered) == 3

data = StringIO("""
  54953.5392583   4.1796191939225127e+02   5.7956622261711652e+02
  54953.5596930   1.9262097401395727e+02   6.0935807794115999e+02
  54953.5801275   5.9470459113319941e+02   6.1677369479803588e+02
  54953.6005623   8.1752961739378622e+02   6.1995181630812556e+02
  54953.6209969   4.5750780589126629e+02   6.2180556616207275e+02
  54953.6414315   6.2782232176861180e+02   6.2298522516003914e+02
  54953.6618661   1.8514406548519878e+02   6.2380191215863124e+02
  54953.6823009   2.7943136812669246e+02   6.2440081595759875e+02
  54953.7027354   2.3804503046531698e+02   6.2485880121563275e+02
""")

def test_filter_2():
    ts = pd.read_csv(data, comment = '#', header = None, delim_whitespace=True)
    ts_filtered,_ = filter(ts)
    print(ts_filtered)
    assert ts_filtered.iat[0, 2] == pytest.approx(-5.4567, 0.001)
