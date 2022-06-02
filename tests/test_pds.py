import pytest

from taco import pds

def test_zero_division():
    with pytest.raises(ZeroDivisionError):
        1 / 0

# def test_pds():
#     ts = pd.Datafram
#     assert pds(ts) == expected
