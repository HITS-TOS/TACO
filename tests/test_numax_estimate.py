import pandas as pd
import pytest
from taco import numax_estimate


def test_numax_001296068():
    """Unit test real data with ID 001296068"""

    pds = pd.read_csv("tests/data/test_numax_4_pds.csv")
    data = pd.read_csv("tests/data/test_numax_4_summary.csv")
    data = numax_estimate(pds, data)
    assert not data["numax0_flag"][0]
    assert data["numax_var"][0] == pytest.approx(43.961062750337454, 1e-6)
    assert data["numax_CWTMexHat"][0] == pytest.approx(61.42446564103499, 1e-6)
    assert data["numax_Morlet"][0] == pytest.approx(59.333183743645144, 1e-6)
    assert data["numax0"][0] == pytest.approx(59.333183743645144, 1e-6)


@pytest.mark.skip(reason="small tests needed to test special cases")
def test_numax_small():
    """Unit test of small data set"""
    pds = pd.read_csv("tests/data/test_numax_3_pds.csv")
    data = pd.read_csv("tests/data/test_numax_3_summary.csv")
    data = numax_estimate(pds, data)
    assert not data["numax0_flag"][0]
    assert data["numax_var"][0] == pytest.approx(43.961062750337454, 1e-6)
    assert data["numax_CWTMexHat"][0] == pytest.approx(61.42446564103499, 1e-6)
