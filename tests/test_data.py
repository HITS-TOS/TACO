import pytest
from taco import data

def test_data():
    d = data.TacoData()
    assert d.modeIDFlag == 0
    assert d.numax == 0.0
    assert d.npeaks == 0
