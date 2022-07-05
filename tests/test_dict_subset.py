import pytest

def test_dict_subset():
    d1 = {'a': 1, 'b': 2}
    d2 = {'a': 1, 'b': 2}
    d3 = {'a': 1}
    d4 = {'a': 3}
    d5 = {'c': 3}
    assert d2.items() <= d1.items()
    assert d3.items() <= d1.items()
    assert d4.keys() <= d1.keys()
    assert not d5.keys() <= d1.keys()
