""" Check behavior of dict subsets """

def test_dict_subset():
    """ Check behavior of dict subsets """
    d_1 = {'a': 1, 'b': 2}
    d_2 = {'a': 1, 'b': 2}
    d_3 = {'a': 1}
    d_4 = {'a': 3}
    d_5 = {'c': 3}
    assert d_2.items() <= d_1.items()
    assert d_3.items() <= d_1.items()
    assert d_4.keys() <= d_1.keys()
    assert not d_5.keys() <= d_1.keys()
