import pandas as pd
import pytest
import yaml
from taco.background_fit import *

yaml_string = """
bins: 300
bkg_model: 'KeplerBg3Comp'
logfile: 'pds_fit.log'
maxsteps: 5000
minsteps: 2000
nwalkers: 50
nwarmup: 1000
output_backg_model: 'out_backg_model.pkl'
output_ofac_pds_bgr: 'ofac_pds_bgr.csv'
output_pds_bgr: 'pds_bgr.csv'
output_quantiles: 'pds_fit_quantiles.csv'
posterior: 'pds_fit_posterior.h5'
save_posteriors: false
"""

def test_settings():

    kwargs = yaml.load(yaml_string, Loader = yaml.Loader)
    settings = Settings(**kwargs)
    assert settings.bins == 300
    assert settings.get_mcmc_settings()['bins'] == 300
    assert settings.get_mcmc_settings()['backend_filename'] == 'pds_fit_posterior.h5'

def test_background_fit():

    pds = pd.DataFrame({"frequency": [1, 2], "power": [1, 0]})
    ofac_pds = pd.DataFrame({"frequency": [0.5, 1, 1.5, 2], "power": [0, 1, 0, 0]})
    data = pd.DataFrame({"nuNyq": [1], "numax0": [1]})
    _, _, data = background_fit(pds, ofac_pds, data, minsteps = 2, maxsteps = 5, nwarmup = 3)
    assert data is None

def test_background_fit_001296068():
    """Unit test real data with ID 001296068"""

    pds = pd.read_csv("tests/data/test_background_fit/pds.csv")
    ofac_pds = pd.read_csv("tests/data/test_background_fit/ofac_pds.csv")
    data = pd.read_csv("tests/data/test_background_fit/summary.csv")

    pds_bgr, ofac_pds_bgr, data = background_fit(pds, ofac_pds, data, seed = 42)

    pds_bgr_reference = pd.read_csv("tests/data/test_background_fit/result/pds_bgr.csv")
    ofac_pds_bgr_reference = pd.read_csv("tests/data/test_background_fit/result/ofac_pds_bgr.csv")
    data_reference = pd.read_csv("tests/data/test_background_fit/result/summary.csv")

    pd.testing.assert_frame_equal(pds_bgr, pds_bgr_reference)
    pd.testing.assert_frame_equal(ofac_pds_bgr, ofac_pds_bgr_reference)
    pd.testing.assert_frame_equal(data, data_reference)
