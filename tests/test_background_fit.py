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

    pds = pd.read_csv("tests/data/test_background_pds.csv")
    ofac_pds = pd.read_csv("tests/data/test_background_ofac_pds.csv")
    data = pd.read_csv("tests/data/test_background_summary.csv")
    _, _, data = background_fit(pds, ofac_pds, data, seed = 42)
    assert data['Hmax'][0] == pytest.approx(10407.676594861878, 1e-6)
    assert data['Bmax'][0] == pytest.approx(3375.215027473658, 1e-6)
    assert data['HBR'][0] == pytest.approx(3.083559568840864, 1e-6)
