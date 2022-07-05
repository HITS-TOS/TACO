import pandas as pd
import pytest
import yaml
from taco import background_fit

yaml_string = """
logfile: 'pds_fit.log'
bkg_model: 'KeplerBg3Comp'
output_pds_pgr: 'pds_bgr.csv'
output_ofac_bgr: 'ofac_pds_bgr.csv'
output_backg_model: 'out_backg_model.pkl'
output_quantiles: 'pds_fit_quantiles.csv'
mcmc:
    bins: 300
    posterior: 'pds_fit_posterior.h5'
    save_posteriors: false
    maxsteps: 5000
    minsteps: 2000
    nwalkers: 50
    nwarmup: 1000
"""

def test_settings():
    kwargs = yaml.load(yaml_string, Loader = yaml.Loader)
    settings = Settings(**kwargs)

def test_background_fit():
    pds = pd.DataFrame({"frequency": [1, 2, 3, 4], "power": [1, 2, 3, 4]})
    background_fit(pds, 1.0, 1.0)
    assert True
