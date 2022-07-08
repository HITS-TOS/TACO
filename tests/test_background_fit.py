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
output_ofac_bgr: 'ofac_pds_bgr.csv'
output_pds_pgr: 'pds_bgr.csv'
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

#@pytest.mark.skip(reason="no way of currently testing this")
def test_background_fit():
    pds = pd.DataFrame({"frequency": [1, 2, 3, 4], "power": [1, 2, 3, 4]})
    Hmax, Bmax, HBR = background_fit(pds, 1.0, 1.0, minsteps = 2, maxsteps = 5, nwarmup = 3)
    assert Hmax == 1.0
    assert Bmax == 1.0
    assert HBR == 1.0