#!/bin/python

import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, '../src')
from taco import filter, pds


def main(argv):

    for directory in [f for f in os.scandir(argv.base_directory) if os.path.isdir(f)]:

        print('Current directory: ', directory.name)
        ts_raw = pd.read_csv(os.path.join(directory, 'raw.dat'), comment = '#', header = None, sep = '\s+')

        ts_filtered, var = filter.filter(ts_raw)
        freq, nyquist = pds.calc_pds(ts_filtered, ofac=2)
        print ('nyquist = ', nyquist)
        # numax_estimate(pds, var, nyquist, filterwidth=0.2)
        # background_fit(bins=300)
        # background_summary()
        # peakFind(snr=1.1, prob=0.0001, minAIC=2)
        # peaksMLE(minAIC=2)
        # peakBagModeId02()
        # peakFind(snr=1.1, prob=0.0001, minAIC=2, removel02=TRUE)
        # peaksMLE(minAIC=2, removel02=TRUE, init=peaksMLE.csv, mixedpeaks=mixedpeaks.csv)
        # peaksMLE(minAIC=2, finalfit=TRUE)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="TACO workflow")
    parser.add_argument('base_directory', default='.',
                        help="Base directory of processable raw data.")
    argv = parser.parse_args()
    
    main(argv)
