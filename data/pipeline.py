#!/bin/python

import argparse
import os
import pathlib

import pandas as pd
from taco import filter, pds


def main(argv):
    """ TACO pipeline """

    for directory in [f for f in os.scandir(argv.input_directory) if os.path.isdir(f)]:

        print('Current directory: ', directory.name)
        pathlib.Path(os.path.join(argv.output_directory, directory.name)).mkdir(exist_ok=True)
        ts_raw = pd.read_csv(os.path.join(directory, 'raw.dat'),
            comment = '#', header = None, sep = '\s+')

        ts_filtered, _ = filter.filter(ts_raw)
        filtered_filename = os.path.join(argv.output_directory, directory.name, 'filtered.cvs')
        print(filtered_filename)
        ts_filtered.to_csv(filtered_filename, index=False)

        _, nyquist = pds.calc_pds(ts_filtered, ofac=2)
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
    parser.add_argument('--input_directory', '-i', default='.',
                        help="Input directory of processable raw data.")
    parser.add_argument('--output_directory', '-o', default='.',
                        help="Output directory for resulting data.")
    argv = parser.parse_args()

    main(argv)
