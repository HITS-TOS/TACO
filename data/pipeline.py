#!/bin/python

""" TACO pipline module """

import argparse
import os
import pathlib
import yaml

import pandas as pd
import taco


def pipeline(argv):
    """ TACO pipeline """

    # Read pipeline settings
    with open(argv.settings_file, 'r', encoding="utf-8") as stream:
        data = yaml.load(stream, Loader=yaml.Loader)

    # Loop over input directories
    for directory in [f for f in os.scandir(argv.input_directory) if os.path.isdir(f)]:

        print('Current directory: ', directory.name)
        pathlib.Path(os.path.join(argv.output_directory, directory.name)).mkdir(exist_ok=True)
        ts_raw = pd.read_csv(os.path.join(directory, 'raw.dat'),
            comment = '#', header = None, delim_whitespace=True)

        # 1) Filter
        ts_filtered, _ = taco.filter(ts_raw, width = data['pipeline'][0]['filter']['width'])
        filtered_filename = os.path.join(argv.output_directory, directory.name, 'filtered.cvs')
        ts_filtered.to_csv(filtered_filename, index=False)

        # 2) PDS
        pds_data, _ = taco.calc_pds(ts_filtered, ofac=2)
        pds_filename = os.path.join(argv.output_directory, directory.name, 'pds.cvs')
        pds_data.to_csv(pds_filename, index=False)

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
    parser = argparse.ArgumentParser(description="TACO pipeline")

    parser.add_argument('--input_directory', '-i', default='.',
                        help="Input directory of processable raw data.")
    parser.add_argument('--output_directory', '-o', default='.',
                        help="Output directory for resulting data.")
    parser.add_argument('--settings-file', '-s', default='pipeline_settings.yaml',
                        help="File with pipeline settings in Yaml.")

    pipeline(parser.parse_args())