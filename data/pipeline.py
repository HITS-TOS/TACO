#!/bin/python

""" TACO pipline module """

import argparse
from pathlib import Path
import yaml

import pandas as pd
import taco


def pipeline(argv):
    """ TACO pipeline """

    # Read pipeline settings
    with open(argv.settings_file, 'r', encoding="utf-8") as stream:
        settings = yaml.load(stream, Loader = yaml.Loader)

    input_files = [f for f in Path(argv.input_directory).iterdir()
        if (f.is_file() and f.suffix == '.dat')]

    if argv.verbose > 0:
        print('Number of input files: ', len(input_files))

    # Loop over input directories
    for input_file in input_files:

        input_name = input_file.stem
        if argv.verbose > 0:
            print('Current input name: ', input_name)

        Path(argv.output_directory, input_name).mkdir(exist_ok = True)
        ts_raw = pd.read_csv(input_file, comment = '#', header = None, delim_whitespace = True)

        # 0) Filter
        filter_settings = settings['pipeline'][0]['filter']
        ts_filtered, variance = taco.filter(ts_raw,
            width = filter_settings['width'],
            remove_gaps = filter_settings['remove_gaps'])
        if 'output' in filter_settings:
            ts_filtered.to_csv(Path(argv.output_directory, input_name, filter_settings['output']), index = False)

        # 1) PDS
        pds_settings = settings['pipeline'][1]['pds']
        pds, nyquist = taco.calc_pds(ts_filtered)
        if 'output' in pds_settings:
            pds.to_csv(Path(argv.output_directory, input_name, pds_settings['output']), index = False)

        # 2) Estimate numax
        # taco.numax_estimate(pds, variance, nyquist,
        #     filterwidth = settings['pipeline'][3]['numax_estimate']['filter_width'])

        # TODO ...
        # background_fit(bins=300)
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
    parser.add_argument('--verbose', '-v', default=0, type=int,
                        help="Print level.")

    pipeline(parser.parse_args())
