#!/bin/python3

""" TACO pipline module """

import argparse
from pathlib import Path
import yaml

import pandas as pd
import taco


def pipeline(argv):
    """ TACO pipeline """

    if not argv.quiet:
        print(" ==========")
        print("    TACO")
        print(" ==========\n")
        print('Print level: ', argv.verbose)

    # Read pipeline settings
    with open(argv.settings_file, 'r', encoding="utf-8") as stream:
        settings = yaml.load(stream, Loader = yaml.Loader)

    input_files = [f for f in Path(argv.input_directory).iterdir()
        if (f.is_file() and f.suffix == '.dat')]

    if not argv.quiet:
        print('Number of input files: ', len(input_files))

    # Loop over input directories
    for input_file in input_files:

        input_name = input_file.stem
        if argv.verbose > 0:
            print('Current input name: ', input_name)

        Path(argv.output_directory, input_name).mkdir(exist_ok = True)
        ts_raw = pd.read_csv(input_file, comment = '#', header = None, delim_whitespace = True)

        # 0) Filter
        ts_filtered, variance = taco.filter(ts_raw, **settings['pipeline'][0]['filter'],
            output_directory = Path(argv.output_directory, input_name))

        # 1) PDS
        pds = taco.calc_pds(ts_filtered, **settings['pipeline'][1]['pds'],
            output_directory = Path(argv.output_directory, input_name))

        # 2) Oversampled PDS
        oversampled_pds = taco.calc_pds(ts_filtered, **settings['pipeline'][2]['oversampled_pds'],
            output_directory = Path(argv.output_directory, input_name))

        # 3) Estimate numax
        nyquist = pds["frequency"].iloc[-1]
        numax0 = taco.numax_estimate(pds, variance, nyquist,
            **settings['pipeline'][3]['numax_estimate'])

        # 4) Background fit
        Hmax, Bmax, HBR = taco.background_fit(pds, numax0, nyquist,
            **settings['pipeline'][4]['background_fit'])
    
        # TODO ...
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
    parser.add_argument('--verbose', '-v', default=0, action='count',
                        help="Print level.")
    parser.add_argument('--quiet', '-q', action='store_true',
                        help="No output")

    pipeline(parser.parse_args())
