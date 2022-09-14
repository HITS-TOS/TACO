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

    if argv.verbose > 1:
        print("settings: ", settings)

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
        ts_filtered, data = taco.filter(ts_raw, **settings['pipeline'][0]['filter'],
            output_directory = Path(argv.output_directory, input_name))

        # 1) PDS
        pds = taco.calc_pds(ts_filtered, **settings['pipeline'][1]['pds'],
            output_directory = Path(argv.output_directory, input_name))

        # 2) Oversampled PDS
        oversampled_pds = taco.calc_pds(ts_filtered, **settings['pipeline'][2]['oversampled_pds'],
            output_directory = Path(argv.output_directory, input_name))

        # 3) Estimate numax
        data["nuNyq"] = pds["frequency"].iloc[-1]
        data = taco.numax_estimate(pds, data,
            **settings['pipeline'][3]['numax_estimate'])

        # 4) Background fit
        pds_bgr, oversampled_pds_bgr, data = taco.background_fit(
            pds, oversampled_pds, data,
            **settings['pipeline'][4]['background_fit'])
    
        # 5) Find peaks
        peaks = taco.peak_find(pds_bgr, oversampled_pds_bgr, data,
            **settings['pipeline'][5]['peak_find'])

        # 6) Peaks MLE
        peaks_mle = taco.peaks_mle(pds_bgr, peaks, data,
            **settings['pipeline'][6]['peaks_mle'])

        # 7) Peaks MLE
        peaks_mle, data = taco.peak_bag_mode_id02(pds_bgr, peaks_mle, data,
            **settings['pipeline'][7]['peak_bag_mode_id02'])

        # 8) Peaks MLE
        mixed_peaks = taco.peak_find(
            pds_bgr, oversampled_pds_bgr, peaks = peaks, removel02 = True,
            **settings['pipeline'][8]['peak_find'])

        # 9) Peaks MLE
        peaks = taco.peaks_mle(
            pds_bgr, peaks_mle, data, removel02 = True, mixed_peaks = mixed_peaks,
            **settings['pipeline'][9]['peaks_mle'])

        # 10) Peaks MLE
        peaks = taco.peaks_mle(pds_bgr, peaks, data, finalfit = True,
            **settings['pipeline'][10]['peaks_mle'])

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="TACO pipeline")
    parser.add_argument('--input_directory', '-i', default='.',
                        help="Input directory of processable raw data (default = '.').")
    parser.add_argument('--output_directory', '-o', default='.',
                        help="Output directory for resulting data (default = '.').")
    parser.add_argument('--settings-file', '-s', default='pipeline_settings.yaml',
                        help="File with pipeline settings in Yaml (default = 'pipeline_settings.yaml').")
    parser.add_argument('--verbose', '-v', default=0, action='count',
                        help="Print level.")
    parser.add_argument('--quiet', '-q', action='store_true',
                        help="No output")

    pipeline(parser.parse_args())
