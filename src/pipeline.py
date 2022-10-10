#!/usr/bin/env python3

""" TACO pipline module """

import argparse
import subprocess
from cmath import nan
from pathlib import Path

import pandas as pd
import yaml
from codetiming import Timer

import taco


def get_kic_id(input_file):
    """ Returns the KIC identifier from raw data file """

    kic = float("nan")
    h = open(input_file, 'r')
    content = h.readlines()

    for line in content:
        if "KIC" in line:
            kic = line.split()[-1]

    return kic


def get_git_revision_short_hash() -> str:
    return subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD']).decode('ascii').strip()


def pipeline(argv):
    """ TACO pipeline """

    t = Timer("total")
    t.start()

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

        # Set Kepler Input Catalogue (KIC) identification number and raw_data filename
        data = pd.DataFrame({"KIC": [get_kic_id(input_file)],
                             "raw_data": [input_name],
                             "git-rev-hash": [get_git_revision_short_hash()]})

        # 0) Filter
        ts_filtered, data = taco.filter(ts_raw, data,
            **settings['pipeline'][0]['filter'],
            output_directory = Path(argv.output_directory, input_name))

        # 1) PDS
        pds = taco.calc_pds(ts_filtered, **settings['pipeline'][1]['pds'],
            output_directory = Path(argv.output_directory, input_name))

        # 2) Oversampled PDS
        oversampled_pds = taco.calc_pds(ts_filtered, **settings['pipeline'][2]['oversampled_pds'],
            output_directory = Path(argv.output_directory, input_name))

        # Set Nyquist frequency
        data["nuNyq"] = pds["frequency"].iloc[-1]

        # 3) Estimate numax
        data = taco.numax_estimate(pds, data,
            **settings['pipeline'][3]['numax_estimate'])

        # 4) Background fit
        pds_bgr, oversampled_pds_bgr, data = taco.background_fit(
            pds, oversampled_pds, data,
            **settings['pipeline'][4]['background_fit'])

        # 5) Find peaks
        if argv.verbose > 0:
            print('5) Find peaks')
        peaks = taco.peak_find(pds_bgr, oversampled_pds_bgr, data,
            **settings['pipeline'][5]['peak_find'])

        # 6) MLE
        if argv.verbose > 0:
            print('5) MLE')
        peaks_mle, data = taco.peaks_mle(pds_bgr, peaks, data,
            **settings['pipeline'][6]['peaks_mle'])

        # 7) Bag mode id02
        peaks_mle, data = taco.peak_bag_mode_id02(pds_bgr, peaks_mle, data)

        # 8) Find mixed peaks
        mixed_peaks = taco.peak_find(
            pds_bgr, oversampled_pds_bgr, data, peaks = peaks_mle, removel02 = True,
            **settings['pipeline'][7]['peak_find'])

        # 9) MLE with mixed peaks
        mixed_peaks, data = taco.peaks_mle(
            pds_bgr, peaks_mle, data, mixed_peaks = mixed_peaks, removel02 = True,
            **settings['pipeline'][8]['peaks_mle'])

        # 10) Final fit
        mixed_peaks, data = taco.peaks_mle(pds_bgr, peaks_mle, data,
            mixed_peaks = mixed_peaks, finalfit = True,
            **settings['pipeline'][9]['peaks_mle'])

        # 11) Bag_period_spacing
        pds_bgr, mixed_peaks, data = taco.peak_bag_period_spacing(pds_bgr, mixed_peaks, data,
            **settings['pipeline'][10]['peak_bag_period_spacing'])

        # Write final results
        data.to_csv(Path(argv.output_directory, input_name, "data.csv"), index = False)
        pds.to_csv(Path(argv.output_directory, input_name, "pds.csv"), index = False)
        pds_bgr.to_csv(Path(argv.output_directory, input_name, "pds_bgr.csv"), index = False)
        peaks_mle.to_csv(Path(argv.output_directory, input_name, "peaks_mle.csv"), index = False)

    t.stop()


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
