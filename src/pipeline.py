#!/usr/bin/env python3

""" TACO pipline module """

import argparse
import subprocess
from pathlib import Path

import pandas as pd
import yaml
from codetiming import Timer

import taco
import csv

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
    result = 'N/A'
    try:
        result = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD']).decode('ascii').strip()
    except Exception:
        print('No git revision hash found.')
    return result


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
    input_files.sort()

    if not argv.quiet:
        print('Number of input files: ', len(input_files))

    #output files
    summaryfile = settings['pipeline'][11]['filenames']['summary']
    pdsfile = settings['pipeline'][11]['filenames']['pds']
    oversampled_pdsfile = settings['pipeline'][11]['filenames']['oversampled_pds']
    background_corrected_pdsfile = settings['pipeline'][11]['filenames']['background_corrected_pds']
    background_corrected_oversampled_pdsfile = settings['pipeline'][11]['filenames']['background_corrected_oversampled_pds']
    resolved_modesfile= settings['pipeline'][11]['filenames']['resolved_modes']
    mle_resolved_modesfile= settings['pipeline'][11]['filenames']['mle_resolved_modes']
    mixed_modesfile= settings['pipeline'][11]['filenames']['mixed_modes']
    mle_mixed_modesfile= settings['pipeline'][11]['filenames']['mle_mixed_modes']
    final_modesfile= settings['pipeline'][11]['filenames']['final_modes']

    # Open csv file necessary for numax-dnu internal flag
    contours = pd.read_csv('~/TACO/contour_90pct_interp.csv')

    # open a csv files to log the stars that have been processed + their flags

    if not argv.start_function:
        path_to_file = 'stars.csv'
    else:
        if argv.start_function == "background":
            path_to_file = 'stars_background.csv'
        if argv.start_function == "resolved_modes":
            path_to_file = 'stars_resolved_modes.csv'
        if argv.start_function == "unresolved_modes":
            path_to_file = 'stars_unresolved_modes.csv'
        if argv.start_function == "final_fit":
            path_to_file = 'stars_final_fit.csv'

    path = Path(argv.output_directory, path_to_file)

    if path.is_file():
      print(f'The file {path_to_file} exists')
    else:
        with open(path, mode = 'w') as fstar:
            writer = csv.writer(fstar, delimiter = ',')
            writer.writerow(['ID', 'flag_cv', 'flag_numax', 'flag_bgr', 'flag_mle_resolved',
                             'flag_numax_dnu', 'flag_02','flag_mle_mixed','flag_mle_final','flag_dP'])
        fstar.close()

    data = pd.read_csv(path, dtype='string')
    stars_done = data['ID'].tolist()
    # Loop over input directories
    for input_file in input_files:
        kic = get_kic_id(input_file)
        print('Current input name: ', kic)
        input_name = str(kic).zfill(9)

        if input_name in stars_done:
          print(f'This star has already been analysed and will not be redone')
        else:
            print('.............................')
            print('       KIC ', kic)
            print('.............................')
            if not argv.start_function:
                Path(argv.output_directory, input_name).mkdir(exist_ok = True)
                ts_raw = pd.read_csv(input_file, comment = '#', header = None, delim_whitespace = True)

            with open(path_to_file, mode = 'a') as fstar:
                writer = csv.writer(fstar, delimiter = ',')
                #set all flags to -1
                flag_cv = -1
                flag_numax = -1
                flag_bgr = -1
                flag_mle_resolved = -1
                flag_02 = -1
                flag_mle_mixed = -1
                flag_mle_final = -1
                flag_dP = -1

                cv_method = settings['pipeline'][12]['cv_method']['cv_method_check']

                if not argv.start_function:
                    # Set Kepler Input Catalogue (KIC) identification number and raw_data filename
                    summary = pd.DataFrame({"KIC": [get_kic_id(input_file)],
                                "raw_data": [input_name],
                                "git-rev-hash": [get_git_revision_short_hash()]})

                    # 0) Filter
                    print('0) Filter lightcurves')
                    ts_filtered, summary = taco.filter(ts_raw, summary,
                        **settings['pipeline'][0]['filter'],
                        output_directory = Path(argv.output_directory, input_name))

                    # 1) PDS
                    print('1) Compute PDS')
                    pds = taco.calc_pds(ts_filtered, **settings['pipeline'][1]['pds'],
                        output_directory = Path(argv.output_directory, input_name))

                    # 2) Oversampled PDS
                    print('2) Compute oversampled PDS')
                    oversampled_pds = taco.calc_pds(ts_filtered, **settings['pipeline'][2]['oversampled_pds'],
                        output_directory = Path(argv.output_directory, input_name))

                    # Set Nyquist frequency
                    summary["nuNyq"] = pds["frequency"].iloc[-1]

                    if cv_method:
                        print('\n3) Estimate number of power excess with CV method')
                        results_cv, pds = taco.cv_method(pds)
                        flag_cv = results_cv['flag_cv'][0]

                        # Add interpolation flag and cv_flag to summary file:
                        summary['interp_pds'] = results_cv['interp_spikes_flag'][0]
                        summary['n_osc'] = flag_cv

                        results_cv.to_csv(Path(argv.output_directory, input_name, "cv_analysis.csv"), index = False)
                        summary.to_csv(Path(argv.output_directory, input_name, summaryfile), index = False)

                    else:
                        print('\n3) CV method deactivated. Continue with numax estimation.')

                    # Only check stars if number of possible power-excess is less than 2:
                    if flag_cv < 2:
                        # 4) Estimate numax
                        print('4) Estimate numax')
                        summary, flag_numax = taco.numax_estimate(pds, summary,
                                            **settings['pipeline'][3]['numax_estimate'])

                        summary.to_csv(Path(argv.output_directory, input_name, summaryfile), index = False)

                    else:
                        print('Found more than one power excess. Stopping analysis!')


                else:
                    pathsummary = Path(argv.output_directory, input_name, summaryfile)
                    pathpds = Path(argv.output_directory, input_name, pdsfile)
                    pathoverpds = Path(argv.output_directory, input_name, oversampled_pdsfile)

                    if not pathsummary.is_file() or not pathpds.is_file() or not pathoverpds.is_file():
                        print('Necessary files not available for background fit, star not analysed')
                        flag_numax = 2
                    else:
                        summary = pd.read_csv(Path(argv.output_directory, input_name, summaryfile), delimiter = ',')
                        flag_numax = summary['initial_numax_flag'][0]
                        if cv_method:
                            flag_cv = summary['flag_cv'][0]
                        pds = pd.read_csv(Path(argv.output_directory, input_name, pdsfile), delimiter = ',')
                        oversampled_pds = pd.read_csv(Path(argv.output_directory, input_name, oversampled_pdsfile), delimiter = ',')

                if  flag_cv < 2 and flag_numax <= 1:
                    if not argv.start_function or argv.start_function == "background":
                        # 5) Background fit
                        print('5) Fit background')
                        pds_bgr, oversampled_pds_bgr, summary, flag_bgr = taco.background_fit(pds, oversampled_pds, summary,
                                                                                        **settings['pipeline'][4]['background_fit'],
                                                                                        output_directory = Path(argv.output_directory, input_name))
                        if flag_bgr == 0:
                            summary.to_csv(Path(argv.output_directory, input_name, summaryfile), index = False)
                            pds_bgr.to_csv(Path(argv.output_directory, input_name, background_corrected_pdsfile), index = False)
                            oversampled_pds_bgr.to_csv(Path(argv.output_directory, input_name, background_corrected_oversampled_pdsfile), index = False)
                    else:
                        pds_bgr = pd.read_csv(Path(argv.output_directory, input_name, background_corrected_pdsfile), delimiter = ',')
                        oversampled_pds_bgr = pd.read_csv(Path(argv.output_directory, input_name, background_corrected_oversampled_pdsfile), delimiter = ',')

                    if flag_bgr <= 0:
                        if not argv.start_function or argv.start_function == "background" or argv.start_function == "resolved_modes":
                            # 6) Find peaks
                            print('6) Find resolved peaks')
                            peaks = taco.peak_find(pds_bgr, oversampled_pds_bgr, summary,
                                                    **settings['pipeline'][5]['peak_find'])
                            peaks.to_csv(Path(argv.output_directory, input_name, resolved_modesfile), index = False)

                            # 7) MLE
                            if (len(peaks.frequency)) >= 1:
                                print('7) MLE fit resolved peaks')
                                peaks_mle, flag_mle_resolved, summary = taco.peaks_mle(pds_bgr, peaks, summary,
                                                                                        **settings['pipeline'][6]['peaks_mle'])
                                summary.to_csv(Path(argv.output_directory, input_name, summaryfile), index = False)
                                peaks_mle.to_csv(Path(argv.output_directory, input_name, mle_resolved_modesfile), index = False)

                                # 8) Bag mode id02
                                if (((len(peaks_mle.frequency)) >= 3) and (flag_mle_resolved == 0.0)):
                                    print('8) Identify 0,2 modes')
                                    peaks_mle, flag_02, flag_numax_dnu, summary = taco.peak_bag_mode_id02(pds_bgr, peaks_mle, summary, contours)
                                    summary.to_csv(Path(argv.output_directory, input_name, summaryfile), index = False)
                                    peaks_mle.to_csv(Path(argv.output_directory, input_name, mle_resolved_modesfile), index = False)
                        else:
                            peaks = pd.read_csv(Path(argv.output_directory, input_name, resolved_modesfile), delimiter = ',')
                            if (len(peaks.frequency)) >= 1:
                                flag_mle_resolved = 0
                                peaks_mle = pd.read_csv(Path(argv.output_directory, input_name, mle_resolved_modesfile), delimiter = ',')
                                if 'l' in peaks_mle:
                                    flag_02 = 0
                                else:
                                    flag_02 = 1
                            else:
                                flag_mle_resolved = 1

                        # 9) Find mixed peaks
                        if flag_02 == 0.0:
                            if not argv.start_function or argv.start_function == "background" or argv.start_function == "resolved_modes" or argv.start_function == "unresolved_modes":
                                print('9) Find mixed peaks')
                                mixed_peaks = taco.peak_find(
                                    pds_bgr, oversampled_pds_bgr, summary, peaks = peaks_mle, removel02 = True,
                                    **settings['pipeline'][7]['peak_find'])
                                summary.to_csv(Path(argv.output_directory, input_name, summaryfile), index = False)
                                mixed_peaks.to_csv(Path(argv.output_directory, input_name, mixed_modesfile), index = False)

                                # 10) MLE with mixed peaks
                                print('10) MLE fit mixed peaks')
                                mixed_peaks, flag_mle_mixed, summary = taco.peaks_mle(
                                    pds_bgr, peaks_mle, summary, mixed_peaks = mixed_peaks, removel02 = True,
                                    **settings['pipeline'][8]['peaks_mle'])
                                summary.to_csv(Path(argv.output_directory, input_name, summaryfile), index = False)
                                mixed_peaks.to_csv(Path(argv.output_directory, input_name, mle_mixed_modesfile), index = False)
                            else:
                                mixed_peaks = pd.read_csv(Path(argv.output_directory, input_name, mle_mixed_modesfile), delimiter = ',')
                                if (len(mixed_peaks.frequency) > 1):
                                    flag_mle_mixed = 0
                                else:
                                    flag_mle_mixed = 1

                            # 11) Final fit
                            if (flag_mle_mixed == 0):
                                if not argv.start_function or argv.start_function == "background" or argv.start_function == "resolved_modes" or argv.start_function == "unresolved_modes" or argv.start_function == "final_fit":
                                    print('11) Final fit all peaks')
                                    all_peaks, flag_mle_final, summary = taco.peaks_mle(pds_bgr, peaks_mle, summary,
                                        mixed_peaks = mixed_peaks, finalfit = True,
                                        **settings['pipeline'][9]['peaks_mle'])
                                    summary.to_csv(Path(argv.output_directory, input_name, summaryfile), index = False)
                                    all_peaks.to_csv(Path(argv.output_directory, input_name, final_modesfile), index = False)

                                # 12) Bag mode id3
                                if (((len(peaks_mle.frequency)) >= 3) and (flag_mle_resolved == 0.0)):
                                    print('12) Identify l=3 modes')
                                    all_peaks, summary = taco.peak_bag_mode_id3(pds_bgr, all_peaks, summary)
                                    summary.to_csv(Path(argv.output_directory, input_name, summaryfile), index = False)
                                    all_peaks.to_csv(Path(argv.output_directory, input_name, final_modesfile), index = False)

                                    # 13) Bag_period_spacing
                                    if (flag_mle_final == 0):
                                        print('13) Find period spacing')
                                        pds_bgr, all_peaks, flag_dP, summary = taco.peak_bag_period_spacing(pds_bgr, all_peaks, summary,
                                            **settings['pipeline'][10]['peak_bag_period_spacing'])
                                        summary.to_csv(Path(argv.output_directory, input_name, summaryfile), index = False)
                                        all_peaks.to_csv(Path(argv.output_directory, input_name, final_modesfile), index = False)
                                        pds_bgr.to_csv(Path(argv.output_directory, input_name, background_corrected_pdsfile), index = False)

                # Write final results
                writer.writerow([input_name, flag_cv, flag_numax, flag_bgr, flag_mle_resolved, flag_numax_dnu, flag_02, flag_mle_mixed, flag_mle_final, flag_dP])
            fstar.close()
    t.stop()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="TACO pipeline")
    parser.add_argument('--input_directory', '-i', default='.',
                        help="Input directory of processable raw data (default = '.').")
    parser.add_argument('--output_directory', '-o', default='.',
                        help="Output directory for resulting data (default = '.').")
    parser.add_argument('--settings-file', '-s', default='pipeline_settings.yaml',
                        help="File with pipeline settings in Yaml (default = 'pipeline_settings.yaml').")
    parser.add_argument('--start_function', '-start', default=0,
                        help="Function to start the analysis with: background, resolved_modes, unresolved_modes, final_fit (default = 'lightcurve')")
    parser.add_argument('--verbose', '-v', default=0, action='count',
                        help="Print level.")
    parser.add_argument('--quiet', '-q', action='store_true',
                        help="No output")

    pipeline(parser.parse_args())