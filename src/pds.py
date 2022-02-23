## Calculating a Lomb-Scargle periodogram.
## The input is assumed to have time-stamps with units of days
## The periodogram is given in units of micro-Hertz

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import lightkurve as lk

def compute_conversion_factor(ts):
    """
    Compute conversion factor from power spectrum of window function
    # THIS SHOULD BE DONE USING LOMB-SCARGLE AND NOT FFT!
    """
    time_steps = np.array(round(ts.time/np.median(np.diff(ts.time))))
    time_steps_idx = (time_steps - time_steps.min() + 1).astype(int)

    new_time = np.linspace(1, max(time_steps_idx), max(time_steps_idx))*np.median(np.diff(ts.time)) + ts.time.iloc[0]
    # Compute window function
    window = np.zeros_like(new_time)
    window[time_steps_idx-1] = 1
    # Compute power spectrum and calibrate properly
    bw = 1/((new_time.max() - new_time.min())*86400.0) * 1e6
    fftpower = np.abs(np.fft.fft(window))**2/len(window)**2 / bw

    return np.sum(fftpower) * bw

def compute_conversion_factor_LOMBSCARGLE(ts):
    lc = lk.LightCurve(time=ts["time"], flux=(ts["flux"]/ts["flux"])).normalize(unit='ppm')
    ps = lc.to_periodogram(normalization="psd")
    print(ps.power.value)
    return np.sum(ps.power.value)

    
def main(argv):

    # Still need to do a check to make sure that timeseries exists
    # Also need to check to make sure summary file exists too

    # Read in summary file
    summary = pd.read_csv(argv.summary)

    # Read in time-series
    ts = pd.read_csv(argv.tseries)

    # Compute Nyquist frequency from time stamps
    nuNyq_d = 1 / (2 * np.median(np.diff(ts["time"])))
    
    # Add change in here so that checks to see if data in ppm or normalized flux already
    # Basically check if standard deviation is very small and mean is close to 1
    if (np.std(ts["flux"]) < 1) and (np.abs(np.mean(ts["flux"]) - 1) < 0.1):
        # Create lightkurve lightcurve from data - flux needs to be normalized about 1
        lc = lk.LightCurve(time=ts["time"], flux=ts["flux"]).normalize(unit='ppm')
    elif (np.std(ts["flux"]) < 1) and (np.abs(np.mean(ts["flux"]) - 0.0) < 0.1):
        # Create lightkurve lightcurve from data - flux needs to be normalized about 1
        lc = lk.LightCurve(time=ts["time"], flux=1+(ts["flux"])).normalize(unit='ppm')
    else:
        # Create lightkurve lightcurve from data - flux needs to be normalized about 1
        lc = lk.LightCurve(time=ts["time"], flux=1+(ts["flux"]/1e6)).normalize(unit='ppm')

    # Compute periodogram with psd normalisation
    ps = lc.to_periodogram(normalization="psd")
    # Length of timeseries in days
    ts_length = ts["time"].max() - ts["time"].min()
    
    # Normalisation factor
    factor = compute_conversion_factor(ts)
    # Normalise power spectrum accounting for "fill"
    pds = pd.DataFrame(data=np.c_[ps.frequency.value, ps.power.value/factor], columns=['frequency', 'power'])

    # Generate an oversampled spectrum
    if ts_length < 365.0:
        ops = lc.to_periodogram(normalization="psd", oversample_factor=5)
    else:
        ops = lc.to_periodogram(normalization="psd", oversample_factor=int(argv.ofac))

    ofac_pds = pd.DataFrame(data=np.c_[ops.frequency.value, ops.power.value/factor], columns=['frequency', 'power'])

    # Add nyquist to summary file
    summary["nuNyq"] = pds["frequency"].iloc[-1]
    
    # Save out the results
    pds.to_csv(argv.output, index=False)
    ofac_pds.to_csv(argv.oversampled_output, index=False)
    summary.to_csv(argv.summary, index=False)

## ======================
## Command-line arguments
## ======================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Make Lomb-Scargle periodogram. We assume the input time is given in days and give the periodogram frequencies in micro-Hertz. The maximum frequency sampled is (1/(2*dt) with dt being the median time sampling interval.")
    parser.add_argument('--tseries', default='filtered.csv',
                        help="File name of the time-series in csv format with two columns: `time` and `flux`.")
    parser.add_argument('--output', default='pds.csv',
                        help="File name of the saved periodogram.")
    parser.add_argument('--oversampled_output', default='ofac_pds.csv',
                        help="File name of the saved oversampled periodogram.")
    parser.add_argument('--ofac', default=2, type=int,
                        help="Oversampling factor to use in oversampled periodogram (default is 2).")
    parser.add_argument('--summary', dest='summary', default='summary.csv',
                        help="Summary file to save some periodogram statistics")
    argv = parser.parse_args()
    
    main(argv)
