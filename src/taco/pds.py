import lightkurve as lk
import numpy as np
import pandas as pd


def compute_conversion_factor(ts):
    """
    Compute conversion factor from power spectrum of window function
    # THIS SHOULD BE DONE USING LOMB-SCARGLE AND NOT FFT!
    """
    time_steps = np.array(round(ts.time / np.median(np.diff(ts.time))))
    time_steps_idx = (time_steps - time_steps.min() + 1).astype(int)

    new_time = (
        np.linspace(1, max(time_steps_idx), max(time_steps_idx))
        * np.median(np.diff(ts.time))
        + ts.time.iloc[0]
    )
    # Compute window function
    window = np.zeros_like(new_time)
    window[time_steps_idx - 1] = 1
    # Compute power spectrum and calibrate properly
    bw = 1 / ((new_time.max() - new_time.min()) * 86400.0) * 1e6
    fftpower = np.abs(np.fft.fft(window)) ** 2 / len(window) ** 2 / bw

    return np.sum(fftpower) * bw


def calc_pds(ts, ofac = 1):
    """
    Calculating a Lomb-Scargle periodogram.

    Parameters:
        ts(pandas.DataFrame):Time series with units of days.
            Columns:
                Name: time, dtype: float[days]
                Name: flux, dtype: float
        ofac(int):Oversampling factor (default is 1).

    Returns:
        pds(pandas.DataFrame):Periodogram
            Columns:
                Name: frequency, dtype: float[micro-Hertz]
                Name: power, dtype: float
        nuNyq(float):Nyquist frequency
    """

    # Add change in here so that checks to see if data in ppm or normalized flux already
    # Basically check if standard deviation is very small and mean is close to 1
    if (np.std(ts["flux"]) < 1) and (np.abs(np.mean(ts["flux"]) - 1) < 0.1):
        # Create lightkurve lightcurve from data - flux needs to be normalized about 1
        lc = lk.LightCurve(time=ts["time"], flux=ts["flux"]).normalize(unit="ppm")
    elif (np.std(ts["flux"]) < 1) and (np.abs(np.mean(ts["flux"]) - 0.0) < 0.1):
        # Create lightkurve lightcurve from data - flux needs to be normalized about 1
        lc = lk.LightCurve(time=ts["time"], flux=1 + (ts["flux"])).normalize(unit="ppm")
    else:
        # Create lightkurve lightcurve from data - flux needs to be normalized about 1
        lc = lk.LightCurve(time=ts["time"], flux=1 + (ts["flux"] / 1e6)).normalize(
            unit="ppm"
        )

    # Compute periodogram with psd normalisation
    ps = lc.to_periodogram(normalization="psd", oversample_factor=ofac)

    # Normalisation factor
    factor = compute_conversion_factor(ts)
    # Normalise power spectrum accounting for "fill"
    pds = pd.DataFrame(
        data=np.c_[ps.frequency.value, ps.power.value / factor],
        columns=["frequency", "power"],
    )

    nuNyq = pds["frequency"].iloc[-1]

    return pds, nuNyq
