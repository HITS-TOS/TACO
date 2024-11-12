import os
from pathlib import Path

import rpy2.robjects as ro
import rpy2.robjects.conversion as cv
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import STAP

def _none2null(none_obj):
    return ro.r("NULL")

def peak_bag_mode_id02(pds, peaks, data, contours):
    """
    Get Δν and label the l=0,2 peaks

    Parameters:
        pds(pandas.DataFrame):Periodogram
            Columns:
                Name: frequency, dtype: float[micro-Hertz]
                Name: power, dtype: float
        peaks(pandas.DataFrame): Identified peaks. It must contain the l=0,2 modes already identified
            Columns:
                Name: frequency, dtype: float[micro-Hertz]
                Name: linewidth, dtype: float
                Name: amplitude, dtype: float
        data(pandas.DataFrame):Summary data
            Columns:
                Name: numax, dtype: float
    Returns:
        peaks(pandas.DataFrame): Identified peaks.
        flag(int):
            Values:
                1: No peaks detected so not proceeding with mode ID
                2: Not enough peaks detected so not proceeding with mode ID
                3: Numax < 5 uHz and is too low for the automated mode identification to be reliable.
                4: Numax > 0.9 * Nyquist frequency. This is too high for the automated mode identification to be reliable.
                5: Not enough radial modes found to go any further. Mode ID not performed and Δν not estimated.
                6: Delta nu and diffence in central radial frequencies not in line
        data(pandas.DataFrame):Summary data
            Columns:
                DeltaNu(float): Large frequency separation [micro-Hertz]
                DeltaNu_sd(float): DeltaNu standard deviation [micro-Hertz]
                dNu02(float): Small frequency separation l=0,2 [micro-Hertz]
                eps_p(float): p-modes phase term
                eps_p_sd(float): eps_p standard deviation
                alpha(float): curvature
                alpha_sd(float): alpha standard deviation
                gamma0(float): Γ_0(ν_max) of Vrard et al. (2018)
                Central_DeltaNu(float): DeltaNu estimate from 3 central radial modes [micro-Hertz]
                Central_DeltaNu_sd(float): Central_DeltaNu standard deviation [micro-Hertz]
                Central_eps_p(float): eps_p estimate from 3 central radial modes [micro-Hertz]
                Central_eps_p(float): Central_eps_p standard deviation [micro-Hertz]
                Central_alpha(float): alpha estimate from 3 central radial modes [micro-Hertz]
                Central_alpha(float): Central_alpha standard deviation [micro-Hertz]



    """

    with open(Path(Path(__file__).parent, 'peak_bag_mode_id02.R'), 'r') as f:
        owd = os.getcwd()
        os.chdir(Path(__file__).parents[2])
        peak_bag_mode_id02 = STAP(f.read(), "peak_bag_mode_id02_r")
        os.chdir(owd)

        none_converter = cv.Converter("None converter")
        none_converter.py2rpy.register(type(None), _none2null)

        with localconverter(ro.default_converter + pandas2ri.converter + none_converter):
            r_pds = ro.conversion.py2rpy(pds)
            r_peaks = ro.conversion.py2rpy(peaks)
            r_data = ro.conversion.py2rpy(data)
            r_contours = ro.conversion.py2rpy(contours)
            result = peak_bag_mode_id02.peak_bag_mode_id02_r(r_pds, r_peaks, r_data, r_contours)

            peaks = ro.conversion.rpy2py(result['peaks'])
            flag = ro.conversion.rpy2py(result['flag'])
            flag_contour = ro.conversion.rpy2py(result['flag_contour'])
            data = ro.conversion.rpy2py(result['data'])

            return peaks, int(flag[0]), int(flag_contour[0]), data
