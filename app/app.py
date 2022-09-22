import site
import sys
site.addsitedir('../')

from bokeh.models.layouts import Column
from numpy.lib.function_base import angle
import streamlit as st
import SessionState

from bokeh.models import ColumnDataSource, Whisker, HoverTool, Span
from bokeh.transform import factor_cmap, factor_mark

import glob
import pandas as pd
import numpy as np
import pathlib
import matplotlib.pyplot as plt

from bokeh.plotting import figure
import app_helpers

from bokeh.palettes import Greys256, Colorblind7

from scipy.interpolate import interp1d

from pathlib import Path

import itertools

def find_stars():
    dir = Path().cwd().parent / "result_reference"
    dirs = dir.rglob('*/summary.csv')
    KICs = ['KIC '+str(i).split('/')[-2].lstrip('0') for i in dirs]
    return KICs

#@st.cache
def load_summary(KIC):
    dir = Path().cwd().parent / "result_reference"
    KIC = KIC.lstrip('KIC ').zfill(9)
    summary = pd.read_csv(str(dir)+'/'+str(KIC)+'/summary.csv')
    return summary

#@st.cache
def load_psd(KIC, background_removed=True):
    dir = Path().cwd().parent / "result_reference"
    KIC = KIC.lstrip('KIC ').zfill(9)
    #print(KIC)
    if background_removed == True:
        psd = pd.read_csv(str(dir)+'/'+str(KIC)+'/pds_bgr.csv')
    else:
        psd = pd.read_csv(str(dir)+'/'+str(KIC)+'/pds.csv')
        #st.write(psd)
    return psd

def load_ts(KIC, filtered=True):
    dir = Path().cwd().parent / "result_reference"
    KIC = KIC.lstrip('KIC ').zfill(9)
    #print(KIC)
    if filtered == True:
        ts = pd.read_csv(str(dir)+'/'+str(KIC)+'/filtered.csv')
    else:
        ts = pd.read_csv(str(dir)+'/'+str(KIC)+'/raw.dat', delimiter=r'\s+', names=['time', 'flux'])


 #   if ts.flux.mean() > 1e3:
    #    ts.flux -= ts.flux.mean()

   # if ts.flux.std() < 1:
  #      ts.flux *= 1e6
    return ts

#@st.cache
def load_peaks(KIC, peakFind=False):
    dir = Path().cwd().parent / "result_reference"
    KIC = KIC.lstrip('KIC ').zfill(9)
    if peakFind == False:
        peaks = pd.read_csv(str(dir)+'/'+str(KIC)+'/peaksMLE.csv')
    else:
        peaks = pd.read_csv(str(dir)+'/'+str(KIC)+'/peaks.csv')

    return peaks

def visualise_timeseries(filtered_ts, unfiltered_ts, summary):
    """
    Visualise the timeseries data
    """
    p = figure(
        title="",
        x_axis_label="Time (days)",
        y_axis_label=r"Flux (ppm)",
        match_aspect=True,
        x_range=(filtered_ts.time.min(), filtered_ts.time.max()),
        #y_range=(psd.power.min()*0.5, psd.power.max()*1.1)
    )
    p.line(filtered_ts.time, filtered_ts.flux, color='black', legend_label=r'Filtered timeseries')

    if st.sidebar.checkbox("Show unfiltered timeseries"):
        p.line(unfiltered_ts.time, unfiltered_ts.flux, color='red', legend_label=r'Unfiltered timeseries')

    p.legend.click_policy="hide"
    st.bokeh_chart(p)



def visualise_psd(psd, summary):

    #if st.sidebar.checkbox('Show power spectrum'):

    # If data is short cadence then don't plot it all
    if len(psd) > 6e4:
        st.write("Data is short cadence, only plotting data below 300 uHz!")
        cond = psd.frequency < 300
        psd = psd.loc[cond, ]
        #psd = psd.iloc[::10, :]

    overplot_fit = False
    if st.sidebar.checkbox("Overplot background fit"):
        overplot_fit = True

    if st.sidebar.checkbox('log x-axis'):
        x_axis_type = 'log'
    else:
        x_axis_type = 'linear'
    if st.sidebar.checkbox('log y-axis'):
        y_axis_type = 'log'
    else:
        y_axis_type = 'linear'

    Initial_numax = False
    if st.sidebar.checkbox('Plot initial numax estimates'):
        Initial_numax = True
    else:
        Initial_numax = False

    p = figure(
        title="",
        x_axis_label="Frequency (μHz)",
        x_axis_type=x_axis_type,
        y_axis_label=r"PSD (ppm^{2}/μHz)",
        y_axis_type=y_axis_type,
        match_aspect=True,
        x_range=(1, psd.frequency.max()),
        y_range=(psd.power.min()*0.5, psd.power.max()*1.1)
    )

    p.line(psd.frequency, psd.power, color="black", legend_label=r'Data')

    if overplot_fit == True:

        theta = summary.loc[:, 'Pn':'sigmaEnv']

        n_comps = np.array([i.startswith("A") for i in list(theta)]).sum()
        n_gauss = np.array([i.startswith("Pg") for i in list(theta)]).sum()

        comps, comp_names, bg_fit, bg_fit_no_osc = app_helpers.bgModel(psd.frequency,
                                                                        theta,
                                                                        summary['nuNyq'].values,
                                                                        n_comps,
                                                                        n_gauss)

        gran_comp = 0
        gauss_comp = 0
        # create a color iterator
        colours = itertools.cycle(Colorblind7)
        for (i, colour) in zip(range(len(comps)), colours):
            if comp_names[i].startswith("gran"):
                if gran_comp == 0:
                    label = "1st Granulation Component"
                elif gran_comp == 1:
                    label = "2nd Granulation Component"
                else:
                    label = f"{i+1}th Granulation Component"
                p.line(psd.frequency, comps[i], line_dash='dashed', color=colour, line_width=3, legend_label=label)
                gran_comp += 1
            if comp_names[i].startswith("gauss"):
                if (gauss_comp == 0) & (n_gauss == 1):
                    label = "Oscillations"
                elif (gauss_comp == 0) & (n_gauss > 1):
                     label = "1st set of Oscillations"
                elif gauss_comp == 1:
                    label = "2nd set of Oscillations"
                else:
                    label = f"{i+1}th set of Oscillations"
                p.line(psd.frequency, comps[i], line_dash='dashed', color=colour, line_width=3, legend_label=label)
                gauss_comp += 1
        #p.line(psd.frequency, comp2, color="green", line_dash='dashed', line_width=3, legend_label=r'2nd Granulation Component')
        #p.line(psd.frequency, comp3, color="magenta", line_dash='dashed', line_width=3, legend_label=r'3rd Granulation Component')
        #p.line(psd.frequency, gauss, color="yellow", line_dash='dashed', line_width=3, legend_label=r'Oscillations')
        # White noise is always the last in the components list
        p.line(psd.frequency, comps[-1], color="blue", line_dash='dashed', line_width=3, legend_label=r'White Noise')
        p.line(psd.frequency, bg_fit_no_osc, color="red", line_dash='dashed', line_width=3, legend_label=r'Full model (no oscillations)')
        p.line(psd.frequency, bg_fit, color="red", line_width=3, legend_label=r'Full model')

    if Initial_numax == True:
        power = np.linspace(psd.power.min()*0.5, psd.power.max()*1.1, 1000)
        p.line(summary.numax_var.values*np.ones_like(power), power, color="turquoise", line_width=3, line_dash='dotted', legend_label=r'numax estimate from variance')
        try:
            p.line(summary.numax_CWTMexHat.values*np.ones_like(power), power, color="violet", line_width=3, line_dash='dotted', legend_label=r'numax estimate from CWT MexHat')
        except:
            p.line(summary.numax_CWTTree.values*np.ones_like(power), power, color="violet", line_width=3, line_dash='dotted', legend_label=r'numax estimate from CWT MexHat')
        finally:
            pass
        p.line(summary.numax_Morlet.values*np.ones_like(power), power, color="goldenrod", line_width=3, line_dash='dotted', legend_label=r'numax estimate from Morlet wavelet')


    p.legend.click_policy="hide"
    if (x_axis_type == 'log') and (y_axis_type == 'log'):
        p.legend.location = 'bottom_left'
    st.bokeh_chart(p)

def visualise_pds_bgr(psd_bgr, peaks, summary, peakFind=None):
    #if st.sidebar.checkbox('Show background-removed power spectrum'):
    psd_bgr = psd_bgr.loc[np.abs(psd_bgr['frequency'].values - summary['numax'].values) < 3*summary['sigmaEnv'].values, ]

    overplot_peaks = st.checkbox("Overplot peaks from peakFind")

    overplot_fit = st.checkbox("Overplot MLE fit")

    overplot_fit_even = st.checkbox("Overplot MLE fit even modes")


    p = figure(
        title="",
        x_axis_label="Frequency (μHz)",
        y_axis_label="SNR",
        match_aspect=True,
        x_range=(psd_bgr.frequency.min(), psd_bgr.frequency.max()),
        y_range=(0, psd_bgr.power.max()*1.1),
        plot_height=400,
        plot_width=600
    )
    #p.min_border_left = 0
    p.line(psd_bgr.frequency, psd_bgr.power, color="black")

    AIC_cut = st.sidebar.slider('AIC cut', min_value=-100.0, max_value=100., value=2.)

    if overplot_peaks:

        if peakFind is not None:
            model_pf = app_helpers.construct_peaksmodel(psd_bgr, peakFind.loc[peakFind.AIC.values > AIC_cut, ])
            p.line(psd_bgr.frequency, model_pf, color="purple", line_width=3)

    if overplot_fit:

        model = app_helpers.construct_peaksmodel(psd_bgr, peaks.loc[peaks.AIC.values > AIC_cut, ])
        p.line(psd_bgr.frequency, model, color="orange", line_width=3)

    if overplot_fit_even:

        model = app_helpers.construct_peaksmodel(psd_bgr, peaks.loc[(peaks.l % 2 == 0) & (peaks.AIC.values > AIC_cut), ])
        p.line(psd_bgr.frequency, model, color="red", line_width=3)



        #model_l02, model_l1 = app_helpers.construct_MLEmodel(psd_bgr, peaks)
        # Add 1 because no background added to individual models
        #full_model = model_l02 + model_l1 + 1
        #if overplot_fit == "Full fit":
        #    p.line(psd_bgr.frequency, full_model, color="red", line_width=3, legend_label="Full MLE fit")
        #elif overplot_fit == "Odd/Even degrees":
            # Have to add background of 1 as not included when individual models created
        #    p.line(psd_bgr.frequency, model_l02+1, color="green", line_width=3, legend_label="l=0/2 (Even)")
        #    p.line(psd_bgr.frequency, model_l1+1, color="red", line_width=3, legend_label="l=1/3 (Odd)")
    if st.sidebar.checkbox("Mode identification"):
        p.circle(peaks.loc[(peaks.l == 0) & (peaks.AIC.values > AIC_cut), 'frequency'], len(peaks.loc[(peaks.l == 0) & (peaks.AIC.values > AIC_cut), 'frequency'])*[psd_bgr.power.max()/2],
                size=10, color="red", alpha=0.5, legend_label='l=0')
        p.square(peaks.loc[(peaks.l == 2) & (peaks.AIC.values > AIC_cut), 'frequency'],
                len(peaks.loc[(peaks.l == 2) & (peaks.AIC.values > AIC_cut), 'frequency'])*[psd_bgr.power.max()/2],
                size=10, color="green", alpha=0.5, legend_label='l=2')
        p.hex(peaks.loc[(peaks.l == 3) & (peaks.AIC.values > AIC_cut), 'frequency'],
                len(peaks.loc[(peaks.l == 3) & (peaks.AIC.values > AIC_cut), 'frequency'])*[psd_bgr.power.max()/2],
                size=10, color="orange", alpha=0.5, legend_label='l=3')
         # Colour by rotation splitting as well
   #     if st.sidebar.checkbox("Include rotational splitting"):
   #         p.triangle(peaks.loc[(peaks.l == 1) & (peaks.m == -1) & (peaks.AIC.values > AIC_cut), 'frequency'],
   #                    len(peaks.loc[(peaks.l == 1) & (peaks.m == -1) & (peaks.AIC.values > AIC_cut), 'frequency'])*[psd_bgr.power.max()/2],
   #                    angle=-np.pi/2, size=10, alpha=0.5, color="blue", legend_label=r'l=1, m=-1')
   #         p.triangle(peaks.loc[(peaks.l == 1) & (peaks.m == 0) & (peaks.AIC.values > AIC_cut), 'frequency'],
   #                    len(peaks.loc[(peaks.l == 1) & (peaks.m == 0) & (peaks.AIC.values > AIC_cut), 'frequency'])*[psd_bgr.power.max()*0.51],
   #                 angle=np.pi, size=10, alpha=0.5, color="blue", legend_label=r'l=1, m=0')
   #         p.triangle(peaks.loc[(peaks.l == 1) & (peaks.m == 1) & (peaks.AIC.values > AIC_cut), 'frequency'],
   #                    len(peaks.loc[(peaks.l == 1) & (peaks.m == 1) & (peaks.AIC.values > AIC_cut), 'frequency'])*[psd_bgr.power.max()/2],
    #                   angle=np.pi/2, size=10, alpha=0.5, color="blue", legend_label=r'l=1, m=+1')
    #    else:
    #        p.triangle(peaks.loc[(peaks.l == 1) & (peaks.AIC.values > AIC_cut), 'frequency'],
    #                    len(peaks.loc[(peaks.l == 1) & (peaks.AIC.values > AIC_cut), 'frequency'])*[psd_bgr.power.max()/2],
    #                    size=10, color="blue", alpha=0.5, legend_label='l=1')
    p.legend.click_policy="hide"
    st.bokeh_chart(p)

    res_plot = st.sidebar.selectbox(
              "plot residuals",
              ["","full fit","even mode fit"])

    if res_plot == "full fit":
        model = app_helpers.construct_peaksmodel(psd_bgr, peaks.loc[peaks.AIC.values > AIC_cut, ])

    if res_plot == "even mode fit":
        model = app_helpers.construct_peaksmodel(psd_bgr, peaks.loc[(peaks.l % 2 == 0) & (peaks.AIC.values > AIC_cut), ])


    if res_plot != "":
        res = psd_bgr.power / model
        p = figure(
            title="",
            x_axis_label="Frequency (μHz)",
            y_axis_label="SNR",
            match_aspect=True,
            x_range=(psd_bgr.frequency.min(), psd_bgr.frequency.max()),
            y_range=(0, res.max()),
            plot_height=400,
            plot_width=600
        )
        p.line(psd_bgr.frequency, psd_bgr.power / model, color="black")

  #      if overplot_fit == "Full fit":
  #          p.line(psd_bgr.frequency, psd_bgr.power.values/full_model, color="black")
  #      elif overplot_fit == "Odd/Even degrees":
  #          p.line(psd_bgr.frequency, psd_bgr.power.values/(model_l02+1), color="black")
  #          p.line(psd_bgr.frequency, model_l1+1, color="red", line_width=3, legend_label="l=1/3 (Odd)")
        st.bokeh_chart(p)




def visualise_echelle(psd_bgr, peaks, summary, session):

    MARKERS = ["circle", "square", "triangle", "hex"]
    MODE_DEGREE = ["l=0", "l=2", "l=1", "l=3"]
    COLOURS = ["red", "green", "blue", "orange"]

    numax = summary['numax'].item()
    orig_DeltaNu = summary['DeltaNu'].item()
    if not np.isfinite(orig_DeltaNu):
        orig_DeltaNu = 0.254 * numax **0.756
    orig_eps = summary['eps_p'].item()
    if not np.isfinite(orig_eps):
        orig_eps = 1.0
    orig_alpha = summary['alpha'].item()
    if not np.isfinite(orig_alpha):
        orig_alpha = 0.0
    orig_d02 = summary['dNu02'].item()
    if not np.isfinite(orig_d02):
        orig_d02 = 0.125 * orig_DeltaNu
    if (orig_eps) > 1.7 or (orig_eps < 0):
        st.write('Epsilon p is outside correct range, something went wrong!')
        orig_eps = 0.0

    # Save original deltanu in case reset needed
    DeltaNu = orig_DeltaNu
    eps_p = orig_eps
    alpha = orig_alpha
    d02 = orig_d02

    # Set up sliders
    reset = st.sidebar.button('Reset parameters')
    if reset:
        session.run_id += 1
    save = st.sidebar.button('Save parameters')
    if save:
        st.write('Doesn\'t do anything at the moment, but will do soon!')
    DeltaNu = st.sidebar.slider('Δν (μHz)', min_value=orig_DeltaNu*0.85, max_value=orig_DeltaNu*1.15, value=orig_DeltaNu, step = orig_DeltaNu/1e3, key=session.run_id, format="%.3f")
    eps_p = st.sidebar.slider('ε', min_value=0.0, max_value=1.7, value=orig_eps, key=session.run_id, format="%.3f")


    AIC_cut = st.sidebar.slider('AIC cut', min_value=-10.0, max_value=20., value=2.0, key=session.run_id)

    psd_bgr = psd_bgr.loc[np.abs(psd_bgr['frequency'].values - summary['numax'].values) < 3*summary['sigmaEnv'].values, ]

    p = figure(
    title="",
    x_axis_label="Frequency modulo Δν (μHz)",
    y_axis_label="Frequency (μHz)",
    match_aspect=True,
    x_range=(0, DeltaNu),
    y_range=(psd_bgr.frequency.min(), psd_bgr.frequency.max()),
    plot_height=400,
    plot_width=600
    )

    dh = DeltaNu

    echelle_offset = 0.2

    # l=0 x and y values
    #'eps_p mod Dnu', eps_p*DeltaNu % DeltaNu
    l0_x = ((peaks.loc[(peaks.l == 0) & (peaks.AIC.values >= AIC_cut), 'frequency'] - eps_p*DeltaNu + echelle_offset*DeltaNu) % DeltaNu).values
    l0_y = peaks.loc[(peaks.l == 0) & (peaks.AIC.values >= AIC_cut), 'frequency'].values
    # l=2 x and y values
    l2_x = ((peaks.loc[(peaks.l == 2) & (peaks.AIC.values >= AIC_cut), 'frequency'] - eps_p*DeltaNu + echelle_offset*DeltaNu) % DeltaNu).values
    l2_y = peaks.loc[(peaks.l == 2) & (peaks.AIC.values >= AIC_cut), 'frequency'].values
    # l=1/3 x and y values
    #'', np.isnan(peaks.l.values)
    l1_x = ((peaks.loc[((peaks.l == 1) | np.isnan(peaks.l.values)) & (peaks.AIC.values >= AIC_cut), 'frequency'] - eps_p*DeltaNu + echelle_offset*DeltaNu) % DeltaNu).values
    l1_y = peaks.loc[((peaks.l == 1) | np.isnan(peaks.l.values)) & (peaks.AIC.values >= AIC_cut), 'frequency'].values
    # l=3 x and y values
    l3_x = ((peaks.loc[(peaks.l == 3) & (peaks.AIC.values >= AIC_cut), 'frequency'] - eps_p*DeltaNu + echelle_offset*DeltaNu) % DeltaNu).values
    l3_y = peaks.loc[(peaks.l == 3) & (peaks.AIC.values >= AIC_cut), 'frequency'].values

    if st.sidebar.checkbox("Overplot theoretical frequencies"):

        alpha = st.sidebar.slider('α', min_value=0.0, max_value=5e-2, step=1e-4, value=orig_alpha, key=session.run_id, format="%.4f")
        d02 = st.sidebar.slider('d02 (μHz)', min_value=0.0, max_value=orig_DeltaNu*0.3, value=orig_d02, key=session.run_id, format="%.4f")

        # Estimate minimum and maximum radial orders
        min_n = np.floor(psd_bgr.frequency.min() / DeltaNu - eps_p).astype(int)
        n_max = (summary.numax.values / DeltaNu) - eps_p
        max_n = np.floor(psd_bgr.frequency.max() / DeltaNu - eps_p).astype(int)

        ridge_n = np.arange(min_n-1, max_n+3, 1)

        l0_ridge_freq = (ridge_n + eps_p + alpha/2 * (ridge_n - n_max)**2) * DeltaNu
        l2_ridge_freq = (ridge_n + eps_p + alpha/2 * (ridge_n - n_max)**2) * DeltaNu - d02

        ridge = (ridge_n*DeltaNu) % DeltaNu
        ridge[np.abs(ridge - DeltaNu) < 1e-6] = 0.0
        l0_ridge = (l0_ridge_freq - eps_p*DeltaNu + echelle_offset*DeltaNu) % DeltaNu
        l2_ridge = (l2_ridge_freq - eps_p*DeltaNu + echelle_offset*DeltaNu) % DeltaNu

        #l0_ridge = (eps_p*DeltaNu) % DeltaNu + (alpha/2 * (ridge_n - n_max)**2 * DeltaNu) % DeltaNu + echelle_offset# + 0.3
        #l2_ridge = (eps_p*DeltaNu) % DeltaNu + (alpha/2 * (ridge_n - n_max)**2 * DeltaNu) % DeltaNu - d02 % DeltaNu  + echelle_offset# + 0.3

        #'ridge', ridge

        ridge_source = ColumnDataSource(data = {
            "red_freq": [l0_ridge, l2_ridge],
            "frequency": [l0_ridge_freq, l2_ridge_freq],
            "l": ["l=0", "l=2"],
            "line_color": ["red", "green"]
            })

        p.multi_line(xs="red_freq", ys="frequency",
               line_width=2,
               line_color="line_color",
               line_dash=[5,5], alpha=0.4, source=ridge_source)
        #p.line(x="l2_ridge", y="l2_ridge_freq", line_width=2, color="green", line_dash=[5,5], alpha=0.4)



    # if st.sidebar.checkbox('\"Prettify Frequencies\"'):

    #     # Want middle of bin, not bin edge
    #     yn_tmp = yn + DeltaNu/2

    #     # Find middle of bin in y that is closes to frequency
    #     for i in range(len(l0_y)):
    #         l0_y[i] = app_helpers.find_nearest(yn_tmp, l0_y[i])
    #     for i in range(len(l2_y)):
    #         l2_y[i] = app_helpers.find_nearest(yn_tmp, l2_y[i])
    #     for i in range(len(l1_y)):
    #         l1_y[i] = app_helpers.find_nearest(yn_tmp, l1_y[i])
    #     for i in range(len(l3_y)):
    #         l3_y[i] = app_helpers.find_nearest(yn_tmp, l3_y[i])

    # Symbol size in echelle defined by mode amplitude, defaults to on
    if st.sidebar.checkbox("Set symbol size by mode amplitude", value=True):

        const = st.sidebar.number_input("Symbol size multiplier.", value=10.0, step=0.5)

        l0_size = const*peaks.loc[(peaks.l == 0) & (peaks.AIC.values >= AIC_cut), 'amplitude']
        l2_size = const*peaks.loc[(peaks.l == 2) & (peaks.AIC.values >= AIC_cut), 'amplitude']
        l1_size = const*peaks.loc[((peaks.l == 1) | np.isnan(peaks.l.values)) & (peaks.AIC.values >= AIC_cut), 'amplitude']
        l3_size = const*peaks.loc[(peaks.l == 3) & (peaks.AIC.values >= AIC_cut), 'amplitude']

    else:
        l0_size = 10.0
        l2_size = 10.0
        l1_size = 10.0
        l3_size = 10.0

    source = ColumnDataSource(data = {
            "red_freq": np.r_[l0_x, l2_x, l1_x, l3_x],
            "frequency": np.r_[l0_y, l2_y, l1_y, l3_y],
            "l": np.r_[np.array(["l=0"]*len(l0_x)),
                       np.array(["l=2"]*len(l2_x)),
                       np.array(["l=1"]*len(l1_x)),
                       np.array(["l=3"]*len(l3_x))],
            "size": np.r_[l0_size, l2_size, l1_size, l3_size]
    })




    p.scatter(x="red_freq", y="frequency", size="size", source=source,
              alpha = 0.5,
              legend = "l",
              marker=factor_mark("l", MARKERS, MODE_DEGREE),
              color=factor_cmap("l", COLOURS, MODE_DEGREE))

    p.add_tools(HoverTool(
                tooltips=[
                    ( 'frequency',   '$data_y μHz'            ),
                    ( 'mode ID', '@l'),
                ],
            ))

    p.legend.title = "Mode degree"
    p.legend.click_policy="hide"

    st.bokeh_chart(p)

def visualise_stretched_echelle(psd_bgr, peaks, summary, session):

    #st.write(psd_bgr)

    from sloscillations import frequencies, mixed_modes_utils
    from src.lib.rotation import rotation_utils

    orig_DPi1 = summary['DeltaPi1'].item()
    orig_numax = summary['numax'].item()
    orig_eps_p = summary['eps_p'].item()
    orig_alpha = summary['alpha'].item()
    #'', summary.columns
    orig_DeltaNu = summary['DeltaNu'].item()
    orig_coupling = summary['coupling'].item()
    orig_eps_g = summary['eps_g'].item()

    orig_d01 = 0.0
    orig_splitting = 0.0

    # Set up sliders
    reset = st.sidebar.button('Reset parameters')
    if reset:
        session.run_id += 1
    save = st.sidebar.button('Save parameters')
    if save:
        st.write('Doesn\'t do anything at the moment, but will do soon!')


    AIC_cut = st.sidebar.slider('AIC cut', min_value=-10.0, max_value=20., value=2.0, key=session.run_id)
    DeltaNu = st.sidebar.slider('Δν (μHz)', min_value=orig_DeltaNu*0.9, max_value=orig_DeltaNu*1.1, value=orig_DeltaNu, key=session.run_id)
    DPi1 = st.sidebar.number_input('ΔΠ₁ (s)', min_value=orig_DPi1*0.5, max_value=orig_DPi1*1.1, value=orig_DPi1, key=session.run_id)
    coupling = st.sidebar.slider('q', min_value=0.0, max_value=0.8, value=orig_coupling, key=session.run_id)
    eps_g = st.sidebar.slider('ε_g', min_value=-1.0, max_value=1.0, value=float(orig_eps_g), key=session.run_id)
    d01 = st.sidebar.slider('δν₀₁ (μHz)', min_value=-orig_DeltaNu/3, max_value=orig_DeltaNu/3, value=orig_d01, key=session.run_id)
    splitting = st.sidebar.slider('δν_rot (μHz)', min_value=0., max_value=1.2, value=orig_splitting, key=session.run_id)

    peaks['x'] = ((peaks['frequency'] % DeltaNu - summary['eps_p'].values) / summary['DeltaNu'].values) % 1
    #st.write(peaks.loc[peaks['l'] == 0, ])
    #st.write([(np.minimum(np.min(peaks.loc[peaks['l'] == 0, 'x']), np.min(peaks.loc[peaks['l'] == 2, 'x'])) - 0.05) % 1,
    #           (np.maximum(np.max(peaks.loc[peaks['l'] == 0, 'x']), np.max(peaks.loc[peaks['l'] == 2, 'x'])) + 0.05) % 1]
    #)

    l1_peaks = rotation_utils.prepare_l1_peaks(peaks, summary, AIC_cut)

    # Create frequencies instance
    freqs = frequencies.Frequencies(frequency=psd_bgr.frequency.values,
                            numax=orig_numax,
                            delta_nu=DeltaNu,
                            epsilon_p=orig_eps_p,
                            alpha=orig_alpha)

    params = {'calc_l0': True, # Compute radial mode properties
            'calc_l2': True, # Compute l=2 mode properties
            'calc_l3': False, # Don't need to calculate l=3 theoretical freqs
            'calc_nom_l1': True, # Compute nominal l=1 p-mode properties
            'calc_mixed': True, # Don't compute mixed modes (as not needed)
            'calc_rot': True, # Don't compute rotation
            'd01': d01,
            'DPi1': DPi1,
            'coupling': coupling,
            'eps_g': eps_g,
            'split_core': splitting,
            'split_env': 0.0,
            'l': 1, # Mixed modes are dipole mixed modes
            }

    # Make computation - in our case this is for the computation of zeta
    freqs(params)
    freqs.generate_tau_values()

    new_peaks_tau = mixed_modes_utils.peaks_stretched_period(l1_peaks.frequency.values,
                                                                psd_bgr.frequency.values,
                                                                freqs.tau)
    new_peaks_zeta = mixed_modes_utils.peaks_stretched_period(l1_peaks.frequency.values,
                                                                psd_bgr.frequency.values,
                                                                freqs.zeta)
    #st.write(freqs.l1_mixed_freqs_p1, freqs.l1_mixed_freqs_n1)
    # Plot stretched echelle
    model_freqs = np.c_[freqs.l1_mixed_freqs, freqs.l1_mixed_freqs_p1, freqs.l1_mixed_freqs_n1]
    model_tau = np.c_[freqs.l1_mixed_tau, freqs.l1_mixed_tau_p1, freqs.l1_mixed_tau_n1]
    real_tau = new_peaks_tau
    real_freqs = l1_peaks.frequency.values
    heights = l1_peaks.amplitude*10
    y_real = (real_tau - DPi1*(1/2 + (-eps_g)))  % DPi1 - DPi1/2
    # The shift is already accounted for in the calculation of the theoretical frequencies, so don't add it in here
    y_theo = (model_tau[:,0] - DPi1/2) % DPi1  - DPi1/2
    y_theo_p1 = (model_tau[:,1] - DPi1/2) % DPi1  - DPi1/2
    y_theo_n1 = (model_tau[:,2] - DPi1/2) % DPi1  - DPi1/2


    # plt.scatter(y_real, real_freqs, s=heights)
    # plt.scatter(y_theo, model_freqs[:,0], marker='x')
    # plt.scatter(y_theo_p1, model_freqs[:,1], marker='x')
    # plt.scatter(y_theo_n1, model_freqs[:,2], marker='x')

    p = figure(
        title="",
        x_axis_label="Stretched Period modulo period spacing (s)",
        y_axis_label="Frequency (μHz)",
        match_aspect=True,
        x_range=(-DPi1/2, DPi1/2),
        y_range=(psd_bgr.frequency.min(),
                    psd_bgr.frequency.max()),
        plot_height=400,
        plot_width=600)


    p.triangle(y_real, real_freqs, size=10*l1_peaks.amplitude,
                color='blue', alpha=0.5)

    p.circle(y_theo,
                model_freqs[:,0],
                #angle=np.pi, size=10,
                #alpha=0.5,
                color="red", legend_label=r'l=1, m=0')

    if splitting > 0.0:

        p.circle(y_theo_n1,
                    model_freqs[:,2],
                    #angle=-np.pi/2, size=10, alpha=0.5,
                    color="cyan",
                    legend_label=r'l=1, m=-1')
        p.circle(y_theo_p1,
                    model_freqs[:,1],
                    #angle=np.pi/2, size=10, alpha=0.5,
                    color="purple", legend_label=r'l=1, m=+1')

    p.legend.click_policy = 'hide'
    st.bokeh_chart(p)

def visualise_reggae(psd_bgr, peaks, summary, session):

    #st.write(psd_bgr)

    from sloscillations import frequencies, mixed_modes_utils
    from src.lib.rotation import rotation_utils

    orig_numax = summary['numax'].item()
    orig_DeltaNu = summary['DeltaNu'].item()
    orig_eps_p = summary['eps_p'].item()
    orig_alpha = summary['alpha'].item()
    orig_d02 = summary['dNu02'].item()

    orig_d01 = 0.0
    #'', summary.columns
    orig_DPi1 = summary['DeltaPi1'].item()
    orig_coupling = summary['coupling'].item()
    orig_eps_g = summary['eps_g'].item()
    orig_splitting = 0.0

    # Set up sliders
    reset = st.sidebar.button('Reset parameters')
    if reset:
        session.run_id += 1
    save = st.sidebar.button('Save parameters')
    if save:
        st.write('Doesn\'t do anything at the moment, but will do soon!')

    AIC_cut = st.sidebar.slider('AIC cut', min_value=-100.0, max_value=100., value=2.)

    overlay_modes = st.sidebar.checkbox("Overlay theoretical mixed mode frequencies as lines")

    mixed_mode_config = st.sidebar.selectbox('Mixed mode configuration',
                                               ['Triplet', 'Doublet', 'Singlet']
                                            )
    nmax = orig_numax / orig_DeltaNu - orig_eps_p

    n_current = st.sidebar.number_input("n", min_value = 0, max_value=int(nmax) + 10,
                                        value=int(nmax), key=session.run_id)

    DeltaNu = st.sidebar.slider('Δν (μHz)', min_value=orig_DeltaNu*0.9, max_value=orig_DeltaNu*1.1, value=orig_DeltaNu, key=session.run_id)
    eps_p = st.sidebar.slider('ε_p', min_value=0.0, max_value=2.0, value = orig_eps_p, key = session.run_id)
    alpha = st.sidebar.slider('α', min_value=0.0, max_value=0.05, value=orig_alpha, key=session.run_id)
    d02 = st.sidebar.slider('δν₀₂ (μHz)', min_value=0.0, max_value=orig_d02*1.5, value=orig_d02, key=session.run_id)
    #d01 = st.sidebar.slider('δν₀₁ (μHz)', min_value=-orig_DeltaNu/3, max_value=orig_DeltaNu/3, value=orig_d01, key=session.run_id)
    d01 = st.sidebar.number_input('δν₀₁ (μHz)', min_value=-orig_DeltaNu/3, max_value=orig_DeltaNu/3, value=orig_d01, key=session.run_id)


    #DPi1 = st.sidebar.slider('ΔΠ₁ (s)', min_value=orig_DPi1*0.9, max_value=orig_DPi1*1.1, value=orig_DPi1, key=session.run_id)
    DPi1 = st.sidebar.number_input('ΔΠ₁ (s)', min_value=0.0, max_value=400.0, value=orig_DPi1, key=session.run_id)
    #coupling = st.sidebar.slider('q', min_value=0.0, max_value=0.8, value=orig_coupling, key=session.run_id)
    coupling = st.sidebar.number_input('q', min_value=0.0, max_value=0.8, value=orig_coupling, key=session.run_id)
    #eps_g = st.sidebar.slider('ε', min_value=-1.0, max_value=1.0, value=float(orig_eps_g), key=session.run_id)
    eps_g = st.sidebar.number_input('ε_g', min_value=-1.0, max_value=1.0, value=float(orig_eps_g), key=session.run_id)
    # splitting = st.sidebar.slider('δνₛ (μHz)', min_value=0., max_value=500.0/1e3, value=orig_splitting, key=session.run_id)
    splitting = st.sidebar.number_input('δν_rot (μHz)', min_value=0., max_value=1.0, value=orig_splitting, key=session.run_id)

    # Create frequencies instance
    freqs = frequencies.Frequencies(frequency=psd_bgr.frequency.values,
                            numax=orig_numax,
                            delta_nu=DeltaNu,
                            epsilon_p=eps_p,
                            alpha=alpha,
                            radial_order_range = [-4, 4])

    params = {'calc_l0': True, # Compute radial mode properties
            'calc_l2': True, # Compute l=2 mode properties
            'calc_l3': False, # Don't need to calculate l=3 theoretical freqs
            'calc_nom_l1': True, # Compute nominal l=1 p-mode properties
            'calc_mixed': True, # Don't compute mixed modes (as not needed)
            'calc_rot': True, # Don't compute rotation
            'd02': d02,
            'd01': d01,
            'DPi1': DPi1,
            'coupling': coupling,
            'eps_g': eps_g,
            'split_core': splitting,
            'split_env': 0.0,
            'l': 1, # Mixed modes are dipole mixed modes
            }

    # Make computation - in our case this is for the computation of zeta
    freqs(params)

    freqs.generate_tau_values()

    #new_peaks_tau = mixed_modes_utils.peaks_stretched_period(l1_peaks.frequency.values,
    #                                                            psd_bgr.frequency.values,
    #                                                            freqs.tau)
    #new_peaks_zeta = mixed_modes_utils.peaks_stretched_period(l1_peaks.frequency.values,
    #                                                            psd_bgr.frequency.values,
    #                                                            freqs.zeta)

    n_current = int(n_current - freqs.n.min())
    #n_current = (freqs.n[np.argmin(np.abs(nmax - freqs.n))] - freqs.n.min()).astype(int)

    #st.write(nmax, freqs.n, n_start)
    lower_lim = freqs.l0_freqs[n_current] - freqs.d02 * 1.5
    upper_lim = freqs.l0_freqs[n_current+1] - freqs.d02 * 0.5
    cond = (psd_bgr.frequency >= lower_lim) & (psd_bgr.frequency <= upper_lim)

    p = figure(
        title="",
        y_axis_label="Signal-to-Noise Spectrum",
        x_axis_label="Frequency (μHz)",
        match_aspect=True,
        x_range=(lower_lim,
                 upper_lim),
        y_range=(0.0,
                    psd_bgr.power[cond].max()),
        plot_height=500,
        plot_width=800)



    p.line(psd_bgr.frequency, psd_bgr.power, color="black")
    model = app_helpers.construct_peaksmodel(psd_bgr, peaks.loc[peaks.AIC.values > AIC_cut, ])
    p.line(psd_bgr.frequency, model, color="orange", line_width=3)

    if st.sidebar.checkbox("Mode identification"):
        p.circle(peaks.loc[(peaks.l == 0) & (peaks.AIC.values > AIC_cut), 'frequency'], len(peaks.loc[(peaks.l == 0) & (peaks.AIC.values > AIC_cut), 'frequency'])*[psd_bgr.power.max()/2],
                size=10, color="red", alpha=0.5, legend_label='l=0')
        p.square(peaks.loc[(peaks.l == 2) & (peaks.AIC.values > AIC_cut), 'frequency'],
                len(peaks.loc[(peaks.l == 2) & (peaks.AIC.values > AIC_cut), 'frequency'])*[psd_bgr.power.max()/2],
                size=10, color="green", alpha=0.5, legend_label='l=2')
        p.hex(peaks.loc[(peaks.l == 3) & (peaks.AIC.values > AIC_cut), 'frequency'],
                len(peaks.loc[(peaks.l == 3) & (peaks.AIC.values > AIC_cut), 'frequency'])*[psd_bgr.power.max()/2],
                size=10, color="orange", alpha=0.5, legend_label='l=3')


    l0_freqs = freqs.l0_freqs[n_current:n_current+2]
    l2_freqs = freqs.l2_freqs[n_current:n_current+2]
    l1_freqs = freqs.l1_nom_freqs[n_current:n_current+2]

    for i in range(len(l0_freqs)):
        p.line([l0_freqs[i]]*2, [0, 300], line_color='red', alpha=0.5,
                line_width=3, legend_label=r'l=0')
        p.line([l2_freqs[i]]*2, [0, 300], line_color='green', alpha=0.5,
                line_width=3, legend_label=r'l=2')
        p.line([l1_freqs[i]]*2, [0, 300], line_color='blue',
                line_dash="dashed", alpha=0.5,
                line_width=3, legend_label=r'nominal dipole mode')

    idx = (freqs.l1_np == int(n_current + (freqs.n.min())))

    l1_mixed_freqs_m0 = freqs.l1_mixed_freqs[idx]
    #st.write(l1_mixed_freqs_m0)
    height = np.percentile(psd_bgr.power[cond], [99.5])
    #'', len(l1_mixed_freqs_m0)
    #st.write([height]*len(l1_mixed_freqs_m0))

    if (mixed_mode_config == "Triplet") or (mixed_mode_config == "Singlet"):
        p.triangle(l1_mixed_freqs_m0, np.ones(len(l1_mixed_freqs_m0))*height,
                size=10, color='blue', angle=-np.pi)

    # for i in range(len(l1_mixed_freqs_m0)):
    #     vline_l1_mixed = Span(location=l1_mixed_freqs_m0[i],
    #                           line_color='blue',
    #                           dimension='height',
    #                           line_width=2,
    #                           line_dash="dashed")
    #     p.renderers.extend([vline_l1_mixed])


    if (splitting > 0.0) and ((mixed_mode_config == "Triplet") or (mixed_mode_config == "Doublet")):
        l1_mixed_freqs_mp1 = freqs.l1_mixed_freqs_p1[idx]
        l1_mixed_freqs_mn1 = freqs.l1_mixed_freqs_n1[idx]
        p.triangle(l1_mixed_freqs_mn1, np.ones(len(l1_mixed_freqs_mn1))*height * 0.9,
                   size=10, color='cyan', angle=-np.pi/2)
        p.triangle(l1_mixed_freqs_mp1, np.ones(len(l1_mixed_freqs_mn1))*height * 0.9,
                   size=10, color='purple', angle=np.pi/2)

        if mixed_mode_config == "Triplet":
            for i in range(len(l1_mixed_freqs_mn1)):
                lines = p.line([l1_mixed_freqs_mn1[i], l1_mixed_freqs_m0[i]],
                            [height * 0.9, height], color='black',
                            line_dash="2 2", alpha=0.5)
                lines.level = "underlay"
                lines = p.line([l1_mixed_freqs_m0[i], l1_mixed_freqs_mp1[i]],
                            [height * 0.9, height][::-1], color='black',
                            line_dash="2 2", alpha=0.5)
                lines.level = "underlay"

    if overlay_modes & ((mixed_mode_config == "Triplet") or (mixed_mode_config == "Singlet")):
        for i in range(len(l1_mixed_freqs_m0)):
            p.line([l1_mixed_freqs_m0[i], l1_mixed_freqs_m0[i]],
                           [0, height], color='blue', alpha=0.5)

    if overlay_modes & (splitting > 0.0) & ((mixed_mode_config == "Triplet") or (mixed_mode_config == "Doublet")):
        for i in range(len(l1_mixed_freqs_mn1)):
            p.line([l1_mixed_freqs_mn1[i], l1_mixed_freqs_mn1[i]],
                        [0, 0.9*height], color='cyan', alpha=0.5)
            p.line([l1_mixed_freqs_mp1[i], l1_mixed_freqs_mp1[i]],
                        [0, 0.9*height], color='purple', alpha=0.5)

    p.legend.click_policy = 'hide'
    st.bokeh_chart(p)


def main():

    #st.title("TACO Explorer App")

    #st.markdown("Explore the stars you have analysed to your hearts content!!!")

    session = SessionState.get(run_id=0)
    # File uploader
    #if uploaded_file is not None:
    KICs = find_stars()
    st.sidebar.header("Please select a star to analyse")
    selected_KIC = st.sidebar.selectbox(
                    "",
                    ["", *KICs]
    )
    if selected_KIC != "":

        #'You selected ', selected_KIC
        st.sidebar.header("Choose a page")
        page = st.sidebar.selectbox("",
                                    ["", "Timeseries", "Background Fit",
                                     "MLE Fit", "Frequency Echelle",
                                     "Stretched Period Echelle",
                                     "Mode matching"])
        summary = load_summary(selected_KIC)
        if page == "Timeseries":
            st.header("Timeseries Data")
            filtered_ts = load_ts(selected_KIC, filtered=True)
            unfiltered_ts = load_ts(selected_KIC, filtered=False)
            visualise_timeseries(filtered_ts, unfiltered_ts, summary)

        elif page == "Background Fit":
            st.header("Background Fit Explorer")

            psd = load_psd(selected_KIC, background_removed=False)
            visualise_psd(psd, summary)

        elif page == "MLE Fit":
            st.header("MLE Fit Explorer")
            st.sidebar.header("MLE Fit Settings")
            psd_bgr = load_psd(selected_KIC, background_removed=True)
            peaksMLE = load_peaks(selected_KIC)
            peakFind = load_peaks(selected_KIC, peakFind=True)
            visualise_pds_bgr(psd_bgr, peaksMLE, summary, peakFind)
            if st.checkbox('Show MLE fit parameters'):
                st.write('MLE parameters')
                st.write(peaksMLE)
                st.write("peakFind parameters")
                st.write(peakFind)

        elif page == "Frequency Echelle":
            st.header("Frequency Echelle Explorer")
            st.sidebar.header("Echelle Parameters")
            psd_bgr = load_psd(selected_KIC, background_removed=True)
            peaksMLE = load_peaks(selected_KIC)
            visualise_echelle(psd_bgr, peaksMLE, summary, session)
            if st.checkbox('Show MLE fit parameters'):
                st.write('MLE parameters')
                st.write(peaksMLE)

        elif page == "Stretched Period Echelle":
            st.header("Stretched Period Echelle Explorer")
            st.sidebar.header("Echelle Parameters")
            psd_bgr = load_psd(selected_KIC, background_removed=True)
            peaksMLE = load_peaks(selected_KIC)
            visualise_stretched_echelle(psd_bgr, peaksMLE, summary, session)

        elif page == "Mode matching":
            #st.header("Reggae!")
            st.sidebar.header("Asteroseismic Parameters")
            psd_bgr = load_psd(selected_KIC, background_removed=True)
            peaksMLE = load_peaks(selected_KIC)
            visualise_reggae(psd_bgr, peaksMLE, summary, session)

if __name__=="__main__":
    main()
