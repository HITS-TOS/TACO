#!/bin/python

import pandas as pd
import taco

for dir in list:
    ts = pd.read_csv('001296068/raw.dat', comment = '#', header = None, sep = '\s+')

    print("filter lightcurve")
    filter_lightcurve(ts)

    print("compute pds")
    calc_pds(ts)

    print("numax estimate")
    numax_estimate()

    print("background fit")
    background_fit(bins=300)

    print("background summary")
    background_summary()

    print("find peaks")
    peakFind(snr=1.1, prob=0.0001, minAIC=2)

    print("MLE optimisation for peaks")
    peaksMLE(minAIC=2)

    print("mode ID 02")
    peakBagModeId02()

    print("find peaks")
    peakFind(snr=1.1, prob=0.0001, minAIC=2, removel02=TRUE)

    print("MLE optimisation for peaks")
    peaksMLE(minAIC=2, removel02=TRUE, init=peaksMLE.csv, mixedpeaks=mixedpeaks.csv)

    print("MLE optimisation for peaks")
    peaksMLE(minAIC=2, finalfit=TRUE)
