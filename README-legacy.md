Automatic Peak-bagging utilities
================================

These instructions were written for the TACO version of october 2021 and was run with python-version 3.8 and R-version 4.1.1

How to install TACO on Windows
==============================
Some steps can be skipped if the software has already been installed previously

1. Download TACO
2. Download Anaconda -> Python (specify adding anaconda to PATH-variable, still possible to add manually later)
3. Install the following packages: 
    - `argparse`
    - `numpy`
    - `pandas`
    - `scipy`
    - `sklearn`
    - `emcee`
    - `lightkurve`
    - `loguru`
    Two methods (in Anaconda Prompt!): `conda install packagename` OR `pip install packagename`
4. Install R (https://cran.r-project.org/)
5. Open the RGui as Administrator and install the following packages by executing `install.packages("packagename")` in the R console:
    - `argparser`
    - `ifultools`
    - `numDeriv`
    - `lomb`
    - `tidyverse`
    - `modelr`
6. Check that all packages are in the R-library (in directory C:\Program Files\R\R-X.X.X\library)
7. Install Git (https://git-scm.com/)
8. Update the path variable: search for 'Environment variables' and select 'Edit environment variables' (open as Admin, otherwise no right to edit system variables!), select `Add` next to system variables.
    Paths to be added (if not already present):
      - `C:\Program Files\Git\cmd`
      - `C:\Program Files\Git\bin\bash.exe`
      - `C:\Program Files\R\R-X.X.X\bin`
      - `C:\Users\Your username\Anaconda3\Scripts`
      - `C:\Users\Your username\Anaconda3\Library\bin`
      - `C:\Users\Your username\Anaconda3`



(Optional) Download Visual Studio Code (https://code.visualstudio.com/)
              No special recommandations during the installation  
           Once installed open VSC and go to Extensions (ctrl+shift+x) and add following extensions if not already installed:
              - `Python`
              - `R`
              - `Excel Viewer` (not required but useful when interpreting the data from .csv-files)




How to install TACO on Ubuntu
=============================

1. Install R (https://cran.r-project.org/) and install following packages by executing `install.packages("packagename")` in the R console (to open console run `R` in terminal):
    - `argparser`
    - `ifultools`
    - `numDeriv`
    - `lomb`
    - `tidyverse`
    - `modelr`
2. Install Anaconda (for Linux installer: https://www.anaconda.com/products/individual#linux ; for documentation: https://docs.anaconda.com/anaconda/install/linux/)
3. Install the following packages: 
    - `argparse`
    - `numpy`
    - `pandas`
    - `scipy`
    - `sklearn`
    - `emcee`
    - `lightkurve`
    - `loguru`
    Two methods (in terminal): `conda install packagename` OR `pip install packagename`
3. Download and unpack TACO



Run TACO on Windows 
===================

Use the shellscript TACO.sh to analyse a star. 

1. Create a specific folder for each new star in the directory `data`.
    Copy the `.dat`-file of the desired star into the folder of the star and rename this file to `raw.dat`.
    All output for this star will be inside this folder.
    Copy `./TACO.sh` into the directory where the dataset of the star is located (you might have to change all python3 to python).

2.  a. Open VSC and start a Git Bash in the terminal window:
      In terminal window (bottom window usually), click on down arrow next to plus sign and select `Git Bash`
      In the bash terminal: go to the folder of your star (where `raw.dat` is located) using `cd xxx/xxx/xxx`
    

    b. Open a Git Bash terminal (right click) in the directory of the star. 

3. Run the command `./TACO.sh` in the Git bash terminal. You might have to run the following in the bash-terminal: `conda init bash` 


Note: By adding or removing additional lines and/or input parameters in `./TACO.sh`, one can finetune the analysis for each dataset.


Run TACO on Ubuntu 
==================

Use the shellscript TACO.sh to analyse a star. 


1. Open a terminal
2. Create a specific folder for each new star in the directory `data`.
    Copy the `.dat`-file of the desired star into the folder of the star and rename this file to `raw.dat`.
    All output for this star will be inside this folder.
    Copy `./TACO.sh` into the directory where the dataset of the star is located (you might have to change all python3 to python).
3. In the terminal: go to the folder of your star (where `raw.dat` is located) using `cd`
   All output for this star will be inside this folder.
4. Copy `./TACO.sh` into the directory where the dataset of the star is located (you might have to change all python3 to python).

Note: By adding or removing additional lines and/or input parameters in `./TACO.sh`, one can finetune the analysis for each dataset.

Description of all the scripts
==============================


Filter light-curves (filter_lightcurve.R)
-----------------------------------------

Filter the lightcurves with a high-pass filter (e.g. 40 days) in the form of a triangular smooth of the provided width.
Additionally it interpolates the single-point missing values.
It also calculates the mean and variance of the light curve and saves them in the summary file.
I assume the light curve is as provided by the APOKASC group which is an ascii table with three columns representing the time (in days), the observed ppm and the error (which is ignored from now on).

```
usage: Rscript filter_lightcurve.R [--] [--help] [--opts OPTS] [--tseries TSERIES] [--output OUTPUT] [--summary SUMMARY] [--width WIDTH] [--remove_gaps REMOVE_GAPS]

Make a triangular filter of the time-series.


flags:
  -h, --help			        show this help message and exit

optional arguments:
  --opts OPTS			      RDS file containing argument values
  --tseries TSERIES			  File name of the time-series in tsv format with the first column giving the time and the second the flux. [default: raw.dat]
  --output OUTPUT			  File name on which to save the filtered time-series [default: filtered.csv]
  --summary SUMMARY			  Summary file to save some time-series statistics [default: summary.csv]
  --width WIDTH			      Width of the triangular filter [default: 40]
  --remove_gaps REMOVE_GAPS   Remove gaps greater than this value (in days). If set to -1, do not remove gaps. [default: -1]
```

Lomb-Scargle periodogram (pds.py)
--------------------------------

Calculate the Lomb-Scargle periodogram of a time series normalized using the spectral window function as described by [Kallinger et al. 2014](http://dx.doi.org/10.1051/0004-6361/201424313). The input lightcurve is assumed to have time stamps given in days, the flux in ppm and the periodogram is given in units of micro-Hertz. The functionality for computing the periodogram comes from the `lightkurve` package.

```
usage: python pds.py [--] [--help] [--tseries TSERIES] [--output OUTPUT] [--oversampled_output OVERSAMPLED_OUTPUT] [--ofac OFAC] [--summary SUMMARY]

Make Lomb-Scargle periodogram. We assume the input time is given in days and give the periodogram frequencies in micro-Hertz.


flags:
  -h, --help			        show this help message and exit

optional arguments:
  --tseries TSERIES			              File name of the time-series in csv format with two columns: `time` and `flux`. [default: filtered.csv]
  --output OUTPUT			              File name of the saved periodogram. [default: pds.csv]
  --oversampled_output OVERSAMPLED_OUTPUT File name on which to save the oversampled periodogram. [default: ofac_pds.csv]
  --ofac OFAC                             Oversampling factor to use in oversampled periodogram [default: 2]
  --summary SUMMARY			              Summary file to save some periodogram statistics [default: summary.csv]
```

Numax estimation (numax_estimate.R)
--------------------------------------

Make three independent possible estimates of the frequency of maximum oscillation power (nu_max) and compares them to get a final guess. The estimates are based on the following

- Variance of the time series: It was shown by [Hekker et al. 2012](http://dx.doi.org/10.1051/0004-6361/201219328) that the variance of the time series has an exponential relation with nu_max.
- Peak detections with a tree-map of the mexican-hat continuous wavelet transform: Using method by [Du et al. 2006](http://bioinformatics.oxfordjournals.org/content/22/17/2059) to detect significant peaks in a data set.
- Morlet wavelet transfroms: The maximum of a Morlet wavelet transform usually occurs near numax.

Updates the summary file with new columns: `numax_var`, `numax_CWTTree`, `numax_Morlet` and `numax0` which are the 3 different nu_max estimations and the final guess respectively.
It also adds a column `numax0_flag`. If true, the three estimates are different. 

```
usage: Rscript numax_estimate.R [--] [--help] [--opts OPTS] [--pds PDS] [--summary SUMMARY]

Make a rough numax estimation using the periodogram and the variance of the time-series.
We assume that the summary file contains the variance of the time-series and the Nyquist
frequency of the periodogram in columns named 'var' and 'nuNyq', respectively.


flags:
  -h, --help			    show this help message and exit

optional arguments:
  -x, --opts OPTS			RDS file containing argument values
  -p, --pds PDS			    File name of the PDS in csv format with headers 'frequency' and 'power' [default: pds.csv]
  -s, --summary SUMMARY		Summary file in csv format with columns named 'var' and 'nuNyq'. [default: summary.csv]
```

Background fitting (background_fit.py)
--------------------------------------

Use an MCMC algorithm (emcee) to make a background estimation for the power-spectrum density normalized according with the above prescription. Note that this script uses some initial guesses that critically depend on the normalization of the power-spectrum density so, if you use a different normalization, this script might need some changes before it's able to work. You need to have installed the `emcee` python package (version 3 or greater).

The background model is similar to the one proposed by [Kallinger et al. 2015](http://www.aanda.org/articles/aa/abs/2014/10/aa24313-14/aa24313-14.html) albeit with some differences.
We use by default a unbinned version of the PDS to do the fit but this can be changed with the command-line arguments.

```
usage: background_fit.py [-h] [--pds PDS] [--summary SUMMARY]                                                                          
                         [--posterior POSTERIOR]                                                     
                         [--logfile LOGFILE] [--nwalkers NWALKERS]                                                                     
                         [--nwarmup NWARMUP] [--minsteps MINSTEPS]                                                                     
                         [--maxsteps MAXSTEPS] [--bins BINS]                                                                           

Make an MCMC background fit.

optional arguments:
  -h, --help            show this help message and exit
  --pds PDS             A csv file with the periodogram on which to make the
                        fit
  --summary SUMMARY     A csv file with a 'numax0' column giving an initial
                        numax estimation
  --save_posteriors      A flag which when given saves out the posterior distributions
                        (i.e. MCMC chains) to the file fiven in the `--posterior` flag.
                        Otherwise if not given then they are not saved.
  --posterior POSTERIOR
                        Destination on which to write the MCMC posterior
                        (collapsed chains)
  --logfile LOGFILE     Location on which to write the MCMC log
  --nwalkers NWALKERS   Number of walkers (chains) for the MCMC fit
  --nwarmup NWARMUP     Number of steps for the MCMC warmup
  --minsteps MINSTEPS   Minimum number of steps for the MCMC estimation
  --maxsteps MAXSTEPS   Maximum number of steps for the whole MCMC run
  --bins BINS           The PDS will be binned by this number of bins for the
                        MCMC. Setting it to 1 will not bin it.
```

By default TACO does *NOT* save out the chains/posterior distributions from the background fit (since the file is ~50-100 MB per star), only the summary statistics. If you explicitly want to save the MCMC chains then please add the `--save_posteriors` flag when running the background fit.


Peak-detections (peakFind.R)
---------------------------

Takes the background-removed power-spectrum density and identifies the relevant oscillations. The oscillations are all modelled as lorentzians. This script estimates their position, width and height which can be later used to construct a prior in a more precise analysis. You need to have installed the R packages: `wmtsa`, `argparser` and `dplyr`. The auxiliary file `peakFind_lib.R` must be present in the same directory as `peakFind.R`.

```
usage: Rscript peakFind.R [--] [--help] [--opts OPTS] [--pds PDS] [--ofac_pds OFAC_PDS][--summary SUMMARY] [--output OUTPUT] [--snr SNR] [--prob PROB] [--maxlwd MAXLWD] [--removel02 REMOVEL02] [--peaks PEAKS] [--mixedpeaks MIXEDPEAKS] [--minAIC MINAIC] [--navg NAVG]

Peak-finding in a background-removed PDS. We only look for peaks in the region numax +- 4*sigmaEnv.


flags:
  -h, --help                    show this help message and exit

optional arguments:
  -x, --opts OPTS               RDS file containing argument values
  -p, --pds PDS                 File name of the background-removed PDS in csv format with headers 'frequency' and 'power' [default: pds_bgr.csv]
  --ofac_pds OFAC_PDS           File name of the background-removed PDS in csv format with headers 'frequency' and 'power' [default: ofac_pds_bgr.csv]
  -s, --summary SUMMARY         File name of the summary file in csv format with at least the headers 'numax' and 'sigmaEnv' [default: summary.csv]
  -o, --output OUTPUT           File name on which to save the results [default: peaks.csv]
  --snr SNR                     Minimum signal-to-noise ratio (on CWT space) for resolved peaks [default: 2]
  --prob PROB                   Minimum (frequentist) probability threshold for unresolved peaks [default: 1e-04]
  -m, --maxlwd MAXLWD           Maximum search linewidth for resolved peaks in the CWT search [default: 1]
  --removel02 REMOVEL02         Whether or not the l02 peaks should be divided out before running the CWT search [default: FALSE]
  --peaks PEAKS               File name of the identified peaks. It must contain the l=0,2 modes already identified [default: peaksMLE.csv]
  --mixedpeaks MIXEDPEAKS     File name of the csv file containing the mixed mode peaks from peak finding [default: mixedpeaks.csv]
  --minAIC MINAIC             Minimum AIC value for a peak to have in order to be considered significant [default: 2]
  --navg NAVG                 Number of power spectra averaged to create current power spectrum [default: 1]
```

Peak optimisations using MLE (peaksMLE.R)
-----------------------------------------

Makes an MLE parameter optimisation taking the peaks found previously with `peakFind.R`

```
usage: Rscript peaksMLE.R [--] [--help] [--opts OPTS] [--pds PDS] [--init INIT] [--summary SUMMARY] [--mixedpeaks MIXEDPEAKS] [--mixedoutput MIXEDOUTPUT] [--output OUTPUT] [--maxlwd MAXLWD] [--minAIC MINAIC] [--removel02 REMOVEL02] [--finalfit FINALFIT] [--navg NAVG]

MLE estimation of the peaks found in a background-removed PDS.


flags:
  -h, --help                    show this help message and exit

optional arguments:
  -x, --opts OPTS               RDS file containing argument values
  -p, --pds PDS                 File name of the background-removed PDS in csv format with headers 'frequency' and 'power' [default: pds_bgr.csv]
  -i, --init INIT               File name of the intial guesses [default: peaks.csv]
  -s, --summary SUMMARY         File name of the csv file with summary values. It must contain numax and sigmaEnv [default: summary.csv]
  --mixedpeaks MIXEDPEAKS       File name of the csv file containing the mixed mode peaks from peak finding [default: mixedpeaks.csv]
  --mixedoutput MIXEDOUTPUT     File name of the csv file to save the mixed mode MLE results to [default: mixedpeaksMLE.csv]
  -o, --output OUTPUT           File name on which to save the results [default: peaksMLE.csv]
  -m, --maxlwd MAXLWD           Maximum linewidth for peaks in the MLE optimization [default: 0.5]
  --minAIC MINAIC               Minimum AIC value for a peak to have in order to be considered significant [default: 2]
  --removel02 REMOVEL02         Whether or not the l02 peaks should be divided out before running the CWT search [default: FALSE]
  --finalfit FINALFIT           Whether or not this is the final MLE optimisation [default: FALSE]
  --navg NAVG                   Number of power spectra averaged to create current power spectrum [default: 1]
```

Spherical degree identification (peakBagModeId02.R)
-------------------------------------------------

Takes the oscillations found and tags them according to their spherical degree. It also calculate the large frequency separation Δν

```
usage: Rscript peakBagModeId02.R [--] [--help] [--opts OPTS] [--peaks PEAKS] [--summary SUMMARY] [--pds PDS] [--output OUTPUT]      
                                                                                                                                  
Mode identification for the peaks found.


flags:
  -h, --help                    show this help message and exit

optional arguments:
  -x, --opts OPTS                       RDS file containing argument values
  -p, --peaks PEAKS                     File name of the identified peaks [default: peaksMLE.csv]
  -s, --summary SUMMARY                 File name of the csv file with summary values. It must contain numax and sigmaEnv [default: summary.csv]
  --pds PDS                     File name of the csv file with the PDS It must contain the columns 'frequency' and 'power [default: pds_bgr.csv]
  -o, --output OUTPUT                   File name of the csv file on which to write the results [default: peaksMLE.csv]
```

Period spacing and coupling determination (peakBagPeriodSpacing.py)
--------------------------------------------------------------------

Computes the period spacing and coupling using the power spectrum of the stretched power spectrum from Vrard et al. (2016). This will also require the `sloscillations` python package for the fast mixed mode calculations.

```
usage: python peakBagPeriodSpacing.py [--peaks PEAKS] [--summary SUMMARY] [--pds PDS] [--full_pds FULLPDS] [--maxiters MAXITERS] [--niters NITERS] [--dpi_only DPIONLY] [--ncores NCORES]

Period spacing and coupling factor calculation

flags:
  -h, --help                    show this help message and exit

optional arguments:
  --peaks PEAKS           File name of the identified peaks [default: peaksMLE.csv]
  --summary SUMMARY       File name of the csv file with summary values (it must contain the parameters numax and sigmaEnv) [default: summary.csv]
  --pds PDS               File name of the csv file with the background-subtracted PDS It must contain the columns "frequency" and "power" [default: pds_bgr.csv]
  --full_pds FULLPDS      File name of the csv file with the original PDS It must contain the columns "frequency" and "power" [default: pds.csv]
  --maxiters MAXITERS     Maximum number of iterations to use in the self-consistent calculation of ΔΠ1 [default: 10]
  --niters NITERS         Number of iterations to repeat the calculation of ΔΠ1, q and ε_g [default: 5]
  --dpi_only DPIONLY      Only infer the period spacing and don't calculate tau or q [default: FALSE]
  --ncores NCORES         Number of cores to use for a parallel calculation [default: 1]
```

Rotational splitting determination (peakBagRotation.py)
--------------------------------------------------------------------

Computes the rotational splitting using the echelle of the stretched power spectrum roughly following Gehan (2018). This will also require the `sloscillations` python package for the fast mixed mode calculations.

```
usage: python peakBagRotation.py [--peaks PEAKS] [--summary SUMMARY] [--pds PDS] [--full_pds FULLPDS] [--output OUTPUT] [--ncores NCORES] [--cut CUT]

Period spacing and coupling factor calculation

flags:
  -h, --help                    show this help message and exit

optional arguments:
  --peaks PEAKS           File name of the identified peaks [default: peaksMLE.csv]
  --summary SUMMARY       File name of the csv file with summary values (it must contain the parameters numax and sigmaEnv) [default: summary.csv]
  --pds PDS               File name of the csv file with the background-subtracted PDS It must contain the columns "frequency" and "power" [default: pds_bgr.csv]
  --full_pds FULLPDS      File name of the csv file with the original PDS It must contain the columns "frequency" and "power" [default: pds.csv]
  --output OUTPUT         File name of the csv file on which to write the results [default: peaksMLE.csv]
  --ncores NCORES         Number of cores to use for a parallel calculation [default: 1]
  --cut CUT               Cut in zeta for rotational splitting finder [default: 0.7]
```




