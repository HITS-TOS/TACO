## Find the relevant solar-like oscillations in a background-removed PDS
## We use a tree-map of the local maxima found by a (mexican-hat) wavelet
## transform in the PDS and also look at the residuals to identify "unresolved"
## oscillations.

library(readr)
library(argparser, quietly=TRUE)
library(lomb, quietly=TRUE)
#library(tcltk)
#library(ggplot2)


## Where are we located?
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

## Get the file names and peak-finding options
p <- arg_parser("Peak-finding in a background-removed PDS. We only look for peaks in the region numax +- 4*sigmaEnv.")
p <- add_argument(p, arg = "--pds", nargs=1, default = "pds_bgr.csv",
                  help = "File name of the background-removed PDS in csv format with headers 'frequency' and 'power'",
                  type = "character", flag = FALSE)
p <- add_argument(p, arg = "--ofac_pds", nargs=1, default = "ofac_pds_bgr.csv",
                  help = "File name of the background-removed PDS in csv format with headers 'frequency' and 'power'",
                  type = "character", flag = FALSE)
p <- add_argument(p, arg = "--summary", nargs=1, default = "summary.csv",
                  help = "File name of the summary file in csv format with at least the headers 'numax' and 'sigmaEnv'",
                  type = "character", flag = FALSE)
p <- add_argument(p, arg = "--output", nargs=1, default = "peaks.csv",
                  help = "File name on which to save the results",
                  type = "character", flag = FALSE)
p <- add_argument(p, arg = "--snr", nargs=1, default = 1.2,
                  help = "Minimum signal-to-noise ratio (on CWT space) for resolved peaks",
                  type = "numeric", flag = FALSE)
p <- add_argument(p, arg = "--prob", nargs=1, default = 0.0001,
                  help = "Minimum (frequentist) probability threshold for unresolved peaks",
                  type = "numeric", flag = FALSE)
# 22/12/19 changed from numeric and default=0.5 to default = NA and type string
# This enables better detection of when argument is given
p <- add_argument(p, arg = "--maxlwd", nargs=1, default = NA,
                  help = "Maximum search linewidth for resolved peaks in the CWT search",
                  flag = FALSE)
p <- add_argument(p, arg = "--removel02", nargs=1, default = FALSE,
                  help = "Whether or not the l02 peaks should be divided out before running the CWT search",
                  type = "boolean", flag = FALSE)
p <- add_argument(p, arg = "--peaks", nargs=1, default = "peaksMLE.csv",
                  help = "File name of the identified peaks. It must contain the l=0,2 modes already identified",
                  type = "character", flag = FALSE)
p <- add_argument(p, arg ="--mixedpeaks", nargs=1, default = "mixedpeaks.csv",
                  help = "File name of the csv file containing the mixed mode peaks from peak finding.")
p <- add_argument(p, arg = "--minAIC", nargs=1, default = 2,
                  help = "Minimum AIC value for a peak to have in order to be considered significant")
p <- add_argument(p, arg = "--navg", nargs=1, default = 1,
                  help = "Number of power spectra averaged to create current power spectrum")
argv <- parse_args(p)

## Check the provided PDS file is valid
if(!file.exists(argv$pds))
    stop(paste("Could not find the background-removed PDS file:", argv$pds))
if(!file.exists(argv$summary))
    stop(paste("Could not find the summary file:", argv$summary))

## We are good to go!
source(file.path(script.basename, "peakFind_lib.R"))
# 22/12/19 adding this to enable estimation of d02 from numax
source(file.path(script.basename, "l02_modes_id.R"))

d.summary <- read_csv(argv$summary,
                      col_types = cols(numax = col_number(),
                                       sigmaEnv = col_number()))

d.pds <-
    read_csv(argv$pds,
             col_types = cols(frequency = col_number(),
                              power     = col_number())) %>%
    filter(frequency > d.summary$numax - 3*d.summary$sigmaEnv,
           frequency < d.summary$numax + 3*d.summary$sigmaEnv)
d.ofac_pds <-
    read_csv(argv$ofac_pds,
             col_types = cols(frequency = col_number(),
                              power     = col_number())) %>%
    filter(frequency > d.summary$numax - 3*d.summary$sigmaEnv,
           frequency < d.summary$numax + 3*d.summary$sigmaEnv)

deltanu <- diff(d.pds$frequency[1:2])
ofac_deltanu <- diff(d.ofac_pds$frequency[1:2])
print(deltanu)
print(ofac_deltanu)

d.pds <- d.ofac_pds
deltanu <- ofac_deltanu
# Need to actually oversample!
# If timeseries is too short then "oversample"



# #new_freq <- seq(from=min(d.pds$frequency), to=max(d.pds$frequency), by=deltanu/2)
# d.new_pds <- as_tibble(approx(d.pds$frequency, d.pds$power, new_freq))
# colnames(d.new_pds) <- c("frequency", "power")
# deltanu <- diff(d.new_pds$frequency[1:2])


if (argv$removel02 == TRUE){
    # Finding mixed modes!

    #new_freq <- seq(from=min(d.pds$frequency), to=max(d.pds$frequency), by=deltanu/2)
    #d.pds <- as_tibble(approx(d.pds$frequency, d.pds$power, new_freq))
    #colnames(d.pds) <- c("frequency", "power")
    #deltanu <- diff(d.pds$frequency[1:2])


    print("Removing l=0,2,3 first")
    d.peaks <- try(read_csv(argv$peaks, col_types = cols()) %>%
                        filter(frequency > d.summary$numax - 3*d.summary$sigmaEnv,
                               frequency < d.summary$numax + 3*d.summary$sigmaEnv))
    if("try-error" %in% class(d.peaks)){
        d.pds_l02_removed <- d.pds

    } else{

        if(nrow(d.peaks) != 0){
            d.l02_peaks <-
                d.peaks %>%
                filter(l == 0 | l == 2 | l == 3)
            d.pds_l02_removed <-
                d.pds %>%
                mutate(power = power / fit_model(pds = ., peaks = d.l02_peaks))
        } else{
            d.pds_l02_removed <- d.new_pds
        }
    } 

    #X11()
    #plot(d.pds, type='l')
    #print(CWTTree[[branch_number]])
    #plot(CWTTree[[branch_number]]$scale, CWTTree[[branch_number]]$extrema)
    #print(peaks_idx)
    #prompt  <- "hit spacebar to close plots"
    #extra   <- "some extra comment"
    #capture <- tk_messageBox(message = prompt, detail = extra) 

    # 01/05/2020 Oversample if resolution is too bad to find "unresolved" modes.
    #deltanu = mean(diff(d.pds$frequency))
    #new_freq = seq(from=min(d.pds$frequency), to=max(d.pds$frequency), by=deltanu/50.)
    ##print(d.pds)
    #d.pds_l02_removed <- data.frame(approx(d.pds$frequency, d.pds$power, new_freq))
    #colnames(d.pds_l02_removed) <- c("frequency", "power")

    # If max linewidth argument not set
    #stop()
    if (is.na(argv$maxlwd)){
        # Since this is for finding mixed modes we add in constraint that linewidth must be less than Gamma0
        if(is.null(d.summary$gamma0)) {
            maxlwd <- deltanu/2 + 0.1*deltanu/2
            print(paste("Gamma0 from summary file is NA therefore setting to slightly larger than bin width: ", maxlwd, "uHz"))
            #stop("Gamma0 from summary file is NA please supply a maximum linewidth")
        } else {
            maxlwd <- 1.5*d.summary$gamma0
            print(paste("Maximum peak linewidth (HWHM) not set, therefore set to Gamma0: ", maxlwd, "uHz"))
        }
    } else{
        maxlwd <- 1.5*as.numeric(argv$maxlwd)
        print(paste("Maximum peak linewidth (HWHM) set, using value ", maxlwd, "uHz"))
    }
    # 23/12/19 Because we are fitting in terms of HWHM the minimum linewidth should be bw/2 not bw!

    d.new_peaks <- peak_find(d.pds_l02_removed, min.snr = argv$snr, p = argv$prob,
                        linewidth.range = c(deltanu/2, maxlwd),
                        find.resolved.only=FALSE, naverages=argv$navg)
    ## Choose only the significant ones

    d.new_peaks <-
        d.new_peaks %>%
        arrange(frequency) %>%
        #mutate(EFP = nrow(d.pds_l02_removed) * EFPp) %>%
        filter(AIC > argv$minAIC)# %>%
        #filter(EFP < 1) # Also added this line to see if this helps trim away bad peaks!
    write_csv(d.new_peaks, file = argv$mixedpeaks)#argv$output)

} else {
    #  Only find resolved peaks
    # If value not given then set it to be equal to 0.1 * dnu or just less than ~d02
    if (is.na(argv$maxlwd)){
        deltanu_est <- DeltaNu_from_numax(d.summary$numax)
        # Set to be < d02 from scaling relation (~0.125 dnu)
        # Divide by 2 because HWHM defined here and want FWHM to be less than ~d02
        maxlwd <-  0.125 * deltanu_est / 2
        if(maxlwd < deltanu/2){
            print("Maximum peak linewidth (HWHM) is less than frequency resolution! Increasing upper bound slightly")
            maxlwd <- deltanu/2 + 0.1*deltanu/2
        }
        print(paste("Maximum peak linewidth (HWHM) not set, therefore taking estimated value ", maxlwd, "uHz"))
    } else{
        maxlwd <- as.numeric(argv$maxlwd)
        print(paste("Maximum peak linewidth (HWHM) set, using value ", maxlwd, "uHz"))
    }

    print(deltanu)

    d.peaks <- peak_find(d.pds, min.snr = argv$snr, p = argv$prob,
                        linewidth.range = c(deltanu/2, maxlwd),
                        find.resolved.only=TRUE, naverages=argv$navg)

    if (is.null(d.peaks)){
        d.peaks <- tibble(frequency=double(),
                              linewidth=double(),
                              height=double(),
                              snr=double(),
                              AIC=double(),
                              amplitude=double(),
                              EFPp=double(),
                              EFP=double())
        write_csv(d.peaks, file=argv$output)
        print("No resolved peaks found!")
    } else{

        ## Choose only the significant ones
        d.peaks <-
            d.peaks %>%
            arrange(frequency) %>%
            #mutate(EFP = nrow(d.pds) * EFPp) %>%
            filter(AIC > argv$minAIC)# %>%
            #filter(EFP < 0.5) # Also added this line to see if this helps trim away bad peaks!

        write_csv(d.peaks, file = argv$output)
    }
}
warnings()
