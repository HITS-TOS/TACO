## MLE optimization for the peaks found by peakFind.R (or another way)
## It takes a csv file with columns named frequency, height and linewidth
## and makes an MLE optimization on a PDS. If linewidth is NA, the modes are
## assumed to be sinc^2 functions, otherwise they are assumed as Lorentzians

library(argparser, quietly = TRUE)
## Where are we located?
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

## Get the file names and options
p <- arg_parser("MLE estimation of the peaks found in a background-removed PDS.")
p <- add_argument(p, arg = "--pds", nargs=1, default = "pds_bgr.csv",
                  help = "File name of the background-removed PDS in csv format with headers 'frequency' and 'power'",
                  type = "character", flag = FALSE)
p <- add_argument(p, arg = "--init", nargs=1, default = "peaks.csv",
                  help = "File name of the intial guesses",
                  type = "character", flag = FALSE)
p <- add_argument(p, arg = "--summary", nargs=1, default = "summary.csv",
                  help = "File name of the csv file with summary values. It must contain numax and sigmaEnv",
                  type = "character", flag = FALSE)
p <- add_argument(p, arg ="--mixedpeaks", nargs=1, default = "mixedpeaks.csv",
                  help = "File name of the csv file containing the mixed mode peaks from peak finding.")
p <- add_argument(p, arg = "--mixedoutput", nargs=1, default = "mixedpeaksMLE.csv",
                  help = "File name of the csv file to save the mixed mode MLE results to.")
p <- add_argument(p, arg = "--output", nargs=1, default = "peaksMLE.csv",
                  help = "File name on which to save the results",
                  type = "character", flag = FALSE)
# 22/12/19 same change as in peakFind.R
p <- add_argument(p, arg = "--maxlwd", nargs=1, default = NA,
                  help = "Maximum linewidth for peaks in the MLE optimization",
                  flag = FALSE)
p <- add_argument(p, arg = "--minAIC", nargs=1, default = 2,
                   help = "Drop peaks that have an AIC lower that minAIC",
                   type = "numeric", flag = FALSE)
## p <- add_argument(p, arg = "--no_sd", nargs=0,
##                   help = "Use this flag to skip calculating the MLE uncertainties for each parameter",
##                   flag = TRUE)
p <- add_argument(p, arg = "--removel02", nargs=1, default = FALSE,
                  help = "Whether or not the l02 peaks should be divided out before running the MLE fit",
                  type = "boolean", flag = FALSE)
p <- add_argument(p, arg = "--finalfit", nargs=1, default = FALSE,
                  help = "Whether or not this is the final MLE optimisation.",
                  type = "boolean", flag = FALSE)
p <- add_argument(p, arg = "--navg", nargs=1, default = 1,
                  help = "Number of power spectra averaged to create current power spectrum")
argv <- parse_args(p)

library(readr, quietly = TRUE)

## Check the provided PDS file is valid
if(!file.exists(argv$pds))
    stop(paste("Could not find the background-removed PDS file:", argv$pds))
if(!file.exists(argv$summary))
    stop(paste("Could not find the summary file:", argv$summary))
if(!file.exists(argv$init))
    stop(paste("Could not find the initial-guesses file:", argv$init))
if(!file.exists(file.path(script.basename, "peakFind_lib.R")))
    stop("Could not find peakFind_lib.R")

source(file.path(script.basename, "peakFind_lib.R"))
# 22/12/19 adding this to enable estimation of d02 from numax
source(file.path(script.basename, "l02_modes_id.R"))
## Read the data
d.pds_bgr <-
    read_csv(argv$pds,
             col_types = cols(frequency = col_number(),
                              power     = col_number()))
d.summary <-
    read_csv(argv$summary,
             col_types = cols(numax    = col_number(),
                              sigmaEnv = col_number()))
d.peaks <-
    read_csv(argv$init,
             col_types = cols(frequency = col_number(),
                              amplitude = col_number(),
                              linewidth = col_number())) %>%
    arrange(frequency)

## Trim to relevant region
d.peaks <-
    d.peaks %>%
    filter(frequency > d.summary$numax - 3*d.summary$sigmaEnv &
           frequency < d.summary$numax + 3*d.summary$sigmaEnv)

d.pds_bgr <-
    d.pds_bgr %>%
    filter(frequency > d.summary$numax - 3*d.summary$sigmaEnv &
           frequency < d.summary$numax + 3*d.summary$sigmaEnv)

deltanu <- diff(d.pds_bgr$frequency[1:2])


# Check to see if this is the final MLE fit
if (argv$finalfit == TRUE) {
    ## Do the calculations

    if(nrow(d.peaks) == 0){
        d.peaks.mle <- tibble(frequency=double(),
                              amplitude=double(),
                              linewidth=double(),
                              frequency_sd=double(),
                              amplitude_sd=double(),
                              linewidth_sd=double(),
                              height=double(),
                              AIC=double())
        write_csv(x = d.peaks.mle, file = argv$output)
        stop("No peaks detected so not proceeding with MLE fit")
    }


    d.peaks.mle <-
        peaks_MLE_final_sd(peaks = d.peaks, pds = d.pds_bgr, 
                           final_fit_factor = 0.1, naverages=argv$navg) %>%
        arrange(frequency) 

    # 31/01/2020 Add n and l ID back in as no frequencies will have been removed
    d.peaks.mle$n <- d.peaks$n
    d.peaks.mle$l <- d.peaks$l

    # Filter now by AIC - this is done to enable line above (so can use d.peaks)
    # rather than having to filter to find right peaks
    d.peaks.mle <- d.peaks.mle %>%
                        filter(AIC > argv$minAIC)

    DeltaNu <- d.summary$DeltaNu
    eps_p <- d.summary$eps_p
    
    ## Write the output
    write_csv(x = d.peaks.mle, file = argv$output)
    write_csv(x = d.summary,   file = argv$summary)

    print(warnings())

} else {
    #   Only find resolved peaks from mixed modes
    if (argv$removel02 == TRUE){
        print("Removing l=0,2,3 first")
        if (nrow(d.peaks) != 0){
            d.l02_peaks <-
                d.peaks %>%
                filter(l == 0 | l == 2 | l == 3)
            d.pds_l02_removed <-
                d.pds_bgr %>%
                mutate(power = power / fit_model(pds = ., peaks = d.l02_peaks))
        } else {
            d.l02_peaks <- d.peaks
            d.pds_l02_removed <- d.pds_bgr
        }
        ## Do the calculations
        d.rest_peaks <-     
                read_csv(argv$mixedpeaks,
                col_types = cols(frequency = col_number(),
                                amplitude = col_number(),
                                linewidth = col_number())) %>%
                arrange(frequency)

        if(nrow(d.rest_peaks) == 0){
            stop("No peaks detected so not proceeding with MLE fit")
        }

        print(d.rest_peaks)

        # If max linewidth argument not set
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
    
        d.peaks.mle <-
            peaks_MLE_sd(peaks = d.rest_peaks, pds = d.pds_l02_removed, maxLWD = maxlwd, naverages=argv$navg) %>%
            arrange(frequency) %>%
            filter(AIC > argv$minAIC)
            
        # Add n,l & m columns so consistent with l=0,2 peaks file
        d.peaks.mle[,c("n","l")] <- NA

        ## Write the output
        write_csv(x = d.peaks.mle, file = argv$mixedoutput)
        ## Concatenate l=0,2 from first peaks file and mixed peaks file and save to peaksMLE.csv
      
        d.peaks.full <- rbind(d.l02_peaks, d.peaks.mle)

        write_csv(x = d.peaks.full, file= argv$output)
        #write_csv(x = d.summary,   file = argv$summary)

        print(warnings())
    } else {

        if(nrow(d.peaks) == 0){
            d.peaks.mle <- data.frame(frequency=double(),
                                    amplitude=double(),
                                    linewidth=double(),
                                    frequency_sd=double(),
                                    amplitude_sd=double(),
                                    linewidth_sd=double(),
                                    height=double(),
                                    AIC=double())
            write_csv(x = d.peaks.mle, file = argv$output)
            stop("No peaks detected so not proceeding with MLE fit")
        }

        # 22/12/19 Add in maxlwd default so consistent with peakfind
        if (is.na(argv$maxlwd)){
            deltanu_est <- DeltaNu_from_numax(d.summary$numax)
            # Set to be < d02 from scaling relation (~0.125 dnu)
            # Divide by 2 because HWHM defined here and want FWHM to be less than ~d02
            maxlwd <- 0.125 * deltanu_est / 2
            print(paste("Maximum peak linewidth (HWHM) not set, therefore taking estimated value ", maxlwd, "uHz"))
            if(maxlwd < deltanu/2){
                maxlwd <- deltanu/2 + 0.1*deltanu/2
                print(paste("Estimated value less than frequency bin width, therefore setting to slightly larger than bin width", maxlwd, "uHz"))
            }
        } else{
            maxlwd <- as.numeric(argv$maxlwd)
            print(paste("Maximum peak linewidth (HWHM) set, using value ", maxlwd, "uHz"))
        }
        ## Do the calculations
        d.peaks.mle <-
            peaks_MLE_sd(peaks = d.peaks, pds = d.pds_bgr, maxLWD = maxlwd, naverages=argv$navg) %>%
            arrange(frequency) %>%
            filter(AIC > argv$minAIC)
            
            if (max(d.peaks.mle$linewidth) >= 0.95*maxlwd){
                print("modes larger than 0.95*maxlwd, redoing the fit")

                # large_peaks <- d.peaks.mle %>%
                #                 filter(linewidth >= 0.9 * maxlwd)


                peaks_split1 <-
                d.peaks.mle %>%
                filter(linewidth >= 0.95*maxlwd) %>%
                mutate(frequency = frequency-0.5*maxlwd) %>%
                mutate(linewidth = maxlwd/3.0)
        
                peaks_split2 <-
                d.peaks.mle %>%
                filter(linewidth >= 0.95*maxlwd) %>%
                mutate(frequency = frequency+0.5*maxlwd) %>%
                mutate(linewidth = maxlwd/3.0)
                    
                peaks_all <-
                d.peaks.mle %>%
                filter(linewidth < 0.9*maxlwd)
                    
                peaks_new <- bind_rows(peaks_split1, peaks_split2, peaks_all)
            
                N <- nrow(d.pds_bgr)
                LL1 <- log_likelihood(d.pds_bgr, fit_model(d.pds_bgr, d.peaks.mle), naverages=1)
                k1 <- 3*nrow(d.peaks.mle) - sum(is.na(d.peaks.mle$linewidth))
                AIC1 <- model_AIC(LL1, k1, N)
                d.peaks.mle2 <-
                peaks_MLE_sd(peaks = peaks_new, pds = d.pds_bgr, maxLWD = maxlwd, naverages=argv$navg) %>%
                    arrange(frequency) %>%
                    filter(AIC > argv$minAIC)
                                         
                LL2 <- log_likelihood(d.pds_bgr, fit_model(d.pds_bgr, d.peaks.mle2), naverages=1)
                k2 <- 3*nrow(d.peaks.mle2) - sum(is.na(d.peaks.mle2$linewidth))
                AIC2 <- model_AIC(LL2, k2, N)

                print("AIC")
                print(AIC1)
                print(AIC2)

                # If second model is better then keep it
                if(AIC2 < AIC1){
                    d.peaks.mle <- d.peaks.mle2
                }

                d.peaks.mle <- d.peaks.mle %>%
                    filter(AIC > argv$minAIC)
            }
        d.summary <-
            d.summary %>%
            mutate(npeaks = nrow(d.peaks.mle %>% filter(AIC > argv$minAIC)))

        ## Write the output
        
        write_csv(x = d.peaks.mle, file = argv$output)
        write_csv(x = d.summary,   file = argv$summary)

        print(warnings())

    }
}
