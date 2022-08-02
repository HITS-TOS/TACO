## Find the relevant solar-like oscillations in a background-removed PDS
## We use a tree-map of the local maxima found by a (mexican-hat) wavelet
## transform in the PDS and also look at the residuals to identify "unresolved"
## oscillations.

library(readr)
library(argparser, quietly = TRUE)
library(lomb, quietly = TRUE)

source("src/peakFind_lib.R", chdir = TRUE)
# 22/12/19 adding this to enable estimation of d02 from numax
source("src/l02_modes_id.R", chdir = TRUE)

peak_find_r <- function(pds, ofac_pds) {

    deltanu <- diff(pds$frequency[1:2])
    ofac_deltanu <- diff(ofac_pds$frequency[1:2])

    pds <- ofac_pds
    deltanu <- ofac_deltanu
    # Need to actually oversample!
    # If timeseries is too short then "oversample"

    if (removel02 == TRUE) {
        # Finding mixed modes!

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
}
