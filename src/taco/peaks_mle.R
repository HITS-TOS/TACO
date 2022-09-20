## MLE optimization for the peaks found by peakFinR (or another way)
## It takes a csv file with columns named frequency, height and linewidth
## and makes an MLE optimization on a PDS. If linewidth is NA, the modes are
## assumed to be sinc^2 functions, otherwise they are assumed as Lorentzians

library(readr, quietly = TRUE)

source("src/peakFind_lib.R", chdir = TRUE)
source("src/l02_modes_id.R", chdir = TRUE)

peaks_mle_r <- function(pds, peaks, data, mixedpeaks, maxlwd,
                        removel02, minAIC, navg, finalfit) {

    peaks <- peaks %>% arrange(frequency)

    ## Trim to relevant region
    peaks <-
        peaks %>%
        filter(frequency > data$numax - 3 * data$sigmaEnv &
               frequency < data$numax + 3 * data$sigmaEnv)

    pds <-
        pds %>%
        filter(frequency > data$numax - 3 * data$sigmaEnv &
               frequency < data$numax + 3 * data$sigmaEnv)

    deltanu <- diff(pds$frequency[1:2])

    if (finalfit == TRUE) {

        if (nrow(peaks) == 0) {
            peaks.mle <- tibble(frequency=double(),
                                amplitude=double(),
                                linewidth=double(),
                                frequency_sd=double(),
                                amplitude_sd=double(),
                                linewidth_sd=double(),
                                height=double(),
                                AIC=double())
            stop("No peaks detected so not proceeding with MLE fit")
            return(peaks.mle)
        }

        peaks.mle <-
            peaks_MLE_final_sd(peaks = peaks, pds = pds,
                            final_fit_factor = 0.1, naverages = navg) %>%
            arrange(frequency)

        # 31/01/2020 Add n and l ID back in as no frequencies will have been removed
        peaks.mle$n <- peaks$n
        peaks.mle$l <- peaks$l

        # Filter now by AIC - this is done to enable line above (so can use peaks)
        # rather than having to filter to find right peaks
        peaks.mle <- peaks.mle %>% filter(AIC > minAIC)

        DeltaNu <- data$DeltaNu
        eps_p <- data$eps_p

    } else {
        #   Only find resolved peaks from mixed modes
        if (removel02 == TRUE) {
            print("Removing l=0,2,3 first")
            if (nrow(peaks) != 0) {
                l02_peaks <-
                    peaks %>%
                    filter(l == 0 | l == 2 | l == 3)
                pds_l02_removed <-
                    pds %>%
                    mutate(power = power / fit_model(pds = ., peaks = l02_peaks))
            } else {
                l02_peaks <- peaks
                pds_l02_removed <- pds
            }

            rest_peaks <- mixedpeaks %>% arrange(frequency)

            if (nrow(rest_peaks) == 0) {
                stop("No peaks detected so not proceeding with MLE fit")
            }

            # If max linewidth argument not set
            if (is.null(maxlwd)) {
                # Since this is for finding mixed modes we add in constraint that linewidth must be less than Gamma0
                if (is.null(data$gamma0)) {
                    maxlwd <- deltanu / 2 + 0.1 * deltanu / 2
                    print(paste("Gamma0 from summary file is NA therefore setting to slightly larger than bin width: ", maxlwd, "uHz"))
                } else {
                    maxlwd <- 1.5 * data$gamma0
                    print(paste("Maximum peak linewidth (HWHM) not set, therefore set to Gamma0: ", maxlwd, "uHz"))
                }
            } else {
                maxlwd <- 1.5 * as.numeric(maxlwd)
                print(paste("Maximum peak linewidth (HWHM) set, using value ", maxlwd, "uHz"))
            }

            peaks.mle <-
                peaks_MLE_sd(peaks = rest_peaks, pds = pds_l02_removed, maxLWD = maxlwd, naverages = navg) %>%
                arrange(frequency) %>%
                filter(AIC > minAIC)

            # Add n,l & m columns so consistent with l=0,2 peaks file
            peaks.mle[,c("n","l")] <- NA

            ## Write the output
            write_csv(x = peaks.mle, file = argv$mixedoutput)
            ## Concatenate l=0,2 from first peaks file and mixed peaks file and save to peaksMLE.csv

            peaks.full <- rbind(l02_peaks, peaks.mle)

        } else {

            if(nrow(peaks) == 0){
                peaks.mle <- data.frame(frequency=double(),
                                        amplitude=double(),
                                        linewidth=double(),
                                        frequency_sd=double(),
                                        amplitude_sd=double(),
                                        linewidth_sd=double(),
                                        height=double(),
                                        AIC=double())
                stop("No peaks detected so not proceeding with MLE fit")
            }

            # 22/12/19 Add in maxlwd default so consistent with peakfind
            if (is.null(maxlwd)) {
                deltanu_est <- DeltaNu_from_numax(data$numax)
                # Set to be < d02 from scaling relation (~0.125 dnu)
                # Divide by 2 because HWHM defined here and want FWHM to be less than ~d02
                maxlwd <- 0.125 * deltanu_est / 2
                print(paste("Maximum peak linewidth (HWHM) not set, therefore taking estimated value ", maxlwd, "uHz"))
                if (maxlwd < deltanu / 2) {
                    maxlwd <- deltanu / 2 + 0.1 * deltanu / 2
                    print(paste("Estimated value less than frequency bin width, therefore setting to slightly larger than bin width", maxlwd, "uHz"))
                }
            } else {
                maxlwd <- as.numeric(maxlwd)
                print(paste("Maximum peak linewidth (HWHM) set, using value ", maxlwd, "uHz"))
            }
            ## Do the calculations
            peaks.mle <-
                peaks_MLE_sd(peaks = peaks, pds = pds, maxLWD = maxlwd, naverages = navg) %>%
                arrange(frequency) %>%
                filter(AIC > minAIC)

                if (max(peaks.mle$linewidth) >= 0.95 * maxlwd) {
                    print("modes larger than 0.95*maxlwd, redoing the fit")

                    peaks_split1 <-
                    peaks.mle %>%
                    filter(linewidth >= 0.95*maxlwd) %>%
                    mutate(frequency = frequency-0.5*maxlwd) %>%
                    mutate(linewidth = maxlwd/3.0)

                    peaks_split2 <-
                    peaks.mle %>%
                    filter(linewidth >= 0.95*maxlwd) %>%
                    mutate(frequency = frequency+0.5*maxlwd) %>%
                    mutate(linewidth = maxlwd/3.0)

                    peaks_all <-
                    peaks.mle %>%
                    filter(linewidth < 0.9*maxlwd)

                    peaks_new <- bind_rows(peaks_split1, peaks_split2, peaks_all)

                    N <- nrow(pds)
                    LL1 <- log_likelihood(pds, fit_model(pds, peaks.mle), naverages=1)
                    k1 <- 3*nrow(peaks.mle) - sum(is.na(peaks.mle$linewidth))
                    AIC1 <- model_AIC(LL1, k1, N)
                    peaks.mle2 <-
                    peaks_MLE_sd(peaks = peaks_new, pds = pds, maxLWD = maxlwd, naverages = navg) %>%
                        arrange(frequency) %>%
                        filter(AIC > minAIC)

                    LL2 <- log_likelihood(pds, fit_model(pds, peaks.mle2), naverages=1)
                    k2 <- 3*nrow(peaks.mle2) - sum(is.na(peaks.mle2$linewidth))
                    AIC2 <- model_AIC(LL2, k2, N)

                    print("AIC")
                    print(AIC1)
                    print(AIC2)

                    # If second model is better then keep it
                    if (AIC2 < AIC1) {
                        peaks.mle <- peaks.mle2
                    }

                    peaks.mle <- peaks.mle %>% filter(AIC > minAIC)
                }
            data <- data %>%
                mutate(npeaks = nrow(peaks.mle %>% filter(AIC > minAIC)))
        }
    }
    return(list(peaks.mle, data))
}
