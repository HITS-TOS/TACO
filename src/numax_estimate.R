## Quick numax estimation from different methods.

library(argparser, quietly = TRUE)
library(readr, quietly = TRUE)
library(tidyr, quietly = TRUE)
library(pracma, quietly = TRUE)

source("src/peakFind_lib.R", chdir = TRUE)
source("src/wavelets.R", chdir = TRUE)

numax_estimate_r <- function(pds, data, filter_width) {
    numpeaks <- 5
    do_estimation <- TRUE

    if (("numax0_flag" %in% names(data)) && ("numax0" %in% names(data))) {
        data <- data %>%
                mutate(numax0_flag = as.factor(numax0_flag))
        if (is.numeric(data$numax0)) {
            if(data$numax0_flag == FALSE) {
                print(paste("Skipping initial numax estimation "))
                do_estimation <- FALSE
            }
        }
    }

    if (do_estimation) {
        if(!("numax0_flag" %in% names(data)))
            data$numax0_flag <- FALSE
        if(data$numax0_flag == FALSE) {
            print(paste("Initial estimation of numax"))
        } else {
            print(paste("Attempting again the estimation of numax"))
        }

        ## Approximate the power spectrum background using a moving median filter in log-space
        ## ===================================================================================

        x0 <- log10(pds$frequency[1])
        # Corrective factor to convert median to mean
        # for chi-squared 2 d.o.f distributed data
        corr_factor <- (8.0 / 9.0)**3

        bkg <- c(1:length(pds$frequency)) * 0
        count <- c(1:length(pds$frequency)) * 0
        while (x0 < log10(pds$frequency[length(pds$frequency)])){
            m = abs(log10(pds$frequency) - x0) < filter_width
            if (length(bkg[m] > 0)){
                bkg[m] <- bkg[m] + median(pds$power[m], na.rm = TRUE) / corr_factor
                count[m] <- count[m] + 1
            }
            x0 <- x0 + 0.5 * filter_width
        }
        bkg <- bkg / count

        pds$power_bgr <- pds$power / bkg
        print(length(pds$power_bgr))
        #pds.SS <- splus2R::signalSeries(data = pds$power_bgr, positions. = pds$frequency)
        
        if (length(pds$power_bgr) < 40000){
            power <- pds$power_bgr
            pos <- pds$frequency
        }
        if (length(pds$power_bgr) > 40000){
            npoints <- 20000
            pos <- seq(1.,max(pds$frequency),length.out = npoints)
            difpos <- 0.5 * ((pos[2])-(pos[1]))
            print(difpos)
            power <- c(1:npoints) * 0
            i <- 1
            while (i < npoints){
                m = abs((pds$frequency) - (pos[i])) < difpos
                if (length(pds$power_bgr[m]) > 0){
                    power[i] <- median(pds$power_bgr[m], na.rm = TRUE)
                }
            i <- i+1
            }
        }
        pds.SS <- splus2R::signalSeries(data = power, positions. = pos)
       
        
        ## Estimating numax from the variance of the time-series
        ## ====================================================

        ## The variance and numax are related by the empirical relation
        ## numax = exp(fit.coef[1])*var^(fit.coef[2])

        fit.coef <- c(13.1981739, -0.7486405)
        data$numax_var <- exp(fit.coef[1]) * data$var^(fit.coef[2])
        
        #if(data$numax_var < 0.9*max(pds$frequency)) {
        #    pds.SS <-
        #        pds %>%
        #        filter(frequency > 0.2*data$numax_var) %>%
        #        splus2R::signalSeries(data       = .$power,
        #                              positions. = .$frequency)
        #}

        ## Estimating numax with a tree-map of the local maxima of a CWT
        ## =============================================================

        ## We take a mexican-hat CWT, find the 5 most significant peaks and
        ## take the median as an estimation of numax.

        max_scale <- sqrt(tail(pds$frequency, 1)) / 2
        print(max_scale)
        pds.CWT.MH <- wavCWT(x = pds.SS, wavelet = "gaussian2",
                            scale.range = c(deltat(pds.SS), max_scale))
        pds.CWTMexHat <- wavCWTTree(pds.CWT.MH)
        pds.Peaks <- cwtPeaks(pds.CWTMexHat)

        if (is.null(pds.Peaks) || dim(pds.Peaks)[1] == 0) {
            data$numax_CWTMexHat <- FALSE
        } else {
            pds.Peaks <- pds.Peaks[order(-pds.Peaks$snr)[1:numpeaks],]
            data$numax_CWTMexHat <- median(pds.Peaks[,"frequency"], na.rm = TRUE)
        }

        ## Estimating numax from a Morlet-CWT
        ## ==================================

        pds.CWT.Morlet <- wavCWT(pds.SS, wavelet = "morlet",
                                scale.range = c(deltat(pds.SS), max_scale))
        pds.CWT.Morlet.M <- Mod(as.matrix(pds.CWT.Morlet))

        ## Removing edge effects
        nuNyq <- tail(attr(pds.CWT.Morlet, "time"), 1)
        time.M <- matrix(rep(attr(pds.CWT.Morlet, "time"),
                            times = attr(pds.CWT.Morlet, "n.scale")),
                        nrow = attr(pds.CWT.Morlet, "n.sample"))
        scale.M <- matrix(rep(attr(pds.CWT.Morlet, "scale"),
                            times = attr(pds.CWT.Morlet, "n.sample")),
                        nrow = attr(pds.CWT.Morlet, "n.scale"))

        aux <- matrix(data = TRUE, nrow = attr(pds.CWT.Morlet, "n.sample"),
                    ncol = attr(pds.CWT.Morlet, "n.scale"))
        aux[t(scale.M) > (time.M - 10)/5] <- FALSE
        aux[t(scale.M) > (nuNyq - time.M)/5] <- FALSE
        aux[t(scale.M) > attr(pds.CWT.Morlet, "scale")[3]] <- FALSE
        pds.CWT.Morlet.M <- pds.CWT.Morlet.M * aux

        ## Supppose numax is the maximum value
        data$numax_Morlet <- attr(pds.CWT.Morlet,
                                    "time")[which(pds.CWT.Morlet.M == max(pds.CWT.Morlet.M),
                                                    arr.ind= TRUE)[1]]

        ## Estimate the final numax from the previous results
        ## ==================================================

        #if (abs(1 - data$numax_CWTMexHat/data$numax_Morlet) < 0.2 |
        #    abs(1 - data$numax_var/data$numax_Morlet) < 0.2)
        #{
        #    data$numax0 <-data$numax_Morlet
        #    data$numax0_flag <- FALSE
        #    flag <- 0
        #} else if (abs(1 - data$numax_CWTMexHat/data$numax_var) < 0.2) {
        #    data$numax0 <- data$numax_CWTMexHat
        #    flag <- 0
        #    data$numax0_flag <- FALSE
        #} else {
        #    # 29/06/2020 Does this even work for MS stars? var estimate will
        #    # always be dodgy!
        #    data$numax0_flag <- TRUE
        #    flag <- 1
        #    #27.10.2021 after a visual check with Nathalie; numax_CWTMexHat seems much better!!!
        #   data$numax0 <- data$numax_CWTMexHat
        #}
       
       
       if (abs(1 - min(data$numax_CWTMexHat,data$numax_Morlet)/max(data$numax_CWTMexHat,data$numax_Morlet)) < 0.5 &
           abs(1 - min(data$numax_var,data$numax_Morlet)/max(data$numax_var,data$numax_Morlet)) < 0.5 &
           abs(1 - min(data$numax_CWTMexHat,data$numax_var)/max(data$numax_CWTMexHat,data$numax_var)) < 0.5){
            data$numax0 <- (data$numax_var+data$numax_Morlet+data$numax_CWTMexHat)/3.0
            data$numax0_sd <- (max(data$numax_var,data$numax_Morlet,data$numax_CWTMexHat) - min(data$numax_var,data$numax_Morlet,data$numax_CWTMexHat))/2.0
            data$numax0_flag <- FALSE
            flag <- 0
        } else if (abs(1 - min(data$numax_CWTMexHat,data$numax_Morlet)/max(data$numax_CWTMexHat,data$numax_Morlet)) < 0.5 ){
            data$numax0 <- (data$numax_Morlet+data$numax_CWTMexHat)/2.0
            data$numax0_sd <- (max(data$numax_Morlet,data$numax_CWTMexHat) - min(data$numax_Morlet,data$numax_CWTMexHat))/2.0
            data$numax0_flag <- FALSE
            flag <- 0
        } else if (abs(1 - min(data$numax_CWTMexHat,data$numax_var)/max(data$numax_CWTMexHat,data$numax_var)) < 0.5 ){
            data$numax0 <- (data$numax_var+data$numax_CWTMexHat)/2.0
            data$numax0_sd <- (max(data$numax_var,data$numax_CWTMexHat) - min(data$numax_var,data$numax_CWTMexHat))/2.0
            data$numax0_flag <- FALSE
            flag <- 0
        } else if (abs(1 - min(data$numax_var,data$numax_Morlet)/max(data$numax_var,data$numax_Morlet)) < 0.5 ){
            data$numax0 <- (data$numax_var+data$numax_Morlet)/2.0
            data$numax0_sd <- (max(data$numax_var,data$numax_Morlet) - min(data$numax_var,data$numax_Morlet))/2.0
            data$numax0_flag <- FALSE
            flag <- 0
        } else {
            if (data$numax_CWTMexHat > 0) {
                data$numax0 <- min(data$numax_var,data$numax_Morlet,data$numax_CWTMexHat)
            } else {
                data$numax0 <- min(data$numax_var,data$numax_Morlet)
            }
            data$numax0_sd <- data$numax0/2.0
            data$numax0_flag <- TRUE
            flag <- 1
        }
    }

    return(list(data,flag))

}
