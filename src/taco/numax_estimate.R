## Quick numax estimation from different methods.

library(argparser, quietly=TRUE)
library(readr, quietly=TRUE)
library(tidyr, quietly=TRUE)

source(file.path(script.basename, "peakFind_lib.R"), echo=FALSE)

if(!file.exists(file.path(script.basename, "wavelets.R"))) {
    stop("Could not find wavelets.R")
}
source(file.path(script.basename, "wavelets.R"))

numax_estimate_r = function(pds, variance, nyquist, filterwidth)
{
    NUMPEAKS <- 5

    do_estimation <- TRUE

    if (("numax0_flag" %in% names(d.summary)) & ("numax0" %in% names(d.summary))) {
        d.summary <-
            d.summary %>%
            mutate(numax0_flag = as.factor(numax0_flag))
        if (is.numeric(d.summary$numax0)) {
            if(d.summary$numax0_flag == FALSE) {
                print(paste("Skipping initial numax estimation for", basename(getwd())))
                do_estimation <- FALSE
            }
        }
    }

    if(do_estimation) {
        if(!("numax0_flag" %in% names(d.summary)))
            d.summary$numax0_flag <- FALSE
        if(d.summary$numax0_flag == FALSE) {
            print(paste("Initial estimation of numax for KIC", basename(getwd())))
        } else {
            print(paste("Attempting again the estimation of numax for KIC", basename(getwd())))
        }

        d.pds <- read_csv(argv$pds,
                        col_types = cols(frequency = col_number(),
                                        power     = col_number()))

        ## Approximate the power spectrum background using a moving median filter in log-space
        ## ===================================================================================

        x0 <- log10(d.pds$frequency[1])
        # Corrective factor to convert median to mean for chi-squared 2 d.o.f distributed data
        corr_factor <- (8.0/9.0)**3
        #
        bkg <- c(1:length(d.pds$frequency))*0
        count <- c(1:length(d.pds$frequency))*0
        filter_width <- argv$filterwidth
        while (x0 < log10(d.pds$frequency[length(d.pds$frequency)])){
            m = abs(log10(d.pds$frequency) - x0) < filter_width
            if (length(bkg[m] > 0)){
                bkg[m] <- bkg[m] + median(d.pds$power[m], na.rm=TRUE) / corr_factor
                count[m] <- count[m] + 1
            }
            x0 <- x0 + 0.5 * filter_width
        }
        bkg <- bkg / count

        d.pds$power_bgr <- d.pds$power / bkg


        #X11()
        #print(res)
        #d.pds$bkg <- bkg
        #h <- ggplot() + 
        #     geom_line(data = d.pds, aes(frequency, power_bgr), color='black')
        #h <- ggplot() + 
        #     geom_line(data=d.pds, aes(frequency, power), color='black') + 
        #     geom_line(data=d.pds, aes(frequency, bkg), color='red') + 
        #     scale_x_continuous(trans='log10') + 
        #     scale_y_continuous(trans='log10') # + coord_cartesian(x)
        #plot(h)
        #prompt  <- "hit spacebar to close plots"
        #extra   <- "some extra comment"
        #capture <- tk_messageBox(message = prompt, detail = extra) 
        #stop()

        d.pds.SS <- splus2R::signalSeries(data = d.pds$power_bgr, positions. = d.pds$frequency)

        ## Estimating numax from the variance of the time-series
        ## ====================================================

        ## The variance and numax are related by the empirical relation
        ## numax = exp(fit.coef[1])*var^(fit.coef[2])

        fit.coef <- c(13.1981739, -0.7486405)
        d.summary$numax_var <- exp(fit.coef[1])*d.summary$var^(fit.coef[2])

        #if(d.summary$numax_var < 0.9*max(d.pds$frequency)) {
        #    d.pds.SS <-
        #        d.pds %>%
        #        filter(frequency > 0.2*d.summary$numax_var) %>%
        #        splus2R::signalSeries(data       = .$power,
        #                              positions. = .$frequency)
        #}

        ## Estimating numax with a tree-map of the local maxima of a CWT
        ## =============================================================

        ## We take a mexican-hat CWT, find the 5 most significant peaks and
        ## take the median as an estimation of numax.

        max_scale <- sqrt(tail(d.pds$frequency, 1))/2
        d.pds.CWT.MH <- wavCWT(x = d.pds.SS, wavelet = "gaussian2",
                            scale.range = c(deltat(d.pds.SS), max_scale))
        d.pds.CWTMexHat <- wavCWTTree(d.pds.CWT.MH)
        d.pds.Peaks <- cwtPeaks(d.pds.CWTMexHat)

        if(dim(d.pds.Peaks)[1] == 0) {
            d.summary$numax_CWTMexHat <- FALSE
        } else {
            d.pds.Peaks <- d.pds.Peaks[order(-d.pds.Peaks$snr)[1:NUMPEAKS],]
            d.summary$numax_CWTMexHat <- median(d.pds.Peaks[,"frequency"], na.rm=TRUE)
        }


        ## Estimating numax from a Morlet-CWT
        ## ==================================

        d.pds.CWT.Morlet <- wavCWT(d.pds.SS, wavelet = "morlet",
                                scale.range = c(deltat(d.pds.SS), max_scale))
        d.pds.CWT.Morlet.M <- Mod(as.matrix(d.pds.CWT.Morlet))

        ## Removing edge effects
        nuNyq <- tail(attr(d.pds.CWT.Morlet, "time"), 1)
        time.M <- matrix(rep(attr(d.pds.CWT.Morlet, "time"),
                            times = attr(d.pds.CWT.Morlet, "n.scale")),
                        nrow = attr(d.pds.CWT.Morlet, "n.sample"))
        scale.M <- matrix(rep(attr(d.pds.CWT.Morlet, "scale"),
                            times = attr(d.pds.CWT.Morlet, "n.sample")),
                        nrow = attr(d.pds.CWT.Morlet, "n.scale"))

        aux <- matrix(data = TRUE, nrow = attr(d.pds.CWT.Morlet, "n.sample"),
                    ncol = attr(d.pds.CWT.Morlet, "n.scale"))
        aux[t(scale.M) > (time.M - 10)/5] <- FALSE
        aux[t(scale.M) > (nuNyq - time.M)/5] <- FALSE
        aux[t(scale.M) > attr(d.pds.CWT.Morlet, "scale")[3]] <- FALSE
        d.pds.CWT.Morlet.M <- d.pds.CWT.Morlet.M * aux

        ## Supppose numax is the maximum value
        d.summary$numax_Morlet <- attr(d.pds.CWT.Morlet,
                                    "time")[which(d.pds.CWT.Morlet.M == max(d.pds.CWT.Morlet.M),
                                                    arr.ind= TRUE)[1]]

        ## Estimate the final numax from the previous results
        ## ==================================================

        if (abs(1 - d.summary$numax_CWTMexHat/d.summary$numax_Morlet) < 0.2 |
            abs(1 - d.summary$numax_var/d.summary$numax_Morlet) < 0.2)
        {
            d.summary$numax0 <-d.summary$numax_Morlet
            d.summary$numax0_flag <- FALSE
        } else if (abs(1 - d.summary$numax_CWTMexHat/d.summary$numax_var) < 0.2) {
            d.summary$numax0 <- d.summary$numax_CWTMexHat
            d.summary$numax0_flag <- FALSE
        } else {
            # 29/06/2020 Does this even work for MS stars? var estimate will
            # always be dodgy!
            d.summary$numax0_flag <- TRUE
            #27.10.2021 after a visual check with Nathalie; numax_CWTMexHat seems much better!!!
            d.summary$numax0 <- d.summary$numax_CWTMexHat
        #   if (d.summary$numax_var < d.summary$nuNyq) {
        #       d.summary$numax0 <- d.summary$numax_var
        #   } else {
        #       d.summary$numax0 <- d.summary$numax_Morlet
        #   }
        }

        ## Save the results
        ## ================

        write_csv(d.summary, argv$summary)

        ## Warning if the estimate is unreliable
        if (d.summary$numax0_flag) {
            write.table(data.frame(), col.names=FALSE,
                        file = file.path(dirname(argv$summary), "NUMAX0_FLAG"))
        }
    }
}