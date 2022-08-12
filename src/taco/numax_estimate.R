## Quick numax estimation from different methods.

library(argparser, quietly = TRUE)
library(readr, quietly = TRUE)
library(tidyr, quietly = TRUE)

source("src/peakFind_lib.R", chdir = TRUE)
source("src/wavelets.R", chdir = TRUE)

numax_estimate_r <- function(pds, data, filterwidth) {
    numpeaks <- 5
    do_estimation <- TRUE

    # if (("numax0_flag" %in% names(d.summary)) && ("numax0" %in% names(d.summary))) {
    #     d.summary <-
    #         d.summary %>%
    #         mutate(numax0_flag = as.factor(numax0_flag))
    #     if (is.numeric(d.summary$numax0)) {
    #         if(d.summary$numax0_flag == FALSE) {
    #             print(paste("Skipping initial numax estimation for", basename(getwd())))
    #             do_estimation <- FALSE
    #         }
    #     }
    # }

    # if (do_estimation) {
    #     if(!("numax0_flag" %in% names(d.summary)))
    #         d.summary$numax0_flag <- FALSE
    #     if(d.summary$numax0_flag == FALSE) {
    #         print(paste("Initial estimation of numax for KIC", basename(getwd())))
    #     } else {
    #         print(paste("Attempting again the estimation of numax for KIC", basename(getwd())))
    #     }

    #     ## Approximate the power spectrum background using a moving median filter in log-space
    #     ## ===================================================================================

    #     x0 <- log10(pds$frequency[1])
    #     # Corrective factor to convert median to mean
    #     # for chi-squared 2 d.o.f distributed data
    #     corr_factor <- (8.0 / 9.0)**3

    #     bkg <- c(1:length(pds$frequency))*0
    #     count <- c(1:length(pds$frequency))*0
    #     filter_width <- argv$filterwidth
    #     while (x0 < log10(pds$frequency[length(pds$frequency)])){
    #         m = abs(log10(pds$frequency) - x0) < filter_width
    #         if (length(bkg[m] > 0)){
    #             bkg[m] <- bkg[m] + median(pds$power[m], na.rm=TRUE) / corr_factor
    #             count[m] <- count[m] + 1
    #         }
    #         x0 <- x0 + 0.5 * filter_width
    #     }
    #     bkg <- bkg / count

    #     pds$power_bgr <- pds$power / bkg


    #     #X11()
    #     #print(res)
    #     #pds$bkg <- bkg
    #     #h <- ggplot() + 
    #     #     geom_line(data = pds, aes(frequency, power_bgr), color='black')
    #     #h <- ggplot() + 
    #     #     geom_line(data=pds, aes(frequency, power), color='black') + 
    #     #     geom_line(data=pds, aes(frequency, bkg), color='red') + 
    #     #     scale_x_continuous(trans='log10') + 
    #     #     scale_y_continuous(trans='log10') # + coord_cartesian(x)
    #     #plot(h)
    #     #prompt  <- "hit spacebar to close plots"
    #     #extra   <- "some extra comment"
    #     #capture <- tk_messageBox(message = prompt, detail = extra) 
    #     #stop()

    #     pds.SS <- splus2R::signalSeries(data = pds$power_bgr, positions. = pds$frequency)

    #     ## Estimating numax from the variance of the time-series
    #     ## ====================================================

    #     ## The variance and numax are related by the empirical relation
    #     ## numax = exp(fit.coef[1])*var^(fit.coef[2])

    #     fit.coef <- c(13.1981739, -0.7486405)
    #     d.summary$numax_var <- exp(fit.coef[1])*d.summary$var^(fit.coef[2])

    #     #if(d.summary$numax_var < 0.9*max(pds$frequency)) {
    #     #    pds.SS <-
    #     #        pds %>%
    #     #        filter(frequency > 0.2*d.summary$numax_var) %>%
    #     #        splus2R::signalSeries(data       = .$power,
    #     #                              positions. = .$frequency)
    #     #}

    #     ## Estimating numax with a tree-map of the local maxima of a CWT
    #     ## =============================================================

    #     ## We take a mexican-hat CWT, find the 5 most significant peaks and
    #     ## take the median as an estimation of numax.

    #     max_scale <- sqrt(tail(pds$frequency, 1))/2
    #     pds.CWT.MH <- wavCWT(x = pds.SS, wavelet = "gaussian2",
    #                         scale.range = c(deltat(pds.SS), max_scale))
    #     pds.CWTMexHat <- wavCWTTree(pds.CWT.MH)
    #     pds.Peaks <- cwtPeaks(pds.CWTMexHat)

    #     if(dim(pds.Peaks)[1] == 0) {
    #         d.summary$numax_CWTMexHat <- FALSE
    #     } else {
    #         pds.Peaks <- pds.Peaks[order(-pds.Peaks$snr)[1:numpeaks],]
    #         d.summary$numax_CWTMexHat <- median(pds.Peaks[,"frequency"], na.rm=TRUE)
    #     }


    #     ## Estimating numax from a Morlet-CWT
    #     ## ==================================

    #     pds.CWT.Morlet <- wavCWT(pds.SS, wavelet = "morlet",
    #                             scale.range = c(deltat(pds.SS), max_scale))
    #     pds.CWT.Morlet.M <- Mod(as.matrix(pds.CWT.Morlet))

    #     ## Removing edge effects
    #     nuNyq <- tail(attr(pds.CWT.Morlet, "time"), 1)
    #     time.M <- matrix(rep(attr(pds.CWT.Morlet, "time"),
    #                         times = attr(pds.CWT.Morlet, "n.scale")),
    #                     nrow = attr(pds.CWT.Morlet, "n.sample"))
    #     scale.M <- matrix(rep(attr(pds.CWT.Morlet, "scale"),
    #                         times = attr(pds.CWT.Morlet, "n.sample")),
    #                     nrow = attr(pds.CWT.Morlet, "n.scale"))

    #     aux <- matrix(data = TRUE, nrow = attr(pds.CWT.Morlet, "n.sample"),
    #                 ncol = attr(pds.CWT.Morlet, "n.scale"))
    #     aux[t(scale.M) > (time.M - 10)/5] <- FALSE
    #     aux[t(scale.M) > (nuNyq - time.M)/5] <- FALSE
    #     aux[t(scale.M) > attr(pds.CWT.Morlet, "scale")[3]] <- FALSE
    #     pds.CWT.Morlet.M <- pds.CWT.Morlet.M * aux

    #     ## Supppose numax is the maximum value
    #     d.summary$numax_Morlet <- attr(pds.CWT.Morlet,
    #                                 "time")[which(pds.CWT.Morlet.M == max(pds.CWT.Morlet.M),
    #                                                 arr.ind= TRUE)[1]]

    #     ## Estimate the final numax from the previous results
    #     ## ==================================================

    #     if (abs(1 - d.summary$numax_CWTMexHat/d.summary$numax_Morlet) < 0.2 |
    #         abs(1 - d.summary$numax_var/d.summary$numax_Morlet) < 0.2)
    #     {
    #         d.summary$numax0 <-d.summary$numax_Morlet
    #         d.summary$numax0_flag <- FALSE
    #     } else if (abs(1 - d.summary$numax_CWTMexHat/d.summary$numax_var) < 0.2) {
    #         d.summary$numax0 <- d.summary$numax_CWTMexHat
    #         d.summary$numax0_flag <- FALSE
    #     } else {
    #         # 29/06/2020 Does this even work for MS stars? var estimate will
    #         # always be dodgy!
    #         d.summary$numax0_flag <- TRUE
    #         #27.10.2021 after a visual check with Nathalie; numax_CWTMexHat seems much better!!!
    #         d.summary$numax0 <- d.summary$numax_CWTMexHat
    #     #   if (d.summary$numax_var < d.summary$nuNyq) {
    #     #       d.summary$numax0 <- d.summary$numax_var
    #     #   } else {
    #     #       d.summary$numax0 <- d.summary$numax_Morlet
    #     #   }
    #     }

    #     ## Save the results
    #     ## ================

    #     write_csv(d.summary, argv$summary)

    #     ## Warning if the estimate is unreliable
    #     if (d.summary$numax0_flag) {
    #         write.table(data.frame(), col.names=FALSE,
    #                     file = file.path(dirname(argv$summary), "NUMAX0_FLAG"))
    #     }
    # }

    data$numax0 = 1.0
    return(data)

}
