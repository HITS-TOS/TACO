library(dplyr, quietly = TRUE)
library(broom, quietly = TRUE)
library(tcltk)
library(ggplot2, quietly = TRUE)

#' For each observed_peak Get the distance squared to the closest theoretical_peak
peaks_closest <- function(observed_peaks, theoretical_peaks) {
    do.call(
        bind_rows,
        Map(f = function(i) {
            peak <- observed_peaks[i,]
            d2 <-
                theoretical_peaks %>%
                mutate(d2 = (frequency - peak$frequency)^2) %>%
                arrange(d2) %>%
                slice(1)
            return(peak %>% mutate(d2 = d2$d2))
        }, seq_len(nrow(observed_peaks))))
}

#' Reference: Hekker et al. 2011
#' http://dx.doi.org/10.1051/0004-6361/201016303
DeltaNu_from_numax <- function(numax) {
    a <- 0.260
    b <- 0.761
    return(a * numax^b)
}

eps_p_from_Dnu <- function(Dnu){
    return(0.634 + 0.546*log10(Dnu))
}

#' Make a linear fit lm(frequency ~ n) and approximate DeltaNu as the slope
#' If .l is given, we select peaks having l == .l
DeltaNu_from_peaks_lm <- function(peaks, .l = NULL, numax = NULL, alpha = TRUE) {
    if(!is.null(.l))
        peaks <- peaks %>% filter(l == .l)

    peaks$weights <- 1/(peaks$frequency_sd^2)

    if((alpha == TRUE) & (!is.null(numax))) {
        # Approximate nmax
        n_max <- numax / mean(diff(peaks$frequency)) - eps_p_from_Dnu(DeltaNu_from_numax(numax))
        peaks$quad <- (peaks$n - n_max)^2 
  
        res <-
            peaks %>%
            # 11/02/2020 added in uncertainties as weights
            lm(frequency ~ n + quad, data = ., weights=weights) %>%
            tidy() %>%
            select(c(estimate, std.error))

        # Remove weights columns
        peaks <- peaks %>% select(-weights)
        # Delta Nu and uncertainty
        Dnu <- list(DeltaNu=res$estimate[2], DeltaNu_sd=res$std.error[2])
        # Compute epsilon_p and error estimate
        c <- list(c=res$estimate[1], c_sd=res$std.error[1])    
        eps <- c$c / Dnu$DeltaNu
        eps_sd <- eps * sqrt((c$c_sd/c$c)^2 + (Dnu$DeltaNu_sd/Dnu$DeltaNu)^2)
        eps <- list(eps=eps, eps_sd=eps_sd)
        # Compute alpha and error estimate
        quad <- list(quad=res$estimate[3], quad_sd=res$std.error[3])
        alpha <- (2*quad$quad)/Dnu$DeltaNu
        alpha_sd <- alpha * sqrt((quad$quad_sd/quad$quad)^2 + (Dnu$DeltaNu_sd/Dnu$DeltaNu)^2)
        alpha <- list(alpha=alpha, alpha_sd=alpha_sd)

#        X11()
#        n <- seq(from = min(peaks$n)-1, to = max(peaks$n)+1, by = 0.1)
#        plot(peaks$n, peaks$frequency)
#        lines(n, (n + eps$eps + (alpha$alpha/2)*(n-n_max)^2)*Dnu$DeltaNu)
#       prompt  <- "hit spacebar to close plots"
#       extra   <- "some extra comment"
#       capture <- tk_messageBox(message = prompt, detail = extra)

        # Return list with results
        res <- list(Dnu=Dnu, eps=eps, alpha=alpha)
#        print(res)
#        stop()
#        return(res)
    } else {
        res <-
            peaks %>%
            # 11/02/2020 added in uncertainties as weights
                lm(frequency ~ n, data = ., weights=weights) %>%
            tidy() %>%
            #filter(term == "n") %>%
            select(c(estimate, std.error))# %>%
            #as.numeric()
        #print(res)
        # Remove weights columns
        peaks <- peaks %>% select(-weights)
        # Delta Nu and uncertainty
        Dnu <- list(DeltaNu=res$estimate[2], DeltaNu_sd=res$std.error[2])
        # Compute epsilon_p and error estimate
        c <- list(c=res$estimate[1], c_sd=res$std.error[1])    
        eps <- c$c / Dnu$DeltaNu
        eps_sd <- eps * sqrt((c$c_sd/c$c)^2 + (Dnu$DeltaNu_sd/Dnu$DeltaNu)^2)
        eps <- list(eps=eps, eps_sd=eps_sd)
        # Return list with results
        res <- list(Dnu=Dnu, eps=eps)
        #print(res)
        return(res)
    }
}

#' Estimate the l=0,2 small frequency separation from DeltaNu
#'
#' @details
#' References:
#' - Bedding et al. 2010 https://doi.org/10.1088%2F2041-8205%2F713%2F2%2Fl176
#' - Corsaro et al. 2012 http://iopscience.iop.org/article/10.1088/0004-637X/757/2/190/meta
#d02_from_DeltaNu <- function(DeltaNu) {
#    c <- 0.125
#    return(c*DeltaNu)
#}
#' - Huber et al. 2010 https://iopscience.iop.org/article/10.1088/0004-637X/723/2/1607/pdf
d02_from_DeltaNu <- function(DeltaNu) {
    c <- 0.121
    c1 <- 0.047
    return(c*DeltaNu + c1)
}

#' Estimate the l=0,1 small frequency separation from DeltaNu
#'
#' @details
#' References:
#' Mosser et al. 2018 (2014) https://www.aanda.org/articles/aa/pdf/2018/10/aa32777-18.pdf
d01_from_DeltaNu <- function(DeltaNu) {
    A <- 0.0553
    B <- -0.036
    return(A + B*log10(DeltaNu))
}

#' alpha parameter from DeltaNu for red giants.
#'
#' Eqs. (17-18) from Mosser et al. (2013)
#'
#' @param n_max nu_max / DeltaNu
#'
#' #' Mosser et al. (2013): http://dx.doi.org/10.1051/0004-6361/201220435
alpha_obs_from_n_max <- function(n_max) {
    if(n_max < 15) {
        return(2 * 0.038 / n_max)
    } else {
        return(2 * 0.57 / n_max^2)
    }
}

#' Calculate the heights of the l=0/2 pair as a function of numax
l02_amplitudes <- function(l0_freqs, l2_freqs, numax, HBR, sigmaEnv, DeltaNu) {
    # Set to total relative visibility for Kepler - for TESS this should be 
    # set to 2.94 and 4.09 for SONG
    vis_tot <- 3.16
    vis_l2 <- 0.58 # for Kepler, 0.46 for TESS and 1.04 for SONG
    amax <- sqrt((HBR * DeltaNu) / vis_tot)
    # Evaluate height envelope at l0_freqs
    l0_amps <- amax# * sqrt(exp(-(l0_freqs - numax)^2/(2*sigmaEnv^2)))
    # Evaluate at l2 freqs and account for relative l=2 visability
    l2_amps <- (amax * sqrt(vis_l2))# * sqrt(exp(-(l2_freqs - numax)^2/(2*sigmaEnv^2)))
    return (list(l0_amps=l0_amps, l2_amps=l2_amps))
}

l02_pattern <- function(pds, DeltaNu, numax, HBR, sigmaEnv, d02, l0_freq, linewidth, include_dipole=TRUE) {
    n_max <- numax / DeltaNu - eps_p_from_Dnu(DeltaNu)
    l0_freqs <-
        # 24/05/2021 - extended n range to -4:4 from -3:3
        l0_from_UP(
            N = -4:4 + l0_freq / DeltaNu  - eps_p_from_Dnu(DeltaNu),
            eps_p = eps_p_from_Dnu(DeltaNu_from_numax(numax)),
            alpha = alpha_obs_from_n_max(n_max),
            n_max = n_max,
            DeltaNu = DeltaNu)
    # 11/02/2020 Added in specific l=0/2 heights
    #amplitudes <- l02_amplitudes(l0_freqs, l0_freqs-d02, numax, HBR, sigmaEnv, DeltaNu)
    # 11/02/2020 This is arbitrary!
    #l2_linewidth <- linewidth * 2
    #l0_heights <- (amplitudes$l0_amps^2) / (pi * linewidth)
    #l2_heights <- (amplitudes$l2_amps^2) / (pi * linewidth)

    if(include_dipole == TRUE){
        d01 <- d01_from_DeltaNu(DeltaNu)
        # 24/05/2021 - Added l=1 into mix as well to better cope with "clean" dipole RGB stars. (SH)
        # 10/06/2021 - widened the linewidth of the dipole modes to make sure to incorporate mixed mode structure (SH)
        res <- fit_model(
            pds = pds,
            peaks = bind_rows(
                tibble(
                    frequency = l0_freqs,
                    linewidth = linewidth,
                    height    = quantile(pds$power, 0.999)),
                tibble(
                    frequency = l0_freqs + (0.5 + d01)*DeltaNu,
                    linewidth = 0.5*d02,
                    height = 0.25*quantile(pds$power, 0.999)),
                tibble(
                    frequency = l0_freqs - d02,
                    linewidth = linewidth,
                    height    = quantile(pds$power, 0.999))
                    ))
    } else {
        res <- fit_model(
            pds = pds,
            peaks = bind_rows(
                tibble(
                    frequency = l0_freqs,
                    linewidth = linewidth,
                    height    = quantile(pds$power, 0.999)),
                tibble(
                    frequency = l0_freqs - d02,
                    linewidth = linewidth,
                    height    = 0.58*quantile(pds$power, 0.999))))
    }
    #X11()
    #pds$res <- res
   # h <- ggplot() +
     #    geom_line(data=pds, aes(frequency, power), color='black') +
    #     geom_line(data=pds, aes(frequency, res), color='red')# + coord_cartesian(x)
   # plot(h)
    #plot(pds$frequency, res)
   # prompt  <- "hit spacebar to close plots"
    #extra   <- "some extra comment"
   # capture <- tk_messageBox(message = prompt, detail = extra)
    #stop()
    return(res)
}

tag_central_l02 <- function(peaks, pds, DeltaNu, d02, numax,
                            HBR, sigmaEnv,
                            linewidth) {
    deltanu <- abs(diff(pds$frequency[1:2]))
    pds <-
        pds %>%
        mutate(
            l02_pattern = l02_pattern(
                pds     = .,
                DeltaNu = DeltaNu,
                numax   = numax,
                HBR     = HBR,
                sigmaEnv = sigmaEnv,
                d02     = d02,
                l0_freq = numax,
                linewidth = linewidth,
                include_dipole = TRUE))
    # 24/05/2021 Extended Np to go slightly above +/- DeltaNu from numax
    # Deals with issue of having l=1 at numax.
    Np <-
        pds %>%
        filter(frequency > numax - 1.1*DeltaNu,
               frequency < numax + 1.1*DeltaNu) %>%
        nrow()
    # TODO: CHECK MAX LAG: EFFECTIVELY NUMAX+/-DELTANU/2, BUT NEEDS TO BE A BIT MORE
    pds.ccf <-
        ccf(x = pds$power,
            y = pds$l02_pattern,
            lag.max = Np/4, plot = FALSE)

    # X11()
    # print(d02)
    # plot(pds$frequency, pds$l02_pattern)
    # lines(pds$frequency, pds$l02_pattern)

    # prompt  <- "hit spacebar to close plots"
    # extra   <- "some extra comment"
    # capture <- tk_messageBox(message = prompt, detail = extra)


#     X11()
#    # # print(d02)
#     plot(numax + deltanu * pds.ccf$lag, pds.ccf$acf )
#     lines(numax + deltanu * pds.ccf$lag, pds.ccf$acf )
#     abline(v=243.5)

#     prompt  <- "hit spacebar to close plots"
#     extra   <- "some extra comment"
#     capture <- tk_messageBox(message = prompt, detail = extra)

    # 12/01/2021
    #pds.ccf$acf <- pds.ccf$acf * exp(-((deltanu * pds.ccf$lag)/DeltaNu)^2/(2*(d02/DeltaNu/2)^2))


    central_l0_expected <-
        numax + deltanu * pds.ccf$lag[which.max(pds.ccf$acf)]
    print(central_l0_expected)
  
    continue_c <- TRUE
    count <- 0
    while(continue_c){
        central_l0 <-
            peaks %>%
            filter(!(is.na(linewidth)) & ((amplitude/mean(amplitude,na.rm = TRUE)) > 1.) & ((amplitude/(mean(amplitude,na.rm = TRUE)+ 3.0*mad(amplitude,na.rm = TRUE))) < 1.)) %>%
            mutate(l0_dist = abs(frequency - central_l0_expected)) %>%
            arrange(l0_dist) %>%
            slice(1) %>%
            select(-l0_dist) %>%
            mutate(l = 0)
    # 31/05/2021 moved the assignment of the order up and included the correct one for both of the l=0 and l=2 mode (SH)
    # 10/08/2020 Changed from round to floor!!!!!
         n_order <- floor(central_l0$frequency/DeltaNu - eps_p_from_Dnu(DeltaNu))
    #n_order <- round(central_l0$frequency/DeltaNu - eps_p_from_Dnu(DeltaNu))
   
        central_l0 <-
            peaks %>%
            filter(!(is.na(linewidth)), !is.na(linewidth_sd), linewidth > 0.6*deltanu,((amplitude/mean(amplitude, na.rm = TRUE)) > 1.),((amplitude/(mean(amplitude,na.rm = TRUE)+ 3.0*mad(amplitude,na.rm = TRUE))) < 1.))  %>%
            mutate(l0_dist = abs(frequency - central_l0_expected)) %>%
            arrange(l0_dist) %>%
            slice(1) %>%
            select(-l0_dist) %>%
            mutate(l = 0) %>%
            mutate(n = n_order)
        print(central_l0)
        print(central_l0$frequency - 1.7*d02)
        print(central_l0$frequency - 0.5*d02)
        
        central_l2 <-
            peaks %>%
            filter(!is.na(linewidth), !is.na(linewidth_sd), linewidth > 0.3*deltanu,
                frequency > central_l0$frequency - 1.7*d02,
                frequency < central_l0$frequency - 0.5*d02) %>%
        # 03/03/2020 change so that relative to expected l=0 frequency,
        # not found peak, in case only l=2 is present.
                #frequency > central_l0_expected - 1.5*d02,
                #frequency < central_l0_expected - 0.5*d02) %>%
                arrange(-amplitude) %>%
                #slice(1) %>%
                mutate(l = 2) %>%
                mutate(n = n_order-1)
       print(nrow(central_l2))
       if(nrow(central_l2) > 0){
           continue_c <- FALSE
       } else if (nrow(central_l2) == 0){
           central_l0_expected <- central_l0_expected + DeltaNu/2.
       }
       if(count == 2) {
           continue_c <- FALSE   #to make sure this is not an endless loop...
       }
       count <- count + 1
    }
    central_l02 <-
        bind_rows(central_l2, central_l0)
        peaks <-
        bind_rows(
            central_l02,
            peaks %>% filter(!(frequency %in% central_l02$frequency))
        ) %>%
        arrange(frequency)
    return(peaks)
}

#' @param sign Whether or not we are searching above or below numax.
tag_l02_pair <- function(peaks, pds, DeltaNu, d02, alpha, search.range, current_radial_order, central_radial_order, sign) {
   
    # Bin width
    deltanu <- abs(diff(pds$frequency[1:2]))
    
    # Define current l=0 and l=2 modes
    current_l0 <-
        peaks %>%
        filter((l==0) & (n == current_radial_order)) %>%
        arrange(-amplitude) %>%
        slice(1)
    current_l2 <-
        peaks %>%
        filter((l==2) & (n == current_radial_order-1)) %>%
        arrange(-amplitude) %>%
        slice(1)
    l2mnfreq <-weighted.mean(current_l2$frequency,current_l2$amplitude,.,na.rm=TRUE)
  #   if(sign < 0){
  #       print("PREVIOUS BIT ===================")
  #       print(current_l0)
  #       print(current_l2)
  #   }
    
    # If have no current_l0 but have l2. This can occur when there is a radial order with a 
    # detected l=2, but no detected l=0
    if((nrow(current_l0) == 0) & (nrow(current_l2) > 0)){
        current_l0 <- current_l2 %>% mutate(frequency = l2mnfreq + d02, l=0)
        current_l0 <- current_l0 %>% mutate(n = current_radial_order, l=0)
    }
    else if((nrow(current_l0) > 0) & (nrow(current_l2) == 0)){
        # If have l0 but no l2
        current_l2 <- current_l0 %>% mutate(frequency = frequency - d02, l=2)
        current_l2 <- current_l2 %>% mutate(n = current_radial_order-1, l=2)
    } 
    else if((nrow(current_l0) == 0) & (nrow(current_l2) == 0)){
        # If have no l=0 or l=2 then just return peaks with no mode ID performed
        # in this radial order
        return (peaks)
    }
    # Combine into single tibble
    current_l02 <- bind_rows(current_l0, current_l2)
    print(current_l02)
    #print(l2mnfreq)
    
    # Estimate d02 from frequencies
    #d02 <- abs(diff(current_l02$frequency))
    #print(alpha)
    #print(DeltaNu)
    if(sign > 0) predictedl0 <- (current_l0$frequency + (DeltaNu * (1.0 + (alpha) * (current_radial_order + 1 - central_radial_order))))
    if(sign < 0) predictedl0 <- (current_l0$frequency - (DeltaNu * (1.0 + (alpha) * (current_radial_order - 1 - central_radial_order))))
    if(sign > 0) predictedl2 <- (l2mnfreq + (DeltaNu * (1.0 + (alpha) * (current_radial_order  - central_radial_order))))
    if(sign < 0) predictedl2 <- (l2mnfreq - (DeltaNu * (1.0 + (alpha) * (current_radial_order  - 2 - central_radial_order))))
    
    print(predictedl2)
    print(predictedl0)
    
    closest_peak_info <- peaks %>%
                            filter(!is.na(linewidth), !is.na(linewidth_sd), linewidth > 0.3*deltanu, # This line might be a bit suspect as will fail when get to shorter datasets! Maybe better to take widest peak e.g.
                                if(sign > 0) frequency > current_l0$frequency + search.range[1]*(DeltaNu * (1.0 + (alpha) * (current_radial_order + 1 - central_radial_order)))
                                else frequency > current_l0$frequency - search.range[2]*(DeltaNu * (1.0 + (alpha) * (current_radial_order - 1 - central_radial_order))),
                                if(sign > 0) frequency < current_l0$frequency + search.range[2]*(DeltaNu * (1.0 + (alpha) * (current_radial_order + 1 - central_radial_order)))
                                else frequency < current_l0$frequency - search.range[1]*(DeltaNu * (1.0 + (alpha) * (current_radial_order - 1 - central_radial_order))))
                                
   if(nrow(closest_peak_info) < 1) {
       closest_peak <- closest_peak_info
       }
   if(nrow(closest_peak_info) >= 1) {
       closest_peak <- closest_peak_info %>%
       filter(((amplitude / max(amplitude)) > 0.55)) #%>%
#        arrange(-amplitude) %>%
 #       slice(1)
       }
#    if(nrow(closest_peak_info) >= 3) {
#        min_dist0 <- abs(closest_peak_info$frequency-predictedl0)
#        min_dist2 <- abs(closest_peak_info$frequency-predictedl2)
 #       if ((min(min(min_dist0) - min(min_dist2))) > 0) {
#            predictedl = predictedl2
#        }
#        else {
 #           predictedl = predictedl0
 #       }
 #       closest_peak <- closest_peak_info %>%
 #                       filter(((amplitude / max(amplitude)) > 0.55)) %>%
 #                       mutate(min_dist = abs(frequency-predictedl)) %>%
 #                       #mutate(min_dist2 = abs(frequency-predictedl2)) %>%
 #                       #mutate(min_dist = min(abs(frequency-predictedl0),abs(frequency-predictedl2)) %>%
#                        arrange(-amplitude) %>%
 #                       slice(1) #%>%
 #                     #  select(-min_dist0) %>%
  #                    #  select(-min_dist2) %>%
  #                    #  select(-amplitude)
  #      }
    
    # Find closest peak
    #closest_peak <- peaks %>%
    # This line might be a bit suspect as will fail when get to shorter datasets! Maybe better to take widest peak e.g.
    #                    filter(!is.na(linewidth), !is.na(linewidth_sd),
    #                            linewidth > 0.6*deltanu,  !is.na(amplitude), !is.na(amplitude_sd),
    #                        if(sign > 0) frequency > current_l0$frequency + search.range[1]*(DeltaNu * (1.0 + (alpha) * (current_radial_order + 1 - central_radial_order)))
    #                        else frequency > current_l0$frequency - search.range[2]*(DeltaNu * (1.0 + (alpha) * (current_radial_order - 1 - central_radial_order))),
    #                        if(sign > 0) frequency < current_l0$frequency + search.range[2]*(DeltaNu * (1.0 + (alpha) * (current_radial_order + 1 - central_radial_order)))
    #                        else frequency < current_l0$frequency - search.range[1]*(DeltaNu * (1.0 + (alpha) * (current_radial_order - 1 - central_radial_order)))
    #                        ) %>%
    #                    #mutate(maxampl = max(amplitude,na.rm=TRUE)) %>%
    #                    arrange(-amplitude) %>%
    #                    slice(1)

    # if(sign > 0){
    #     print("2 NEXT BIT ===================")
    #     print(current_l0)
    #     print(current_l2)
    #     print(closest_peak)
    # }
    # If we have a closest peak
    
    if(nrow(closest_peak) >= 1){
        # Compute distances between predicted l=0,2 frequencies and closest peaks
        dist_l0 <- (closest_peak$frequency - ifelse(sign > 0,
                                                    current_l0$frequency + (DeltaNu * (1.0 + (alpha) * (current_radial_order + 1 - central_radial_order))),
                                                    current_l0$frequency - (DeltaNu * (1.0 + (alpha) * (current_radial_order - 1 - central_radial_order)))
                                                    )) / DeltaNu
        dist_l2 <- (closest_peak$frequency - ifelse(sign > 0,
                                                    current_l2$frequency + (DeltaNu * (1.0 + (alpha) * (current_radial_order - central_radial_order))),
                                                    current_l2$frequency - (DeltaNu * (1.0 + (alpha) * (current_radial_order - 2 - central_radial_order)))
                                                    )) / DeltaNu
      
        # 31/05/2021 changed the values a bit to incorporate more room for curvature
     #   if( DeltaNu < 4){
     #       bounds0 <- -0.05
     #       bounds1 <- 0.15
     #   } else{
     #       bounds0 <- -0.05
     #       bounds1 <- 0.15
     #   }
        
        
        
        # If peak is closer to l=0 than l=2, assign l=0
        if((min(abs(dist_l0)) < min(abs(dist_l2))) & nrow(closest_peak) > 1){ #}& (dist_l0 > bounds0) & (dist_l0 < bounds1)){
            l0_peak <- closest_peak[which.min(abs(dist_l0)),] %>%
                        mutate(l=0, n=ifelse(sign > 0, current_l0$n+1, current_l0$n-1))
            l2_peak <- NULL
        # Otherise assign l=2
        } else if((min(abs(dist_l2)) < min(abs(dist_l0))) & nrow(closest_peak) > 1){#} & (dist_l2 > bounds0) & (dist_l2 < bounds1)) {
            l2_peak <- closest_peak[which.min(abs(dist_l2)),] %>%
                        mutate(l=2, n=ifelse(sign > 0, current_l0$n, current_l0$n-2))
            l0_peak <- NULL
        }else if(nrow(closest_peak) == 1){#} & (dist_l2 > bounds0) & (dist_l2 < bounds1)) {
            l0_peak <- closest_peak[which.min(abs(dist_l0)),] %>%
                        mutate(l=0, n=ifelse(sign > 0, current_l0$n+1, current_l0$n-1))
            l2_peak <- NULL
        } else {
            return(peaks)
        }
    # If don't find a close peak then return 
    } else{
        return(peaks)
    }
    
    # if(sign > 0){
    #     print("New Peaks")
    #     print(current_l0)
    #     print(current_l2)
    #     print(l0_peak)
    #     print(l2_peak)
    # }

    # If peak was an l=2 then we don't have a currently assigned l=0 peak.
    # So look for a peak ~0.5-1.5*d02 above l=2 peak position.
    if(is.null(l0_peak)){
        closest_peak_info <-
            peaks %>%
            filter(!is.na(linewidth),linewidth > 0.3*deltanu,
                frequency < l2_peak$frequency + 1.5*d02,
                frequency > l2_peak$frequency + 0.5*d02)
         #   if(nrow(closest_peak_info) <= 1) {
                l0_peak <- closest_peak_info %>%
                # 03/03/2020
                #frequency > l0_loc - 1.5*d02,
                #frequency < l0_loc - 0.5*d02) %>%
                arrange(-amplitude) %>%
                slice(1) %>%
                mutate(l = 0, n = ifelse(sign > 0, current_l0$n+1, current_l0$n - 1))
    } else if(is.null(l2_peak)){
        closest_peak_info <-
            peaks %>%
            filter(!is.na(linewidth),linewidth > 0.3*deltanu,
                    #24/05/2021 decreased limits slightly
                frequency > l0_peak$frequency - 1.5*d02,
                frequency < l0_peak$frequency - 0.5*d02) #,
       #         amplitude / max(amplitude) > 0.1)
      #      if(nrow(closest_peak_info) <= 1) {
                l2_peak <- closest_peak_info %>%
                arrange(-amplitude) %>%
       #          slice(1) %>%
                 mutate(l = 2, n = ifelse(sign > 0, current_l0$n, current_l0$n-2))
       #     }
       #    if(nrow(closest_peak_info) >= 2) {
       #     l2_peak <- closest_peak_info %>%
        #                    filter(((amplitude / max(amplitude)) > 0.25)) %>%
        #                    mutate(min_dist2 = abs(frequency-predictedl2)) %>%
        #                    arrange(min_dist2) %>%
         #                   slice(1) %>%
         #                   mutate(l = 2, n = ifelse(sign > 0, current_l0$n, current_l0$n-2)) #%>%
         #                   select(-min_dist2)
        #    }
    }

    # if(sign > 0){
    #     print("Newer Peaks")
    #     print(l0_peak$frequency)
    #     print(current_l0$frequency)
    #     print(l2_peak)
    #     stop()
    # }

 #   if(nrow(l0_peak) > 0)
 #       if(l0_peak$frequency == current_l0$frequency)
 #           return(peaks)
 #   if(nrow(l2_peak) > 0)
 #       if (l2_peak$frequency == current_l2$frequency)
 #           return(peaks)

    l02_pair <-
        bind_rows(l0_peak, l2_peak)
   
    res <-
        bind_rows(
            l02_pair,
            peaks %>% filter(!(frequency %in% l02_pair$frequency))
        ) %>%
        arrange(frequency)
    return(res)
}

#' Tag the l=0 and l=2 modes from the set of detected peaks.
#'
#' @param peaks A tibble containing the detected peaks.
#' @param pds A tibble containing the power spectrum.
#' @param DeltaNu The first guess delta nu which the tagging uses first.
#' @param d02 The small frequency separation to use in the tagging.
#' @param numax The frequency of maximum oscillation power.
#' @param HBR The height-background ratio at numax.
#' @param sigmaEnv The width of the oscillation envelope.
#' @param nuNyq The nyquist frequency of the data.
#' @param search.range Look for the next l=0,2 modes using this range away from Δν.
tag_l02_peaks <- function(peaks, pds, DeltaNu, d02, alpha, numax, HBR, sigmaEnv, nuNyq, search.range) {
    peaks <-
        peaks %>%
        tag_central_l02(
            pds = pds, DeltaNu = DeltaNu, d02 = d02, numax = numax,
            HBR = HBR, sigmaEnv = sigmaEnv,
            linewidth = search_l02_linewidth(numax, nuNyq))
    n.l02 <-
        peaks %>%
        filter(l == 0 | l == 2) %>%
        nrow()
    
    continue_p <- TRUE
    DeltaNu_lower <- DeltaNu
    DeltaNu_upper <- DeltaNu

    #print("Starting DNU")

    l0 <- peaks %>%
            filter(l == 0) %>%
            arrange(-amplitude) %>%
            slice(1)
    l2 <- peaks %>%
            filter(l == 2) #%>%
            #arrange(-amplitude) %>%
            #slice(1)
    l2mnfreq <-weighted.mean(l2$frequency,l2$amplitude,.,na.rm=TRUE)
    #print(l2mnfreq)
            
    print("CENTRAL PEAKS")
    print(l0)
    print(l2)
    print("=======================")
            
    if ((nrow(l2) > 0) & (nrow(l0) > 0)) {
        d02 <- l0$frequency-l2mnfreq
    }
    central_l0 <- max(l0$frequency)
    radial_order_high <- max(l0$n)
    radial_order_low <- max(l0$n)
    central_radial_order <- max(l0$n)
    # This is designed to improve estimate of delta nu upper and lower.
    # Accounts for the fact that if central l0 is above numax then will 
    # calculated delta nu lower wrong (have one too many radial orders),
    # and vice versa for delta nu upper.
    # if(l0$frequency > numax){
    #     peak_lower_max_radial_order <- max(l0$n) - 1
    #     peak_upper_max_radial_order <- max(l0$n)
    # } else{
    #     peak_lower_max_radial_order <- max(l0$n)
    #     peak_upper_max_radial_order <- max(l0$n) + 1
    # }

    first_attempt <- TRUE

    while(continue_p) {
        # Tag l=0,2 pair above and below current radial order using DeltaNu estimated from previously tagged
        # modes.
        # 13/05/2021
        # Bear in mind that if no l=0 is returned then it will stall because it will add DeltaNu onto the last l=0 it found and won't move anywhere. So try to
        # account for this in the tagging functions.
        # Need to give radial order as well so that it knows where it should be if no l=0 is found.
        print("CURRENT RADIAL ORDER")
        print(radial_order_high)
        print(radial_order_low)
        #print("CENTRAL RADIAL ORDER")
        #print(central_radial_order)
        peaks <-
            peaks %>%
            tag_l02_pair(pds = pds, DeltaNu = DeltaNu, d02 = d02, alpha = alpha, search.range = search.range, current_radial_order=radial_order_high, central_radial_order=central_radial_order, sign=+1) %>%
            tag_l02_pair(pds = pds, DeltaNu = DeltaNu, d02 = d02, alpha = alpha, search.range = search.range, current_radial_order=radial_order_low, central_radial_order=central_radial_order, sign=-1)
        
        radial_order_high <- radial_order_high + 1
        radial_order_low <- radial_order_low - 1
        # Take delta nu as average of differences of l=0 below and above central l=0 peak
        #peaks_lower <- peaks %>%
        #                filter((l==0) & (frequency <= central_l0))
       # if(nrow(peaks_lower) > 1){
       #     DeltaNu_lower <- sum(diff(peaks_lower$frequency)) / (max(peaks_lower$n) - min(peaks_lower$n))
            #print("DeltaNu")
            #print(DeltaNu)
            #stop()
       # }
       
        
        #peaks_upper <- peaks %>%
        #                filter((l==0) & (frequency >= central_l0))
        #if(nrow(peaks_upper) > 1){
        #    DeltaNu_upper <- sum(diff(peaks_upper$frequency)) / (max(peaks_upper$n) - max(l0$n))
            #print("DeltaNu")
            #print(DeltaNu)
            #print("DeltaNu_lower")
            #print(DeltaNu_lower)
            #stop()
        #}
        #print("lower")
        #print(peaks_lower)
        #print(DeltaNu_lower)
        #print("upper")
        #print(peaks_upper)
        #print(DeltaNu_upper)
        #stop()
        # Compute DeltaNu from currently detected radial modes. - THIS IS A BAD IDEA AS CURVATURE MESSES WITH THIS!                        
        #res <- DeltaNu_from_peaks_lm(peaks, .l = 0)
        #DeltaNu <- res$Dnu$DeltaNu
        
        # Set of tagged peaks, included the most recently tagged pairs.
        n.l02.new <-
            peaks %>%
            filter(l == 0 | l == 2) %>%
            nrow()
#        print(peaks %>%
#            filter(l == 0 | l == 2))
        #print("NUMBER OF NEW TAGGED MODES")
        #print(n.l02.new)
        #print("NUMBER OF OLD TAGGED MODES")
        #print(n.l02)
        # If no new modes are added then stop.
        if(n.l02.new == n.l02){
            if(first_attempt == FALSE){
                #print("EXITING")
                continue_p <- FALSE
            }
            first_attempt <- FALSE
            #print("FIRST ATTEMPT FAILURE")
        }
        # If new modes have been tagged then make a note and keep going.
        n.l02 <-
            peaks %>%
            filter(l == 0 | l == 2) %>%
            nrow()
    }
    #print("DONE TAGGING")
    return(peaks)
}

#' Frequency estimations from the universal pattern for l=0 modes
#'
#' We can estimate the l=0 frequencies of a given radial order (N)
#' using the large frequency separation (DeltaNu), the asymptotic
#' offset (eps_p) and n_max = numax/DeltaNu.
#'
#' @param N integer radial order of the l=0 mode
#' @param eps_p Asymptotic offset
#' @param alpha_ magnitude of the second order term (curvature in the
#' echelle diagram)
#' @param n_max numax/DeltaNu
#' @param DeltaNu large frequency separation
#'
#' @details
#' Reference?
l0_from_UP <- function(N, eps_p, alpha, n_max, DeltaNu) {
    return((N + eps_p + ((alpha/2.0) * (N - n_max)^2)) * DeltaNu)
}
l2_from_UP <- function(N, eps_p, alpha, n_max, DeltaNu, d02) {
    return(((N + eps_p + 1.0 + ((alpha/2.0) * (N - n_max)^2)) * DeltaNu) - d02)
}
search_l02_linewidth <- function(numax, nuNyq) {
    if(numax < 30) {
        Q1 <- 0.02
        Q2 <- 0.075
        n1 <- 3
        n2 <- 30
    } else {
        Q1 <- 0.075
        Q2 <- 0.035
        n1 <- 30
        n2 <- nuNyq
    }
    a <- log(Q2/Q1) / log(n2/n1)
    return(Q1 * (numax/n1)^a)
}

#' Fit for the observed DeltaNu taking into account the curvature.
#'
#' MLE fit of the observed frequencies with the predicted frequencies
#' as a function of DeltaNu, epsilonp and alpha for a given spherical
#' degree. We obtain the best combination for the observed frequencies
#'
#' @param peaks A data_frame with at least 3 columns: frequency, l
#' and n. We only use the l=0 modes.
#' @param numax A single number with the numax value
#' @param DeltaNu0 Initial guess for DeltaNu
#' @param weight.by.amplitudes Should the square differences by weighted
#' by the peak amplitudes?
#'
#' @return A list of numbers with elements `DeltaNu`, `epsilonp` and `alpha`
#'
#' @details
#' This fits Eq. (9) of Mosser et al. (2013) for the observed parameters.
#'
#' Initial values for:
#' - DeltaNu: Provided DeltaNu0
#' - epsilonp: set as 1.25 as representative value taken from
#' Mosser et al. (2013)
#' - alpha: From Mosser et al. (2013) Eq. (17). It uses
#' the `alpha_obs_from_n_max` function
#'
#' Mosser et al. (2013): http://dx.doi.org/10.1051/0004-6361/201220435
DeltaNu_l0_fit <- function(peaks, numax, DeltaNu0,alpha0,
                          weight.by.amplitudes = FALSE,
                          return_res = FALSE) {
    #alpha0 <- alpha_obs_from_n_max(n_max = (numax/DeltaNu0))
    #print(peaks)
    #print(5*alpha0)
    l0_peaks <-
        peaks %>%
        arrange(frequency) %>%
        filter(l==0)
    if(nrow(l0_peaks) == 0)
        stop(paste("'peaks' does not have l=0 modes"))
    res <-
        optim(
            # I use theta = (DeltaNu, epsilonp, alpha)
            par = c(DeltaNu0, eps_p_from_Dnu(DeltaNu0), alpha0),   # initial values
            fn = function(theta) {
                n_max <- numax/theta[1] - theta[2] # Using updated expression from Mosser et al. (2018)
                pks <-
                    l0_peaks %>%
                    mutate(predFreq =
                               l0_from_UP(
                                   N       = .$n,
                                   eps_p   = theta[2],
                                   alpha   = theta[3],
                                   n_max   = n_max,
                                   DeltaNu = theta[1]),
                           diffsq = (frequency - predFreq)^2/(frequency_sd^2))
                if(weight.by.amplitudes) {
                    pks <-
                        pks %>% mutate(diffsq = amplitude * diffsq)
                }
                res <-
                    pks %>%
                    summarise(sqdiff_sum = sum(diffsq)) %>%
                    as.numeric()
                return(res)
            },
            control = list(
                parscale = c(1e-2, 1e-2, 1e-3)
            ),
            method = "L-BFGS-B",
            ## Limits on:  DeltaNu, epsilonp,      alpha)
            lower = c(0.8*DeltaNu0,      0.4, 0.1*alpha0),
            upper = c(1.2*DeltaNu0,      2.8,   10*alpha0),
            hessian = TRUE
        )
    if(return_res == TRUE){
        return(res)
    }
    sd <- sqrt(diag(solve(res$hessian)))

    #if(res$convergence != 0) stop("Error (DeltaNu_l_fit): optim didn't converge")
    return(
        list(
            DeltaNu  = res$par[1],
            DeltaNu_sd = sd[1],
            eps_p = res$par[2],
            eps_p_sd = sd[2],
            alpha    = res$par[3],
            alpha_sd = sd[3],
            message  = res$message))
}

DeltaNu_l2_fit <- function(peaks, numax, DeltaNu0, alpha0, eps_p0, d020,
                          weight.by.amplitudes = FALSE,
                          return_res = FALSE) {
    #alpha0 <- alpha_obs_from_n_max(n_max = (numax/DeltaNu0))
    #print(5*alpha0)
    print("d020")
    print(d020)
    l2_peaks <-
        peaks %>%
        arrange(frequency) %>%
        filter(l==2)
    
    if(nrow(l2_peaks) == 0)
        print("'peaks' does not have l=2 modes")
        list(
             d02      = 0.0,
             d02_sd   = 0.0,
             message  = '')
             
    if(nrow(l2_peaks) > 0) {
        tmp_l2_0 <- l2_peaks %>%
            filter(n==min(l2_peaks$n))
        l2mnfreq <- weighted.mean(tmp_l2_0$frequency,tmp_l2_0$amplitude,.,na.rm=TRUE)
        l2totampl <- sum(tmp_l2_0$amplitude)
        l2_unique <- tmp_l2_0 %>%
            mutate(frequency = l2mnfreq, amplitude = l2totampl) %>%
            slice(1)

        for (i in unique(l2_peaks$n)){
            tmp_l2 <- l2_peaks %>%
                filter(n==i)
            l2mnfreq <- weighted.mean(tmp_l2$frequency,tmp_l2$amplitude,.,na.rm=TRUE)
            l2totampl <- sum(tmp_l2$amplitude)
            tmp_l2 <- tmp_l2 %>%
                mutate(frequency = l2mnfreq, amplitude = l2totampl) %>%
                slice(1)
       
            if (i > min(l2_peaks$n)){
                l2_unique <-
                    bind_rows(l2_unique,tmp_l2)
             }
        }
    
        l2_peaks <- l2_unique
        res2 <-
            optim(
                # I use theta = (DeltaNu, epsilonp, alpha)
                par = c(d020),   # initial values
                fn = function(theta) {
                    n_max <- numax/DeltaNu0 - eps_p0 # Using updated expression from Mosser et al. (2018)
                    pks <-
                        l2_peaks %>%
                        mutate(predFreq =
                                  l2_from_UP(
                                  N       = .$n,
                                  eps_p   = eps_p0,
                                  alpha   = alpha0,
                                  n_max   = n_max,
                                  DeltaNu = DeltaNu0,
                                  d02     = theta[1]),
                           diffsq = (frequency - predFreq)^2)
                    
                    if(weight.by.amplitudes) {
                        pks <-
                            pks %>% mutate(diffsq = amplitude * diffsq)
                    }
                   
                    res2 <-
                        pks %>%
                        summarise(sqdiff_sum = sum(diffsq)) %>%
                        as.numeric()
                return(res2)
            },
            control = list(
                    parscale = c(1e-3)
                ),
                method = "L-BFGS-B",
                ## Limits on:  DeltaNu, epsilonp,      alpha)
                lower = c(0.5*d020),
                upper = c(2.0*d020),
                hessian = TRUE
                )
                if(return_res == TRUE){
                    return(res2)
            }
        sd <- sqrt(diag(solve(res2$hessian)))
        #if(res$convergence != 0) stop("Error (DeltaNu_l_fit): optim didn't converge")
    }
   
    return(
        list(
            d02      = res2$par[1],
            d02_sd   = sd[1],
            message  = res2$message))
}

radial_order_shift_from_l0 <- function(peaks, numax, DeltaNu, epsilonp, alphaobs) {
    n_max <- numax/DeltaNu - epsilonp
    shifts <-
        do.call(
            bind_rows,
            Map(function(i) {
                peaks %>%
                    filter(l == 0) %>%
                    mutate(n = n + i) %>%
                    mutate(predFreq =
                               l0_from_UP(
                                   N = .$n,
                                   eps_p = epsilonp,
                                   alpha_ = alphaobs,
                                   n_max = n_max,
                                   DeltaNu = DeltaNu
                               )) %>%
                    mutate(diffsq = (frequency - predFreq)^2) %>%
                    summarise(dsq = sum(diffsq)) %>%
                    bind_cols(nshift = i)
            }, -5:5))
    return(shifts$nshift[which.min(shifts$dsq)])
}

#' Frequency estimations from the universal pattern for l=0 modes
#'
#' We can estimate the l=0 frequencies of a given radial order (N)
#' using the large frequency separation (DeltaNu) and the frequency
#' of maximum oscillation power (numax). The asymptotic
#' offset (eps_p) and second order term magnitude (alpha) are
#' estimated from Dnu.
#'
#' @param numax frequency of maximum oscillation power
#' @param Dnu large frequency separation
#' @param num_freqs number of frequencies around numax to return
l0_theoretical <- function(numax, DeltaNu, num_freqs, eps_p = NULL, alpha_ = NULL) {
    if(is.null(eps_p))
        eps_p <- eps_p_from_Dnu(DeltaNu)
    if(is.null(alpha_))
        alpha_ <- alpha_obs_from_n_max(n_max = numax / DeltaNu - eps_p)
    n_max <- (numax / DeltaNu) - eps_p
    N <- seq(from = floor(n_max - (num_freqs / 2)),
             length.out = num_freqs,
             by   = 1)
    return(l0_from_UP(
        N = N,
        eps_p = eps_p,
        alpha_ = alpha_,
        n_max = n_max,
        DeltaNu = DeltaNu))
}

central_d02 <- function(peaks, numax) {
    n.central <-
        peaks %>%
        filter(l == 0) %>%
        arrange(abs(frequency - numax)) %>%
        slice(1) %>%
        select(n) %>%
        as.numeric()
    
    l0_central <-
        peaks %>%
        filter(n == n.central & l == 0)
    
    l2_central <-
            peaks %>%
            filter(n == n.central-1 & l == 2)
    
    l02_central <- bind_rows(l0_central, l2_central)
    
    d02 <- ifelse(nrow(l02_central) == 2, abs(diff(l02_central$frequency)), NA)
    return(d02)
}
