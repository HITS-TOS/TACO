## Get Δν and label the l=0,2 peaks

library(lomb, quietly = TRUE, warn.conflicts = FALSE)
library(argparser, quietly = TRUE, warn.conflicts = FALSE)
library(readr, quietly = TRUE, warn.conflicts = FALSE)
library(magrittr, quietly = TRUE, warn.conflicts = FALSE)
options(tibble.width = Inf)

#https://stackoverflow.com/questions/38351820/negation-of-in-in-r
`%notin%` <- Negate(`%in%`)

source("src/peakFind_lib.R", chdir = TRUE)
source("src/l02_modes_id.R", chdir = TRUE)

peak_bag_mode_id02_r <- function(pds, peaks, data) {

    ## Arbitrary parameters that I use
    MAX.ERROR <- 0.03    # Discard taggings where the square difference between expected and predicted frequencies is greater than MAX.ERROR*Δν
    search.range <- c(0.75, 1.25)  # Look for the next l=0,2 modes using this range away from Δν

    ## Read the data
    d.summary <- read_csv(argv$summary, col_types = cols())

    d.peaks <-
        read_csv(argv$peaks, col_types = cols()) %>%
        filter(frequency > d.summary$numax - 3*d.summary$sigmaEnv,
            frequency < d.summary$numax + 3*d.summary$sigmaEnv) %>%
        select(-one_of("n", "l", "m")) %>%
        filter(AIC > 0)
    d.pds <-
        read_csv(argv$pds, col_types = cols()) %>%
        filter(frequency > d.summary$numax - 3*d.summary$sigmaEnv,
            frequency < d.summary$numax + 3*d.summary$sigmaEnv)
    deltanu <- d.pds$frequency[2] - d.pds$frequency[1]

    ## Check to see if there are any peaks in file, if not then stop
    if(nrow(d.peaks)  == 0){
        d.peaks <- tibble(frequency=double(),
                            amplitude=double(),
                            linewidth=double(),
                            frequency_sd=double(),
                            amplitude_sd=double(),
                            linewidth_sd=double(),
                            height=double(),
                            AIC=double(),
                            n=integer(),
                            l=integer())
        write_csv(x = d.peaks, file = argv$peaks)
        # Central Δν
        d.summary <-
            d.summary %>%
            mutate(DeltaNu = NaN,
                DeltaNu_sd = NaN,
                dNu02     = NaN,
                eps_p     = NaN,
                eps_p_sd  = NaN,
                alpha     = NaN,
                alpha_sd  = NaN,
                Central_DeltaNu = NaN,
                Central_DeltaNu_sd = NaN,
                Central_eps_p     = NaN,
                Central_eps_p_sd  = NaN,
                Central_alpha     = NaN,
                Central_alpha_sd  = NaN,
                gamma0    = NULL,
                modeIDFlag = 1)
        write_csv(x = d.summary, file = argv$summary)
        stop("No peaks detected so not proceeding with mode ID")
    } else if(nrow(d.peaks) < 3){
        d.peaks$l <- NA
        d.summary <-
            d.summary %>%
            mutate(DeltaNu = NaN,
                DeltaNu_sd = NaN,
                dNu02     = NaN,
                eps_p     = NaN,
                eps_p_sd  = NaN,
                alpha     = NaN,
                alpha_sd  = NaN,
                Central_DeltaNu = NaN,
                Central_DeltaNu_sd = NaN,
                Central_eps_p     = NaN,
                Central_eps_p_sd  = NaN,
                Central_alpha     = NaN,
                Central_alpha_sd  = NaN,
                gamma0    = NULL,
                modeIDFlag = 1)
        write_csv(x = d.peaks, file = argv$peaks)
        write_csv(x = d.summary, file = argv$summary)
        stop("Not enough peaks detected so not proceeding with mode ID")
    }

    if(d.summary$numax < 5) {
        d.summary <-
            d.summary %>%
            mutate(DeltaNu = NaN,
                DeltaNu_sd = NaN,
                dNu02     = NaN,
                eps_p     = NaN,
                eps_p_sd  = NaN,
                alpha     = NaN,
                alpha_sd  = NaN,
                Central_DeltaNu = NaN,
                Central_DeltaNu_sd = NaN,
                Central_eps_p     = NaN,
                Central_eps_p_sd  = NaN,
                Central_alpha     = NaN,
                Central_alpha_sd  = NaN,
                gamma0    = NULL,
                modeIDFlag = 2)
        write_csv(x = d.summary, file = argv$summary)
        stop("Numax < 10uHz and is too low for the automated mode identification to be reliable.")
    }

    ## Expected Δν from ν_max
    d.Dnu <- DeltaNu_from_numax(d.summary$numax)

    ## Refine it using the PS of the PS. We look at Δν and Δν/2 and take the highest peak
    d.psxps1 <- lsp(times = d.pds$frequency, x = d.pds$power, 
                    from = 0.3*d.Dnu, to = 0.7*d.Dnu,
                    type = "period", plot = FALSE, ofac = 10)
    d.psxps2 <- lsp(times = d.pds$frequency, x = d.pds$power, 
                    from = 0.7*d.Dnu, to = 1.3*d.Dnu,
                    type = "period", plot = FALSE, ofac = 10)
    d.Dnu <- ifelse(d.psxps1$peak > d.psxps2$peak, d.psxps1$peak.at[1] * 2, d.psxps2$peak.at[1])

    print(paste0("Initial Dnu from PSxPS ", d.Dnu))

    # Expected δν_02 from Δν
    d.d02 <- d02_from_DeltaNu(d.Dnu)
    #print(d.d02)
    ## Refine the δν_02 with the PSxPS
    d.psxps <- lsp(times = d.pds$frequency, x = d.pds$power, 
                from = 0.7*d.d02, to = 1.3*d.d02,
                type = "period", plot = FALSE, ofac = 10)
    d.d02 <- d.psxps$peak.at[1]
    #print(d.d02)

    # Expected alpha from Δν
    n_max <- d.summary$numax / d.Dnu - eps_p_from_Dnu(d.Dnu)
    d.alpha <- 1.5*alpha_obs_from_n_max(n_max)    #overestimating alpha a bit to avoid issues at extreme frequencies
    print(paste0("Expected alpha from nmax ",d.alpha))

    ## Mode identification for l=0,2 modes
    d.peaks <-
        d.peaks %>%
        tag_l02_peaks(pds = d.pds, DeltaNu = d.Dnu, d02 = d.d02, alpha=d.alpha, numax = d.summary$numax,
                    HBR = NULL, sigmaEnv = d.summary$sigmaEnv,
                    nuNyq = d.summary$nuNyq,
                    search.range = search.range)

    d.peaks.l0 <- d.peaks %>% filter(l==0)
    dnu_est <- d.Dnu

    # Pull out l=0 and l=2 modes
    d.peaks.l0 <- d.peaks %>% filter(l==0)
    d.peaks.l2 <- d.peaks %>% filter(l==2)

    # If there are fewer than two modes exit mode ID
    if (nrow(d.peaks.l0) < 3) {

        d.peaks <- d.peaks %>%
                    mutate(n = NaN,
                        l = NaN)

        d.summary <-
                d.summary %>%
                mutate(DeltaNu = NaN,
                        DeltaNu_sd = NaN,
                        dNu02     = NaN,
                        eps_p     = NaN,
                        eps_p_sd  = NaN,
                        alpha     = NaN,
                        alpha_sd  = NaN,
                        Central_DeltaNu = NaN,
                        Central_DeltaNu_sd = NaN,
                        Central_eps_p     = NaN,
                        Central_eps_p_sd  = NaN,
                        Central_alpha     = NaN,
                        Central_alpha_sd  = NaN,
                        gamma0    = NaN,
                        modeIDFlag = 3)
        
        write_csv(x = d.peaks, file = argv$peaks)
        write_csv(x = d.summary, file = argv$summary)
        stop("Not enough radial modes found to go any further. Mode ID not performed and Δν not estimated")
    }

    # Fit through frequencies to estimate delta nu, epsilon and aplha
    res <- DeltaNu_l0_fit(
                peaks = d.peaks %>%
                    filter(l==0) %>%
                    drop_na() %>% # in case have nan in frequency_sd
                    arrange(frequency),
                    numax = d.summary$numax,
                    DeltaNu0 = dnu_est,
                    alpha0 = d.alpha)

    # Extract Delta Nu and eps_p values
    d.Dnu <- res$DeltaNu
    d.Dnu_sd <- res$DeltaNu_sd
    d.eps_p <- res$eps_p
    d.alpha <- res$alpha
    d.alpha_sd <- res$alpha_sd

    # For consistency with epsilon from Kallinger et al. (2012)
    if( d.eps_p < 0){
        print("Epsilon p is negative! There may be a problem with the mode ID.")
    } else if( (d.eps_p < 0.5) & (d.eps_p > 0) & (d.Dnu > 3) ){
        d.eps_p <- d.eps_p + 1
        d.peaks$n <- d.peaks$n - 1
        print("Epsilon p was too small, corrected")
    } else if( (d.eps_p > 1.8)){
        d.eps_p <- d.eps_p - 1
        d.peaks$n <- d.peaks$n + 1
        print("Epsilon p was too large, corrected")
    }
    # Extract uncertainties
    d.eps_p_sd <- res$eps_p_sd

    peaks = d.peaks %>%
        filter(l==2)
    #print(peaks)
    # fit for d02, while keeping delta nu, epsilon and alpha fixed
    res2 <- DeltaNu_l2_fit(
                peaks = d.peaks %>%
                    filter(l==2) %>%
                # drop_na() %>% # in case have nan in frequency_sd
                    arrange(frequency),
                numax = d.summary$numax,
                DeltaNu0 = d.Dnu,
                alpha0 = d.alpha,
                eps_p0 = d.eps_p,
                d020 = d.d02)

    # extract d02 values
    d.d02 <- res2$d02
    d.d02_sd <- res2$d02_sd

    print(paste0("Initial fit to radial modes gives dnu: ", round(d.Dnu, 2), ", eps: ", round(d.eps_p, 2), ", alpha: ", round(d.alpha, 4), ", d02: ", round(d.d02, 4)))

    ## Make the tagging again with the new Δν
    #d.peaks <-
    #    d.peaks %>%
    #    select(-n, -l) %>%
    #    tag_l02_peaks(pds = d.pds, DeltaNu = d.Dnu, d02 = d.d02, alpha=d.alpha, numax = d.summary$numax,
    #                  HBR = NULL, sigmaEnv = d.summary$sigmaEnv,
    #                  nuNyq = d.summary$nuNyq,
    #                  search.range = search.range)
    d.peaks.l0 <- d.peaks %>% filter(l==0)

    ## Get the central small frequency separation (δν_02) if possible
    d.d02_est <- central_d02(peaks = d.peaks, numax = d.summary$numax)

    ## Get the central small frequency separation (δν_02) if possible
    d.d02.central <- central_d02(peaks = d.peaks, numax = d.summary$numax)

    # 31/01/2020 Tag l=3 as wide modes around where expected
    d.peaks$x <- (d.peaks$frequency / d.Dnu - d.eps_p) %% 1
    # l=3 occur at l=0 + deltanu/2 - 0.280 according to Mosser et al. (2010)
    l3 <- d.peaks %>%
            filter((x > 0.15) & (x < 0.26))
    l3['n'] = floor((l3$frequency / d.Dnu) - d.eps_p)
    print("Tagging any possible l=3 modes...")
    if (nrow(l3) > 0){

        
        
        tmp_l0 <- d.peaks %>%
                    filter(l == 0)

        for (i in unique(l3$n)){
            # Take widest l=3 candidate
            tmp_l3 <- l3 %>%
                        filter(n == i) %>%
                        arrange(-linewidth) %>%
                        slice(1)
        #   print(tmp_l3$frequency)
        #   print(tmp_l3$linewidth)
        #   print(tmp_l3$amplitude)
        #   print(tmp_l3$l)
        #   print(deltanu)
            # TODO: Need to make sure checking against l=0 with right radial order!
            closest_l0_amp = tmp_l0 %>% 
                                filter(n == i) %>% 
                                select(amplitude)
            closest_l0_width = tmp_l0 %>% 
                                filter(n == i) %>% 
                                select(linewidth)        

            # Make sure there is a nearest l=0 before doing this
            if(is.na(tmp_l3$l) & !is.na(tmp_l3$linewidth)){
            if((nrow(closest_l0_amp) > 0) & (nrow(closest_l0_width) > 0) & tmp_l3$linewidth > 2*deltanu){
            #   if(is.na(tmp_l3$linewidth)){
            #       if((closest_l0_amp * sqrt(0.05)*0.9 < tmp_l3$amplitude) & (closest_l0_amp * sqrt(0.05)*1.8 > tmp_l3$amplitude)){
                        d.peaks$l[d.peaks$frequency == tmp_l3$frequency] <- 3
                        # Subtract one from nearest l=0 radial order to ensure correct
                        d.peaks$n[d.peaks$frequency == tmp_l3$frequency] <- tmp_l3$n - 1
            #      }
            #  } else {
            #      if(((closest_l0_amp * sqrt(0.05)*0.9 < tmp_l3$amplitude) & (closest_l0_amp * sqrt(0.05)*1.8 > tmp_l3$amplitude)) | (closest_l0_width <= tmp_l3$linewidth)){
            #          d.peaks$l[d.peaks$frequency == tmp_l3$frequency] <- 3
                        # Subtract one from nearest l=0 radial order to ensure correct
            #          d.peaks$n[d.peaks$frequency == tmp_l3$frequency] <- tmp_l3$n - 1
            #       }
            #    }
            }
            }

        }

    }

    # Drop x column as no longer needed
    d.peaks <- d.peaks %>%
                select(-x)
    #print(paste0("ESTIMATE DELTA NU: ", dnu_est))
    # Fit through frequencies to estimate delta nu, epsilon and aplha
    res <- DeltaNu_l0_fit(        
                peaks = d.peaks %>%
                    filter(l==0) %>%
                    drop_na() %>% # in case have nan in frequency_sd
                    arrange(frequency),
                    numax = d.summary$numax,
                    DeltaNu0 = dnu_est,
                    alpha0 = d.alpha)

    # Extract Delta Nu and eps_p values
    d.Dnu <- res$DeltaNu
    d.Dnu_sd <- res$DeltaNu_sd
    d.eps_p <- res$eps_p

    # For consistency with epsilon from Kallinger et al. (2012)
    if( d.eps_p < 0){
        print("Epsilon p is negative! There may be a problem with the mode ID.")
    } else if( (d.eps_p < 0.5) & (d.eps_p > 0) & (d.Dnu > 3) ){
        d.eps_p <- d.eps_p - 1
    } else if( (d.eps_p > 1.6)) {
        print("Epsilon p too large, altering radial orders and rerunning fit")
        d.peaks$n <- d.peaks$n + 1
        # Fit through frequencies to estimate delta nu, epsilon and aplha
        res <- DeltaNu_l0_fit(        
                    peaks = d.peaks %>%
                        filter(l==0) %>%
                        drop_na() %>% # in case have nan in frequency_sd
                        arrange(frequency),
                        numax = d.summary$numax,
                        DeltaNu0 = dnu_est,
                        alpha0 = d.alpha)

        # Extract Delta Nu and eps_p values
        d.Dnu <- res$DeltaNu
        d.Dnu_sd <- res$DeltaNu_sd
        d.eps_p <- res$eps_p
    }
    # Extract uncertainties
    d.eps_p_sd <- res$eps_p_sd
    d.alpha <- res$alpha
    d.alpha_sd <- res$alpha_sd

    # fit for d02, while keeping delta nu, epsilon and alpha fixed
    res2 <- DeltaNu_l2_fit(
                peaks = d.peaks %>%
                    filter(l==2) %>%
                    #drop_na() %>% # in case have nan in frequency_sd
                    arrange(frequency),
                numax = d.summary$numax,
                DeltaNu0 = d.Dnu,
                alpha0 = d.alpha,
                eps_p0 = d.eps_p,
                d020 = d.d02)

    # extract d02 values
    d.d02 <- res2$d02
    d.d02_sd <- res2$d02_sd

    print(paste0("Final dnu: ", format(round(d.Dnu, 3), nsmall=3), "+/- ", format(round(d.Dnu_sd, 3), nsmall=3), " uHz, and eps: ", format(round(d.eps_p, 3), nsmall=3), "+/- ", format(round(d.eps_p_sd, 3), nsmall=3), ", alpha: ", format(round(d.alpha, 4),nsmall=3),"+/-",format(round(d.alpha_sd, 4), nsmall=4), ", d02: ", format(round(d.d02, 4))))

    ## Estimate central Δν from a linear fit through the 3 l=0 peaks closest to numax
    central_res <- DeltaNu_l0_fit(        
                    peaks = d.peaks %>%
                        filter(l==0) %>%
                        mutate(absdiff = abs(frequency - d.summary$numax)) %>%
                        arrange(absdiff) %>%
                        slice(1:3) %>%
                        arrange(frequency),
                        numax = d.summary$numax,
                        DeltaNu0 = d.Dnu,
                        alpha0 = d.alpha)

    # Extract Delta Nu and eps_p values
    d.central_Dnu <- central_res$DeltaNu
    d.central_Dnu_sd <- central_res$DeltaNu_sd
    d.central_eps_p <- central_res$eps_p
    # For consistency with epsilon from Kallinger et al. (2012)
    if( d.central_eps_p < 0){
        print("Central epsilon p is negative! There may be a problem with the mode ID.")
    } else if( (d.central_eps_p > 0) & (d.central_eps_p < 0.5) & (d.Dnu > 3) ){
        d.central_eps_p <- d.central_eps_p + 1
        d.peaks$n <- d.peaks$n + 1
    }
    d.central_eps_p_sd <- central_res$eps_p_sd
    d.central_alpha <- central_res$alpha
    d.central_alpha_sd <- central_res$alpha_sd

    print(paste0("Central dnu: ", format(round(d.central_Dnu, 3), nsmall=3), "+/- ", format(round(d.central_Dnu_sd, 3), nsmall=3), " uHz, and eps: ", format(round(d.central_eps_p, 3), nsmall=3), "+/- ", format(round(d.central_eps_p_sd, 3), nsmall=3), ", alpha: ", format(round(d.central_alpha, 4),nsmall=3),"+/-",format(round(d.central_alpha_sd, 4), nsmall=4)))


    # Calculate the Γ_0(ν_max) of Vrard et al. 2018 http://dx.doi.org/10.1051/0004-6361/201832477
    d.gamma0 <-
        d.peaks %>%
        filter(l == 0) %>%
        arrange(abs(frequency - d.summary$numax)) %>%
        slice(1:3) %>%
        summarise(gamma0 = sum(amplitude * linewidth) / sum(amplitude)) %>%
        as.numeric()

    # Central Δν
    d.summary <-
        d.summary %>%
        mutate(DeltaNu = d.Dnu,
            DeltaNu_sd = d.Dnu_sd,
            dNu02     = d.d02,
            eps_p     = d.eps_p,
            eps_p_sd  = d.eps_p_sd,
            alpha     = d.alpha,
            alpha_sd  = d.alpha_sd,
            Central_DeltaNu = d.central_Dnu,
            Central_DeltaNu_sd = d.central_Dnu_sd,
            Central_eps_p     = d.central_eps_p,
            Central_eps_p_sd  = d.central_eps_p_sd,
            Central_alpha     = d.central_alpha,
            Central_alpha_sd  = d.central_alpha_sd,
            gamma0    = d.gamma0,
            modeIDFlag = 0)

    return(list(peaks, data))

}
