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
    MAX.ERROR <- 0.03 # Discard taggings where the square difference
                      # between expected and predicted frequencies
                      # is greater than MAX.ERROR*Δν
    search.range <- c(0.75, 1.25) # Look for the next l=0,2 modes
                                  # using this range away from Δν

    peaks <- peaks %>%
        filter(frequency > data$numax - 3 * data$sigmaEnv,
               frequency < data$numax + 3 * data$sigmaEnv) %>%
        select(-one_of("n", "l", "m")) %>%
        filter(AIC > 0)

    pds <- pds %>%
        filter(frequency > data$numax - 3 * data$sigmaEnv,
               frequency < data$numax + 3 * data$sigmaEnv)

    deltanu <- pds$frequency[2] - pds$frequency[1]
    
    flag <- 0

    ## Check to see if there are any peaks in file, if not then stop
    if (nrow(peaks) == 0) {
        peaks <- tibble(frequency = double(),
                        amplitude = double(),
                        linewidth = double(),
                        frequency_sd = double(),
                        amplitude_sd = double(),
                        linewidth_sd = double(),
                        height = double(),
                        AIC = double(),
                        n = integer(),
                        l = integer())
        # Central Δν
        data <- data %>%
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
        flag <- 1
        print("No peaks detected so not proceeding with mode ID")
        return(list(peaks, flag, data))
    } else if (nrow(peaks) < 3) {
        peaks$l <- NA
        data <- data %>%
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
        flag <- 1
        print("Not enough peaks detected so not proceeding with mode ID")
        return(list(peaks, flag, data))
    }

    if (data$numax < 5) {
        data <- data %>%
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
        flag <- 1
        print("Numax < 10uHz and is too low for the automated mode identification to be reliable.")
        return(list(peaks, flag, data))
    }

    ## Expected Δν from ν_max
    Dnu <- DeltaNu_from_numax(data$numax)

    ## Refine it using the PS of the PS. We look at Δν and Δν/2 and take the highest peak
    psxps1 <- lsp(times = pds$frequency, x = pds$power,
                    from = 0.3 * Dnu, to = 0.7 * Dnu,
                    type = "period", plot = FALSE, ofac = 10)
    psxps2 <- lsp(times = pds$frequency, x = pds$power,
                    from = 0.7 * Dnu, to = 1.3 * Dnu,
                    type = "period", plot = FALSE, ofac = 10)
    Dnu <- ifelse(psxps1$peak > psxps2$peak, psxps1$peak.at[1] * 2, psxps2$peak.at[1])

    print(paste0("Initial Dnu from PSxPS ", Dnu))

    # Expected δν_02 from Δν
    d02 <- d02_from_DeltaNu(Dnu)

    ## Refine the δν_02 with the PSxPS
    psxps <- lsp(times = pds$frequency, x = pds$power,
                 from = 0.7 * d02, to = 1.3 * d02,
                type = "period", plot = FALSE, ofac = 10)
    d02 <- psxps$peak.at[1]

    # Expected alpha from Δν
    n_max <- data$numax / Dnu - eps_p_from_Dnu(Dnu)
    alpha <- 1.5 * alpha_obs_from_n_max(n_max) # overestimating alpha a bit to avoid issues at extreme frequencies
    print(paste0("Expected alpha from nmax ", alpha))

    ## Mode identification for l=0,2 modes
    peaks <- peaks %>%
        tag_l02_peaks(pds = pds, DeltaNu = Dnu, d02 = d02,
            alpha = alpha, numax = data$numax,
            HBR = NULL, sigmaEnv = data$sigmaEnv,
            nuNyq = data$nuNyq,
            search.range = search.range)

    peaks.l0 <- peaks %>% filter(l == 0)
    dnu_est <- Dnu

    # Pull out l=0 and l=2 modes
    peaks.l0 <- peaks %>% filter(l == 0)
    peaks.l2 <- peaks %>% filter(l == 2)

    # If there are fewer than two modes exit mode ID
    if (nrow(peaks.l0) < 3) {

        peaks <- peaks %>%
                    mutate(n = NaN,
                        l = NaN)

        data <- data %>%
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
        flag <- 1
        print("Not enough radial modes found to go any further. Mode ID not performed and Δν not estimated")
        return(list(peaks, flag, data))
    }

    # Fit through frequencies to estimate delta nu, epsilon and alpha
    res <- DeltaNu_l0_fit(
                peaks = peaks %>%
            filter(l == 0) %>%
                    drop_na() %>% # in case have nan in frequency_sd
                    arrange(frequency),
                    numax = data$numax,
                    DeltaNu0 = dnu_est,
                    alpha0 = alpha)

    # Extract Delta Nu and eps_p values
    Dnu <- res$DeltaNu
    Dnu_sd <- res$DeltaNu_sd
    eps_p <- res$eps_p
    alpha <- res$alpha
    alpha_sd <- res$alpha_sd

    # For consistency with epsilon from Kallinger et al. (2012)
    if (eps_p < 0) {
        print("Epsilon p is negative! There may be a problem with the mode ID.")
    } else if ((eps_p < 0.5) & (eps_p > 0) & (Dnu > 3)) {
        eps_p <- eps_p + 1
        peaks$n <- peaks$n - 1
        print("Epsilon p was too small, corrected")
    } else if ((eps_p > 1.8)) {
        eps_p <- eps_p - 1
        peaks$n <- peaks$n + 1
        print("Epsilon p was too large, corrected")
    }
    # Extract uncertainties
    eps_p_sd <- res$eps_p_sd

    # fit for d02, while keeping delta nu, epsilon and alpha fixed
    res2 <- DeltaNu_l2_fit(
                peaks = peaks %>%
                filter(l == 2) %>%
                    arrange(frequency),
                numax = data$numax,
                DeltaNu0 = Dnu,
                alpha0 = alpha,
                eps_p0 = eps_p,
                d020 = d02)

    # extract d02 values
    d02 <- res2$d02
    d02_sd <- res2$d02_sd

    print(paste0("Initial fit to radial modes gives dnu: ", round(Dnu, 2),
          ", eps: ", round(eps_p, 2),
          ", alpha: ", round(alpha, 4),
          ", d02: ", round(d02, 4)))

    ## Make the tagging again with the new Δν
    peaks.l0 <- peaks %>% filter(l == 0)

    ## Get the central small frequency separation (δν_02) if possible
    d02_est <- central_d02(peaks = peaks, numax = data$numax)

    ## Get the central small frequency separation (δν_02) if possible
    d02.central <- central_d02(peaks = peaks, numax = data$numax)

    # 31/01/2020 Tag l=3 as wide modes around where expected
    peaks$x <- (peaks$frequency / Dnu - eps_p) %% 1
    # l=3 occur at l=0 + deltanu/2 - 0.280 according to Mosser et al. (2010)
    l3 <- peaks %>% filter((x > 0.15) && (x < 0.26))
    l3['n'] <- floor((l3$frequency / Dnu) - eps_p)

    print("Tagging any possible l=3 modes...")
    if (nrow(l3) > 0) {

        tmp_l0 <- peaks %>% filter(l == 0)

        for (i in unique(l3$n)) {
            # Take widest l=3 candidate
            tmp_l3 <- l3 %>%
                        filter(n == i) %>%
                        arrange(-linewidth) %>%
                        slice(1)

            # TODO: Need to make sure checking against l=0 with right radial order!
            closest_l0_amp = tmp_l0 %>%
                                filter(n == i) %>%
                                select(amplitude)
            closest_l0_width = tmp_l0 %>%
                                filter(n == i) %>%
                                select(linewidth)

            # Make sure there is a nearest l=0 before doing this
            if (is.na(tmp_l3$l) && !is.na(tmp_l3$linewidth)) {
                if ((nrow(closest_l0_amp) > 0) && (nrow(closest_l0_width) > 0) && tmp_l3$linewidth > 2 * deltanu) {
                        peaks$l[peaks$frequency == tmp_l3$frequency] <- 3
                        # Subtract one from nearest l=0 radial order to ensure correct
                        peaks$n[peaks$frequency == tmp_l3$frequency] <- tmp_l3$n - 1
            }
            }
        }
    }

    # Drop x column as no longer needed
    peaks <- peaks %>% select(-x)

    # Fit through frequencies to estimate delta nu, epsilon and alpha
    res <- DeltaNu_l0_fit(
                peaks = peaks %>%
                    filter(l == 0) %>%
                    drop_na() %>% # in case have nan in frequency_sd
                    arrange(frequency),
                    numax = data$numax,
                    DeltaNu0 = dnu_est,
                    alpha0 = alpha)

    # Extract Delta Nu and eps_p values
    Dnu <- res$DeltaNu
    Dnu_sd <- res$DeltaNu_sd
    eps_p <- res$eps_p

    # For consistency with epsilon from Kallinger et al. (2012)
    if (eps_p < 0) {
        print("Epsilon p is negative! There may be a problem with the mode ID.")
    } else if ((eps_p < 0.5) && (eps_p > 0) && (Dnu > 3)) {
        eps_p <- eps_p - 1
    } else if ((eps_p > 1.6)) {
        print("Epsilon p too large, altering radial orders and rerunning fit")
        peaks$n <- peaks$n + 1
        # Fit through frequencies to estimate delta nu, epsilon and alpha
        res <- DeltaNu_l0_fit(
                    peaks = peaks %>%
                        filter(l == 0) %>%
                        drop_na() %>% # in case have nan in frequency_sd
                        arrange(frequency),
                        numax = data$numax,
                        DeltaNu0 = dnu_est,
                        alpha0 = alpha)

        # Extract Delta Nu and eps_p values
        Dnu <- res$DeltaNu
        Dnu_sd <- res$DeltaNu_sd
        eps_p <- res$eps_p
    }
    # Extract uncertainties
    eps_p_sd <- res$eps_p_sd
    alpha <- res$alpha
    alpha_sd <- res$alpha_sd

    # fit for d02, while keeping delta nu, epsilon and alpha fixed
    res2 <- DeltaNu_l2_fit(
                peaks = peaks %>%
                    filter(l == 2) %>%
                    arrange(frequency),
                numax = data$numax,
                DeltaNu0 = Dnu,
                alpha0 = alpha,
                eps_p0 = eps_p,
                d020 = d02)

    # extract d02 values
    d02 <- res2$d02
    d02_sd <- res2$d02_sd

    print(paste0("Final dnu: ", format(round(Dnu, 3), nsmall = 3),
                 "+/- ", format(round(Dnu_sd, 3), nsmall = 3),
                 " uHz, and eps: ", format(round(eps_p, 3), nsmall = 3),
                 "+/- ", format(round(eps_p_sd, 3), nsmall = 3),
                 ", alpha: ", format(round(alpha, 4),nsmall = 3),
                 "+/-",format(round(alpha_sd, 4), nsmall = 4),
                 ", d02: ", format(round(d02, 4))))

    ## Estimate central Δν from a linear fit through the 3 l=0 peaks closest to numax
    central_res <- DeltaNu_l0_fit(
                    peaks = peaks %>%
                        filter(l == 0) %>%
                        mutate(absdiff = abs(frequency - data$numax)) %>%
                        arrange(absdiff) %>%
                        slice(1:3) %>%
                        arrange(frequency),
                        numax = data$numax,
                        DeltaNu0 = Dnu,
                        alpha0 = alpha)

    # Extract Delta Nu and eps_p values
    central_Dnu <- central_res$DeltaNu
    central_Dnu_sd <- central_res$DeltaNu_sd
    central_eps_p <- central_res$eps_p
    # For consistency with epsilon from Kallinger et al. (2012)
    if (central_eps_p < 0) {
        print("Central epsilon p is negative! There may be a problem with the mode ID.")
    } else if ((central_eps_p > 0) && (central_eps_p < 0.5) && (Dnu > 3)) {
        central_eps_p <- central_eps_p + 1
        peaks$n <- peaks$n + 1
    }
    central_eps_p_sd <- central_res$eps_p_sd
    central_alpha <- central_res$alpha
    central_alpha_sd <- central_res$alpha_sd

    print(paste0("Central dnu: ", format(round(central_Dnu, 3), nsmall = 3),
                 "+/- ", format(round(central_Dnu_sd, 3), nsmall = 3),
                 " uHz, and eps: ", format(round(central_eps_p, 3), nsmall = 3),
                 "+/- ", format(round(central_eps_p_sd, 3), nsmall = 3),
                 ", alpha: ", format(round(central_alpha, 4),nsmall = 3),
                 "+/-", format(round(central_alpha_sd, 4), nsmall = 4)))

    # Calculate the Γ_0(ν_max) of Vrard et al. 2018 http://dx.doi.org/10.1051/0004-6361/201832477
    gamma0 <-
        peaks %>%
        filter(l == 0) %>%
        arrange(abs(frequency - data$numax)) %>%
        slice(1:3) %>%
        summarise(gamma0 = sum(amplitude * linewidth) / sum(amplitude)) %>%
        as.numeric()

    # Central Δν
    data <-
        data %>%
        mutate(
            DeltaNu    = Dnu,
            DeltaNu_sd = Dnu_sd,
            dNu02     = d02,
            eps_p     = eps_p,
            eps_p_sd  = eps_p_sd,
            alpha     = alpha,
            alpha_sd  = alpha_sd,
            Central_DeltaNu = central_Dnu,
            Central_DeltaNu_sd = central_Dnu_sd,
            Central_eps_p     = central_eps_p,
            Central_eps_p_sd  = central_eps_p_sd,
            Central_alpha     = central_alpha,
            Central_alpha_sd  = central_alpha_sd,
            gamma0    = gamma0,
            modeIDFlag = 0)

    return(list(peaks, flag, data))

}
