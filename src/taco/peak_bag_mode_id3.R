## Label the l=3 peaks

library(lomb, quietly = TRUE, warn.conflicts = FALSE)
library(argparser, quietly = TRUE, warn.conflicts = FALSE)
library(readr, quietly = TRUE, warn.conflicts = FALSE)
library(magrittr, quietly = TRUE, warn.conflicts = FALSE)
options(tibble.width = Inf)

#https://stackoverflow.com/questions/38351820/negation-of-in-in-r
`%notin%` <- Negate(`%in%`)

source("src/peakFind_lib.R", chdir = TRUE)
source("src/l02_modes_id.R", chdir = TRUE)

peak_bag_mode_id3_r <- function(pds, peaks, data) {

    ## Arbitrary parameters that I use
    MAX.ERROR <- 0.05 # Discard taggings where the square difference
                      # between expected and predicted frequencies
                      # is greater than MAX.ERROR*Δν
    search.range <- c(0.75, 1.25) # Look for the next l=0,2 modes
                                  # using this range away from Δν

    peaks <- peaks %>%
        filter(frequency > data$numax - 3 * data$sigmaEnv,
               frequency < data$numax + 3 * data$sigmaEnv) %>%
        #select(-one_of("n", "l", "m")) %>%
        filter(AIC > 0)

    pds <- pds %>%
        filter(frequency > data$numax - 3 * data$sigmaEnv,
               frequency < data$numax + 3 * data$sigmaEnv)

    deltanu <- pds$frequency[2] - pds$frequency[1]

    # 31/01/2020 Tag l=3 as wide modes around where expected
    peaks$x <- (peaks$frequency / data$DeltaNu - data$eps_p) %% 1 #include alpha
    #print(peaks)
    # l=3 occur at l=0 + deltanu/2 - 0.280 according to Mosser et al. (2010)
    l3 <- peaks %>% filter((x > 0.13))
    l3 <- l3 %>% filter((x < 0.26))   # was selecting from peaks
    l3['n'] <- floor((l3$frequency / data$DeltaNu) - data$eps_p)
    #print(l3)
    if (nrow(l3) > 0) {
        print("Tagging any possible l=3 modes...")
        #print(l3)
        tmp_l0 <- peaks %>% filter(l == 0)

        for (i in unique(l3$n)) {
            # TODO: Need to make sure checking against l=0 with right radial order!
            closest_l0_amp = tmp_l0 %>%
                                filter(n == i) %>%
                                select(amplitude)
            closest_l0_width = tmp_l0 %>%
                                filter(n == i) %>%
                                select(linewidth)

            # Take widest l=3 candidate, can be multiple
            if (nrow(closest_l0_width)> 0){
                tmp_l3 <- l3 %>%
                        filter(n == i) %>%
                        #arrange(-linewidth) %>%
                        #slice(1)
                        filter(linewidth > 0.3 * closest_l0_width)

            # Make sure there is a nearest l=0 before doing this
                if (is.na(tmp_l3$l) && !is.na(tmp_l3$linewidth)) {
                    if ((nrow(closest_l0_amp) > 0) && (nrow(closest_l0_width) > 0) && tmp_l3$linewidth > deltanu) {
                        peaks$l[peaks$frequency == tmp_l3$frequency] <- 3
                        # Subtract one from nearest l=0 radial order to ensure correct
                        peaks$n[peaks$frequency == tmp_l3$frequency] <- tmp_l3$n - 1
                    }
                }
            }
        }
    }

    # Drop x column as no longer needed
    peaks <- peaks %>% select(-x)


    return(list(peaks=peaks, data=data))

}
