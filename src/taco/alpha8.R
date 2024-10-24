## Get Δν and label the l=0,2 peaks

library(lomb, quietly = TRUE, warn.conflicts=FALSE)
library(argparser, quietly = TRUE, warn.conflicts=FALSE)
library(readr, quietly = TRUE, warn.conflicts=FALSE)
library(magrittr, quietly = TRUE, warn.conflicts=FALSE)
#library(tcltk)
#library(ggplot2)
options(tibble.width = Inf)

#https://stackoverflow.com/questions/38351820/negation-of-in-in-r
`%notin%` <- Negate(`%in%`)

#source("~/TACO/src/peakFind_lib.R", chdir = TRUE)
source("~/TACO/src/l02_modes_id.R", chdir = TRUE)


## Where are we located?
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

## Get the file names and options
p <- arg_parser("Degree identification for radial and quadrupole (l=0,2) modes in the peaks found. Also estimates Δν and δν_02")
p <- add_argument(p, arg = "--peaks", nargs=1, default = "peaksMLE.csv",
                  help = "File name of the identified peaks",
                  type = "character", flag = FALSE)
p <- add_argument(p, arg = "--summary", nargs=1, default = "summary.csv",
                  help = "File name of the csv file with summary values. It must contain numax and sigmaEnv",
                  type = "character", flag = FALSE)
p <- add_argument(p, arg = "--output", nargs=1, default = "alpha6.csv",
                  help = "File name of the identified peaks",
                  type = "character", flag = FALSE)
argv <- parse_args(p)

## Check the provided PDS file is valid
if(!file.exists(argv$summary))
    stop(paste("Could not find the summary file:", argv$summary))
if(!file.exists(argv$peaks))
    stop(paste("Could not find the file with the peaks:", argv$peaks))


## Arbitrary parameters that I use
MAX.ERROR <- 0.03    # Discard taggings where the square difference between expected and predicted frequencies is greater than MAX.ERROR*Δν
search.range <- c(0.75, 1.25)  # Look for the next l=0,2 modes using this range away from Δν

## Read the data
d.summary <- read_csv(argv$summary, col_types = cols())

d.peaks <-
    read_csv(argv$peaks, col_types = cols()) %>%
    filter(frequency > d.summary$numax - 3*d.summary$sigmaEnv,
           frequency < d.summary$numax + 3*d.summary$sigmaEnv) %>%
    #select(-one_of("n", "l")) %>%
    filter(AIC > 0)

d.peaks.l0 <- d.peaks %>% filter(l==0)
# Fit through frequencies to estimate delta nu, epsilon and aplha
alpha <- alpha_obs_from_DeltaNu(d.summary$DeltaNu)
res <- DeltaNu_l0_fit(
            peaks = d.peaks %>%
                filter(l==0) %>%
                mutate(absdiff = abs(frequency - d.summary$numax)) %>%
                arrange(absdiff) %>%
                slice(1:8) %>%
                arrange(frequency),
                numax = d.summary$numax,
                DeltaNu0 = d.summary$DeltaNu,
                alpha0 = d.summary$alpha)

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

print(paste0("Initial fit to radial modes gives dnu: ", round(d.Dnu, 2), ", eps: ", round(d.eps_p, 2), ", alpha: ", round(d.alpha, 4)))

d.summary <-
    d.summary %>%
    mutate(DeltaNu = d.Dnu,
           DeltaNu_sd = d.Dnu_sd,
           eps_p     = d.eps_p,
           eps_p_sd  = d.eps_p_sd,
           alpha     = d.alpha,
           alpha_sd  = d.alpha_sd)

write_csv(x = d.summary,   file = argv$output)
