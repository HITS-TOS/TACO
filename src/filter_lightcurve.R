## Filtering the lightcurves with a triangular smooth
## I use two rectangular smooths with half the provided width.
## Additionally it interpolates the single-point missing values
## It also calculates the mean and variance of the light curve
## and saves them in the summary file.

library(argparser, quietly=TRUE)
library(readr,     quietly = TRUE, warn.conflicts = FALSE)
library(dplyr,     quietly = TRUE, warn.conflicts = FALSE)
library(tcltk)

## Where are we located?

initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

## Get the file names and options
p <- arg_parser("Make a triangular filter of the time-series.")
p <- add_argument(p, arg = "--tseries", nargs=1, default = "raw.dat",
                  help = "File name of the time-series in tsv format with the first column giving the time and the second the flux.",
                  type = "character", flag = FALSE)
p <- add_argument(p, arg = "--output", nargs=1, default = "filtered.csv",
                  help = "File name on which to save the filtered time-series",
                  type = "character", flag = FALSE)
p <- add_argument(p, arg = "--summary", nargs=1, default = "summary.csv",
                  help = "Summary file to save some time-series statistics",
                  type = "character", flag = FALSE)
p <- add_argument(p, arg = "--width", nargs=1, default = 40,
                  help = "Width of the triangular filter",
                  type = "numeric", flag = FALSE)
p <- add_argument(p, arg = "--remove_gaps", nargs=1, default = -1,
                  help = "Remove gaps greater than this value (in days). If set to -1, do not remove gaps.",
                  type = "numeric", flag = FALSE)
argv <- parse_args(p)
#argv <- list(tseries = "raw.dat", output = "filtered.csv", summary = "summary.csv", width = 40, "remove_gaps" = -1)

## Check the provided tseries file is valid
if(!file.exists(argv$tseries))
    stop(paste("Could not find the time-series file:", argv$tseries))

## Try reading the summary file, if it exists
if (file.exists(argv$summary)) {
    d.summary <- read_csv(argv$summary, col_types = cols())
} else {
    d.summary <- tibble(KIC = NA, raw_data=NA)
}

## Read the lightcurve
d.lc <- read.table(argv$tseries, 
                   header = FALSE, stringsAsFactors = FALSE) %>%
    as_tibble() %>%
    rename(time = V1, flux = V2)

d.lc <-
    d.lc %>%
    filter(is.finite(flux))

## Fill the single-point gaps with the mean of the two adjacent points
d.lc <-
    d.lc %>%
    mutate(dt = time - lag(time))



d.deltat <- median(d.lc$dt,  na.rm=TRUE)

d.lc <-
    d.lc %>%
    mutate(dt_int = round(dt / d.deltat))

d.lc <-
    bind_rows(
        lapply(X = which(d.lc$dt_int == 2),
               FUN = function(idx) {
                   tibble(
                       time = d.lc$time[idx - 1] + d.deltat,
                       flux = (d.lc$flux[idx - 1] + d.lc$flux[idx]) / 2
                   ) %>% mutate(time_raw = time)
               })) %>%
    bind_rows(d.lc) %>%
    arrange(time) %>%
    select(-dt, -dt_int)

d.lc <-
    d.lc %>%
    mutate(time_raw = time) %>%
    mutate(dt = time - lag(time))

## Remove the large gaps
if(argv$remove_gaps != -1){
   d.gaps.idx <- which(d.lc$dt > argv$remove_gaps)

    for(idx in d.gaps.idx) {
        d.gap <- d.lc$dt[idx]
        d.time_mod <- d.lc$time[idx-1]
        d.lc <-
            d.lc %>%
            mutate(time = if_else(
                    time > d.time_mod,
                    time - d.gap + d.deltat,
                    time))               
    } 
}



## Make the filter
d.smooth <- ksmooth(x = d.lc$time, 
                    y = d.lc$flux,
                    x.points = d.lc$time,
                    kernel = "box", bandwidth = argv$width/2)
d.smooth <- ksmooth(x = d.smooth$x, 
                    y = d.smooth$y,
                    x.points = d.lc$time,
                    kernel = "box", bandwidth = argv$width/2)
d.filtered <- tibble(
    time     = d.lc$time,
    time_raw = d.lc$time_raw,
    flux     =  d.lc$flux - d.smooth$y)

## Some useful quantities
d.summary <-
    d.summary %>%
    mutate(
        mean        = mean(d.filtered$flux),
        var         = var(d.filtered$flux),
        start_date  = min(d.filtered$time),
        end_date    = max(d.filtered$time),
        fill_factor = nrow(d.filtered) *
            median(diff(d.filtered$time)) /
            diff(range(d.filtered$time)))


## Save the results
write_csv(d.filtered, file = argv$output)
write_csv(d.summary,  file = argv$summary)
