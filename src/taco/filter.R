library(readr, quietly = TRUE, warn.conflicts = FALSE)
library(dplyr, quietly = TRUE, warn.conflicts = FALSE)

filter = function(lc){
    lc <-
    lc %>%
    filter(is.finite(flux))
}

lc <- read.table("../../data/001296068/raw.dat", header = FALSE, stringsAsFactors = FALSE) %>%
      as_tibble() %>%
      rename(time = V1, flux = V2)

print(lc)
