library(dplyr, quietly = TRUE, warn.conflicts = FALSE)

filter_r = function(lc)
{
    lc <-
        lc %>%
        filter(is.finite(flux))
    return(lc)
}
