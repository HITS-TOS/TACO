library(dplyr, quietly = TRUE, warn.conflicts = FALSE)

filter_r <- function(lc, data, width, remove_gaps) {

    names(lc)[1] <- "time"
    names(lc)[2] <- "flux"

    lc <- lc %>% filter(is.finite(flux))

    # Fill the single-point gaps with the mean of the two adjacent points
    lc <- lc %>% mutate(dt = time - lag(time))

    deltat <- median(lc$dt, na.rm = TRUE)

    lc <- lc %>% mutate(dt_int = round(dt / deltat))

    lc <-
        bind_rows(
            lapply(X = which(lc$dt_int == 2),
                FUN = function(idx) {
                    tibble(
                        time = lc$time[idx - 1] + deltat,
                        flux = (lc$flux[idx - 1] + lc$flux[idx]) / 2
                    ) %>% mutate(time_raw = time)
                })) %>%
        bind_rows(lc) %>%
        arrange(time) %>%
        select(-dt, -dt_int)

    lc <-
        lc %>%
        mutate(time_raw = time) %>%
        mutate(dt = time - lag(time))

    # Remove the large gaps
    if (remove_gaps != -1) {
        gaps.idx <- which(lc$dt > remove_gaps)

        for (idx in gaps.idx) {
            gap <- d.lc$dt[idx]
            time_mod <- d.lc$time[idx - 1]
            lc <-
                lc %>%
                mutate(time = if_else(
                        time > time_mod,
                        time - gap + deltat,
                        time))
        }
    }

    # Make the filter
    smooth <- ksmooth(x = lc$time,
                      y = lc$flux,
                      x.points = lc$time,
                      kernel = "box", bandwidth = width / 2)
    smooth <- ksmooth(x = smooth$x,
                      y = smooth$y,
                      x.points = lc$time,
                      kernel = "box", bandwidth = width / 2)
    filtered <- tibble(
        time     = lc$time,
        time_raw = lc$time_raw,
        flux     = lc$flux - smooth$y)

    # Some useful quantities
    data$mean        <- mean(filtered$flux)
    data$var         <- var(filtered$flux)
    data$start_date  <- min(filtered$time)
    data$end_date    <- max(filtered$time)
    data$fill_factor <- nrow(filtered) *
                        median(diff(filtered$time)) /
                        diff(range(filtered$time))

    return(list(filtered=filtered, data=data))

}
