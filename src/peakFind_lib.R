## Automatic peak-finding with a CWT-based pattern matching algorithm
## using a tree-map of the local maxima in CWT-space.
## The main output is a data_frame with the peaks found. These
## can be lorentzians or sinc functions. If it's a sinc function
## the value for 'linewidth' will be NA

#library(wmtsa, quietly=TRUE)
library(tibble, quietly=TRUE, warn.conflicts=FALSE)
library(tidyr, quietly=TRUE, warn.conflicts=FALSE)
library(readr, quietly=TRUE, warn.conflicts=FALSE)
library(purrr, quietly=TRUE, warn.conflicts=FALSE)
library(dplyr, quietly=TRUE, warn.conflicts=FALSE)
library(modelr, quietly=TRUE)
#library(numDeriv)
#library(tcltk)

initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

if(is_empty(script.basename)) {
    source("wavelets.R")
} else {
    source(file.path(script.basename, "wavelets.R"))
}

if(!is.null(sys.frame(1)$ofile)) {
  EFP_model <- readRDS(file.path(dirname(sys.frame(1)$ofile), "false-alarm-model.rds"))
  False_negatives_model <- readRDS(file.path(dirname(sys.frame(1)$ofile), "false-positives-model.rds"))
}

GAMMA_PROP    <- 4  ## how many HWHM to consider around each peak while making comparisons

sinc <- function (x) {
  y <- x
  z <- x == 0
  y[z] <- 1
  y[!z] <- sin(x[!z])/(x[!z])
  return(y)
}

Lor_model <- function(pds, peak) {
  ## PDS with a lorentzian 'peak' added to it
  return(peak$height/(1 + ((pds$frequency - peak$frequency)/peak$linewidth)^2))
}

sinc2_model <- function(pds, peak) {
  ## PDS with an unresolved sinc^2 'peak' added to it
  deltanu <- diff(pds$frequency[1:2])
  return(peak$height * (sinc(pi * (pds$frequency - peak$frequency) / deltanu))^2)
}

Lor_HWHM <- function(scale) {
  ## Transform CWT scale to a lorentzian's HWHM
  return(exp(-0.9022856) * (scale^0.9568811))
}

CWT_scale_from_HWHM <- function(gk) {
    ## Inverse of Lor_HWHM
    return((exp(0.9022856) * gk)^(1/0.9568811))
}

Lor_I <- function(extrema, scale, deltanu) {
  ## Transform CWT extrema and scale to a lorentzian's height
  return(exp(0.7558125) * ((extrema * sqrt(deltanu/scale))^0.9528125))
}

peak_area <- function(peak, deltanu) {
  if (nrow(peak) == 0) return(NULL)
  if (is.na(peak$linewidth)) {
    return(pi * peak$height * deltanu)
  } else {
    return(pi * peak$height * peak$linewidth)
  }
}

peak_amplitude <- function(peak, deltanu) {
  if (nrow(peak) == 0) return(NULL)
  return(sqrt(peak_area(peak, deltanu)))
}

peak_height <- function(peak, deltanu) {
  if(is.na(peak$linewidth)) {
    return(peak$amplitude^2 / (pi*deltanu))
  } else {
    return(peak$amplitude^2 / (pi*peak$linewidth))
  }
}

peaks_with_amplitude <- function(peaks, deltanu) {
  peaks <- peaks %>% select_if(names(peaks) != "amplitude")
  res <- do.call(rbind,
                 Map(function (i) {
                   return(cbind(peaks[i,],
                                data.frame(amplitude = peak_amplitude(peaks[i,], deltanu))))
                 },1:nrow(peaks)))
  return(as_tibble(res))
}

peaks_with_height <- function(peaks, deltanu) {
  peaks <- peaks %>% select_if(names(peaks) != "height")
  res <- do.call(rbind,
                 Map(function (i) {
                   return(cbind(peaks[i,],
                                data.frame(height = peak_height(peaks[i,], deltanu))))
                 },1:nrow(peaks)))
  return(as_tibble(res))
}

fit_model <- function(pds, peaks = NULL) {
  ## Take 'peaks' and return a pds model supposing they are all from
  ## lorentzian profiles or sinc^2
  N <- nrow(pds)
  m <- nrow(peaks)
  res <- vector(mode = "numeric", length = N)
  res <- res + 1
  if (is.null(peaks)) return(res)
  if (m == 0) return(res)
  for (i in 1:m) {
    if (is.na(peaks$linewidth[i]) | is.na(peaks$linewidth[i])) {
      res <- res + sinc2_model(pds, peaks[i,])
    } else {
      res <- res + Lor_model(pds, peaks[i,])
    }
  }
  return(res)
}

add_peaks_to_model <- function(pds, model, peaks) {
  N <- nrow(pds)
  m <- nrow(peaks)
  if (m == 0) return(model)
  for (i in 1:m) {
    if (is.na(peaks$linewidth[i]) | is.na(peaks$linewidth[i])) {
      model <- model + sinc2_model(pds, peaks[i,])
    } else {
      model <- model + Lor_model(pds, peaks[i,])
    }
  }
  return(model)
}

log_likelihood <- function(pds, model, naverages=1) {
  return(-naverages*sum(log(model) + pds$power/model, na.rm=TRUE))
}

## First derivatives of the Lorentzian peaks
partial_M_freq <- function(nu, peak, deltanu = NULL) {
  nuj <- peak$frequency
  Aj  <- peak$amplitude
  gj  <- peak$linewidth
  .e1 <- nu - nuj
  .e3 <- (.e1/gj)^2 + 1
  .e5 <- gj * pi * .e3
  .e6 <- .e5^2
  .e8 <- Aj^2 * pi
  .e9 <- 2 * (.e8 * .e1/(gj * .e6))
  return(.e9)
}

partial_M_ampl <- function(nu, peak, deltanu = NULL) {
  nuj <- peak$frequency
  Aj  <- peak$amplitude
  gj  <- peak$linewidth
  .e1 <- nu - nuj
  .e3 <- (.e1/gj)^2 + 1
  .e5 <- gj * pi * .e3
  return(2 * (Aj/.e5))
}

partial_M_linw <- function(nu, peak, deltanu = NULL) {
  nuj <- peak$frequency
  Aj  <- peak$amplitude
  gj  <- peak$linewidth
  .e1 <- nu - nuj
  .e3 <- (.e1/gj)^2 + 1
  .e5 <- gj * pi * .e3
  .e6 <- .e5^2
  .e8 <- Aj^2 * pi
  return(-(.e8 * (.e3 - 2 * (.e1^2/gj^2))/.e6))
}

## First derivatives of the sinc^2 peaks
partial_M_amplsc <- function(nu, peak, deltanu) {
  nuj <- peak$frequency
  Aj  <- peak$amplitude
  return(2 * Aj * (sinc((nu - nuj)/deltanu))^2 / (deltanu *  pi))
}

partial_M_freqsc <- function(nu, peak, deltanu) {
  nuj <- peak$frequency
  Aj  <- peak$amplitude
  .e1 <- nu - nuj
  .e2 <- .e1/deltanu
  .e3 <- sin(.e2)
  .e5 <- sinc(.e1/deltanu)
  .e7 <- Aj^2
  .e8 <- cos(.e2)
  .e9 <- pi * .e1^2
  res <- 2 * (.e7 * (.e5 - .e8) * .e3/.e9)
  res[near(.e9, 0)] <- 0
  return(res)
}

partial_LL_theta <- function(peak, pds, M, deltanu, partial_M_theta_fn) {
  ## First derivative of the negative log-likelihood
  return(
    (1 - pds$power / M) *
      partial_M_theta_fn(nu = pds$frequency, peak = peak, deltanu = deltanu) /
      M
  )
}

LLgr <- function(peaks, pds, other_peaks = NULL) {
  ## Gradient of the negative log-likelihood function
  deltanu <- diff(pds$frequency[1:2])
  peaks.res <- peaks %>% filter(!is.na(peaks$linewidth))
  peaks.unr <- peaks %>% filter( is.na(peaks$linewidth))
  n.res <- nrow(peaks.res)
  n.unr <- nrow(peaks.unr)
  N <- nrow(pds)
  M <- fit_model(pds = pds, peaks = bind_rows(peaks, other_peaks))
  res_freq <- c()
  res_ampl <- c()
  res_linw <- c()
  unr_freq <- c()
  unr_ampl <- c()
  if (n.res > 0) {
    for(i in 1:n.res) {
      peak.res <- peaks.res[i,]
      res_freq <- c(res_freq, sum(partial_LL_theta(peak.res, pds, M, deltanu, partial_M_freq)))
      res_ampl <- c(res_ampl, sum(partial_LL_theta(peak.res, pds, M, deltanu, partial_M_ampl)))
      res_linw <- c(res_linw, sum(partial_LL_theta(peak.res, pds, M, deltanu, partial_M_linw)))
    }
  }
  if (n.unr > 0) {
    for(i in 1:n.unr) {
      peak.unr <- peaks.unr[i,]
      unr_freq <- c(unr_freq, sum(partial_LL_theta(peak.unr, pds, M, deltanu, partial_M_freqsc)))
      unr_ampl <- c(unr_ampl, sum(partial_LL_theta(peak.unr, pds, M, deltanu, partial_M_amplsc)))
    }
  }
  return(c(res_freq, res_ampl, res_linw, unr_freq, unr_ampl))
}

which_local_max <- function(x) {
  return(which(diff(sign(diff(x)))==-2)+1)
}

branch_peaks <- function(CWTTree, branch_number) {
  ## Take branch_number from CWTTree and get all the local maxima.
  ## Return a data.frame with the results
  if (!is(CWTTree, "wavCWTTree"))
    stop("Input must be an object of class wavCWTTree")

  #X11()
  ##print(CWTTree[[branch_number]])
  #plot(CWTTree[[branch_number]]$scale, CWTTree[[branch_number]]$extrema)


  peaks_idx <- which_local_max(CWTTree[[branch_number]]$extrema)
  #print("REMEMBER TO CHANGE THIS BIT BACK!")
  #peaks_idx <- which(CWTTree[[branch_number]]$extrema == max(CWTTree[[branch_number]]$extrema))
  if (length(peaks_idx) == 0) {
    peaks <- NULL
  } else {
    peaks <- data.frame(CWTTree[[branch_number]])[peaks_idx,]
    names(peaks)[names(peaks) == "time"] <- "frequency"
    for (i in 1:length(peaks_idx)) {
      peaks[i, "snr"] <- peaks[i, "extrema"] / 2
      peaks[i, "branch"] <- branch_number
    }
    peaks <- peaks[,c("frequency", "scale", "extrema", "snr", "branch")]
  }
  #print(peaks_idx)
  #print(peaks)
  #prompt  <- "hit spacebar to close plots"
  #extra   <- "some extra comment"
  #capture <- tk_messageBox(message = prompt, detail = extra)

  return(peaks)
}

cwtPeaks <- function(CWTTree, min.snr = 3,
                     scale.range = attr(CWTTree, "scale")[c(3, length(attr(CWTTree, "scale")))],
                     var.maxlw = FALSE) {
  ## A 'peak' is a local maxima (in scale) from a branch
  ## found in CWTTree. Here I just call brach_peaks for all branches.
  if (!is(CWTTree, "wavCWTTree"))
    stop("Input must be an object of class wavCWTTree")
  peaks <- do.call(rbind, Map(function (i) {
    return(branch_peaks(CWTTree, i))
  }, 1:length(CWTTree)))

  peaks <- peaks[peaks$snr >= min.snr,]

  if (var.maxlw==TRUE) {
    # Updated on 06.05.2024.
    # According to Vrard+2016, radial modes linewidth is expected to increase with frequency
    # This is implemented to change maxlw used to find peaks.

    # 1. Find peaks around numax, within 1.1*Deltanu
    numax <- 10^(((Lor_HWHM(scale.range[2])/1.5) + 0.08792122)/0.16423368)
    estimated_dnu <- 1.1 * DeltaNu_from_numax(numax)
    central_peaks <- peaks[peaks$frequency >= numax - (estimated_dnu/2) & peaks$frequency <= numax + (estimated_dnu/2),]

    # 2. Peak with maximum scale in this range is going to be used as reference ("central mode")
    central <- central_peaks[which.max(central_peaks$scale), ]

    # 3. Assign estimated radial orders to all the peaks with respect to central one.
    if (length(central$frequency) > 0) {
      peaks$n_frac <- (peaks$frequency - central$frequency) / estimated_dnu

      # 4. Correct maxlw estimation
      good_peaks <- peaks[NULL,]

      for (i in 1:length(peaks$scale)){
        peak <- peaks[i,]
        if (peak$n_frac>0){
          slope = 0.249
        } else {
          slope = 0.0
        }

        if (peak$scale >= scale.range[1] & peak$scale <= (scale.range[2]*(1 + slope*peak$n_frac))){
          peak$maxlw <- (scale.range[2]*(1 + slope*peak$n_frac))
          good_peaks <- rbind(good_peaks, peak)
        }
      }
    } else {
      good_peaks <- peaks[peaks$scale >= scale.range[1] & peaks$scale <= scale.range[2],] #[fixed maxlw]
    }

  } else {
    good_peaks <- peaks[peaks$scale >= scale.range[1] & peaks$scale <= scale.range[2],] #[fixed maxlw]
  }


  deltanu <- diff(attr(CWTTree, "time")[1:2])
  if (!is.null(good_peaks)) {
    good_peaks$height <- Lor_I(good_peaks$extrema, good_peaks$scale, deltanu)
    good_peaks$linewidth <- Lor_HWHM(good_peaks$scale)
  }

  return(good_peaks)
}

which_all_peaks_under <- function(peaks, p) {
  ## Look at all peaks laying within plus or minus GAMMA_PROP*FWHM of a
  ## lorentzian at p. Return one peak per branch such that it has
  ## the greatest scale for all the peaks at each branch.
  peaks$idx <- 1:nrow(peaks)
  peaks <- peaks[order(peaks$scale, decreasing = TRUE),]
  s <- p$scale
  w <- GAMMA_PROP * Lor_HWHM(s)
  f <- p$frequency
  peaks <- peaks[peaks$frequency  > f - w & peaks$frequency  < f + w & peaks$scale < s,]
  #peaks_branch <- peaks[duplicated(peaks$branch),]
  #if (nrow(peaks_branch) > 1){
  #    peaks <- rbind(peaks,peaks_branch)
  #}
  #peaks <- peaks[!duplicated(peaks$branch),]
  return(peaks[, "idx"])
}

which_peaks_under <- function(peaks, p) {
  ## Look at all the branches below GAMMA_PROP of a lorentzian's width
  ## at p. Return the index of one peak per branch such that those
  ## peaks don't have any other point *above* them (in the sense
  ## of them being within GAMMA_PROP of a lorentzian's width
  ## of another peak also *below* p).
  all_p_idx <- which_all_peaks_under(peaks, p)
  #peaks <- peaks[order(peaks$scale, decreasing=TRUE),]
  if (length(all_p_idx) <= 1) {
    return(all_p_idx)
  } else {
    exclude_idx <- c()
    for(i in 1:length(all_p_idx)) {
      exclude_idx <- append(exclude_idx,
                            which_all_peaks_under(peaks, peaks[all_p_idx[i],]))
    }
    exclude_idx <- intersect(exclude_idx, all_p_idx)
    return(setdiff(all_p_idx, exclude_idx))
  }
}

model_BIC <- function(logLikelihood, k, n) {
  ## Bayesian Information Criterion for a model with a given logLikelihood,
  ## k degrees of freedom and n data points
  return(-2*logLikelihood + k*log(n))
}

model_AIC <- function(logLikelihood, k, n) {
  ## Akaike information criterion
  return(2*(k - logLikelihood)) # + 2*k*(k+1)/(n - k - 1))  AICc?
}

relative_log_likelihood <- function(pds, peaks1, peaks2 = NULL, other_peaks = NULL, use.AIC = TRUE, naverages=1) {
  ## Compare the model from peaks1 with peaks2 using the Akaike information criteria (default) or
  ## just the log-likelihood. Return the difference of both measurements. If negative, this
  ## means that 'peaks1' is a better model. Positive means 'peaks2' is a better model.
  if(!is.null(peaks2))
      if(nrow(peaks2) == 0)
          peaks2 <- NULL
  LL1 <- log_likelihood(pds, fit_model(pds, bind_rows(peaks1, other_peaks)), naverages=naverages)
  LL2 <- log_likelihood(pds, fit_model(pds, bind_rows(peaks2, other_peaks)), naverages=naverages)
  if (use.AIC) {
    N <- nrow(pds)                                      ## Number of points
    ## sinc^2 functions have 1 less dof than lorentzians
    k1 <- 3*nrow(peaks1) - sum(is.na(peaks1$linewidth))
    AIC1 <- model_AIC(LL1, k1, N)
    if (is.null(peaks2)) {
      AIC2 <- model_AIC(LL2, 0, N)
    } else {
      k2 <- 3*nrow(peaks2) - sum(is.na(peaks2$linewidth))
      AIC2 <- model_AIC(LL2, k2, N)
    }

    return((AIC1 - AIC2)/2)
  } else {
    return(LL2 - LL1)
  }
}

optimise_peaks <- function(pds, peaks, other_peaks = NULL,
                           min.linewidths = NULL, max.linewidths = NULL, naverages=1) {
  ## Make an MLE optimization of all (lorentzian) 'peaks' under 'pds'
  ## assuming a model of rbind(peaks, other_peaks)
  ## The optimisatino is done in (frequency, amplitude, linewidth)
  N <- nrow(peaks)
  deltanu <- diff(pds$frequency[1:2])

  #print(paste("Optimise_peaks bw: ", deltanu))

  peaks <- peaks %>% peaks_with_amplitude(deltanu)
  if (is.null(max.linewidths)) max.linewidths <- rep(diff(range(pds$frequency)), times = N)
  if (is.null(min.linewidths)) min.linewidths <- rep(deltanu/2, times = N)
  #print(paste("MIN LW: ", min.linewidths, " BW/2: ", deltanu/2))
  #stop()
  #max.linewidths <- 0.2
  pre_model <- fit_model(pds, other_peaks)
  peaks_mle <- optim(
      par = c(peaks$frequency, peaks$amplitude, peaks$linewidth),  # Start values
      fn = function(theta) {
          peaks <-
              data.frame(frequency = theta[1:N],
                         amplitude = theta[(N+1):(2*N)],
                         linewidth = theta[(2*N+1):(3*N)]) %>%
              peaks_with_height(deltanu)
          return(
              -log_likelihood(pds, add_peaks_to_model(pds = pds, model = pre_model, peaks = peaks), naverages=naverages))
      },
      gr = function(theta) {
          peaks <-data.frame(frequency = theta[1:N],
                             amplitude = theta[(N+1):(2*N)],
                             linewidth = theta[(2*N+1):(3*N)])
          return(LLgr(peaks = peaks, pds = pds, other_peaks = other_peaks))},
      method = "L-BFGS-B",
      ## Limits on:  frequency, amplitude, linewidth
      lower = c(rep(min(pds$frequency), times = N),
                rep(sqrt(pi*deltanu/2), times = N),
                min.linewidths),
      upper = c(rep(max(pds$frequency), times = N),
                rep(sqrt(pi*max(pds$power)*diff(range(pds$frequency))), times = N),
                max.linewidths),
      hessian = FALSE)
  #print("MLE STUFF")
  #print(peaks_mle)
  res <- tibble(
      frequency = peaks_mle$par[1:N],
      amplitude = peaks_mle$par[(N+1):(2*N)],
      linewidth = peaks_mle$par[(2*N+1):(3*N)])
  res$snr <- peaks$snr
  return(res %>% peaks_with_height(deltanu))
}

optimise_sinc2 <- function(pds, peak, other_peaks, naverages=1) {
  ## Takes a model (given by other_peaks) for the pds and adds to it a new
  ## peak shaped as a sinc2 function. Optimises its location and height
  ## to give the best likelihood possible.
  deltanu <- diff(pds$frequency[1:2])
  l <- which.min(abs(pds$frequency - peak$frequency))
  ## pds.local <- pds[(l-5):(l+5),]
  pds.local <- pds %>% filter(!is.na(frequency))
  pre_model <- fit_model(pds.local, other_peaks)

  peak_mle <- optim(
    par = c(peak$frequency, peak$height),  # Start values
    fn = function(theta) {
      f  <- theta[1]
      A0 <- theta[2]
      return(-log_likelihood(pds.local,
                             add_peaks_to_model(pds.local,
                                                pre_model,
                                                tibble(frequency = f,
                                                      #06/07/20 removed
                                                      #amplitude and inserted
                                                      #height instead
                                                           #amplitude = A0,
                                                           linewidth = NA,
                                                           height = A0)),
                              naverages=naverages))
    },
    method = "L-BFGS-B",
    ## Limits on frequency and height
    lower = c(min(pds.local$frequency), sqrt(pi*deltanu/2)),
    upper = c(max(pds.local$frequency), 2*sqrt(pi*max(pds.local$power))),
    hessian = FALSE,
    control = list(parscale = c(deltanu, sqrt(pi*deltanu))))

  return(tibble(frequency = peak_mle$par[1],
                    height = peak_mle$par[2],
                    linewidth = NA, snr = NA) %>% peaks_with_amplitude(deltanu))
}

best_fit_along_peak <- function(pds, peaks, p, use.AIC = TRUE, naverages=1) {
  ## Look at all the peaks below p, propose a model based on them and
  ## compare it with the model stemming from p and return the best
  ## model for that region.
  pw <- GAMMA_PROP * Lor_HWHM(p$scale)
  deltanu <- diff(pds$frequency[1:2])
  #print("p")
  #print(p)
  #print("pw")
  #print(pw)
  pds.local <- pds[pds$frequency > p$frequency - pw & pds$frequency < p$frequency + pw,]
  peaks_under_p <- peaks[which_peaks_under(peaks, p),]
  #print("peaks_under_p")
  #print(nrow(peaks_under_p))
  p <- optimise_peaks(pds.local, p, naverages=naverages)
  if (nrow(peaks_under_p) == 0) {
    return(p)
  } else {
    #print("peaks_under_p > 0")
    # Andres' original
    alt_peaks <- do.call(rbind,
                         Map(function(peak_idx) {
                           peak <- peaks_under_p[peak_idx,]
                           return(best_fit_along_peak(pds, peaks, peak, use.AIC = use.AIC, naverages=naverages))
                         }, 1:nrow(peaks_under_p))) %>%
      peaks_with_amplitude(deltanu)
    #print(alt_peaks)
    alt_peaks <- optimise_peaks(pds.local, alt_peaks, naverages=naverages)
    #print(alt_peaks)
    rel.log.likelihood <- relative_log_likelihood(pds.local, p, alt_peaks, use.AIC = use.AIC, naverages=naverages)
    #print(rel.log.likelihood)
    if (rel.log.likelihood < 0) {
      return(p)
    } else {
      return(alt_peaks)
    }

  }
}

all_branches_below_peak <- function(peaks, p) {
  ## All branches with a peak within GAMMA_PROP*FHHM of p.
  pw <- GAMMA_PROP * Lor_HWHM(p$scale)
  branches <- peaks[peaks$frequency > p$frequency - pw & peaks$frequency < p$frequency + pw, "branch"]
  return(unique(branches))
}

#original only one peak per branch used, this however can exclude some maxima per branch
find_best_peaks <- function(pds, peaks, use.AIC = TRUE, naverages=1) {
  n_branches <- length(unique(peaks$branch))
  branches_searched <- c()
  best_peaks <- peaks[NULL,]
  remaining_peaks <- peaks

  while (length(branches_searched) < n_branches) {
    current_peak <- tail(remaining_peaks[order(remaining_peaks$scale),], n=1)
    current_branch <- peaks %>%
                      filter(peaks$branch == current_peak$branch)
    if (nrow(current_branch) > 1){
      current_peak <- rbind(colMeans(current_branch),current_branch)
      current_peak <- current_peak %>% slice(1)
    }
    new_peaks  <- best_fit_along_peak(pds, peaks, current_peak, use.AIC = use.AIC, naverages=naverages)
    best_peaks <- rbind(best_peaks, new_peaks)
    branches_searched <- union(branches_searched,
                               all_branches_below_peak(peaks, current_peak))
    remaining_peaks <- peaks[!(peaks$branch %in% branches_searched),]

  }
  return(best_peaks)
}

peak_AIC_diff <- function(pds, peak, other_peaks, naverages=1) {
  all_peaks <- rbind(peak, other_peaks)
  m0 <- model_AIC(log_likelihood(pds, fit_model(pds, other_peaks), naverages=naverages),
                  k=3*nrow(other_peaks) - sum(is.na(other_peaks$linewidth)), n=nrow(pds))
  m1 <- model_AIC(log_likelihood(pds, fit_model(pds, all_peaks), naverages=naverages),
                  k=3*nrow(all_peaks) - sum(is.na(all_peaks$linewidth)), n=nrow(pds))
  return(m0-m1)
}

peaks_with_AIC <- function(peaks, pds, naverages=1) {
  peaks <- peaks %>% select_if(names(peaks) != "AIC")
  res <- do.call(rbind,
                 Map(function (i) {
                   return(cbind(peaks[i,],
                                data.frame(AIC = peak_AIC_diff(pds, peaks[i,], peaks[-i,], naverages=naverages))))
                 },1:nrow(peaks)))
  return(as_tibble(res))
}

peaks_with_low_AIC <- function(peaks, pds, minAIC= 2, naverages=1) {
  peaks['AIC1'] <- NA
  #peaks_low <- peaks %>% filter(AIC < minAIC)
  #peaks_low <- peaks_low %>% arrange(desc(AIC))
  peaks_high <- peaks %>% filter(AIC >= minAIC)
  peaks_resolved <- peaks %>% filter(linewidth > 0) %>%
                           arrange(desc(AIC))
  peaks_unresolved <- peaks %>% filter(is.na(linewidth)) %>%
                             arrange(desc(AIC))
  peaks <- bind_rows(peaks_resolved,peaks_unresolved)
  #res <- do.call(rbind,
  #               Map(function (i) {
  for(i in 1:nrow(peaks)){
    if (peaks[i,9] < minAIC){
        AIC1 = peak_AIC_diff(pds, peaks[i,], peaks_high, naverages=naverages)
        peaks[i,12] <- AIC1
        if (AIC1 >= minAIC){
            peaks_high <- rbind(peaks[i,],peaks_high)
            peaks_high <- peaks_high %>% arrange(frequency)
        }
    }
 }
 return(as_tibble(peaks))
}

peaks_AIC_addvalues <- function(peaks, pds, minAIC = 2, naverages=1) {
  peaks['AIC1'] <- NA
  maxAIC <- max(peaks$AIC, na.rm = TRUE)
  #peaks_low <- peaks %>% filter(AIC < minAIC)
  #peaks_low <- peaks_low %>% arrange(desc(AIC))
  #peaks_high <- peaks %>% filter(AIC >= minAIC)
  peaks_high <- peaks %>% filter(AIC >= maxAIC)
  peaks_resolved <- peaks %>% filter(linewidth > 0) %>%
                           arrange(desc(AIC))
  peaks_unresolved <- peaks %>% filter(is.na(linewidth)) %>%
                             arrange(desc(AIC))
  peaks <- bind_rows(peaks_resolved,peaks_unresolved)
  #res <- do.call(rbind,
  #               Map(function (i) {
  for(i in 1:nrow(peaks)){
    if (peaks[i,9] < maxAIC){
        AIC1 = peak_AIC_diff(pds, peaks[i,], peaks_high, naverages=naverages)
        peaks[i,12] <- AIC1
        if (AIC1 >= minAIC){
            peaks_high <- rbind(peaks[i,],peaks_high)
            peaks_high <- peaks_high %>% arrange(frequency)
        }
    }
 }
 return(as_tibble(peaks))
}

peaks_with_EFPp <- function(peaks, dw) {
  peaks <- peaks %>% select_if(names(peaks) != "EFPp")
  EFPp <-
    peaks %>%
    mutate(linewidth = linewidth / dw) %>%
    peaks_with_amplitude(1) %>%
    add_predictions(EFP_model) %>%
    #mutate(EFPp = exp(pred)) %>%
    mutate(EFPp = ifelse(amplitude > 10 | snr > 6,
                         0,
                         exp(pred))) %>%
    select(EFPp)
  return(bind_cols(peaks, EFPp))
}

lorentzian_peaks <- function(pds, min.snr = 3, linewidth.range = NULL,
                             pds.CWTTree = NULL, use.AIC = TRUE, naverages=1,
                             var.maxlw = FALSE) {
  ## Uses the tree-map of CWT values to identify lorentzian-shaped peaks
  ## Returns a data.frame with columns: frequency, linewidht, height
  if (is.null(pds.CWTTree)) {
    pds.SS      <- splus2R::signalSeries(data = pds$power, positions. = pds$frequency)
    pds.CWT     <- wavCWT(pds.SS, scale.range = CWT_scale_from_HWHM(linewidth.range))

    pds.CWTTree <- try(expr = wavCWTTree(pds.CWT), silent = TRUE)
    #print("PLOTTING")

    if(class(pds.CWTTree) == "try-error") return(NULL)
    # X11()
    # plot(pds.CWT, series=TRUE)
    # #plot(pds.CWT)#, series=TRUE)
    # #plot(pds.CWTTree, extrema=TRUE, add=TRUE)
    # prompt  <- "hit spacebar to close plots"
    # extra   <- "some extra comment"
    # capture <- tk_messageBox(message = prompt, detail = extra)
  } else if (!is(pds.CWTTree, "wavCWTTree")) {
    stop("pds.CWTTree must be of class wavCWTTree")
  }

  if (is.null(linewidth.range)) {
    scale.range <- attr(pds.CWTTree, "scale")[c(1, length(attr(pds.CWTTree, "scale")))]
  } else if (linewidth.range[1] > linewidth.range[2]) {
    stop("linewidth.range[1] should be greater than linewidth.range[2]")
  } else {
    scale.range <- CWT_scale_from_HWHM(linewidth.range)
  }
  #pds.Peaks   <- cwtPeaks(pds.CWTTree, min.snr = min.snr, scale.range = 0.8*scale.range)
  pds.Peaks   <- cwtPeaks(pds.CWTTree, min.snr = min.snr, scale.range = scale.range, var.maxlw = var.maxlw)

  pds.peakFind <- as_tibble(find_best_peaks(pds, pds.Peaks, use.AIC = use.AIC, naverages=naverages))

  #if(nrow(pds.peakFind) == 0) return(NULL)
  # 20/04/2020 Changed from !is.null to if doesn't have 0 rows
  #if(!is.null(pds.peakFind)) {
  if(nrow(pds.peakFind) != 0){
    ## Calculate a lorentzian parameters and subselect what we need
    pds.peakFind <-
      pds.peakFind %>%
      as_tibble() %>%
      select(frequency, linewidth, height, snr) %>%
      peaks_with_amplitude(diff(pds$frequency[1:2])) %>%
      peaks_with_AIC(pds, naverages=naverages)
  } else {
    pds.peakFind <- tibble(frequency=double(),
                              amplitude=double(),
                              linewidth=double(),
                              snr=double(),
                              height=double())
  }
  return(pds.peakFind %>%
                    select(frequency, linewidth, height, snr))
}

unresolved_candidates <- function(pds, peaks, p) {
  ## Looks at the residuals of 'pds' assuming a model of 'peaks' and
  ## identifies points with a probability lower than 'p' of occuring
  ## due to noise. Returns a list with frequencies for such points
  ## ordered by descending power in the residuals heights.
  residuals <- pds$power / fit_model(pds, peaks)
  min.height <- log(1/p)
  res <- data.frame(frequency = pds$frequency[residuals > min.height],
                    power     = residuals[residuals > min.height])
  res <- res[order(-res$power),]
  return(res$frequency)
}

narrow_peak <- function(pds, frequency, other_peaks, naverages=1) {
  ## Fits a lorentzian and a sinc^2 around "frequency"
  ## on the pds. Returns the best fit or NULL depending on
  ## the AIC values computed.
  j <- which(pds$frequency == frequency)
  H <- pds[[j, "power"]] - 1
  if (length(H) == 0) stop("The 'frequency' value must be present in 'pds'")
  minj <- max(1, j-10)
  maxj <- min(nrow(pds), j+10)
  pds.local <-
      pds[minj:maxj,] %>%
      filter(!is.na(frequency))
  deltanu <- diff(pds.local$frequency[1:2])
  ## Lorentzian function
  p.lorentzian <- tibble(frequency = frequency, height = H, linewidth = deltanu, snr = NA)
  p.lorentzian <- optimise_peaks(pds   = pds.local, peaks = p.lorentzian, other_peaks = other_peaks,
                                 min.linewidths = c(deltanu/2), naverages=naverages)
  ## Sinc^2 function
  p.sinc2 <-
    tibble(frequency = frequency, height = H, linewidth = NA, snr = NA) %>%
    peaks_with_amplitude(deltanu)
  p.sinc2 <- optimise_sinc2(pds = pds.local, peak = p.sinc2, other_peaks = other_peaks)

  ## Make the comparisons
  if (relative_log_likelihood(pds.local, p.lorentzian, p.sinc2, other_peaks, naverages=naverages) < 0) {
    p <- p.lorentzian
  } else {
    p <- p.sinc2
  }
  if (relative_log_likelihood(pds.local, p, other_peaks = other_peaks, naverages=naverages) < 0) {
    return(p)
  } else {
    return(NULL)
  }
}

add_unresolved_peaks <- function(peaks, pds, p, naverages=1) {
  u_candidates <- unresolved_candidates(pds, peaks, p)
  if(length(u_candidates) != 0) {
    for (i in 1:length(u_candidates)) {
      new_peak <- narrow_peak(pds, u_candidates[i], peaks, naverages=naverages)
      peaks <- bind_rows(peaks, new_peak)
    }
  }

  return(peaks)
}

fit_all_as_narrow <- function(pds, p, naverages=1) {
    ## Fit all the peaks inside pds as 'narrow' peaks
    add_more <- TRUE
    peaks <- tibble()
    while(add_more) {
        unr_locations <-
            pds %>%
            mutate(power = power / fit_model(pds = pds, peaks = peaks)) %>%
            filter(frequency %in% unresolved_candidates(pds = pds, peaks = peaks, p = p)) %>%
            arrange(-power)
        if(nrow(unr_locations) == 0) {
            add_more <- FALSE
        } else {
            new_peak <- narrow_peak(pds = pds,
                                    frequency = unr_locations$frequency[1],
                                    other_peaks = NULL, naverages=naverages)
            if(is.null(new_peak)) {
                add_more <- FALSE
            } else {
                peaks <- bind_rows(peaks, new_peak)
            }
        }
    }
    return(peaks)
}

check_if_unresolved <- function(peak, pds, min.snr, p, use.AIC = TRUE, naverages=1) {
    ## The function 'lorentzian_peaks' has trouble when the peaks
    ## are too narrow. Here we take a 'peak' and try to find if there
    ## is at least an unresolved peak that better describes the pds
    ## region close to 'peak'
    if(nrow(peak) == 0) return(NULL)
    if(nrow(peak) > 1) stop("'peak' must be a single row data_frame")
    if(is.na(peak$linewidth)) return(peak)   # 'peak' is already unresolved...
    lwdFactor <- 0.5        # Check for peaks narrower than peak$linewidht * lwdFactor
    pw <- GAMMA_PROP * peak$linewidth
    deltanu <- diff(pds$frequency[1:2])
    pds.local <-
        pds %>%
        filter(frequency > peak$frequency - pw,
               frequency < peak$frequency + pw)
    ## 'peak' is too narrow to continue
    if (deltanu >= lwdFactor*peak$linewidth)
        return(peak)

    ## Try to find a set of resolved peaks that are narrower
    ## than peak$linewidth
    p1 <- lorentzian_peaks(pds = pds.local, min.snr = min.snr,
                           linewidth.range = c(deltanu/2, lwdFactor*peak$linewidth),
                           naverages=naverages)
    ## If we do not have any resolved peaks narrower than peak, maybe there
    ## are 1 or more unresolved peaks that better describe the PDS close to it
    if (is.null(p1)) {
        p2 <- fit_all_as_narrow(pds = pds.local, p = p, naverages=naverages)
        if (relative_log_likelihood(pds = pds.local,
                                    peaks1 = peak, peaks2 = p2,
                                    use.AIC = use.AIC, naverages=naverages) > 2) {
            return(p2)
        } else {
            return(peak)
        }
    }
    ## Check if any of these peaks can be better described by unresolved peaks
    p1 <-
        do.call(
            bind_rows,
            Map(function(i) {
                return(check_if_unresolved(peak = p1[i,], pds = pds, min.snr = min.snr, p = p, naverages=naverages))
            }, 1:nrow(p1)))

    ## If we at least one resolved peak maybe there are unresolved peaks nearby
    p2 <- fit_all_as_narrow(
        pds = pds.local %>% mutate(power = power / fit_model(pds = ., peaks = p1)), p = p, naverages=naverages)
    if (relative_log_likelihood(pds = pds.local,
                                peaks1 = peak, peaks2 = bind_rows(p1,p2),
                                use.AIC = use.AIC, naverages=naverages) > 2) {
        return(bind_rows(p1,p2))
    } else {
        return(peak)
    }

    ## If we made it all the way here, something went wrong...
    stop("Something went wrong in 'check_if_unresolved'")
}

check_if_peak_unresolved <- function(peak, pds, naverages=1) {
    ## The function 'lorentzian_peaks' has trouble when the peaks
    ## are too narrow. Here we take a 'peak' and try to see if an unresolved
    ## peak fits better
    deltanu <- diff(pds$frequency[1:2])
    # 29/06/2020 changed from return(NULL) to return(peak)

    if(nrow(peak) == 0) return(peak)
    if(nrow(peak) > 1) stop("'peak' must be a single row data_frame")
    if(is.na(peak$linewidth)) return(peak)   # 'peak' is already unresolved...
    if(peak$linewidth > 1 * deltanu) return(peak) # Peak resolved and not tested
    p <- narrow_peak_check(pds, peak, naverages=1)

    #stop()
    if(is.null(p)){
      return(peak)
    } else {
      return(p)
    }
}

narrow_peak_check <- function(pds, peak, naverages=1) {
  ## Fits a lorentzian and a sinc^2 around "frequency"
  ## on the pds. Returns the best fit or NULL depending on
  ## the AIC values computed.
  pw <- GAMMA_PROP * peak$linewidth

  pds.local <-
      pds %>%
      filter(frequency > peak$frequency - pw,
              frequency < peak$frequency + pw)

  deltanu <- diff(pds.local$frequency[1:2])
  ## Lorentzian function
  p.lorentzian <- tibble(frequency = peak$frequency,
                         height = peak$height,
                         linewidth = peak$linewidth,
                         snr = NA) %>%
                         # 06/07/20 changed peaks_with_height to
                         #peaks_with_amplitude
                         #peaks_with_height(deltanu)
                         peaks_with_amplitude(deltanu)
  p.lorentzian <- optimise_peaks(pds   = pds.local, peaks = p.lorentzian, other_peaks = NULL,
                                 min.linewidths = c(deltanu/2), naverages=naverages)

  ## Sinc^2 function
  p.sinc2 <- tibble(frequency = peak$frequency,
                    height = peak$height,
                    linewidth = NA,
                    snr = NA) %>%
                        # peaks_with_height(deltanu)
                    # %>%
              peaks_with_amplitude(deltanu)
  p.sinc2 <- optimise_sinc2(pds = pds.local, peak = p.sinc2, other_peaks = NULL)

  ## Make the comparisons

  if (relative_log_likelihood(pds.local, p.lorentzian, p.sinc2, other_peaks = NULL, naverages=naverages) < 0) {
    p <- p.lorentzian
  } else {
    p <- p.sinc2
  }

  return(p)

  #if (relative_log_likelihood(pds.local, p, other_peaks = NULL, naverages=naverages) < 0) {
  #  return(p)
  #} else {
  #  return(NULL)
  #}
}

peak_find <- function(pds, min.snr = 3, p=0.0001, linewidth.range = NULL, use.AIC = TRUE,
                      find.first.set=TRUE, naverages=1, var.maxlw = FALSE) {
    ## Strategy:
    ## 1. Find the resolved peaks (p1)
    ## 2. Normalise by p1 to get the unidentified "narrow" peaks (p2)
    ## 3. Normalise by p2 (pds2) and find again the resolved peaks (p3)
    ## 4. Peak identification is now p = p2 + p3
    ## 5. Check if there are "narrow" peaks confounding the p idenfitications
    ## 19/11/19 Added extra keyword find.resolved.only which allows us to only find resolved peaks, ignoring those that are unresolved
    deltanu <- diff(pds$frequency[1:2])

    if(find.first.set == TRUE){
        print("Finding all peaks")
        peaks <-
            lorentzian_peaks(pds = pds, min.snr = min.snr,
                             linewidth.range = linewidth.range,
                             pds.CWTTree = NULL, use.AIC = use.AIC, naverages=naverages,
                             var.maxlw = var.maxlw)%>%
                             add_unresolved_peaks(pds=pds, p = p, naverages=naverages)

        if(is.null(peaks) || nrow(peaks) == 0) return(NULL)
        # AIC computed in lorentzian peaks and again here, why?!?!

        #peaks <-
        #    do.call(bind_rows,
        #            Map(function (i) {
        #                return(check_if_peak_unresolved(peak = peaks[i,], pds = pds))
        #            }, 1:nrow(peaks)))

        #stop()

        peaks <- peaks %>%
                       peaks_with_AIC(pds, naverages=naverages)
        #stop()

    } else {
      peaks <-
          lorentzian_peaks(pds = pds, min.snr = min.snr,
                          linewidth.range = linewidth.range,
                          pds.CWTTree = NULL, use.AIC = use.AIC, naverages=naverages,
                          var.maxlw = var.maxlw) %>%
                          add_unresolved_peaks(pds=pds, p = p, naverages=naverages) %>%
                          peaks_with_amplitude(deltanu) %>%
                          peaks_with_AIC(pds, naverages=naverages)
          #filter(AIC > 2)

      peaks <-
              do.call(bind_rows,
                      Map(function (i) {
                          return(check_if_peak_unresolved(peak = peaks[i,], pds = pds))
                      }, 1:nrow(peaks)))

    }
    if(is.null(peaks) || nrow(peaks) == 0) return(NULL)
    #peaks <-
    #    do.call(bind_rows,
    #            Map(function (i) {
    #                return(check_if_unresolved(peak = peaks[i,], pds = pds,
    #                                           min.snr = min.snr, p = p))
    #            }, 1:nrow(peaks)))

    # AIC calculated AGAIN! WHY?!?!?!?!
    peaks <-
        peaks %>%
        arrange(frequency) %>%
        peaks_with_AIC(pds = pds, naverages=naverages) %>%
        #filter(AIC > 2) %>%
        peaks_with_amplitude(deltanu)# %>%
        #peaks_with_EFPp(deltanu) %>%
        #mutate(EFP = EFPp*nrow(pds))
    #stop()
    return(peaks)
}

peaks_MLE2 <- function(peaks, pds, maxLWD, other_peaks = NULL, hessian = FALSE, final_fit = FALSE, final_fit_factor=0.1, naverages=1) {
  ## MLE fit of 'peaks' using (frequency, amplitude, linewidth) as free parameters.
  ## If 'hessian' is TRUE it will also compute the uncertainties.
  ## If 'final fit' is TRUE then will change upper and lower bounds on MLE fit to be
  ## tighter to reflect the fact that this is the final fit.
  ## final_fit_factor controls how wider the upper and lower bounds of the final_fit are
  ## as a function of the given parameter
  #if(max(peaks$linewidth, na.rm = TRUE) >= maxLWD)
  #  stop("One or more of the peaks has a linewidth larger than or equal to maxLWD")
  #else{
  #  print("All peaks have linewidth less than maxLWD, continuing ...")
  #}
  peaks <- peaks %>% arrange(frequency)
  N <- nrow(peaks)
  deltanu <- diff(pds$frequency[1:2])
  peaks.res <- peaks %>% filter(!is.na(peaks$linewidth))
  peaks.unr <- peaks %>% filter( is.na(peaks$linewidth))
  n.res <- nrow(peaks.res)
  n.unr <- nrow(peaks.unr)
  maxPDS <- max(pds$power)

  ## Limits on the frequencies for the peaks
  if (final_fit == TRUE){
    # Start with +/- 10%
    peaks.freq.low <- peaks$frequency - final_fit_factor*peaks$frequency
    peaks.freq.upp <- peaks$frequency + final_fit_factor*peaks$frequency
  } else {
    # Not final fit so set up upper and lower limits as usual
    if (N == 1) {
      peaks.freq.low <- peaks$frequency - peaks$linewidth
      peaks.freq.upp <- peaks$frequency + peaks$linewidth
    } else {
      # 20/01/2020 Maybe change from 0 to min(frequency), as 0 will almost never be in range as cutting down to +/- 3 sigma_env
      # 31/01/2020 Changed 0.1*diff(frequency) to 0.5*diff(frequency) to ensure frequency parameter allowed to move enough
      # during fitting.
      peaks.freq.low <- c(max(0,peaks$frequency[1] - 0.5*diff(peaks$frequency)[1]),
                          peaks$frequency[2:N] - 0.5*(diff(peaks$frequency)/2))

      peaks.freq.upp <- c(peaks$frequency[1:N-1] + 0.5*(diff(peaks$frequency)/2),
                          min(peaks$frequency[N] + 0.5*diff(peaks$frequency)[N-1], max(pds$frequency)))
    }

  }

  # Set up rest of limits outside optimisation function for clarity
  # Amplitude lower and upper limits - resolved
  ifelse(final_fit == FALSE, peaks.res.amp.low <- rep(sqrt(pi*deltanu), times = n.res),
                             peaks.res.amp.low <- peaks.res$amplitude - final_fit_factor * peaks.res$amplitude)
  ifelse(final_fit == FALSE, peaks.res.amp.upp <- rep(sqrt(pi*maxPDS*maxLWD), times = n.res),
                             peaks.res.amp.upp <- peaks.res$amplitude + final_fit_factor * peaks.res$amplitude)
  # Linewidth lower and upper limits -resolved
  ifelse(final_fit == FALSE, peaks.res.lwd.low <- rep(deltanu/2, times = n.res),
                             ifelse(peaks.res$linewidth - final_fit_factor * peaks.res$linewidth < deltanu,
                                    peaks.res.lwd.low <- deltanu/2,
                                    peaks.res.lwd.low <- peaks.res$linewidth - final_fit_factor * peaks.res$linewidth))
  ifelse(final_fit == FALSE, peaks.res.lwd.upp <- rep(maxLWD, times = n.res),
                             peaks.res.lwd.upp <- peaks.res$linewidth + final_fit_factor * peaks.res$linewidth)

  # Amplitude lower and upper limits - unresolved
  ifelse(final_fit == FALSE, peaks.unr.amp.low <- rep(sqrt(pi*deltanu), times = n.unr),
                             peaks.unr.amp.low <- peaks.unr$amplitude - final_fit_factor * peaks.unr$amplitude)
  ifelse(final_fit == FALSE, peaks.unr.amp.upp <- rep(pi*maxPDS*maxLWD, times = n.unr),
                             peaks.unr.amp.upp <- peaks.unr$amplitude + final_fit_factor * peaks.unr$amplitude)

  #print(peaks.res.amp.low)
  #print(peaks.res.amp.upp)
  #quit()
  peaks.mle <- optim(
    ## Start values
    par = c(peaks.res$frequency, peaks.res$amplitude, peaks.res$linewidth,
            peaks.unr$frequency, peaks.unr$amplitude),
    ## Function to optimise
    fn = function(theta) {
      if(n.res > 0) {
        pks.res <- data.frame(frequency = theta[1:n.res],
                              amplitude = theta[(1+n.res):(2*n.res)],
                              linewidth = theta[(2*n.res + 1):(3*n.res)])
      } else {
        pks.res <- NULL
      }
      if (n.unr > 0) {
        pks.unr  <- data.frame(frequency = theta[(3*n.res + 1):(3*n.res + n.unr)],
                               amplitude = theta[(3*n.res + n.unr + 1):(3*n.res + 2 * n.unr)])
        pks.unr$linewidth <- NA
      } else {
        pks.unr <- NULL
      }
      pks <- bind_rows(pks.res, pks.unr)
      pks <- pks %>% peaks_with_height(deltanu)
      return(-log_likelihood(pds, fit_model(pds, bind_rows(pks, other_peaks)), naverages=naverages))
    },
    gr = function(theta) {
     if(n.res > 0) {
       pks.res <- data.frame(frequency = theta[1:n.res],
                             amplitude = theta[(1+n.res):(2*n.res)],
                             linewidth = theta[(2*n.res + 1):(3*n.res)])
     } else {
       pks.res <- NULL
     }
     if (n.unr > 0) {
       pks.unr  <- data.frame(frequency = theta[(3*n.res + 1):(3*n.res + n.unr)],
                              amplitude = theta[(3*n.res + n.unr + 1):(3*n.res + 2 * n.unr)])
       pks.unr$linewidth <- NA
     } else {
       pks.unr <- NULL
     }
     pks <- bind_rows(pks.res, pks.unr)
     pks <- pks %>% peaks_with_height(deltanu)
     return(LLgr(peaks = pks, pds = pds, other_peaks = other_peaks))
    },
    method = "L-BFGS-B",
    ## Limits on the parameters
    lower = c(
      ## Resolved
      peaks.freq.low[!is.na(peaks$linewidth)],
      # 31/01/2020 Added in ifelse, so if final_fit is FALSE then set limits as before
      # otherwise +/- 10%
      peaks.res.amp.low,
      #rep(sqrt(pi*deltanu), times = n.res),
      # 31/01/2020 Added if else as above with ifelse making sure that linewidth bound doesn't go below bin width for resolved peaks
      #ifelse(final_fit == FALSE, rep(deltanu, times = n.res),
      #                           ifelse(peaks[!is.na(peaks$linewidth)]$linewidth - final_fit_factor*peaks[!is.na(peaks$linewidth)]$linewidth < deltanu,
      #                                  deltanu,
      #                                  peaks[!is.na(peaks$linewidth)]$linewidth - final_fit_factor*peaks[!is.na(peaks$linewidth)]$linewidth)),
      peaks.res.lwd.low,
      #rep(deltanu, times = n.res),
      ## Unresolved
      peaks.freq.low[is.na(peaks$linewidth)],
      # 16/01/2020 Why is times=n.res not n.unr?
      # 31/01/2020 Changed times=n.res to times=n.unr and added final_fit ifelse statement.
      #ifelse(final_fit == FALSE, rep(sqrt(pi*deltanu), times = n.unr),
      #                           peaks[is.na(peaks$linewidth)]$amplitude - final_fit_factor*peaks[is.na(peaks$linewidth)]$amplitude)),
      peaks.unr.amp.low),
      #rep(sqrt(pi*deltanu), times = n.unr)),
    upper = c(
      ## Resolved
      peaks.freq.upp[!is.na(peaks$linewidth)],
      # 31/01/2020 Adding ifelse statement so that limits set according to initial parameters if final_fit is TRUE
      #ifelse(final_fit == FALSE, rep(sqrt(pi*maxPDS*maxLWD), times = n.res),
      #                           peaks[!is.na(peaks$linewidth)]$amplitude + final_fit_factor*peaks[!is.na(peaks$linewidth)]$amplitude),
      peaks.res.amp.upp,
      #rep(sqrt(pi*maxPDS*maxLWD), times = n.res),
      # 31/01/2020 Adding ifelse statement for linewidths
      #ifelse(final_fit == FALSE, rep(maxLWD, times = n.res),
      #                           peaks[!is.na(peaks$linewidth)]$linewidth + final_fit_factor*peaks[!is.na(peaks$linewidth)]$linewidth),
      peaks.res.lwd.upp,
      #rep(maxLWD, times = n.res),
      ## Unresolved
      peaks.freq.upp[is.na(peaks$linewidth)],
      # 31/01/2020 Adding ifelse statement for amplitudes and changing times=n.res to times=n.unr
      # Should change to sqrt(pi*maxPDS*deltanu), as shouldn't be constraining max amplitude using resolved max linewidth!
      #ifelse(final_fit == FALSE, rep(pi*maxPDS*maxLWD, times = n.unr),
      #                           peaks[is.na(peaks$linewidth)]$amplitude + final_fit_factor*peaks[is.na(peaks$linewidth)]$amplitude)),
      peaks.unr.amp.upp),
      #rep(pi*maxPDS*maxLWD, times = n.unr)),
    control = list(
      parscale = c(
        ## Resolved
        rep(deltanu,          times = n.res),
        rep(sqrt(pi*deltanu), times = n.res),
        rep(deltanu,          times = n.res),
        ## Unresolved
        rep(deltanu,          times = n.unr),
        rep(sqrt(pi*deltanu), times = n.unr))),
    hessian = hessian)
  pars <- peaks.mle$par
  if (n.res > 0) {
    pks.res  <- data.frame(frequency = pars[1:n.res],
                           amplitude = pars[(1+n.res):(2*n.res)],
                           linewidth = pars[(2*n.res + 1):(3*n.res)])
  } else {
    pks.res <- NULL
  }
  if (n.unr > 0) {
    pks.unr  <- data.frame(frequency = pars[(3*n.res + 1):(3*n.res + n.unr)],
                           amplitude = pars[(3*n.res + n.unr + 1):(3*n.res + 2 * n.unr)])
    pks.unr$linewidth <- NA
  } else {
    pks.unr <- NULL
  }
  pks <-
    bind_rows(pks.res, pks.unr) %>%
    peaks_with_height(deltanu) %>%
    peaks_with_AIC(pds, naverages=naverages)
  return(list(fit = peaks.mle, peaks = pks))
}

peak_res_sd <- function(peak, pds, maxLWD, other_peaks = NULL, final_fit = FALSE, final_factor = 0.1, naverages=1) {
  ## Uncertainty calculation for a resolved 'peak' using
  ## (frequency, amplitude, linewidth) as free parameters.
  deltanu <- diff(pds$frequency[1:2])

  # Set up upper and lower bounds
  ifelse(final_fit == FALSE, freq.low <- peak$frequency - peak$linewidth,
                             freq.low <- peak$frequency - final_fit_factor*peak$frequency)
  ifelse(final_fit == FALSE, freq.upp <- peak$frequency + peak$linewidth,
                             freq.upp <- peak$frequency + final_fit_factor*peak$frequency)

  ifelse(final_fit == FALSE, amp.low <- sqrt(pi*deltanu/2),
                             amp.low <- peak$amplitude - final_fit_factor * peak$amplitude)
  ifelse(final_fit == FALSE, amp.upp <- peak$amplitude * 2,
                             amp.upp <- peak$amplitude + final_fit_factor * peak$amplitude)
  # Linewidth lower and upper limits -resolved
  #27/07/2021 lowered lower limit to under resolved/unresolved boundary
  # This is because the hard cut caused negative values in inverse of hessian matrix
  # due to the hard cut.
  ifelse(final_fit == FALSE, lwd.low <- deltanu/4,
                             ifelse(peak$linewidth - final_fit_factor * peak$linewidth < deltanu,
                                    #27/06/2021 lowered lower limit to under resolved/unresolved boundary
                                    # This is because the hard cut caused negative values in inverse of hessian matrix
                                    # due to the hard cut.
                                    lwd.low <- deltanu/4,
                                    lwd.low <- peak$linewidth - final_fit_factor * peak$linewidth))
  ifelse(final_fit == FALSE, lwd.upp <- maxLWD,
                             lwd.upp <- peaks$linewidth + final_fit_factor * peak$linewidth)

  peak.optim <- optim(
    par = c(peak$frequency, peak$amplitude, peak$linewidth),
    ## Function to optimise
    fn = function(theta) {
      pk <- data.frame(frequency = theta[1],
                       amplitude = theta[2],
                       linewidth = theta[3])
      pk <- pk %>% peaks_with_height(deltanu)
      return(-log_likelihood(pds, fit_model(pds, bind_rows(pk, other_peaks)), naverages=naverages))
    },
    #gr = function(theta) {
    #    pk <- data.frame(frequency = theta[1],
    #                     amplitude = theta[2],
    #                     linewidth = theta[3])
    #    pk <- pk %>% peaks_with_height(deltanu)
    #    return(LLgr(peaks = pk, pds = pds, other_peaks = other_peaks))
    #},
    control = list(parscale = c(deltanu, sqrt(pi*deltanu), deltanu),
                   maxit=500),
    hessian = TRUE,
    method = "L-BFGS-B",
    ## Limits on the parameters
    lower = c(freq.low,
              amp.low,
              lwd.low),
    upper = c(freq.upp,
              amp.upp,
              lwd.upp)
  )
  peak.hess <- peak.optim$hessian
  #print(peak.hess)
  # numHessian <- hessian(func =
  #                       function(theta) {
  #                                     pk <- data.frame(frequency = theta[1],
  #                                                     amplitude = theta[2],
  #                                                     linewidth = theta[3])
  #                                     pk <- pk %>% peaks_with_height(deltanu)
  #                                     return(-log_likelihood(pds, fit_model(pds, bind_rows(pk, other_peaks)), naverages=naverages))
  #                                   },
  #                       peak.optim$par)
  peak.pars <- peak.optim$par
  sd = c(-99,-99,-99)
  try(sd <- sqrt(diag(solve(peak.hess))),TRUE)
  res <- data.frame(
      frequency = peak.pars[1], amplitude = peak.pars[2], linewidth = peak.pars[3],
      frequency_sd = sd[1], amplitude_sd = sd[2], linewidth_sd = sd[3], convergence=peak.optim$convergence)
  return(res)
}

peak_unr_sd <- function(peak, pds, maxLWD, other_peaks = NULL, final_fit = FALSE, final_factor = 0.1, naverages=1) {
    ## Uncertainty calculation for 'peak' using (frequency, amplitude) as
    ## free parameters.
    deltanu <- diff(pds$frequency[1:2])

    # Set up upper and lower bounds
    ifelse(final_fit == FALSE, freq.low <- peak$frequency - deltanu,
                               freq.low <- peak$frequency - final_fit_factor*peak$frequency)
    ifelse(final_fit == FALSE, freq.upp <- peak$frequency + deltanu,
                               freq.upp <- peak$frequency + final_fit_factor*peak$frequency)
    ifelse(final_fit == FALSE, amp.low <- sqrt(pi*deltanu),
                               amp.low <- peak$amplitude - final_fit_factor * peak$amplitude)
    ifelse(final_fit == FALSE, amp.upp <- peak$amplitude*2,
                               amp.upp <- peak$amplitude + final_fit_factor * peak$amplitude)

    peak.optim <- optim(
        par = c(peak$frequency, peak$amplitude),
        ## Function to optimise
        fn = function(theta) {
            pk <- data.frame(frequency = theta[1],
                             amplitude = theta[2],
                             linewidth = NA)
            pk <- peaks_with_height(pk, deltanu)
            return(-log_likelihood(pds, fit_model(pds, bind_rows(pk, other_peaks)), naverages=naverages))
        },
        #gr = function(theta) {
        #    pk <- data.frame(frequency = theta[1],
        #                     amplitude = theta[2],
        #                     linewidth = NA)
        #    pk <- peaks_with_height(pk, deltanu)
        #    return(LLgr(peaks = pk, pds = pds, other_peaks = other_peaks))
        #},
        control = list(parscale = c(deltanu, sqrt(pi*deltanu))),
        method = "L-BFGS-B",
        ## Limits on the parameters
        lower = c(freq.low,
                  amp.low,
                  deltanu),
        upper = c(freq.upp,
                  amp.upp,
                  peak$linewidth * 2),
        hessian = TRUE)
    peak.hess <- peak.optim$hessian
    peak.pars <- peak.optim$par
    sd <- sqrt(diag(solve(peak.hess)))
    res <- data.frame(
        frequency = peak.pars[1], amplitude = peak.pars[2], linewidth = NA,
        frequency_sd = sd[1], amplitude_sd = sd[2], linewidth_sd = NA, convergence=peak.optim$convergence)
    #stop()
    return(res)
}

peaks_with_MLE_sd <- function(peaks, pds, maxLWD, naverages=1) {
  res <- do.call(bind_rows,
                 Map(function(i) {
                   p <- peaks[i,]
                   others <- peaks[-i,]
                   fn <- if (is.na(p$linewidth)) peak_unr_sd else peak_res_sd
                   res <- fn(peak = p, pds = pds, other_peaks = others, maxLWD = maxLWD, naverages=naverages)
                   return(res)
                 }, seq_len(nrow(peaks))))
  return(res)
}

peaks_with_MLE_sd2 <- function(peaks, pds, maxLWD, naverages=1) {
  res <- do.call(bind_rows,
                 Map(function(i) {
                   p <- peaks[i,]
                   others <- peaks[-i,]
                   MLEfit <- peaks_MLE2(peaks = p, pds = pds, maxLWD = maxLWD,
                                        other_peaks = others, hessian = TRUE, naverages=naverages)
                   p.new <- MLEfit$peaks
                   sd <- sqrt(diag(solve(MLEfit$fit$hessian)))
                   if (is.na(p$linewidth)) {
                     sd.res <- tibble(frequency_sd = sd[1], amplitude_sd = sd[2], linewidth_sd = NA)
                   } else {
                     sd.res <- tibble(frequency_sd = sd[1], amplitude_sd = sd[2], linewidth = sd[3])
                   }
                   return(cbind(p, sd.res))
                 }, seq_len(nrow(peaks))))
  return(res)
}

peaks_MLE_sd <- function(peaks, pds, maxLWD = 0.5, naverages=1) {
  deltanu <- diff(pds$frequency[1:2])
  ## Optimise all the peaks simultaneously
  mle.fit <- peaks_MLE2(peaks, pds, maxLWD, other_peaks = NULL, naverages=naverages)
  ## Get the error estimates one at a time
  print("Getting error estimates")
  peaks <- peaks_with_MLE_sd(peaks = mle.fit$peaks, pds = pds, maxLWD = maxLWD, naverages=naverages)
  ## TODO  if a resolved peak gets a linewidth~deltanu during the optimisation
  ## we should consider using a sinc^2
  peaks <-
    peaks %>%
    peaks_with_height(deltanu) %>%
    peaks_with_AIC(pds, naverages=naverages) %>%
    arrange(frequency)
  return(peaks)
}

## Convert the peaks from the background-normalised PDS to the original PDS

sLor <- function(nu, A, b, c) {
  return(A / (1 + (nu / b)^c))
}

bg_component <- function(nu, A, b, c, nuNyq) {
  return(sinc(pi * nu / (2 * nuNyq))^2 * sLor(nu, A, b, c))
}

bg_model <- function(nu, theta, nuNyq) {
  bg <- theta$Pn
  bg <- bg + bg_component(nu, theta$A1, theta$b1, 4, nuNyq)
  bg <- bg + bg_component(nu, theta$A2, theta$b2, 4, nuNyq)
  bg <- bg + bg_component(nu, theta$A3, theta$b3, 4, nuNyq)
  return(bg)
}

peak_rescale_height <- function(peak, theta, nuNyq) {
  nu <- peak$frequency
  r <- bg_model(nu, theta, nuNyq) / sinc(pi * nu / (2 * nuNyq))^2
  peak$height <- peak$height * r
  if("height_sd" %in% names(peak))
    peak$height_sd <- peak$height_sd * r
  return(peak)
}

peaks_rescale_height <- function(peaks, theta, nuNyq, deltanu) {
  ## Convert the height and amplitude of 'peaks' from the background-normalised PDS to
  ## the not normalised PDS.
  ## Theta is a data frame with the background parameters (Pn, A1, b1, A2, b2, A3, b3)
  res <- do.call(rbind,
                 Map(function (i) {
                   p <- peak_rescale_height(peaks[i,], theta, nuNyq)
                   return(right_join(peaks[i,], p, by = names(peaks)))
                 },1:nrow(peaks)))
  res <- peaks_with_amplitude(res, deltanu)
  return(as_tibble(res))
}

peaks_rescale_height <- function(peaks, theta, nuNyq, deltanu) {
  ## Convert the height and amplitude of 'peaks' from the background-normalised PDS to
  ## the not normalised PDS.
  ## Theta is a data frame with the background parameters (Pn, A1, b1, A2, b2, A3, b3)
  res <- do.call(rbind,
                 Map(function (i) {
                   return(peak_rescale_height(peaks[i,], theta, nuNyq))
                 },1:nrow(peaks)))
  res <- peaks_with_amplitude(res, deltanu)
  return(as_tibble(res))
}

peak_rescale_amplitude <- function(peak, theta, nuNyq) {
  nu <- peak$frequency
  r <- bg_model(nu, theta, nuNyq) / sinc(pi * nu / (2 * nuNyq))^2
  peak$amplitude <- sqrt(r) * peak$amplitude
  if("amplitude_sd" %in% names(peak))
    peak$amplitude_sd <- sqrt(r) * peak$amplitude_sd
  return(peak)
}

peaks_rescale_amplitude <- function(peaks, theta, nuNyq, deltanu) {
  ## Convert the height and amplitude of 'peaks' from the background-normalised PDS to
  ## the not normalised PDS.
  ## Theta is a data frame with the background parameters (Pn, A1, b1, A2, b2, A3, b3)
  res <- do.call(rbind,
                 Map(function (i) {
                   return(peak_rescale_amplitude(peaks[i,], theta, nuNyq))
                 },1:nrow(peaks)))
  res <- peaks_with_height(res, deltanu)
  return(as_tibble(res))
}

false_negative_probability <- function(amplitude, deltanu) {
  a <- amplitude / sqrt(deltanu)
  return(unname(exp(predict(False_negatives_model,
                            newdata = tibble(amplitude = a)))))
}

# 31/01/2020
# Adding functions for final MLE fit where first all peaks are fitted together with
# starting guesses near fit values
peaks_MLE_final_sd <- function(peaks, pds, final_fit_factor=0.1, naverages=1) {
  deltanu <- diff(pds$frequency[1:2])
  # maxLWD isn't used but needs to be given to function so set to NA
  maxLWD <- 10.0
  ## Optimise all the peaks simultaneously with error estimates
  mle.fit <- peaks_MLE2(peaks, pds, maxLWD, hessian = FALSE, final_fit = TRUE, final_fit_factor=final_fit_factor, naverages=naverages)
  ## Get the error estimates one at a time
  print("Getting error estimates")
  peaks <- peaks_with_MLE_sd(peaks = mle.fit$peaks, pds = pds, maxLWD = maxLWD, naverages=naverages)
  ## TODO  if a resolved peak gets a linewidth~deltanu during the optimisation
  ## we should consider using a sinc^2
  peaks <-
    peaks %>%
    peaks_with_height(deltanu) %>%
    peaks_with_AIC(pds, naverages=naverages) %>%
    arrange(frequency)
  return(peaks)
}
