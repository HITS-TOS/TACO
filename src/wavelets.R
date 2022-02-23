###
# wavCWT (constructor)
###

library(ifultools)
#library(splus2R)

# From line 259 of Swrappers.R from splus2R
"ifelse1" <- function(test, x, y, ...)
{
  # if(test) return x, else return y.
  # Like ifelse(), except that test is length 1, and x or y
  # is returned as is (whatever its length).
  #     ifelse1(test, x, y)
  # is equivalent to
  #     if(test){x} else {y}.
  # This is particularly useful for assignment;
  #     answer = ifelse1(test, x, y)
  # is equivalent to
  #     if(test) answer = x else answer = y
  #
  # If more than three arguments are supplied, then y should be
  # a second test;
  #     ifelse1(test1, x, test2, u, v)
  # is equivalent to
  #     if(test){x} else if(test2) {y} else {v}
  # This may be iterated; there should be an odd number of arguments.
  if(test) x else if(missing(..1))
    y
  else ifelse1(y, ...)
}

# lowerCase from line 319 of Swrappers.R from splus2R
"lowerCase" <- function(x)
  casefold(x, upper=FALSE)

wavCWT <- function(x, scale.range=deltat(x) * c(1, length(x)), n.scale=100,
                        wavelet="gaussian2", shift=5, variance=1){

  # check input argument types and lengths
  checkVectorType(scale.range,"numeric")
  checkScalarType(n.scale,"integer")
  checkScalarType(wavelet,"character")
  checkScalarType(shift,"numeric")
  checkScalarType(variance,"numeric")
  checkRange(n.scale, c(1,Inf))

  # obtain series name
  series.name <- deparse(substitute(x))

  if (length(scale.range) != 2)
    stop("scale.range must be a two-element numeric vector")
  if (variance <= 0)
    stop("variance input must be positive")

  # obtain sampling interval
  sampling.interval <- deltat(x)

  # create a vector of log2-spaced scales over the specified range of scale
  octave <- logb(scale.range, 2)
  scale  <- ifelse1(n.scale > 1, 2^c(octave[1] + seq(0, n.scale - 2) * diff(octave) /
    (floor(n.scale) - 1), octave[2]), scale.range[1])

  # project the scale vector onto a uniform grid of sampling interval width
  scale <- unique(round(scale / sampling.interval) * sampling.interval)
  n.scale <- length(scale)

  # check scale range
  if (min(scale) + .Machine$double.eps < sampling.interval)
    stop("Minimum scale must be greater than or equal to sampling interval ",
      "of the time series")

  # map the time vector
  if (inherits(x, "signalSeries"))
  {
    times <- as(x@positions,"numeric")
    x <- x@data
  }
  else
  {
  
    times <- time(x)
    x <- as.vector(x)
  }

  # ensure that x is a double vector
  storage.mode(x) <- "double"

  # map the wavelet filter and corresponding argument
  gauss1 <- c("gaussian1", "gauss1")
  gauss2 <- c("gaussian2", "gauss2", "mexican hat", "sombrero")
  supported.wavelets <- c("haar", gauss1, gauss2, "morlet")
  wavelet <- match.arg(lowerCase(wavelet), supported.wavelets)

  # map the filter type to MUTILS type
  # 4: gaussian1
  # 5: gaussian2, sombrero, mexican hat
  # 6: morlet
  # 7: haar
  filter <- mutilsFilterTypeContinuous(wavelet)

  if (filter == 4)
  {
    filter.arg <- sqrt(variance)
    wavelet    <- "gaussian1"
  }
  else if (filter == 5)
  {
    filter.arg <- sqrt(variance)
    wavelet    <- "gaussian2"
  }
  else if (filter == 6)
  {
    filter.arg <- shift
    wavelet    <- "morlet"
  }
  else if (filter == 7)
  {
    filter.arg <- 0.0
    wavelet    <- "haar"

    # in the case of the Haar, the Euler-Macluarin approximation
    # to the DFT of the wavelet filter is not defined for non-integer
    # multiples of the sampling interval. therefore, we coerce the
    # scales accordingly in this case
    scale <- sampling.interval * unique(round(scale / sampling.interval))
  }
  else
  {
    stop("Unsupported filter type")
  }
  
  # calculate the CWT
  z <- itCall("RS_wavelets_transform_continuous_wavelet",
    as.numeric(x),
    as.numeric(sampling.interval),
    as.integer(filter),
    as.numeric(filter.arg),
    as.numeric(scale))
    #
    #COPY=rep(FALSE, 5),
    #CLASSES=c("numeric", "numeric", "integer", "numeric", "numeric"),
    #PACKAGE="ifultools")

  # if the impulse response of the waveleyt filter is purely real
  # then transform CWT coefficients to purely real as well
  if (wavelet != "morlet")
    z <- Re(z)

  # assign attributes
  attr(z, "scale")       <- scale
  attr(z, "time")        <- as.vector(times)
  attr(z, "wavelet")     <- wavelet
  attr(z, "series")      <- x
  attr(z, "sampling.interval") <- sampling.interval
  attr(z, "series.name") <- series.name
  attr(z, "n.sample")    <- length(x)
  attr(z, "n.scale")     <- n.scale
  attr(z, "filter.arg")  <- filter.arg

  oldClass(z) <- "wavCWT"

  z
}

###
# wavCWTTree
###

"wavCWTTree" <- function(x, n.octave.min=1, tolerance=0.0, type="maxima")
{
  # define local functions
  "WTMM" <- function(x, tolerance=NULL, type="maxima"){

    if (!is(x,"wavCWT"))
      stop("Input object must be of class wavCWT")

    # obtain attributes
    x.attr   <- attributes(x)
    times    <- x.attr$time
    scales   <- x.attr$scale
    n.sample <- x.attr$n.sample
    series   <- x.attr$series

    # check tolerances
    if (is.null(tolerance)){

      # use Donoho and Johnstone universal thresholding
      # tol=sqrt(2 * sigma^2 * log(N)) where sigma^2
      # is the variance of the noise, and N is the length
      # of the original time series. Since we do not
      # know sigma^2 a priori, we use the median absolute deviation
      # the level 1 DWT coefficients using the Haar filter as an
      # approxiation. Finally, we divide by the sqrt(scale) to make
      # the threshold scale dependent

      # noise.variance <- (median(abs((diff(series)[seq(2,n.sample-1,by=2)] / sqrt(2))))/ 0.6745)^2
      # vprint(noise.variance)
      # tolerance <- sqrt(2 * noise.variance * log(n.sample)) / sqrt(scales)
      tolerance <- mad(Mod(x[,1])) / scales
    }

    if (length(tolerance) < length(scales))
      tolerance <- tolerance[1] / sqrt(scales)

    wtmmz <- itCall("RS_wavelets_transform_continuous_wavelet_modulus_maxima",
      as.matrix(x)+0i, tolerance, mutilsTransformPeakType(type))
    #
      #CLASSES=c("matrix","numeric","integer"),
      #COPY=rep(FALSE,3),
      #PACKAGE="ifultools")

    z <- matrix(0, nrow=nrow(x), ncol=ncol(x))
    z[matrix(unlist(wtmmz),ncol=2)+1] <- 1

    z
  }

  "wtmmBranches" <- function(wtmm, extrema.mask, times, scales, span.min=5, gap.max=3, skip=NULL, sampling.interval=1)
  {
    # define normalized scales
    scales <- as.integer(scales / sampling.interval)

    n.scale <- ncol(extrema.mask)
    n.sample <- nrow(extrema.mask)

    if (is.null(scales))
      scales <- 1:n.scale

    iwtmm <- which(extrema.mask[, n.scale] > 0)

    # scan from large to small scale
    iscale <- seq(n.scale-1,1,-1)

    tree <- as.list(iwtmm)
    names(tree) <- iwtmm
    peakStatus <- as.list(rep(0, length(iwtmm)))
    names(peakStatus) <- iwtmm
    orphanRidgeList <- NULL
    orphanRidgeName <- NULL

    n.level <- length(iscale)

    for (j in seq(n.level)){

      iscale.j <- iscale[j]
      scale.j <- scales[iscale.j]

      if (length(iwtmm) == 0){
          iwtmm <- which(extrema.mask[, iscale.j] > 0)
          next
      }

      span <- scale.j * 2 + 1

      if (span < span.min)
        span <- span.min

      remove.j <- selPeak.j <- NULL

      # loop through each point in the current scale's WTMM
      for (k in seq(along=iwtmm)){

        # define search range in the time index
        itime <- iwtmm[k]
        itime.start <- itime - span
        if (itime.start < 1)
          itime.start <- 1
        itime.end <- itime + span
        if (itime.end > n.sample)
          itime.end <- n.sample
        itime.candidates <- which(extrema.mask[itime.start:itime.end, iscale.j] > 0) + itime.start - 1

        if (length(itime.candidates) == 0){

            status.k <- peakStatus[[as.character(itime)]]

            if (status.k > gap.max & scale.j >= 2){
              temp            <- tree[[as.character(itime)]]
              orphanRidgeList <- c(orphanRidgeList, list(temp[1:(length(temp) - status.k)]))
              orphanRidgeName <- c(orphanRidgeName, paste(iscale.j + status.k + 1, itime, sep="_"))
              remove.j        <- c(remove.j, as.character(itime))
              next
            }
            else {
              itime.candidates <- itime
              peakStatus[[as.character(itime)]] <- status.k + 1
            }
        }
        else {
          peakStatus[[as.character(itime)]] <- 0
          if (length(itime.candidates) >= 2)
            itime.candidates <- itime.candidates[which.min(abs(itime.candidates - itime))]

        }

        tree[[as.character(itime)]] <- c(tree[[as.character(itime)]], itime.candidates)
        selPeak.j <- c(selPeak.j, itime.candidates)
      }

      if (length(remove.j) > 0){
        bad.tree   <- which(is.element(names(tree), remove.j))
        tree       <- tree[-bad.tree]
        peakStatus <- peakStatus[-bad.tree]
      }

      dupPeak.j <- unique(selPeak.j[duplicated(selPeak.j)])

      if (length(dupPeak.j) > 0){

          bad.tree <- NULL

          for (dupPeak.jk in dupPeak.j){
            selInd          <- which(selPeak.j == dupPeak.jk)
            selLen          <- sapply(tree[selInd], length)
            bad.tree.jk     <- which.max(selLen)
            bad.tree        <- c(bad.tree, selInd[-bad.tree.jk])
            orphanRidgeList <- c(orphanRidgeList, tree[bad.tree.jk])
            orphanRidgeName <- c(orphanRidgeName, paste(iscale.j, selPeak.j[bad.tree.jk], sep="_"))
          }

          selPeak.j  <- selPeak.j[-bad.tree]
          tree       <- tree[-bad.tree]
          peakStatus <- peakStatus[-bad.tree]
      }
      names(tree) <- selPeak.j
      names(peakStatus) <- selPeak.j

      if (scale.j >= 2){
        maxInd.next      <- which(extrema.mask[, iscale.j] > 0)
        unSelPeak.j      <- maxInd.next[!is.element(maxInd.next, selPeak.j)]
        newPeak.j        <- as.list(unSelPeak.j)
        names(newPeak.j) <- unSelPeak.j
        tree             <- c(tree, newPeak.j)
        iwtmm            <- c(selPeak.j, unSelPeak.j)
        newPeakStatus    <- as.list(rep(0, length(newPeak.j)))
        names(newPeakStatus) <- newPeak.j
        peakStatus       <- c(peakStatus, newPeakStatus)
      }
      else {
        iwtmm <- selPeak.j
      }
    }

    names(tree) <- paste(1, names(tree), sep="_")
    names(orphanRidgeList) <- orphanRidgeName
    tree <- c(tree, orphanRidgeList)
    tree <- lapply(tree, rev)
    tree <- tree[unique(names(tree))]

    tree <- lapply(seq(along=tree), function(i, tree, iscale.min, times, scales, wtmm){
      itime <- tree[[i]]
      iscale <- seq(iscale.min[i], length=length(itime))
      list(itime=itime, iscale=iscale, time=times[itime], scale=scales[iscale], extrema=wtmm[cbind(itime,iscale)])
    },
    tree=tree,
    iscale.min=as.integer(gsub("_.*","",names(tree))),
    times=times,
    scales=scales*sampling.interval,
    wtmm=wtmm)

    # remove any redundant branches
    iflat <- lapply(tree, function(x, nr) (x$iscale-1)*nr + x$itime, nr=nrow(wtmm))

    flatset <- iflat[[1]]
    bad <- NULL

    for (i in seq(2,length(iflat))){

       if (any(is.element(iflat[[i]], flatset)))
         bad <- c(bad, i)
       else
         flatset <- c(flatset, iflat[[i]])
    }

    if (length(bad) > 0)
      tree <- tree[-bad]

    tree
  }

  # obtain attributes
  x.attr   <- attributes(x)
  times    <- x.attr$time
  scales   <- x.attr$scale
  n.sample <- x.attr$n.sample
  sampling.interval <- x.attr$sampling.interval
  border.times <- range(times) + sampling.interval * c(1,-1)

  # locate the the extrema in the CWT matrix
  extrema.mask <- WTMM(x, tolerance=tolerance, type=type)

  if (!identical(dim(x),dim(extrema.mask)))
    stop("Input WTMM dimenions do not match those of the input CWT matrix")

  # develop WTMM tree
  z <- wtmmBranches(ifelse1(is.complex(x), Mod(as.matrix(x)), as.matrix(x)), extrema.mask, times, scales, sampling.interval=sampling.interval)

  # define minimum number of points needed per branch
  n.scale  <- length(scales)
  n.octave <- log2(max(scales) / min(scales))
  n.voice  <- (n.scale - 1) / n.octave
  n.scale.min <- as.integer(n.voice * n.octave.min)

  good <- which(unlist(lapply(z,function(x, n.scale.min) length(x[[1]]) > n.scale.min, n.scale.min=n.scale.min)))
  z <- z[good]
  endtime <- unlist(lapply(z,function(x,iscale) x$itime[iscale], iscale=which.min(scales)))
  isort <- order(endtime)
  z <- z[isort]

  names(z) <- seq(z)

  attr(z, "iendtime")     <- endtime[isort]
  attr(z, "endtime")      <- times[endtime[isort]]
  attr(z, "time")         <- times
  attr(z, "scale")        <- scales
  attr(z, "extrema.mask") <- extrema.mask
  attr(z, "noise")        <- x[,1]
  attr(z, "branch.hist")  <- colSums(extrema.mask*abs(x))
  attr(z, "wavelet")      <- attr(x,"wavelet")
  attr(z, "filter.arg")   <- attr(x,"filter.arg")
  attr(z, "series.name")  <- attr(x,"series.name")
  attr(z, "series")       <- attr(x,"series")
  attr(z, "sampling.interval") <- attr(x,"sampling.interval")

  oldClass(z) <- "wavCWTTree"

  z
}

###
# [.wavCWTTree
###

"[.wavCWTTree" <- function(x, i, ..., time=NULL, range=NULL)
{
  ax    <- attributes(x)
  times <- ax$endtime

  if (!missing(range)){

    if (length(range) == 2){

      i <- which(times >= range[1] & times <= range[2])
    }
    else{
      i <- seq(length(x))
    }
  }
  else if(!missing(time)){

    i <- sort(time)

    min.scale <- median(diff(times))

    itime <- NULL

    for (j in i){

      itime <- c(itime, which(times >= j - min.scale & times <= j + min.scale))
    }

    i <- itime
  }

  i <- i[i <= length(x) & i >= 1]

  z <- oldUnclass(x)[i]

  attributes(z) <- c(attributes(z),
    ax[ setdiff(names(ax), c("names", "dim", "dimnames")) ])

  z
}

###
# plot.wavCWTTree
###

"plot.wavCWTTree" <- function(x, add=FALSE, pch="o",  label=TRUE, log.="y",
  xlab=NULL, ylab=NULL, extrema=FALSE, zoom=NULL,
  fit=FALSE, models= c("lm", "lmsreg", "ltsreg"),
  cex=0.8, col.skip=max(1,round(256/length(x))),  ...)
{
  branches <- names(x)

  if (!length(x))
    stop("No branches found in input object")

  if (fit){

    n.branch <- min(length(x), 4)

    nrow <- 1
    ncol <- 1

    if (n.branch > 1)
      nrow <- 2
    if (n.branch > 2)
      ncol <- 2

    frame()
    gap <- 0.11
    old.plt <- splitplot(nrow,ncol,1,gap=gap)
    on.exit(par(old.plt))

    lwd <- c(2, 1, 2, rep(1,10))
    lty <- c(1, 1, 4, seq(5,length=10))

    # initialize arguments
    slope <- vector(mode="numeric", length=length(models))

    for (i in seq(n.branch)){

      if (i > 1)
        splitplot(nrow,ncol,i,gap=gap)

      x.branch <- log(x[[i]]$scale)
      y.branch <- log(x[[i]]$extrema)

      if (x.branch[length(x.branch)] < x.branch[1]){
        x.branch <- rev(x.branch)
        y.branch <- rev(y.branch)
      }

      breaks <- linearSegmentation(x.branch, y.branch)

      if (is.null(breaks))
        breaks <- length(x.branch)

      plot(x.branch, y.branch, xlab="log(scale)", ylab="log(|EXTREMA|)", type="b",
        pch="o", cex=cex)

      abline(v=x.branch[breaks], lty=2, xpd=FALSE)

      cut      <- seq(breaks[1]-1)
      x.branch <- x.branch[cut]
      y.branch <- y.branch[cut]

      # the lmsreg() and ltsreg() models will return an
      # error if less than 6 points are fit, so restrict
      # the models used accodingly
      if (length(x.branch) > 5)
        admissible.models <- models
      else
        admissible.models <- "lm"

      datalist <- list(x.branch=x.branch, y.branch=y.branch)

      for (imodel in seq(admissible.models)){

        eval(parse(text=paste("fit <- coef(", admissible.models[imodel],
          "(y.branch ~ x.branch,data=datalist))", sep="")))
        abline(fit, lty=lty[imodel], lwd=lwd[imodel], xpd=FALSE)
        slope[imodel] <- fit[2]
      }

      imodel <- seq(admissible.models)

      key.cex <- 0.7

      mdlName <- upperCase(admissible.models)
	  legend("bottomright",
             paste(format(mdlName), "\t: slope=", format(round(slope[imodel], 4))),
    	     lty=lty[imodel],
    	     lwd=lwd[imodel],
    	     cex=key.cex)

      mtext(paste("Branch", branches[i]), adj=1, cex=0.85, line=0.5)
    }

    return(invisible(NULL))
  }

  if (!add)
    frame()

  x.attr <- attributes(x)

  if (extrema){

    pch <- 20
    col <- "blue"

    mask <- which(x.attr$extrema.mask == 1, arr.ind=TRUE)

    data <- scaleZoom(x.attr$time[mask[,1]], x.attr$scale[mask[,2]], zoom=zoom, logxy=log., xy.linked=TRUE,
      xlab="Time", ylab="Scale")
    if (!add)
      plot(data$x, data$y, col=col, pch=pch, xlab=data$xlab, ylab=data$ylab, cex=cex, ...)
    else
      points(data$x, data$y, col=col, pch=pch, cex=cex, ...)

    return(invisible(NULL))
  }

  if (!add){

    times  <- range(x.attr$time)
    scales <- range(x.attr$scale)

    if(is.element(log., "y")){
      scales <- logb(scales, base=2)
      if (is.null(ylab))
        ylab <- "log2(scale)"
    }
    else if (is.null(ylab))
      ylab <- "Scale"

    if(is.element(log., "x")){
      times <- logb(times, base=2)
      if (is.null(xlab))
        xlab <- "log2(time)"
    }
    else if (is.null(xlab))
      xlab <- "Time"

    plot(times, scales, type="n", xlab=xlab, ylab=ylab)
  }

  col <- c("black","red","blue","green","deeppink","orange","cyan","magenta","violet","navy","purple",
     "yellowgreen","orangered","goldenrod","cornflowerblue","plum","steelblue","tomato","pink")

  for (i in seq(along=x)){

    times  <- x[[i]]$time
    scales <- x[[i]]$scale

    if(is.element(log., "y"))
      scales <- logb(scales, base=2)
    if(is.element(log., "x"))
      times <- logb(times, base=2)

    icol <- ((i-1) %% length(col)) + 1

    if (label){

      imaxscale <- order(scales)[length(scales)]
      lines(times[-imaxscale], scales[-imaxscale], col=col[icol], pch=pch, cex=0.5, type="b", ...)
      text(times[imaxscale], scales[imaxscale], as.character(branches[i]), cex=1, col=col[icol])
    }
    else
      lines(times, scales, col=col[icol], lty=1, pch=pch, ...)
  }

  invisible(NULL)
}


###
# print.summary.wavCWTTree
###

"print.summary.wavCWTTree" <- function(x, digits=max(2, .Options$digits - 4), ...)
{
  oldOptions <- options(digits = digits)
  on.exit(options(oldOptions))
  NextMethod("print")
  invisible(x)

#  print(round(data.frame(oldUnclass(x)), digits), ...)
#  invisible(x)
}

###
# print.wavCWTTree
###

"print.wavCWTTree" <- function(x, justify="left", sep=":", ...)
{
  # obtain attributes
  xatt       <- attributes(x)
  times      <- xatt$time
  name       <- xatt$series.name
  scale      <- xatt$scale
  n.scale    <- length(scale)
  series     <- xatt$series
  sampling.interval <- xatt$sampling.interval
  n.branch   <- length(x)
  filter.arg <- xatt$filter.arg

  # pretty print strings
  waveletstr <- "Mexican Hat (Gaussian, second derivative)"
  filtcat1 <- "Wavelet variance"
  filtval1 <- filter.arg^2

  if (is(series, "signalSeries")){
    units.time <- series@units.position
    if (length(units.time) > 0)
      filtval1 <- paste(filtval1, " (", units.time, ")", sep="")
  }

  scale.range <- range(scale)

  main <- paste("Continuous Wavelet Transform Tree of", name)

  z <- list(
    "Wavelet"=waveletstr,
    "Wavelet variance"=filtval1,
    "Length of series"=length(series),
    "Sampling interval"=sampling.interval,
    "Number of scales"=n.scale,
    "Range of scales"=paste(scale.range[1], "to", scale.range[2]),
    "Number of branches"=n.branch
  )

  prettyPrintList(z, header=main, sep=sep, justify=justify, ...)

  invisible(x)
}

###
# summary.wavCWTTree
###

"summary.wavCWTTree" <- function(object, ...)
{
  tmpargs <- lapply(object, function(x){

    scale <- x$scale
    ext <- x$extrema

    c(x$time[which.min(scale)], length(x$time), log2(max(scale) / min(scale)), range(ext), mean(ext), sd(ext), var(ext), mad(ext))})
    
  z <- data.frame(do.call("rbind", tmpargs))

  names(z) <- c("End Time", "Length", "Octaves", "Min", "Max", "Mean", "SD", "Var", "MAD")

  oldClass(z) <- c("summary.wavCWTTree","data.frame")

  z
}