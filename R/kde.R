#' Fit a variable bandwidth kernel density estimate to a vector of data
#'
#' The bandwidth of each kernel is calculated using neighboring points.  We
#' then sample the density on a fixed grid and construct interpolating functions
#' for the density and CDF.
#'
#' The \code{bwmin} parameter allows a minimum allowable bandwidth to be set.  It's
#' a good idea to set this to something like the smallest meaningful difference in
#' measurement values.  When in doubt, the smallest reporting unit of the measurement
#' apparatus is a good choice.  The default value is 1e-6, which is meant to be
#' effectively zero, while protecting against a divide-by-zero in the event that
#' the bandwidth calculation produces zero as a result.
#'
#' TODO: Right now the quantile function is not an exact inverse of the CDF because
#' there is no way to guarantee that the two spline functions interpolate the same
#' way.  Need to fix.
#'
#' @param x Vector of samples for which to construct density.
#' @param npt Size of the neighborhood to use for estimating bandwidth. N
#' data points on each side of the a point will be used in the estimate. Default
#' is 1/20 the total size of the dataset.
#' @param bwfac Bandwidth adjustment factor.  Calculated bandwidths will be
#' uniformly scaled by this factor.  Higher values produce smoother estimates.
#' @param bwmin Minimum allowable bandwidth for an individual kernel.
#' @param ngrid Number of grid points to use for the interpolation function.
#' Default is 256.
#' @export
kde <- function(x, n=NA, bwfac=1, bwmin=1e-6, ngrid=256)
{
  if(is.na(n)) {
    n <- max(min(150, round(length(x)/20)), 4)
  }
  x <- x[!is.na(x)]
  x <- sort(x)

  ## Estimate the kernel width for the ith kernel
  estwidth <- function(i) {
    ## basic calculation (see calcbw)
    ilo <- pmax(1, i-n)
    ihi <- pmin(length(x), i+n)
    xp <- x[ilo:ihi]
    sig <- bwfac * calcbw(xp)

    ## Low-density backstop: require that the bandwidth at each sample be at
    ## least as large as the distance to both neighbors.
    iplus <- pmin(length(x), i+1)
    iminus <- pmax(1, i-1)
    sig <- max(c(sig, x[iplus]-x[i], x[i]-x[iminus]))

    ## backstop for high-density regions: flat minimum on bandwidth.  This will
    ## generally be set to the precision of the measurement. E.g., if a variable
    ## is measured to the nearest integer, bwmin will be 1.
    if(sig < bwmin) {
      sig <- bwmin
    }
    epan_bnd(sig)
  }

  laminv <- sapply(seq_along(x), estwidth)
  lam <- 1/laminv

  ## Set up the grid for estimating the density
  xlo <- min(x - laminv)
  xhi <- max(x + laminv)

  xgrid <- seq(xlo, xhi, length.out = ngrid)
  dengrid <- rep(0, ngrid)

  ## Loop over kernels; accumulate density at each grid point
  for(i in seq_along(lam)) {
    u <- xgrid - x[i]
    dengrid <- dengrid + epan(u, lam[i])
  }
  nsamp <- length(x)
  dengrid <- dengrid / nsamp

  ## We have to deal with any zeros in the density, so we
  ## set a floor on the density.  The floor is roughly equal to the probability
  ## that a number of draws equal to the sample we have would have a X% chance to
  ## fail to land any draws in the zero region.  Right now X is set to 0.75.  We
  ## may make it adjustable in the future, if it looks reasonable to do so.
  xtol <- 1e-16
  pzero <- 0.75
  nzero <- sum(dengrid < xtol)
  denmin <- (1-pzero^(1.0/nsamp)) / (xhi-xlo) * nzero / (ngrid-1)
  message('denmin = ', denmin)
  dengrid[dengrid < denmin] <- denmin

  ## Do our interpolation in log density so that we can ensure that the result
  ## is always positive.
  ldengrid <- log(dengrid)

  logdenfun <- splinefun(xgrid, ldengrid, method='natural')

  ## Calculate the probability mass in the tails.  We need this to make a correction
  ## to the the density function
  lslope <- logdenfun(xgrid[1]) - logdenfun(xgrid[1]-1)
  ltail <- dengrid[1] / lslope
  rslope <- logdenfun(xgrid[ngrid]) - logdenfun(xgrid[ngrid]+1)
  rtail <- dengrid[ngrid] / rslope
  tottail <- ltail + rtail
  message('Left tail: ', ltail, '  Right tail: ', rtail, '  Total: ', tottail)

  ## density function adjusted for tails
  tailfac <- 1.0/(1+tottail)
  logtailfac <- log(tailfac)
  denfun <- function(x, log=FALSE) {
    if(log)  {
      logtailfac + logdenfun(x)
    }
    else {
      tailfac * exp(logdenfun(x))
    }
  }

  ## We also need a function to compute the CDF.  We can't just integrate the
  ## density because the sampled version isn't well behaved, and we don't want
  ## to store the entire list of kernels.
  cdfgrid <- rep(0, ngrid)
  for(i in seq_along(cdfgrid)) {
    u <- xgrid[i] - x
    nplus <- sum(u >= laminv)   # number of kernels for which xgrid is above the kernel's support
    midkerns <- which(u < laminv & u > -laminv)   # kernels for which xgrid falls within the kernel's support
    if(length(midkerns) > 0) {
      nmid <- sum(epan_int(u[midkerns], lam[midkerns]))
    }
    else {
      nmid <- 0
    }
    cdfgrid[i] <- ltail + (nplus + nmid)/nsamp
  }
  ## As with the density, we need a correction here for the tails.  This should be
  ## equal to the tailfac calculated above, but check it anyhow
  tailfac2 <- 1.0 / (cdfgrid[ngrid] + rtail)
  stopifnot(isTRUE(all.equal(tailfac, tailfac2)))
  cdfgrid <- cdfgrid * tailfac2

  pqfuns <- mkpqfun(xgrid, cdfgrid, lslope, ltail, rslope, rtail)

  structure(
    list(dfun = denfun, pfun = pqfuns[[1]], qfun = pqfuns[[2]],
         xgmin = xgrid[1], xgmax = xgrid[ngrid]),
    class = c('varkde', 'list')
  )
}

#' Calculate the KDE bandwidth for a set of observations
#'
#' This calculation is loosely based on the rule of thumb implemented in
#' \code{\link[stats]{bw.nrd}}.
#'
#' This version of the bandwidth takes advantage of the fact that
#' the input vector is already sorted, so the quartile points can be read out
#' of the array directly.  We also dispense with some of the subtleties of
#' calculating empirical quantiles, since this is not a precise calculation.
#'
#' Finally, the OG version of the calculation uses the minimum of the standard
#' deviation of the vector and the IQR estimate of the standard deviation.  The
#' two values are usually quite close, and we are likely to apply a fixed minimum
#' to the bandwidth anyhow, so we save a lot of time by omitting the sd calculation.
#'
#' @param x A sorted numeric vector
#' @noRd
calcbw <- function(x)
{
  n <- length(x)
  iq <- round(c(0.25 * n, 0.75 * n))
  q <- x[iq]
  ## This is the nrd formula.  The 0.7412898 is approximately 1/1.3490, the latter
  ## being the expectation value of the IQR for a standard normal distribution.
  h <- (q[2] - q[1]) * 0.7412898
  1.06 * h * n^(-0.2)
}
