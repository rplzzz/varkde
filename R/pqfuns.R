#' Make CDF and quantile functions for an interpolation grid
#'
#' @param xgrid Grid x coordinates
#' @param cdfgrid Grid CDF values
#' @param lslope Exponential constant for the left tail
#' @param ltail Integrated probability mass for the left tail
#' @param rslope Exponential constant for the right tail
#' @param rtail Integrated probability mass for the right tail
#' @return A list containing the CDF function and the quantile function
#' @noRd
mkpqfun <- function(xgrid, cdfgrid, lslope, ltail, rslope, rtail)
{
  ngrid <- length(xgrid)

  ## Monotonic splines are used for fitting within the grid
  psfun <- splinefun(xgrid, cdfgrid, method='hyman')

  pfun <- function(x, log.p=FALSE) {
    cdfvals <- rep(0, length(x))
    left <- x < xgrid[1]
    right <- x > xgrid[ngrid]
    mid <- !left & !right

    if(any(left)) {
      xl <- x[left]
      if(log.p) {
        cdfvals[left] <- log(ltail) - lslope * (xgrid[1] - xl)
      }
      else {
        cdfvals[left] <- ltail * exp(-lslope * (xgrid[1] - xl))
      }
    }
    if(any(right)) {
      xr <- x[right]
      tf <- -rtail * exp(-rslope * (xr-xgrid[ngrid]))  # tail fraction
      if(log.p) {
        cdfvals[right] <- log1p(tf)
      }
      else {
        cdfvals[right] <- 1 + tf
      }
    }
    if(any(mid)) {
      xm <- x[mid]
      cdfmid <- psfun(xm)
      if(log.p) {
        cdfvals[mid] <- log(cdfmid)
      }
      else {
        cdfvals[mid] <- cdfmid
      }
    }

    cdfvals
  }

  ## To make the quantile function, eliminate duplicated values in the CDF
  keep <- !duplicated(cdfgrid)
  cdfg <- cdfgrid[keep]
  xg <- xgrid[keep]
  ng <- length(xg)
  qsfun <- splinefun(cdfg, xg, method='hyman')

  qfun <- function(p, log.p = FALSE) {
    if(log.p) {
      lp <- p
      p <- exp(p)
    }
    else {
      lp <- log(p)
    }

    left <- p < cdfg[1]
    right <- p > cdfg[ng]
    mid <- !left & !right

    qvals <- rep(0, length(p))

    if(any(left)) {
      lpl <- lp[left]
      qvals[left] <- xg[1] - (log(ltail) - lpl) / lslope
    }
    if(any(right)) {
      lpr <- log1p(-p[right])
      qvals[right] <- xg[ng] + (log(rtail) - lpr) / rslope
    }
    if(any(mid)) {
      qvals[mid] <- qsfun(p[mid])
    }
    qvals
  }

  list(pfun, qfun)
}
