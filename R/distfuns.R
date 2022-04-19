#' Kernel Estimates of Empirical Distributions
#'
#' Density, CDF, quantile, and RNG functions for empirical distributions estimated
#' with variable-bandwidth kernels.
#'
#' In the RNG functions for built-in distributions \code{n} can be passed as a
#' vector, in which case the length of the vector is the number of deviates to
#' generate.  The main purpose of this functionality seems to be to generate subtle
#' and hard to debug errors; therefore, it is \emph{not} supported here.  Passing
#' a vector for \code{n} will generate a warning, and the first value will be used.
#'
#' @section Caution:
#'
#' For all of these functions we evaluate the kernel density estimate on a grid
#' of x (or p) values and construct an interpolating function.  Because of this
#' it isn't guaranteed that the mathematical relationships between the functions
#' will necessarily hold exactly.  That is, \code{pkde} will not necessarily equal the
#' integral of \code{dkde}, and \code{qkde} may not be an exact inverse of \code{qkde}.
#'
#' @param x Vector of quantiles
#' @param p Vector of probabilities
#' @param n Number of random values to generate
#' @param obj A \code{varkde} object returned by \code{\link{kde}}
#' @param log Flag indicating that log density should be returned
#' @param log.p Flag indicating that DF probabilities are logs
#' @name varkernel
NULL

#' @describeIn varkernel Density function
#' @export
dvarkde <- function(x, obj, log=FALSE)
{
  obj$dfun(x, log)
}

#' @describeIn varkernel Cumulative distribution function
#' @export
pvarkde <- function(x, obj, log.p=FALSE)
{
  obj$pfun(x, log.p)
}

#' @describeIn varkernel Quantile function
#' @export
qvarkde <- function(x, obj, log.p=FALSE)
{
  obj$qfun(x, log.p)
}

#' @describeIn varkernel Random deviate generator
#' @export
rvarkde <- function(n, obj)
{
  if(length(n) > 1) {
    warning('rvarkde: Vector passed for n, where single value is expected. Using n=', n[1])
    n <- n[1]
  }
  u <- runif(n)
  qvarkde(u, obj)
}
