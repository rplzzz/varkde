#' Compute the Epanechnikov (quadratic) kernel
#'
#' @param x Values at which to evaluate the kernel. These should be differences
#' from the center of the kernel.
#' @param lam Kernel scale factor
#' @noRd
epan <- function(x, lam)
{
  bnd<- 1/lam

  ifelse(abs(x) < bnd,
         lam * (1-(x*lam)^2),
         0)
}

#' @describeIn epan Compute the support boundary for the Epanechnikov kernel
#' @noRd
epan_bnd <- function(bw)
{
  sqrt(5)*bw
}

#' @describeIn epan Compute the integral of an Epanechnikov kernel.
#'
#' The \code{epan_int} function is only valid for x/lam in [-1,1].  This condition
#' is not checked in the function.
#' @noRd
epan_int <- function(x, lam)
{
  ul <- x*lam
  0.5 + 0.75 * ul - 0.25 * ul^3
}
