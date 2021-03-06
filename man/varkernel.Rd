% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distfuns.R
\name{varkernel}
\alias{varkernel}
\alias{dvarkde}
\alias{pvarkde}
\alias{qvarkde}
\alias{rvarkde}
\title{Kernel Estimates of Empirical Distributions}
\usage{
dvarkde(x, obj, log = FALSE)

pvarkde(x, obj, log.p = FALSE)

qvarkde(x, obj, log.p = FALSE)

rvarkde(n, obj)
}
\arguments{
\item{x}{Vector of quantiles}

\item{obj}{A \code{varkde} object returned by \code{\link{kde}}}

\item{log}{Flag indicating that log density should be returned}

\item{log.p}{Flag indicating that DF probabilities are logs}

\item{n}{Number of random values to generate}

\item{p}{Vector of probabilities}
}
\description{
Density, CDF, quantile, and RNG functions for empirical distributions estimated
with variable-bandwidth kernels.
}
\details{
In the RNG functions for built-in distributions \code{n} can be passed as a
vector, in which case the length of the vector is the number of deviates to
generate.  The main purpose of this functionality seems to be to generate subtle
and hard to debug errors; therefore, it is \emph{not} supported here.  Passing
a vector for \code{n} will generate a warning, and the first value will be used.
}
\section{Functions}{
\itemize{
\item \code{dvarkde}: Density function

\item \code{pvarkde}: Cumulative distribution function

\item \code{qvarkde}: Quantile function

\item \code{rvarkde}: Random deviate generator
}}

\section{Caution}{


For all of these functions we evaluate the kernel density estimate on a grid
of x (or p) values and construct an interpolating function.  Because of this
it isn't guaranteed that the mathematical relationships between the functions
will necessarily hold exactly.  That is, \code{pkde} will not necessarily equal the
integral of \code{dkde}, and \code{qkde} may not be an exact inverse of \code{qkde}.
}

