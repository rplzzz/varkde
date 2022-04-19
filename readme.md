# varkde: Variable Bandwidth Kernel Density Estimate

## Development notes

This package came about as a result of some work I was doing on
Bayesian classifiers.  I was using kernel density estimates to produce
marginal density functions for observations. The process worked
reasonably well for values where the density was relatively high, but
tended to gyrate wildly for values where data from one or both of the
classes was sparse.  The reason why is that the kernel bandwidth
estimate based on the entire dataset tended to produce bandwidths that
were too small in the sparse regions.

The solution here was conceptually inspired by the `muhaz` package,
which uses a variable bandwidth kernel for precisely this reason.
However, `muhaz` computes hazard functions, and I couldn't find a
comparable package for density functions.  Backing out a density
function by multiplying the kernel hazard function by an empirical
survival function seemed unlikely to work well, because of the
discontinuities in the empirical survival function, so I cobbled
together my own.

At this point I've made quite a few tweaks related to using the KDEs
in Bayesian classifiers, to the point that I'm not sure it's actually
useful for general purpose use.  I'm keeping it separate for now
because it makes development on the classifier a little cleaner.
However, it could get merged into the classifier code base at any
time.

## Limitations

The most important limitation, by far, is that `pvarkde` and `qvarkde`
(the cdf and quantile functions) are not exact inverses of one
another.  This is a result of calculating exact integrals over the
kernel functions on a grid and using spline interpolation to get
values in between grid points.  To solve this problem we need to
either downgrade to linear interpolation, or we need to find an
invertible spline.  The latter is not so easy to do, considering that
these functions also need to be monotonic.  I guess a third option
would be to implement one of the functions as an interpolation and the
other as a bisection solver, so as to guarantee that we get an exact
inverse.  

So far, the crude aproximations for the cdf and quantile functions
have been good enough for our purposes.  We do most of our work in
terms of the cdf, so it's important to get that right.  Mostly we just
use the quantile function to generate more analyst-friendly axes for
plots, so high accuracy isn't a big priority.  

The second major limitation is the behavior in the tails.  Past the
last kernel on each end of the distribution, we extrapolate by holding
the logarithmic derivative constant at its value at the end of the
spline.  This has the effect of creating exponential tails on both
sides, which is a reasonable guess, but both the initial level and the
logarithmic slope are a little arbitrary, depending on exactly where
we cut off the last kernel.  In fact, the way we place the end grid
points (at the lowest/highest extent of any kernel), we are pretty
much guaranteeing that the last grid point will be at the background
density, which is itself a little arbitrary.  Among other things it
depends on the number of sample points.  So, if we take log density
ratios in those tails, we are taking the ratio of two exponentials
with rather arbitrary parameters, and we could get all sorts of
strange behavior.

One thing that mitigates this limitation is that for the classifier we
only use the density functions between the 1e-3 and 1-1e-3 quantiles.
Since we normally have tens of thousands of points even in the event
group, that puts us well inside sampling grid.  One edge case is that
the edge of one grid (generally the event group) could be well outside
the edge of the other.  In those cases we would be comparing in-grid
values in one distribution to extrapolated values in the other.
Though this will correctly indicate that the probability density of
the group that is still within its grid is much higher than that of
the group that is out of its grid, the absolute value of the
log-density ratio could be arbitrarily large.  This is potentially
problematic when we sum the log-density ratio for that one variable
with the other terms in the EWS.  

The last limitation (that I can think of so far) is hard boundaries.
We follow the practice of most KDE schemes and simply ignore them.
Because of this, the tails of the distributions extend past the
boundaries.  In EWS usage this probably isn't a big deal.  We will
never have inputs that are on the wrong side of the boundary, so the
overlap has no practical effect, but we probably understate the
density a bit on the right side.  As long as the densities for the two
categories are equally understated, there is no effect, but of course
there is no guarantee that they will be.
