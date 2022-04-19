#### Plot the log-density ratio for a pair of event and non-event kernel density
#### estimates

plot_denratio <- function(denstructs, varname, qrange = c(5e-4, 1-5e-4))
{
  kdetarg <- denstructs$targ
  kdentarg <- denstructs$ntarg

  trng <- qvarkde(qrange, kdetarg)
  ntrng <- qvarkde(qrange, kdentarg)
  xmin <- max(0, min(trng[1], ntrng[1]))
  if(varname == 'spo2') {
    xmax <- 100
  }
  else {
    xmax <- max(trng[2], ntrng[2])
  }

  logratio <- function(x) {dvarkde(x, kdetarg, log=TRUE) - dvarkde(x, kdentarg, log=TRUE)}

  curve(logratio, xmin, xmax, xlab=varname, ylim=c(-5,10))
}
