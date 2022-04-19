library(varkde)

datadir <- here::here('analysis')
testdata <- readRDS(file.path(datadir, 'ews-testdata.rds'))

bwmin <- c(
  age = 1.0,
  temp = 0.2,
  pulse = 1.0,
  resp = 2.0,
  spo2 = 1.0,
  sysbp = 1.0,
  diabp = 1.0,
  white_cell = 0.01,
  BUN = 1.0,
  creatinine = 0.1,
  bicarbonate = 1.0,
  sodium = 1.0,
  bilirubin = 0.1,
  hematocrit = 0.1
)

varnames <- names(bwmin)

procvar <- function(varname) {
  message('Running ', varname)
  vtarg <- testdata$target[,varname]
  vntarg <- testdata$nontarget[,varname]

  kdetarg <- kde(vtarg, bwfac=5, bwmin=bwmin[varname])
  kdentarg <- kde(vntarg, bwfac=5, bwmin=bwmin[varname])

  logratio <- function(x) {
    dvarkde(x, kdetarg, log=TRUE) - dvarkde(x, kdentarg, log=TRUE)
  }

  #pltfile <- file.path(datadir, paste0('logratio-', varname, '.png'))
  trng <- quantile(vtarg, c(0.0005, 0.9995), na.rm=TRUE)
  ntrng <- quantile(vntarg, c(0.0005, 0.9995), na.rm=TRUE)
  xmin <- max(0, min(trng[1], ntrng[1]))
  if(varname == 'spo2') {
    xmax <- 100
  }
  else {
    xmax <- max(trng[2], ntrng[2])
  }

  curve(logratio, xmin, xmax, xlab=varname, ylim=c(-5,10))

  list(targ=kdetarg, ntarg=kdentarg)
}

vardens <- lapply(varnames, procvar)
names(vardens) <- varnames
