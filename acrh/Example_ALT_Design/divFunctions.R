lnIsTrue_observed <- function(y, m1, m2, s1, s2) {
  #
  lpdf1 <- dlnorm(y, m1, s1, log = TRUE)
  #lpdf2 <- dweibull(y, 1/s2, exp(m2), log = TRUE)
  lpdf2 <- log(dweibull(y, 1/s2, exp(m2)) + 1e-12)
  pdf1 <- exp(lpdf1)
  #
  val <- pdf1*(lpdf1 - lpdf2)
  return(val)
}

lnIsTrue_censored <- function(y, m1, m2, s1, s2) {
  #
  lcdf1 <- log(1 - plnorm(y, m1, s1) + 1e-12)
  lcdf2 <- log(1 - pweibull(y, 1/s2, exp(m2)) + 1e-12)
  cdf1 <- exp(lcdf1)
  #
  val <- cdf1*(lcdf1 - lcdf2)
  return(val)
}


wbIsTrue_observed <- function(y, m1, m2, s1, s2) {
  #
  lpdf1 <- dweibull(y, 1/s1, exp(m1), log = TRUE)
  lpdf2 <- dlnorm(y, m2, s2, log = TRUE)
  pdf1 <- exp(lpdf1)
  #
  val <- pdf1*(lpdf1 - lpdf2)
  return(val)
}


wbIsTrue_censored <- function(y, m1, m2, s1, s2) {
  #
  lcdf1 <- log(1 - pweibull(y, 1/s1, exp(m1)) + 1e-12)
  lcdf2 <- log(1 - plnorm(y, m2, s2) + 1e-12)
  cdf1 <- exp(lcdf1)
  #
  val <- cdf1*(lcdf1 - lcdf2)
  return(val)
}

kldiv_censored_lnIsTrue_tc_5000 <- function(xt, xr, st, sr) {
  tc <- 5000
  intVec <- rep(0, length(xt))
  #print(cbind(m1 = xt, m2 = xr, s1 = st, s2 = sr))
  for (i in 1:length(xt)) {
    #print(c(m1 = xt[i], m2 = xr[i], s1 = st[i], s2 = sr[i]))
    intg_part <- integrate(lnIsTrue_observed, 0, tc,
                           m1 = xt[i], m2 = xr[i], s1 = st[i], s2 = sr[i],
                           subdivisions = 100, stop.on.error = FALSE)$value
    cens_part <- lnIsTrue_censored(tc, m1 = xt[i], m2 = xr[i], s1 = st[i], s2 = sr[i])
    intVec[i] <- intg_part + cens_part
  }
  return(intVec)
}

kldiv_censored_wbIsTrue_tc_5000  <- function(xt, xr, st, sr) {
  tc <- 5000
  intVec <- rep(0, length(xt))
  #print(cbind(m1 = xt, m2 = xr, s1 = st, s2 = sr))
  for (i in 1:length(xt)) {
    #print(c(m1 = xt[i], m2 = xr[i], s1 = st[i], s2 = sr[i]))
    intg_part <- integrate(wbIsTrue_observed, 0, tc,
                           m1 = xt[i], m2 = xr[i], s1 = st[i], s2 = sr[i],
                           subdivisions = 100, stop.on.error = FALSE)$value
    cens_part <- wbIsTrue_censored(tc, m1 = xt[i], m2 = xr[i], s1 = st[i], s2 = sr[i])
    intVec[i] <- intg_part + cens_part
  }
  return(intVec)
}

kldiv_censored_lnIsTrue_tc_1000 <- function(xt, xr, st, sr) {
  tc <- 1000
  intVec <- rep(0, length(xt))
  for (i in 1:length(xt)) {
    #print(c(m1 = xt[i], m2 = xr[i], s1 = st[i], s2 = sr[i]))
    intg_part <- integrate(lnIsTrue_observed, 0, tc,
                           m1 = xt[i], m2 = xr[i], s1 = st[i], s2 = sr[i],
                           subdivisions = 100, stop.on.error = FALSE)$value
    cens_part <- lnIsTrue_censored(tc, m1 = xt[i], m2 = xr[i], s1 = st[i], s2 = sr[i])
    intVec[i] <- intg_part + cens_part
  }
  return(intVec)
}

kldiv_censored_wbIsTrue_tc_1000  <- function(xt, xr, st, sr) {
  tc <- 1000
  intVec <- rep(0, length(xt))
  for (i in 1:length(xt)) {
    #print(c(m1 = xt[i], m2 = xr[i], s1 = st[i], s2 = sr[i]))
    intg_part <- integrate(wbIsTrue_observed, 0, tc,
                           m1 = xt[i], m2 = xr[i], s1 = st[i], s2 = sr[i],
                           subdivisions = 100, stop.on.error = FALSE)$value
    cens_part <- wbIsTrue_censored(tc, m1 = xt[i], m2 = xr[i], s1 = st[i], s2 = sr[i])
    intVec[i] <- intg_part + cens_part
  }
  return(intVec)
}




