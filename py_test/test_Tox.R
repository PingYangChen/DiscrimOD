library(DiscrimOD)
library(Rcpp); library(RcppArmadillo)

# R Function
m5_par <- c(4.282, 835.571, 0.739, 3.515)
MODEL_INFO <- list(
  # m5
  list(model = function(x, p) p[1]*(p[3] - (p[3] - 1)*exp(-(x/p[2])^p[4])),
       para = m5_par),
  # m4
  list(model = function(x, p) p[1]*(p[3] - (p[3] - 1)*exp(-(x/p[2]))),
       paraLower = c(0, 0, 0),
       paraUpper = c(20, 5000, 1)),
  # m3
  list(model = function(x, p) p[1]*exp(-(x/p[2])^p[3]),
       paraLower = c(0, 0, 1),
       paraUpper = c(20, 5000, 15)),
  # m2
  list(model = function(x, p) p[1]*exp(-(x/p[2])),
       paraLower = c(0, 0),
       paraUpper = c(20, 5000)),
  # m1
  list(model = function(x, p) rep(p[1], length(x)),
       paraLower = c(0),
       paraUpper = c(20))
)
DISTANCE <- function(xt, xr) (xt - xr)^2

#
ALG_INFO <- getAlgInfo(nSwarm = c(16, 128),
                       maxIter = c(100, 200),
                       LBFGS_RETRY = rep(1, 2))

dsLower <- 0
dsUpper <- 1250

TWO_M_INFO <- list(
  list(TWO_M = list(MODEL_INFO[[1]], MODEL_INFO[[2]]), nSupp = 3),
  list(TWO_M = list(MODEL_INFO[[1]], MODEL_INFO[[3]]), nSupp = 4),
  list(TWO_M = list(MODEL_INFO[[1]], MODEL_INFO[[4]]), nSupp = 3),
  list(TWO_M = list(MODEL_INFO[[1]], MODEL_INFO[[5]]), nSupp = 2)
)


outALL <- NULL
outAll <- lapply(1:length(TWO_M_INFO), function(CaseID) {
  #CaseID <- 1
  out <- DiscrimOD(TWO_M_INFO[[CaseID]]$TWO_M, DISTANCE,
                   TWO_M_INFO[[CaseID]]$nSupp, dsLower, dsUpper,
                   crit_type = "pair_fixed_true",
                   MaxMinStdVals = NULL, IF_INNER_LBFGS = TRUE,
                   ALG_INFO = ALG_INFO, seed = NULL, verbose = TRUE)

  round(out$BESTDESIGN, 3)
  out$BESTVAL

  eqv <- equivalence(PSO_RESULT = out, MODEL_INFO = TWO_M_INFO[[CaseID]]$TWO_M,
                     DISTANCE = DISTANCE, ngrid = 100,
                     dsLower = dsLower, dsUpper = dsUpper, MaxMinStdVals = NULL, ALG_INFO = ALG_INFO)
  plot(eqv$Grid_1, eqv$DirDeriv, type = "l", col = "blue"); abline(h = 0);
  points(out$BESTDESIGN[,1], rep(0, nrow(out$BESTDESIGN)), pch = 16)
  DESIGN1 <- out$BESTDESIGN
  #DESIGN1 <- cbind(c(0, 468.186, 1064.179), c(.249, .498, .253))
  designCriterion(DESIGN1, TWO_M_INFO[[CaseID]]$TWO_M, DISTANCE, dsLower, dsUpper,
                  IF_INNER_LBFGS = TRUE,
                  crit_type = "pair_fixed_true", MaxMinStdVals = NULL, ALG_INFO = ALG_INFO)
  list(res = out, eqv = eqv)
})

CaseID <- 3
plot(outAll[[CaseID]]$eqv$Grid_1, outAll[[CaseID]]$eqv$DirDeriv, type = "l", col = "blue"); abline(h = 0);
points(outAll[[CaseID]]$res$BESTDESIGN[,1], rep(0, nrow(outAll[[CaseID]]$res$BESTDESIGN)), pch = 19)



vals_MM <- sapply(1:length(outAll), function(i) outAll[[i]]$res$BESTVAL)
nSupp_MM <- 4
out <- DiscrimOD(MODEL_INFO, DISTANCE,
                 nSupp_MM, dsLower, dsUpper, crit_type = "maxmin_fixed_true",
                 IF_INNER_LBFGS = TRUE,
                 MaxMinStdVals = vals_MM, ALG_INFO, seed = NULL, verbose = TRUE)
round(out$BESTDESIGN, 3)

eqv <- equivalence(PSO_RESULT = out, MODEL_INFO = MODEL_INFO,
                   DISTANCE = DISTANCE, ngrid = 100, crit_type = "maxmin_fixed_true",
                   dsLower = dsLower, dsUpper = dsUpper, MaxMinStdVals = vals_MM, 
                   ALG_INFO = ALG_INFO)

plot(eqv$Grid_1, eqv$DirDeriv, type = "l", col = "blue"); abline(h = 0);
points(out$BESTDESIGN[,1], rep(0, nrow(out$BESTDESIGN)), pch = 16)

DESIGN1 <- out$BESTDESIGN
designCriterion(DESIGN1, MODEL_INFO, DISTANCE, dsLower, dsUpper,
                IF_INNER_LBFGS = TRUE, crit_type = "maxmin_fixed_true", 
                MaxMinStdVals = vals_MM, ALG_INFO = ALG_INFO)

