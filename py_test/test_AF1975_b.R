library(DiscrimOD)
library(Rcpp); library(RcppArmadillo); library(inline); #library(assertthat)

#path <- "D:/Ping_Yang/Google Drive/PYChen_Statistics_NCKU/Researches/2015 min-max optimal discriminating designs/DiscrimOD"
#cppPath <- file.path(path, "src")
#rPath <- file.path(path, "R")
#sourceCpp(file.path(cppPath, "cppPSO.cpp"))
#source(file.path(rPath, "DiscrimOD_assist.R"))
#source(file.path(rPath, "DiscrimOD_tools.R"))
#source(file.path(rPath, "DiscrimOD.R"))


AF_para_m1 <- c(4.5, -1.5, -2)
# R Function
MODEL_INFO <- list(
  list(model = function(x, p) p[1] + p[2]*exp(x) + p[3]*exp(-x), para = AF_para_m1),
  list(model = function(x, p) p[1] + p[2]*x + p[3]*x^2,
       paraLower = rep(-10, 3),
       paraUpper = rep(10, 3),
       paraInit = runif(3, -10,10)),
  list(model = function(x, p) p[1] + p[2]*sin(0.5*pi*x) + p[3]*cos(0.5*pi*x) + p[4]*sin(pi*x),
       paraLower = rep(-10, 4),
       paraUpper = rep(10, 4),
       paraInit = c(0,0,0,0))
)
DISTANCE <- function(xt, xr) (xt - xr)^2


# C++ Funciton
m1.inc <- 'arma::rowvec m1(SEXP xx, SEXP pp){
arma::rowvec x = Rcpp::as<arma::rowvec>(xx);
arma::rowvec p = Rcpp::as<arma::rowvec>(pp);
arma::rowvec val(x.n_elem);
val = p(0) + p(1)*arma::exp(x) + p(2)*arma::exp(-x);
return val;
}'
m1.body <- '
typedef arma::rowvec (*funcPtr)(SEXP, SEXP);
return(XPtr<funcPtr>(new funcPtr(&m1)));
'

m2.inc <- 'arma::rowvec m2(SEXP xx, SEXP pp){
arma::rowvec x = Rcpp::as<arma::rowvec>(xx);
arma::rowvec p = Rcpp::as<arma::rowvec>(pp);
arma::rowvec val(x.n_elem);
val = p(0) + p(1)*x + p(2)*(x%x);
return val;
}'
m2.body <- '
typedef arma::rowvec (*funcPtr)(SEXP, SEXP);
return(XPtr<funcPtr>(new funcPtr(&m2)));
'

m3.inc <- 'arma::rowvec m3(SEXP xx, SEXP pp){
arma::rowvec x = Rcpp::as<arma::rowvec>(xx);
arma::rowvec p = Rcpp::as<arma::rowvec>(pp);
arma::rowvec val(x.n_elem);
val = p(0) + p(1)*arma::sin(0.5*M_PI*x) + p(2)*arma::cos(0.5*M_PI*x) + p(3)*arma::sin(M_PI*x);
return val;
}'
m3.body <- '
typedef arma::rowvec (*funcPtr)(SEXP, SEXP);
return(XPtr<funcPtr>(new funcPtr(&m3)));
'

m1_Cpp <- cxxfunction(signature(), body = m1.body, inc = m1.inc,
                      plugin = "RcppArmadillo")
m2_Cpp <- cxxfunction(signature(), body = m2.body, inc = m2.inc,
                      plugin = "RcppArmadillo")
m3_Cpp <- cxxfunction(signature(), body = m3.body, inc = m3.inc,
                      plugin = "RcppArmadillo")

MODEL_INFO_Cpp <- list(
  list(model = m1_Cpp(), para = AF_para_m1),
  list(model = m2_Cpp(), paraLower = rep(-10, 3),
       paraUpper = rep(10, 3),
       paraInit = c(0,0,0)),
  list(model = m3_Cpp(), paraLower = rep(-10, 4),
       paraUpper = rep(10, 4),
       paraInit = c(0,0,0,0))
)

dist.inc <- 'arma::rowvec t_optimal(SEXP xt, SEXP xr){
arma::rowvec val_t = Rcpp::as<arma::rowvec>(xt);
arma::rowvec val_r = Rcpp::as<arma::rowvec>(xr);
arma::rowvec val(val_t.n_elem);
val = (val_t - val_r)%(val_t - val_r);
return val;
}'

dist.body <- '
typedef arma::rowvec (*funcPtr)(SEXP, SEXP);
return(XPtr<funcPtr>(new funcPtr(&t_optimal)));
'

DISTANCE_Cpp <- cxxfunction(signature(), body = dist.body, inc = dist.inc,
                            plugin = "RcppArmadillo")

#
ALG_INFO <- getAlgInfo(nSwarm = 32, maxIter = 200, typePSO = 2,
                       LBFGS_RETRY = 1,
                       FVAL_EPS = 0, GRAD_EPS = 1e-6,
                       LINESEARCH_MAX = 1e6)

dsLower <- -1
dsUpper <- 1

TWO_M_INFO_Cpp <- list(
  list(TWO_M = list(MODEL_INFO_Cpp[[1]], MODEL_INFO_Cpp[[2]]), nSupp = 4),
  list(TWO_M = list(MODEL_INFO_Cpp[[1]], MODEL_INFO_Cpp[[3]]), nSupp = 5)
)

outALL <- vector("list", length(TWO_M_INFO_Cpp))
for (CaseID in 1:length(TWO_M_INFO_Cpp)) {

  #sourceCpp(file.path(cppPath, "cppPSO.cpp"))
  #CaseID <- 1
  out <- DiscrimOD(TWO_M_INFO_Cpp[[CaseID]]$TWO_M, DISTANCE_Cpp(),
                   TWO_M_INFO_Cpp[[CaseID]]$nSupp, dsLower, dsUpper,
                   crit_type = "pair_fixed_true",
                   MaxMinStdVals = NULL, ALG_INFO, seed = NULL, verbose = TRUE)
  round(out$BESTDESIGN, 3)

  eqv <- equivalence(PSO_RESULT = out, MODEL_INFO = TWO_M_INFO_Cpp[[CaseID]]$TWO_M,
                     DISTANCE = DISTANCE_Cpp(), ngrid = 100, crit_type = "pair_fixed_true",
                     dsLower = dsLower, dsUpper = dsUpper, MaxMinStdVals = NULL, ALG_INFO = ALG_INFO)
  plot(eqv$Grid_1, eqv$DirDeriv, type = "l", col = "blue"); abline(h = 0);
  points(out$BESTDESIGN[,1], rep(0, nrow(out$BESTDESIGN)), pch = 19)
  DESIGN1 <- out$BESTDESIGN
  #DESIGN1 <- cbind(c(0, 468.186, 1064.179), c(.249, .498, .253))
  designCriterion(DESIGN1, TWO_M_INFO_Cpp[[CaseID]]$TWO_M, DISTANCE_Cpp(), dsLower, dsUpper,
                  crit_type = "pair_fixed_true", MaxMinStdVals = NULL, ALG_INFO)

  outALL[[CaseID]] <- list(res = out, eqv = eqv)
}


CaseID <- 1
outALL[[CaseID]]$res
plot(outALL[[CaseID]]$eqv$Grid_1, outALL[[CaseID]]$eqv$DirDeriv, type = "l", col = "blue"); abline(h = 0);
points(outALL[[CaseID]]$res$BESTDESIGN[,1], rep(0, nrow(outALL[[CaseID]]$res$BESTDESIGN)), pch = 19)


vals_MM <- sapply(1:length(outALL), function(i) outALL[[i]]$res$BESTVAL)
#vals_MM <- c(0.1,0.1)
nSupp_MM <- 5
out_MM <- DiscrimOD(MODEL_INFO_Cpp, DISTANCE_Cpp(),
                    nSupp_MM, dsLower, dsUpper, crit_type = "maxmin_fixed_true",
                    MaxMinStdVals = vals_MM, ALG_INFO, seed = NULL, verbose = TRUE)
round(out_MM$BESTDESIGN, 3)

PARA_SET <- designCriterion(out_MM$BESTDESIGN, MODEL_INFO_Cpp, DISTANCE_Cpp(), dsLower, dsUpper,
                crit_type = "maxmin_fixed_true", MaxMinStdVals = vals_MM, ALG_INFO)$theta2

eqv <- equivalence(PSO_RESULT = out_MM, MODEL_INFO = MODEL_INFO_Cpp, DISTANCE = DISTANCE_Cpp(),
                   ngrid = 100, crit_type = "maxmin_fixed_true",
                   dsLower = dsLower, dsUpper = dsUpper,
                   MaxMinStdVals = vals_MM, ALG_INFO = ALG_INFO)

plot(eqv$Grid_1, eqv$DirDeriv, type = "l", col = "blue"); abline(h = 0);
points(out_MM$BESTDESIGN[,1], rep(0, nrow(out_MM$BESTDESIGN)), pch = 19)


