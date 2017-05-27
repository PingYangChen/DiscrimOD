library(DiscrimOD)
library(Rcpp); library(RcppArmadillo); library(inline)

AF_para_m1 <-c (4.5, -1.5, -2)
# R Function
MODEL_INFO <- list(
  list(model = function(x, p) p[1] + p[2]*exp(x) + p[3]*exp(-x), para = AF_para_m1),
  list(model = function(x, p) p[1] + p[2]*x + p[3]*x^2,
       paraLower = rep(-10, 3),
       paraUpper = rep(10, 3),
       paraInit = runif(3, -10,10))
)
DISTANCE <- function(xt, xr) (xt - xr)^2


# C++ Funciton
m1.inc <- 'arma::rowvec m1(SEXP xx, SEXP pp){
  arma::rowvec x = Rcpp::as<arma::rowvec>(xx);
  arma::rowvec p = Rcpp::as<arma::rowvec>(pp);
  return p(0) + p(1)*arma::exp(x) + p(2)*arma::exp(-x);
}'
m1.body <- '
  typedef arma::rowvec (*funcPtr)(SEXP, SEXP);
  return(XPtr<funcPtr>(new funcPtr(&m1)));
'

m2.inc <- 'arma::rowvec m2(SEXP xx, SEXP pp){
  arma::rowvec x = Rcpp::as<arma::rowvec>(xx);
  arma::rowvec p = Rcpp::as<arma::rowvec>(pp);
  return p(0) + p(1)*x + p(2)*(x%x);
}'
m2.body <- '
  typedef arma::rowvec (*funcPtr)(SEXP, SEXP);
  return(XPtr<funcPtr>(new funcPtr(&m2)));
'

m1_Cpp <- cxxfunction(signature(), body = m1.body, inc = m1.inc,
                      plugin = "RcppArmadillo")
m2_Cpp <- cxxfunction(signature(), body = m2.body, inc = m2.inc,
                      plugin = "RcppArmadillo")

MODEL_INFO_Cpp <- list(
  list(model = m1_Cpp(), para = AF_para_m1),
  list(model = m2_Cpp(), paraLower = rep(-10, 3),
                         paraUpper = rep(10, 3),
                         paraInit = c(0,0,0))
)

dist.inc <- 'arma::rowvec t_optimal(SEXP xt, SEXP xr){
  arma::rowvec val_t = Rcpp::as<arma::rowvec>(xt);
  arma::rowvec val_r = Rcpp::as<arma::rowvec>(xr);
  return (val_t - val_r)%(val_t - val_r);
}'

dist.body <- '
  typedef arma::rowvec (*funcPtr)(SEXP, SEXP);
  return(XPtr<funcPtr>(new funcPtr(&t_optimal)));
'

DISTANCE_Cpp <- cxxfunction(signature(), body = dist.body, inc = dist.inc,
                            plugin = "RcppArmadillo")

#
ALG_INFO <- getAlgInfo(nSwarm = 32, maxIter = 200, typePSO = 0,
                       LBFGS_RETRY = 3,
                       FVAL_EPS = 0, GRAD_EPS = 1e-5,
                       LINESEARCH_MAX = 1, LINESEARCH_ARMIJO = 1e-4)

nSupp <- 4
dsLower <- -1
dsUpper <- 1
MaxMinStdVals = NULL; seed = NULL; verbose = TRUE

out <- DiscrimOD(MODEL_INFO_Cpp, DISTANCE_Cpp(), nSupp, dsLower, dsUpper,
                 crit_type = "pair_fixed_true",
                 MaxMinStdVals, ALG_INFO = ALG_INFO, seed, verbose)
round(out$BESTDESIGN, 3)

eqv <- equivalence(PSO_RESULT = out, MODEL_INFO = MODEL_INFO_Cpp, DISTANCE = DISTANCE_Cpp(),
                   ngrid = 100, crit_type = "pair_fixed_true",
                   dsLower = dsLower, dsUpper = dsUpper,
                   MaxMinStdVals = NULL, ALG_INFO = ALG_INFO)

plot(eqv$Grid_1, eqv$DirDeriv, type = "l", col = "blue"); abline(h = 0);
points(out$BESTDESIGN[,1], rep(0, nrow(out$BESTDESIGN)), pch = 19)

nn <- 100
aa <- numeric(nn)
system.time(
for (i in 1:nn) {
DESIGN1 <- out$BESTDESIGN
#uu <- runif(4); DESIGN1 <- cbind(runif(4, -1, 1), uu/sum(uu))
DESIGN1 <- cbind(c(-1,-.669,.144,.957), c(.253,.428,.247,.072))
aa[i] <-
  designCriterion(DESIGN1, MODEL_INFO_Cpp, DISTANCE_Cpp(), dsLower, dsUpper,
                  crit_type = "pair_fixed_true",
                  MaxMinStdVals = NULL, ALG_INFO)$cri_val
})[3]/nn

max(aa)
sum(aa > 1)

