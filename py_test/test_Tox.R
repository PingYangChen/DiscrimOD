library(DiscrimOD)
library(Rcpp); library(RcppArmadillo); library(inline)

# R Function
m5_par <- c(4.282, 835.571, 0.739, 3.515)
MODEL_INFO <- list(
  # m5
  list(model = function(x, p) p[1]*(p[3] - (p[3] - 1)*exp(-(x/p[2])^p[4])),
       para = m5_par),
  # m4
  list(model = function(x, p) p[1]*(p[3] - (p[3] - 1)*exp(-(x/p[2]))),
       paraLower = c(0, 800, 0),
       paraUpper = c(10, 5000, 1),
       paraInit = m5_par[1:3]),
  # m3
  list(model = function(x, p) p[1]*exp(-(x/p[2])^p[3]),
       paraLower = c(0, 800, 1),
       paraUpper = c(10, 5000, 15),
       paraInit = m5_par[c(1, 2, 4)]),
  # m2
  list(model = function(x, p) p[1]*exp(-(x/p[2])),
       paraLower = c(0, 800),
       paraUpper = c(10, 5000),
       paraInit = m5_par[1:2]),
  # m1
  list(model = function(x, p) p[1],
       paraLower = c(0),
       paraUpper = c(10),
       paraInit = m5_par[1])
)
DISTANCE <- function(xt, xr) (xt - xr)^2


# C++ Funciton
tox5.inc <- 'arma::rowvec tox5(SEXP xx, SEXP pp){
  arma::rowvec x = Rcpp::as<arma::rowvec>(xx);
  arma::rowvec p = Rcpp::as<arma::rowvec>(pp);
  arma::rowvec val(x.n_elem);
  val = p(0)*(p(2) - (p(2) - 1.0)*arma::exp((-1.0)*arma::pow(x/p(1), p(3))));
  return val;
}'
tox5.body <- '
  typedef arma::rowvec (*funcPtr)(SEXP, SEXP);
  return(XPtr<funcPtr>(new funcPtr(&tox5)));
'

tox4.inc <- 'arma::rowvec tox4(SEXP xx, SEXP pp){
  arma::rowvec x = Rcpp::as<arma::rowvec>(xx);
  arma::rowvec p = Rcpp::as<arma::rowvec>(pp);
  arma::rowvec val(x.n_elem);
  val = p(0)*(p(2) - (p(2) - 1.0)*arma::exp((-1.0)*(x/p(1))));
  return val;
}'
tox4.body <- '
  typedef arma::rowvec (*funcPtr)(SEXP, SEXP);
  return(XPtr<funcPtr>(new funcPtr(&tox4)));
'

tox3.inc <- 'arma::rowvec tox3(SEXP xx, SEXP pp){
  arma::rowvec x = Rcpp::as<arma::rowvec>(xx);
  arma::rowvec p = Rcpp::as<arma::rowvec>(pp);
  arma::rowvec val(x.n_elem);
  val = p(0)*arma::exp((-1.0)*arma::pow(x/p(1), p(2)));
  return val;
}'
tox3.body <- '
  typedef arma::rowvec (*funcPtr)(SEXP, SEXP);
  return(XPtr<funcPtr>(new funcPtr(&tox3)));
'

tox2.inc <- 'arma::rowvec tox2(SEXP xx, SEXP pp){
  arma::rowvec x = Rcpp::as<arma::rowvec>(xx);
  arma::rowvec p = Rcpp::as<arma::rowvec>(pp);
  arma::rowvec val(x.n_elem);
  val = p(0)*arma::exp((-1.0)*(x/p(1)));
  return val;
}'
tox2.body <- '
  typedef arma::rowvec (*funcPtr)(SEXP, SEXP);
  return(XPtr<funcPtr>(new funcPtr(&tox2)));
'

tox1.inc <- 'arma::rowvec tox1(SEXP xx, SEXP pp){
  arma::rowvec x = Rcpp::as<arma::rowvec>(xx);
  arma::rowvec p = Rcpp::as<arma::rowvec>(pp);
  arma::rowvec val(x.n_elem);
  val.fill(p(0));
  return val;
}'
tox1.body <- '
  typedef arma::rowvec (*funcPtr)(SEXP, SEXP);
  return(XPtr<funcPtr>(new funcPtr(&tox1)));
'

tox5_Cpp <- cxxfunction(signature(), body = tox5.body, inc = tox5.inc,
                        plugin = "RcppArmadillo")
tox4_Cpp <- cxxfunction(signature(), body = tox4.body, inc = tox4.inc,
                        plugin = "RcppArmadillo")
tox3_Cpp <- cxxfunction(signature(), body = tox3.body, inc = tox3.inc,
                        plugin = "RcppArmadillo")
tox2_Cpp <- cxxfunction(signature(), body = tox2.body, inc = tox2.inc,
                        plugin = "RcppArmadillo")
tox1_Cpp <- cxxfunction(signature(), body = tox1.body, inc = tox1.inc,
                        plugin = "RcppArmadillo")

MODEL_INFO_Cpp <- list(
  # m5
  list(model = tox5_Cpp(),
       para = m5_par),
  # m4
  list(model = tox4_Cpp(),
       paraLower = c(0, 800, 0),
       paraUpper = c(10, 5000, 1),
       paraInit = m5_par[1:3]),
  # m3
  list(model = tox3_Cpp(),
       paraLower = c(0, 800, 1),
       paraUpper = c(10, 5000, 15),
       paraInit = m5_par[c(1, 2, 4)]),
  # m2
  list(model = tox2_Cpp(),
       paraLower = c(0, 800),
       paraUpper = c(10, 5000),
       paraInit = m5_par[1:2]),
  # m1
  list(model = tox1_Cpp(),
       paraLower = c(0),
       paraUpper = c(10),
       paraInit = m5_par[1])
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
ALG_INFO <- getAlgInfo(nSwarm = 32, maxIter = 200, typePSO = 0,
                       LBFGS_RETRY = 2,
                       FVAL_EPS = 0, GRAD_EPS = 1e-6, LINESEARCH_MAX = 1e9)

dsLower <- 0
dsUpper <- 1250

#SUB_MODEL_INFO_Cpp <- list(MODEL_INFO_Cpp[[1]], MODEL_INFO_Cpp[[2]])
#SUB_MODEL_INFO <- list(MODEL_INFO[[1]], MODEL_INFO[[2]])

TWO_M_INFO_Cpp <- list(
  list(TWO_M = list(MODEL_INFO_Cpp[[1]], MODEL_INFO_Cpp[[2]]), nSupp = 3),
  list(TWO_M = list(MODEL_INFO_Cpp[[1]], MODEL_INFO_Cpp[[3]]), nSupp = 4),
  list(TWO_M = list(MODEL_INFO_Cpp[[1]], MODEL_INFO_Cpp[[4]]), nSupp = 3),
  list(TWO_M = list(MODEL_INFO_Cpp[[1]], MODEL_INFO_Cpp[[5]]), nSupp = 2)
)

outALL <- NULL
outAll <- lapply(1:length(TWO_M_INFO_Cpp), function(CaseID) {
  CaseID <- 3
  out <- DiscrimOD(TWO_M_INFO_Cpp[[CaseID]]$TWO_M, DISTANCE_Cpp(),
                   TWO_M_INFO_Cpp[[CaseID]]$nSupp, dsLower, dsUpper,
                   crit_type = "pair_fixed_true",
                   MaxMinStdVals = NULL, ALG_INFO, seed = NULL, verbose = TRUE)
  round(out$BESTDESIGN, 3)

  eqv <- equivalence(PSO_RESULT = out, MODEL_INFO = TWO_M_INFO_Cpp[[CaseID]]$TWO_M,
                     DISTANCE = DISTANCE_Cpp(), ngrid = 100,
                     dsLower = dsLower, dsUpper = dsUpper, MaxMinStdVals = NULL, ALG_INFO = ALG_INFO)
  plot(eqv$Grid_1, eqv$DirDeriv, type = "l", col = "blue"); abline(h = 0);
  points(out$BESTDESIGN[,1], rep(0, nrow(out$BESTDESIGN)), pch = 19)
  DESIGN1 <- out$BESTDESIGN
  #DESIGN1 <- cbind(c(0, 468.186, 1064.179), c(.249, .498, .253))
  designCriterion(DESIGN1, TWO_M_INFO_Cpp[[CaseID]]$TWO_M, DISTANCE_Cpp(), dsLower, dsUpper,
                  MaxMinStdVals = NULL, ALG_INFO)
  list(res = out, eqv = eqv)
})

CaseID <- 2
plot(outAll[[CaseID]]$eqv$Grid_1, outAll[[CaseID]]$eqv$DirDeriv, type = "l", col = "blue"); abline(h = 0);
points(outAll[[CaseID]]$res$BESTDESIGN[,1], rep(0, nrow(outAll[[CaseID]]$res$BESTDESIGN)), pch = 19)



vals_MM <- sapply(1:length(outAll), function(i) outAll[[i]]$res$BESTVAL)
nSupp_MM <- 4
out <- DiscrimOD(MODEL_INFO_Cpp, DISTANCE_Cpp(),
                 nSupp_MM, dsLower, dsUpper, crit_type = "maxmin_fixed_true",
                 MaxMinStdVals = vals_MM, ALG_INFO, seed = NULL, verbose = TRUE)
round(out$BESTDESIGN, 3)




