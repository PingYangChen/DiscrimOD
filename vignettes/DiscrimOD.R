## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(cache = 2, cache.lazy = FALSE, tidy = FALSE, warning = FALSE)
library(knitcitations)

## ---- eval=FALSE, echo=TRUE----------------------------------------------
#  install.packages(c("devtools", "Rcpp", "RcppArmadillo"))
#  devtools::install_github("PingYangChen/DiscrimOD")

## ---- cache = FALSE, eval = FALSE, echo = TRUE---------------------------
#  # xt is the vector of mean values of the true model at each support point
#  # xr is the vector of mean values of the rival model at each support point
#  DISTANCE <- function(xt, xr) "write the R codes for the distance measure between 'xt' and 'xr'"

## ---- cache = FALSE, eval = FALSE, echo = TRUE---------------------------
#  # heteroscedastic Gaussian case
#  heter_norm <- function(xt, xr) {
#    var_t <- xt^2; var_r <- xr^2
#    (var_t + (xt - xr)^2)/var_r - log(var_t/var_r)
#  }

## ----call_empty, cache = FALSE, echo = TRUE------------------------------
emptyModelList(N_model = 2)

## ---- cache = FALSE, eval = FALSE, echo = TRUE---------------------------
#  DiscrimOD(MODEL_INFO, DISTANCE, nSupp, dsLower, dsUpper,
#            crit_type = "pair_fixed_true", MaxMinStdVals = NULL,
#  					PSO_INFO = NULL, LBFGS_INFO = NULL, seed = NULL, verbose = TRUE)
#  
#  equivalence(ngrid = 100, DESIGN = NULL, PSO_RESULT = NULL, MODEL_INFO, DISTANCE, dsLower, dsUpper,
#              crit_type = "pair_fixed_true", MaxMinStdVals = NULL,
#  						PSO_INFO = NULL, LBFGS_INFO = NULL, ALPHA_PSO_INFO = NULL)
#  
#  designCriterion(DESIGN1, MODEL_INFO, DISTANCE, dsLower, dsUpper,
#                  crit_type = "pair_fixed_true", MaxMinStdVals = NULL,
#  								PSO_INFO = NULL, LBFGS_INFO = NULL)

## ----TOX_MODEL-----------------------------------------------------------
tox5_par <- c(4.282, 835.571, 0.739, 3.515)
tox_info <- list(
  # tox_5
  list(model = function(x, p) p[1]*(p[3] - (p[3] - 1)*exp(-(x/p[2])^p[4])),
       para = tox5_par),
  # tox_4
  list(model = function(x, p) p[1]*(p[3] - (p[3] - 1)*exp(-(x/p[2]))),
       paraLower = c(0, 0, 0),
       paraUpper = c(20, 5000, 1)),
  # tox_3
  list(model = function(x, p) p[1]*exp(-(x/p[2])^p[3]),
       paraLower = c(0, 0, 1),
       paraUpper = c(20, 5000, 15)),
  # tox_2
  list(model = function(x, p) p[1]*exp(-(x/p[2])),
       paraLower = c(0, 0),
       paraUpper = c(20, 5000)),
  # tox_1
  list(model = function(x, p) rep(p[1], length(x)),
       paraLower = c(0),
       paraUpper = c(20))
)

## ----TOX_PAIR------------------------------------------------------------
# Number of support points
two_tox_n <- c(3, 4, 3, 2)
# Cases of pairwise discrimination
two_tox_info <- list(
  list(tox_info[[1]], tox_info[[2]]), # tox_5 vs tox_4
  list(tox_info[[1]], tox_info[[3]]), # tox_5 vs tox_3
  list(tox_info[[1]], tox_info[[4]]), # tox_5 vs tox_2
  list(tox_info[[1]], tox_info[[5]])  # tox_5 vs tox_1
)
#
tox_res_pair <- vector("list", 4)
for (CaseID in 1:4) {
  res <- DiscrimOD(MODEL_INFO = two_tox_info[[CaseID]], DISTANCE = sq_diff, 
                   nSupp = two_tox_n[CaseID], dsLower = 0, dsUpper = 1250, 
                   crit_type = "pair_fixed_true", 
                   PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO, seed = 100, verbose = FALSE)
  
  eqv <- equivalence(ngrid = 100, PSO_RESULT = res, 
                     MODEL_INFO = two_tox_info[[CaseID]], DISTANCE = sq_diff,
                     dsLower = 0, dsUpper = 1250,
                     crit_type = "pair_fixed_true", 
                     PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)
  
  tox_res_pair[[CaseID]] <- list(res = res, eqv = eqv)
}

## ---- cache=F, eval=F, echo=T--------------------------------------------
#  # R Code
#  l2_diff <- function(xt, xr) (xt - xr)^2
#  # C++ Code
#  l2_diff_cpp <- cppFunction('Rcpp::NumericVector l2_diff(SEXP xt, SEXP xr) {
#    arma::rowvec eta_t = Rcpp::as<arma::rowvec>(xt);
#    arma::rowvec eta_r = Rcpp::as<arma::rowvec>(xr);
#    arma::rowvec div = (eta_t - eta_r) % (eta_t - eta_r);
#    return Rcpp::wrap(div);
#  }', depends = "RcppArmadillo")

## ---- cache=F, eval=F, echo=T--------------------------------------------
#  # R Code
#  heter_norm <- function(xt, xr) {
#    var_t <- xt^2
#    var_r <- xr^2
#    (var_t + (xt - xr)^2)/var_r - log(var_t/var_r)
#  }
#  # C++ Code
#  heter_norm_cpp <- cppFunction('Rcpp::NumericVector heter_norm(SEXP xt, SEXP xr) {
#    arma::rowvec eta_t = Rcpp::as<arma::rowvec>(xt);
#    arma::rowvec eta_r = Rcpp::as<arma::rowvec>(xr);
#    arma::rowvec var_t = eta_t % eta_t;
#    arma::rowvec var_r = eta_r % eta_r;
#    arma::rowvec div = (var_t + arma::pow(eta_t - eta_r, 2))/var_r - arma::log(var_t/var_r);
#    return Rcpp::wrap(div);
#  }', depends = "RcppArmadillo")

## ---- cache=F, eval=F, echo=T--------------------------------------------
#  # R Code
#  log_norm_A <- function(xt, xr) {
#    vSQ <- 1.0
#    s_t <- log(1.0 + (vSQ/(xt^2)))
#    s_r <- log(1.0 + (vSQ/(xr^2)))
#    mu_t <- log(xt) - s_t
#    mu_r <- log(xr) - s_r
#    0.5*log(s_r/s_t) - (s_r - s_t - (mu_r - mu_t)^2)/(2*s_r)
#  }
#  # C++ Code
#  log_norm_A_cpp <- cppFunction('Rcpp::NumericVector log_norm_A(SEXP xt, SEXP xr) {
#    double vSQ = 1.0;
#    arma::rowvec eta_t = Rcpp::as<arma::rowvec>(xt);
#    arma::rowvec eta_r = Rcpp::as<arma::rowvec>(xr);
#    arma::rowvec s_t = arma::log(1.0 + (vSQ/(eta_t % eta_t)));
#    arma::rowvec s_r = arma::log(1.0 + (vSQ/(eta_r % eta_r)));
#    arma::rowvec mu_t = arma::log(eta_t) - s_t;
#    arma::rowvec mu_r = arma::log(eta_r) - s_r;
#    arma::rowvec div = 0.5*arma::log(s_r/s_t) -
#      (s_r - s_t - (mu_r - mu_t) % (mu_r - mu_t))/(2*s_r);
#    return Rcpp::wrap(div);
#  }', depends = "RcppArmadillo")

## ---- cache=F, eval=F, echo=T--------------------------------------------
#  # R Code
#  log_norm_B <- function(xt, xr) {
#    sigsq <- 1.0
#    var_t <- (exp(sigsq) - 1.0)*(xt^2)
#    var_r <- (exp(sigsq) - 1.0)*(xr^2)
#    mu_t <- log(xt) - 0.5*log(1.0 + (var_t/(xt^2)))
#    mu_r <- log(xr) - 0.5*log(1.0 + (var_r/(xr^2)))
#    ((mu_r - mu_t)^2)/(2*sigsq)
#  }
#  # C++ Code
#  log_norm_B_cpp <- cppFunction('Rcpp::NumericVector log_norm_B(SEXP xt, SEXP xr) {
#    double sigsq = 1.0;
#    arma::rowvec eta_t = Rcpp::as<arma::rowvec>(xt);
#    arma::rowvec eta_r = Rcpp::as<arma::rowvec>(xr);
#    arma::rowvec var_t = (std::exp(sigsq) - 1.0)*(eta_t % eta_t);
#    arma::rowvec var_r = (std::exp(sigsq) - 1.0)*(eta_r % eta_r);
#    arma::rowvec mu_t = arma::log(eta_t) - 0.5*arma::log(1.0 + (var_t/(eta_t % eta_t)));
#    arma::rowvec mu_r = arma::log(eta_r) - 0.5*arma::log(1.0 + (var_r/(eta_r % eta_r)));
#    arma::rowvec div = ((mu_r - mu_t) % (mu_r - mu_t))/(2*sigsq);
#    return Rcpp::wrap(div);
#  }', depends = "RcppArmadillo")

## ---- cache=F, eval=F, echo=T--------------------------------------------
#  # R Code
#  log_norm_C <- function(xt, xr) {
#    c <- 1; d <- 1
#    var_t <- d*exp(c*xt)
#    var_r <- d*exp(c*xr)
#    s_t <- log(1.0 + (var_t/(xt^2)))
#    s_r <- log(1.0 + (var_r/(xr^2)))
#    mu_t <- log(xt) - s_t
#    mu_r <- log(xr) - s_r
#    0.5*log(s_r/s_t) - (s_r - s_t - (mu_r - mu_t)*(mu_r - mu_t))/(2*s_r)
#  }
#  # C++ Code
#  log_norm_C_cpp <- cppFunction('Rcpp::NumericVector log_norm_C(SEXP xt, SEXP xr) {
#    double c = 1.0; double d = 1.0;
#    arma::rowvec eta_t = Rcpp::as<arma::rowvec>(xt);
#    arma::rowvec eta_r = Rcpp::as<arma::rowvec>(xr);
#    arma::rowvec var_t = d * arma::exp(c*eta_t);
#    arma::rowvec var_r = d * arma::exp(c*eta_r);
#    arma::rowvec s_t = arma::log(1.0 + (var_t/(eta_t % eta_t)));
#    arma::rowvec s_r = arma::log(1.0 + (var_r/(eta_r % eta_r)));
#    arma::rowvec mu_t = arma::log(eta_t) - s_t;
#    arma::rowvec mu_r = arma::log(eta_r) - s_r;
#    arma::rowvec div = 0.5*arma::log(s_r/s_t) -
#      (s_r - s_t - (mu_r - mu_t) % (mu_r - mu_t))/(2*s_r);
#    return Rcpp::wrap(div);
#  }', depends = "RcppArmadillo")

## ---- cache=F, eval=F, echo=T--------------------------------------------
#  # R Code
#  logit_diff <- function(xt, xr) {
#    exp_t <- exp(xt)
#    exp_r <- exp(xr)
#    mu_t <- exp_t/(1 + exp_t)
#    mu_r <- exp_r/(1 + exp_r)
#    mu_t*(log(mu_t) - log(mu_r)) + (1 - mu_t)*(log(1.0 - mu_t) - log(1.0 - mu_r))
#  }
#  # C++ Code
#  logit_diff_cpp <- cppFunction('Rcpp::NumericVector logit_diff(SEXP xt, SEXP xr) {
#    arma::rowvec eta_t = Rcpp::as<arma::rowvec>(xt);
#    arma::rowvec eta_r = Rcpp::as<arma::rowvec>(xr);
#    arma::rowvec exp_t = arma::exp(eta_t);
#    arma::rowvec exp_r = arma::exp(eta_r);
#    arma::rowvec mu_t = exp_t/(1.0 + exp_t);
#    arma::rowvec mu_r = exp_r/(1.0 + exp_r);
#    arma::rowvec div = mu_t % (arma::log(mu_t) - arma::log(mu_r)) +
#      (1.0 - mu_t) % (arma::log(1.0 - mu_t) - arma::log(1.0 - mu_r));
#    return Rcpp::wrap(div);
#  }', depends = "RcppArmadillo")

## ---- cache=F, eval=F, echo=T--------------------------------------------
#  # R Code
#  gamma_diff <- function(xt, xr) log(xr/xt) + (xt - xr)/xr
#  # C++ Code
#  gamma_diff_cpp <- cppFunction('Rcpp::NumericVector gamma_diff(SEXP xt, SEXP xr) {
#    arma::rowvec eta_t = Rcpp::as<arma::rowvec>(xt);
#    arma::rowvec eta_r = Rcpp::as<arma::rowvec>(xr);
#    arma::rowvec div = arma::log(eta_r/eta_t) + (eta_t - eta_r)/eta_r;
#    return Rcpp::wrap(div);
#  }', depends = "RcppArmadillo")

