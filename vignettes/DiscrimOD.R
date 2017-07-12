## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(cache = 2, cache.lazy = FALSE, tidy = FALSE, warning = FALSE)
library(knitcitations)
set.seed(round(runif(1, 1, 1e4)))

## ---- eval=FALSE, echo=TRUE----------------------------------------------
#  install.packages(c("devtools", "Rcpp", "RcppArmadillo"))
#  devtools::install_github("PingYangChen/DiscrimOD")

## ---- cache = FALSE, eval = FALSE, echo = TRUE---------------------------
#  # xt is the vector of mean values of the true model at each support point
#  # xr is the vector of mean values of the rival model at each support point
#  DISTANCE <- function(xt, xr)
#    "write the R codes for the distance measure between 'xt' and 'xr'"

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
#            PSO_INFO = NULL, LBFGS_INFO = NULL, seed = NULL, verbose = TRUE)
#  
#  equivalence(ngrid = 100, DESIGN = NULL, PSO_RESULT = NULL, MODEL_INFO, DISTANCE,
#              dsLower, dsUpper,
#              crit_type = "pair_fixed_true", MaxMinStdVals = NULL,
#              PSO_INFO = NULL, LBFGS_INFO = NULL, ALPHA_PSO_INFO = NULL)
#  
#  designCriterion(DESIGN1, MODEL_INFO, DISTANCE, dsLower, dsUpper,
#                  crit_type = "pair_fixed_true", MaxMinStdVals = NULL,
#                  PSO_INFO = NULL, LBFGS_INFO = NULL)

## ----MM_MODEL------------------------------------------------------------
# Two types of Michaelise-Menten models
enzyme2 <- function(x, p) p[1]*x/(p[2] + x) + p[3]*x  # Modified Michaelis-Menten
enzyme1 <- function(x, p) p[1]*x/(p[2] + x)           # Michaelise-Menten
# Specify the nominal value for the true model, 'linlogi4'.
para_mmm_2 <- c(1, 1, 1)
# Create model list
model_mmm <- list(
  list(model = enzyme2, para = para_mmm_2),
  list(model = enzyme1, paraLower = c(1e-4, 1e-4), paraUpper = c(20, 20))
)

## ----LOGNORM_B_GAMMA-----------------------------------------------------
# Log-normal model
log_norm_B <- function(xt, xr) {
  sigsq <- 1 # nuisance parameter
  var_t <- (exp(sigsq) - 1.0)*(xt^2) # variance (true mdoel)
  var_r <- (exp(sigsq) - 1.0)*(xr^2) # variance (rival model)
  mu_t <- log(xt) - 0.5*log(1.0 + (var_t/(xt^2))) # mean (true model)
  mu_r <- log(xr) - 0.5*log(1.0 + (var_r/(xr^2))) # mean (rival model)
  ((mu_r - mu_t)^2)/(2*sigsq) # KL-divergence
}
# Gamma model
gamma_diff <- function(xt, xr) log(xr/xt) + (xt - xr)/xr # KL-divergence

## ----MM_RES_LOGNORM_B----------------------------------------------------
# Run for finding KL-optimal design under lognormal assumption
mmm_res_lognB <- DiscrimOD(MODEL_INFO = model_mmm, DISTANCE = log_norm_B, 
                           nSupp = 3, dsLower = 0.1, dsUpper = 5, 
                           crit_type = "pair_fixed_true", 
                           PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO, 
                           seed = NULL, verbose = FALSE)

## ----MM_DEISGN_LOGNORM_B-------------------------------------------------
round(mmm_res_lognB$BESTDESIGN, 3) 

## ----MM_RES_GAMMA--------------------------------------------------------
# Run for finding KL-optimal design under lognormal assumption
mmm_res_gamma <- DiscrimOD(MODEL_INFO = model_mmm, DISTANCE = gamma_diff, 
                           nSupp = 3, dsLower = 0.1, dsUpper = 5, 
                           crit_type = "pair_fixed_true", 
                           PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO, 
                           seed = NULL, verbose = FALSE)

## ----MM_DEISGN_GAMMA-----------------------------------------------------
round(mmm_res_gamma$BESTDESIGN, 3) 

## ----MM_EQV--------------------------------------------------------------
# Check equivalence theorem for KL-optimal design under lognormal error assumption
mmm_eqv_lognB <- equivalence(ngrid = 100, PSO_RESULT = mmm_res_lognB, 
                             MODEL_INFO = model_mmm, DISTANCE = log_norm_B, 
                             dsLower = 0.1, dsUpper = 5,  
                             crit_type = "pair_fixed_true", 
                             PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)
# Check equivalence theorem for KL-optimal design under gamma error assumption
mmm_eqv_gamma <- equivalence(ngrid = 100, PSO_RESULT = mmm_res_gamma, 
                             MODEL_INFO = model_mmm, DISTANCE = gamma_diff, 
                             dsLower = 0.1, dsUpper = 5,  
                             crit_type = "pair_fixed_true", 
                             PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

## ----LOGIT_MODEL---------------------------------------------------------
# Tommasi et al. (2016)
linlogi4 <- function(x, p) p[1] + p[2]*x + p[3]*(x^2)  # Logistic 4
linlogi3 <- function(x, p) x*(p[1] + p[2]*x)           # Logistic 3
linlogi2 <- function(x, p) p[1] + p[2]*x               # Logistic 2
linlogi1 <- function(x, p) p[1]*x                      # Logistic 1
# Specify the nominal value for the true model, 'linlogi4'.
para_logit_4 <- c(1, 1, 1)
# Create model list
model_logit <- list(
  list(model = linlogi4, para = para_logit_4),
  list(model = linlogi3, paraLower = c(-10, -10), paraUpper = c(10, 10)),
  list(model = linlogi2, paraLower = c(-10, -10), paraUpper = c(10, 10)),
  list(model = linlogi1, paraLower = -10, paraUpper = 10)
)

## ----LOGIT---------------------------------------------------------------
logit_diff <- function(xt, xr) {
  exp_t <- exp(xt)
  exp_r <- exp(xr)
  mu_t <- exp_t/(1 + exp_t)
  mu_r <- exp_r/(1 + exp_r)
  # Return KL-divergence
  mu_t*(log(mu_t) - log(mu_r)) + (1 - mu_t)*(log(1 - mu_t) - log(1 - mu_r))
}

## ----LOGIT_PAIR----------------------------------------------------------
# Cases of pairwise discrimination
two_model_logit <- list(
  list(model_logit[[1]], model_logit[[2]]),
  list(model_logit[[1]], model_logit[[3]]),
  list(model_logit[[1]], model_logit[[4]])
)
# Number of support points 
two_logit_nSupp <- c(3, 3, 3)
# Run for each case
logit_res_pair <- vector("list", 3)
for (CaseID in 1:3) {
  res <- DiscrimOD(MODEL_INFO = two_model_logit[[CaseID]], DISTANCE = logit_diff, 
                   nSupp = two_logit_nSupp[CaseID], dsLower = 0, dsUpper = 1, 
                   crit_type = "pair_fixed_true", 
                   PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO,
                   seed = NULL, verbose = FALSE) 
  
  eqv <- equivalence(ngrid = 100, PSO_RESULT = res, 
                     MODEL_INFO = two_model_logit[[CaseID]], DISTANCE = logit_diff,
                     dsLower = 0, dsUpper = 1, 
                     crit_type = "pair_fixed_true", 
                     PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)
  
  logit_res_pair[[CaseID]] <- list(res = res, eqv = eqv)
}

## ----LOGIT_RES_MM--------------------------------------------------------
# Get criterion values from pairwise discrimination results
logit_vals_pair <- sapply(1:length(logit_res_pair), function(i) {
  logit_res_pair[[i]]$res$BESTVAL
})
# Run for finding max-min KL-optimal design 
logit_res_mxmn <- DiscrimOD(MODEL_INFO = model_logit, DISTANCE = logit_diff, 
                            nSupp = 3, dsLower = 0, dsUpper = 1, 
                            crit_type = "maxmin_fixed_true", 
                            MaxMinStdVals = logit_vals_pair,
                            PSO_INFO = PSO_INFO_MAXMIN, LBFGS_INFO = LBFGS_INFO, 
                            seed = NULL, verbose = FALSE)

## ----LOGIT_DESIGN_MM-----------------------------------------------------
round(logit_res_mxmn$BESTDESIGN, 3) 

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

## ----LOGNORM_B_TOX-------------------------------------------------------
log_norm_B <- function(xt, xr) {
  sigsq <- 0.01 # nuisance parameter
  var_t <- (exp(sigsq) - 1.0)*(xt^2) # variance (true mdoel)
  var_r <- (exp(sigsq) - 1.0)*(xr^2) # variance (rival model)
  mu_t <- log(xt) - 0.5*log(1.0 + (var_t/(xt^2))) # mean (true model)
  mu_r <- log(xr) - 0.5*log(1.0 + (var_r/(xr^2))) # mean (rival model)
  ((mu_r - mu_t)^2)/(2*sigsq) # KL-divergence
}

## ----TOX_PAIR------------------------------------------------------------
# Cases of pairwise discrimination
two_tox_info <- list(
  list(tox_info[[1]], tox_info[[2]]), # tox_5 vs tox_4
  list(tox_info[[1]], tox_info[[3]]), # tox_5 vs tox_3
  list(tox_info[[1]], tox_info[[4]]), # tox_5 vs tox_2
  list(tox_info[[1]], tox_info[[5]])  # tox_5 vs tox_1 
)
# Number of support points 
two_tox_nSupp <- c(3, 4, 3, 2)
# Run for each case
tox_res_pair <- vector("list", 4)
for (CaseID in 1:4) {
  res <- DiscrimOD(MODEL_INFO = two_tox_info[[CaseID]], DISTANCE = log_norm_B, 
                   nSupp = two_tox_nSupp[CaseID], dsLower = 0, dsUpper = 1250, 
                   crit_type = "pair_fixed_true", 
                   PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO,
                   seed = NULL, verbose = FALSE) 
  
  eqv <- equivalence(ngrid = 100, PSO_RESULT = res, 
                     MODEL_INFO = two_tox_info[[CaseID]], DISTANCE = log_norm_B,
                     dsLower = 0, dsUpper = 1250, 
                     crit_type = "pair_fixed_true", 
                     PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)
  
  tox_res_pair[[CaseID]] <- list(res = res, eqv = eqv)
}

## ----TOX_RES_MM----------------------------------------------------------
# Get criterion values from pairwise discrimination results
tox_vals_pair <- sapply(1:length(tox_res_pair), function(i) tox_res_pair[[i]]$res$BESTVAL)
# Run for finding max-min KL-optimal design 
tox_res_maxmin <- DiscrimOD(MODEL_INFO = tox_info, DISTANCE = log_norm_B, 
                            nSupp = 4, dsLower = 0, dsUpper = 1250, 
                            crit_type = "maxmin_fixed_true", 
                            MaxMinStdVals = tox_vals_pair,
                            PSO_INFO = PSO_INFO_MAXMIN, LBFGS_INFO = LBFGS_INFO, 
                            seed = NULL, verbose = FALSE)

## ----TOX_DESIGN_MM-------------------------------------------------------
round(tox_res_maxmin$BESTDESIGN, 3) 

## ----TOX_VAL_MM----------------------------------------------------------
tox_res_maxmin$BESTVAL 

## ---- cache=F, eval=F, echo=T--------------------------------------------
#  getPSOInfo(nSwarm, maxIter, checkConv = 0, freeRun = 0.25, tol = 1e-06,
#             c1 = 2.05, c2 = 2.05, w0 = 0.95, w1 = 0.2, w_var = 0.8, vk = 4)

## ----cache=F,eval=F,echo=T-----------------------------------------------
#  getLBFGSInfo(IF_INNER_LBFGS = TRUE, LBFGS_RETRY = 1, LBFGS_MAXIT = 0,
#               LBFGS_LM = 6, FVAL_EPS = 0, GRAD_EPS = 1e-05,
#               LINESEARCH_MAXTRIAL = 20, LINESEARCH_MAX = 1e+20,
#               LINESEARCH_MIN = 1e-20, LINESEARCH_ARMIJO = 1e-04,
#               LINESEARCH_WOLFE = 0.9)

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

