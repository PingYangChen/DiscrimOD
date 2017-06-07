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

