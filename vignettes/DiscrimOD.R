## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(cache = 2, cache.lazy = FALSE, tidy = FALSE, warning = FALSE)
library(knitcitations)

## ---- eval=FALSE, echo=TRUE----------------------------------------------
#  install.packages(c("devtools", "Rcpp", "RcppArmadillo"))
#  devtools::install_github("PingYangChen/DiscrimOD")

## ----LOAD_PKG, cache=TRUE,eval=TRUE,echo=TRUE----------------------------
library(DiscrimOD)

## ----ALG_SETTING, cache=T------------------------------------------------
# PSO basic settings
PSO_INFO <- getPSOInfo(nSwarm = 32, maxIter = 200)
# L-BFGS basic settings
LBFGS_INFO <- getLBFGSInfo(LBFGS_RETRY = 2)

## ----ALG_SETTING_MAXMIN, cache=T-----------------------------------------
PSO_INFO_MAXMIN <- getPSOInfo(nSwarm = 32, maxIter = 400)

## ----ALG_SETTING_NESTEDPSO, cache=T--------------------------------------
# Set NestedPSO options. The length of setting indicates the number of loops
NESTEDPSO_INFO <- getPSOInfo(nSwarm = c(32, 32), maxIter = c(200, 100))
# Turn off L-BFGS implementation for the inner optimization loop
LBFGS_NOTRUN <- getLBFGSInfo(IF_INNER_LBFGS = FALSE)

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

## ----AF1975b_MODEL, cache=T----------------------------------------------
af1 <- function(x, p) p[1] + p[2]*exp(x) + p[3]*exp(-x)
af2 <- function(x, p) p[1] + p[2]*x + p[3]*x^2
af3 <- function(x, p) p[1] + p[2]*sin(0.5*pi*x) + p[3]*cos(0.5*pi*x) + p[4]*sin(pi*x)

## ----AF1975b_CASE_PAIR, cache=T------------------------------------------
AF_para_af1 <- c(4.5, -1.5, -2)
# The first pair: \eta_1 vs \eta_2
af_info_12 <- list(
  # The first list should be the true model and the specified nominal values
  list(model = af1, para = AF_para_af1),
  # Then the rival models are listed accordingly. We also need to specify the model space.
  list(model = af2, paraLower = rep(-10, 3), paraUpper = rep(10, 3))
)
# The second pair: \eta_1 vs \eta_3
af_info_13 <- list(
  list(model = af1, para = AF_para_af1),
  list(model = af3, paraLower = rep(-10, 4), paraUpper = rep(10, 4))
)

## ----T_OPTIMAL, cache=T--------------------------------------------------
# xt is the mean values of the true model
# xr is the mean values of the rival model
sq_diff <- function(xt, xr) (xt - xr)^2

## ----AF1975b_RESULT_PAIR_12, cache=T-------------------------------------
# Run PSO-QN Algorithm for discriminating \eta_1 and \eta_2
af_res_12 <- DiscrimOD(MODEL_INFO = af_info_12, DISTANCE = sq_diff,
                       nSupp = 4, dsLower = -1, dsUpper = 1, 
                       crit_type = "pair_fixed_true",
                       PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO, 
                       seed = NULL, verbose = FALSE)

## ----AF1975b_DESIGN_PAIR_12, cache=T-------------------------------------
round(af_res_12$BESTDESIGN, 3) 

## ----AF1975b_VAL_PAIR_12, cache=T----------------------------------------
af_res_12$BESTVAL 

## ----AF1975b_CPU_PAIR_12, cache=T----------------------------------------
af_res_12$CPUTIME 

## ----AF1975b_EQUV_PAIR_12, cache=T, fig.align='center', fig.height = 4, fig.width = 4, fig.cap='Directional derivative function of T-optimal design in the case \'af_res_12\'.'----
af_eqv_12 <- equivalence(ngrid = 100, PSO_RESULT = af_res_12, 
                         MODEL_INFO = af_info_12, DISTANCE = sq_diff, 
                         dsLower = -1, dsUpper = 1,
                         crit_type = "pair_fixed_true",
                         PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)
# Draw the directional derivative curve
plot(af_eqv_12$Grid_1, af_eqv_12$DirDeriv, type = "l", col = "blue", 
     main = "af_res_12", xlab = "x", ylab = "Directional Derivative"); abline(h = 0)
points(af_res_12$BESTDESIGN[,1], rep(0, nrow(af_res_12$BESTDESIGN)), pch = 16)

## ----AF1975b_RESULT_PAIR_13, cache=T-------------------------------------
# Run PSO-QN Algorithm for discriminating \eta_1 and \eta_3
af_res_13 <- DiscrimOD(MODEL_INFO = af_info_13, DISTANCE = sq_diff, 
                       nSupp = 5, dsLower = -1.0, dsUpper = 1.0, 
                       crit_type = "pair_fixed_true",
                       PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO, 
                       seed = NULL, verbose = FALSE)

## ----AF1975b_DESIGN_PAIR_13, cache=T-------------------------------------
round(af_res_13$BESTDESIGN, 3)

## ----AF1975b_VAL_PAIR_13, cache=T----------------------------------------
af_res_13$BESTVAL

## ----AF1975b_EQUV_PAIR_13, cache=T, fig.align='center', fig.height = 4, fig.width = 4, fig.cap='FIGURE. Directional derivative function of T-optimal design in the case \'af_res_13\'.'----
af_eqv_13 <- equivalence(ngrid = 100, PSO_RESULT = af_res_13, 
                         MODEL_INFO = af_info_13, DISTANCE = sq_diff, 
                         dsLower = -1, dsUpper = 1, 
                         crit_type = "pair_fixed_true",
                         PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)
# Draw the directional derivative curve
plot(af_eqv_13$Grid_1, af_eqv_13$DirDeriv, type = "l", col = "blue", 
     main = "af_res_13", xlab = "x", ylab = "Directional Derivative"); abline(h = 0)
points(af_res_13$BESTDESIGN[,1], rep(0, nrow(af_res_13$BESTDESIGN)), pch = 16)

## ----AF1975b_CASE_MAXMIN, cache=T----------------------------------------
af_info_maxmin <- list(af_info_12[[1]], af_info_12[[2]], af_info_13[[2]])
# Define the vector of optimal criterion values for efficiency computations
af_vals_pair <- c(af_res_12$BESTVAL, af_res_13$BESTVAL)

## ----AF1975b_RESULT_MAXMIN, cache=T--------------------------------------
af_res_maxmin <- DiscrimOD(MODEL_INFO = af_info_maxmin, DISTANCE = sq_diff, 
                           nSupp = 5, dsLower = -1, dsUpper = 1, 
                           crit_type = "maxmin_fixed_true", 
                           MaxMinStdVals = af_vals_pair,
                           PSO_INFO = PSO_INFO_MAXMIN, LBFGS_INFO = LBFGS_INFO, 
                           seed = NULL, verbose = FALSE)

## ----AF1975b_DESIGN_MAXMIN, cache=T--------------------------------------
round(af_res_maxmin$BESTDESIGN, 3)

## ----AF1975b_VAL_MAXMIN, cache=T-----------------------------------------
af_res_maxmin$BESTVAL

## ----AF1975b_CPU_MAXMIN, cache=T-----------------------------------------
af_res_maxmin$CPUTIME

## ----AF1975b_EQUV_MAXMIN, cache=T, fig.align='center', fig.height = 4, fig.width = 4, fig.cap='FIGURE. Directional derivative function of max-min T-optimal design in the case af_res_maxmin.'----
af_eqv_maxmin <- equivalence(ngrid = 100, PSO_RESULT = af_res_maxmin, 
                             MODEL_INFO = af_info_maxmin, DISTANCE = sq_diff, 
                             dsLower = -1, dsUpper = 1, 
                             crit_type = "maxmin_fixed_true", 
                             MaxMinStdVals = af_vals_pair, 
                             PSO_INFO = PSO_INFO_MAXMIN, LBFGS_INFO = LBFGS_INFO)
# The weight of efficiency values
af_eqv_maxmin$alpha
# Draw the directional derivative curve
plot(af_eqv_maxmin$Grid_1, af_eqv_maxmin$DirDeriv, type = "l", col = "blue", 
     main = "af_res_maxmin", xlab = "x", ylab = "Directional Derivative"); abline(h = 0)
points(af_res_maxmin$BESTDESIGN[,1], rep(0, nrow(af_res_maxmin$BESTDESIGN)), pch = 16)

## ----MM_MODEL------------------------------------------------------------
# 
enzyme2 <- function(x, p) p[1]*x/(p[2] + x) + p[3]*x  # Modified Michaelis-Menten
enzyme1 <- function(x, p) p[1]*x/(p[2] + x)           # Michaelise-Menten
# Specify the nominal value for the true model, 'linlogi4'.
para_mmm_2 <- c(1, 1, 1)
# Create model list
model_mmm <- list(
  list(model = enzyme2, para = para_mmm_2),
  list(model = enzyme1, paraLower = c(-20, -20), paraUpper = c(20, 20))
)

## ----LOGNORM_B-----------------------------------------------------------
log_norm_B <- function(xt, xr) {
  sigsq <- 1 # nuisance parameter
  var_t <- (exp(sigsq) - 1.0)*(xt^2) # variance (true mdoel)
  var_r <- (exp(sigsq) - 1.0)*(xr^2) # variance (rival model)
  mu_t <- log(xt) - 0.5*log(1.0 + (var_t/(xt^2))) # mean (true model)
  mu_r <- log(xr) - 0.5*log(1.0 + (var_r/(xr^2))) # mean (rival model)
  ((mu_r - mu_t)^2)/(2*sigsq) # KL-divergence
}

## ----MM_RES_LOGNORM_B----------------------------------------------------
# Run for finding max-min KL-optimal design 
mmm_res_lognormB <- DiscrimOD(MODEL_INFO = model_mmm, DISTANCE = log_norm_B, 
                              nSupp = 3, dsLower = 0.1, dsUpper = 5, 
                              crit_type = "pair_fixed_true", 
                              PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO, 
                              seed = NULL, verbose = FALSE)

## ----MM_EQV_LOGNORM_B, cache=T, fig.align='center', fig.height = 4, fig.width = 4, fig.cap='FIGURE. Directional derivative function of KL-optimal design in the case \'mmm_res_lognormB\'.'----
mmm_eqv_lognormB <- equivalence(ngrid = 100, PSO_RESULT = mmm_res_lognormB, 
                                MODEL_INFO = model_mmm, DISTANCE = log_norm_B, 
                                dsLower = 0.1, dsUpper = 5,  
                                crit_type = "pair_fixed_true", 
                                PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)
# Draw the directional derivative curve
plot(mmm_eqv_lognormB$Grid_1, mmm_eqv_lognormB$DirDeriv, type = "l", col = "blue", 
     main = "mmm_res_lognormB", xlab = "x", ylab = "Directional Derivative"); abline(h = 0)
points(mmm_res_lognormB$BESTDESIGN[,1], rep(0, nrow(mmm_res_lognormB$BESTDESIGN)), pch = 16)

## ----MM_DEISGN_LOGNORM_B-------------------------------------------------
round(mmm_res_lognormB$BESTDESIGN, 3) 

## ----MM_VAL_LOGNORM_B----------------------------------------------------
mmm_res_lognormB$BESTVAL 

## ----MM_CPU_LOGNORM_B----------------------------------------------------
mmm_res_lognormB$CPUTIME 

## ----LOGIT_MODEL---------------------------------------------------------
# Tommasi et al. (2016)
linlogi4 <- function(x, p) p[1] + p[2]*x + p[3]*(x^2)  # Logistic 4
linlogi3 <- function(x, p) x*(p[1] + p[2]*x)           #Logistic 3
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
logit_vals_pair <- sapply(1:length(logit_res_pair), function(i) logit_res_pair[[i]]$res$BESTVAL)
# Run for finding max-min KL-optimal design 
logit_res_maxmin <- DiscrimOD(MODEL_INFO = model_logit, DISTANCE = logit_diff, 
                              nSupp = 4, dsLower = 0, dsUpper = 1, 
                              crit_type = "maxmin_fixed_true", 
                              MaxMinStdVals = logit_vals_pair,
                              PSO_INFO = PSO_INFO_MAXMIN, LBFGS_INFO = LBFGS_INFO, 
                              seed = NULL, verbose = FALSE)

## ----LOGIT_EQV_MM, cache=T, fig.align='center', fig.height = 4, fig.width = 4, fig.cap='FIGURE. Directional derivative function of max-min KL-optimal design in the case \'logit_res_maxmin\'.'----
logit_eqv_maxmin <- equivalence(ngrid = 100, PSO_RESULT = logit_res_maxmin, 
                                MODEL_INFO = model_logit, DISTANCE = logit_diff, 
                                dsLower = 0, dsUpper = 1,  
                                crit_type = "maxmin_fixed_true", 
                                MaxMinStdVals = logit_vals_pair, 
                                PSO_INFO = PSO_INFO_MAXMIN, LBFGS_INFO = LBFGS_INFO)
# The weight of efficiency values
logit_eqv_maxmin$alpha
# Draw the directional derivative curve
plot(logit_eqv_maxmin$Grid_1, logit_eqv_maxmin$DirDeriv, type = "l", col = "blue", 
     main = "logit_res_maxmin", xlab = "x", ylab = "Directional Derivative"); abline(h = 0)
points(logit_res_maxmin$BESTDESIGN[,1], rep(0, nrow(logit_res_maxmin$BESTDESIGN)), pch = 16)

## ----LOGIT_DESIGN_MM-----------------------------------------------------
round(logit_res_maxmin$BESTDESIGN, 3) 

## ----LOGIT_VAL_MM--------------------------------------------------------
logit_res_maxmin$BESTVAL 

## ----LOGIT_CPU_MM--------------------------------------------------------
logit_res_maxmin$CPUTIME 

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

