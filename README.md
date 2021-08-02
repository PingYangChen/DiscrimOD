
DiscrimOD
=================
Using Hybridized PSO algorithm to search for the optimal approximate discrimination design

* Chen, R. B., Chen, P. Y., Hsu, C. L., and Wong, W. K. (2020). Hybrid algorithms for generating optimal designs for discriminating multiple nonlinear models under various error distributional assumptions. PloS one, 15(10), e0239864.

Installation
------------
Please install the latest development version from github with

``` r
install.packages("devtools")
devtools::install_github("PingYangChen/DiscrimOD")
```

If you encounter a bug, please file a reproducible example on [github](https://github.com/PingYangChen/DiscrimOD/issues).

Examples
--------
``` r
# Atkinson and Fedorov (1975a): T-optimal
# Two R functions of competing models are given by
af1 <- function(x, p) p[1] + p[2]*exp(x) + p[3]*exp(-x)
af2 <- function(x, p) p[1] + p[2]*x + p[3]*x^2

# Set the model information
# The nominla value in 'm1' is 4.5, -1.5, -2.0
# For 'af2', we set the parameter space to be [-10, 10]^3 and
# the initial guess (for LBFGS) of the rival model parameter is zero vector
AF_para_af1 <- c(4.5, -1.5, -2)
af_info_12 <- list(
  # The first list should be the true model and the specified nominal values
  list(model = af1, para = AF_para_af1),
  # Then the rival models are listed accordingly. We also need to specify the model space.
  list(model = af2, paraLower = rep(-10, 3), paraUpper = rep(10, 3))
)
# Define the R function for the distance measure in T-optimal criterion
# xt is the mean values of the true model
# xr is the mean values of the rival model
sq_diff <- function(xt, xr) (xt - xr)^2

# Initialize PSO and BFGS options
PSO_INFO <- getPSOInfo(nSwarm = 32, maxIter = 100)
LBFGS_INFO <- getLBFGSInfo(LBFGS_RETRY = 2)

# Find T-optimal design for models af1 and af2
af_res_12 <- DiscrimOD(MODEL_INFO = af_info_12, DISTANCE = sq_diff,
  nSupp = 4, dsLower = -1, dsUpper = 1, crit_type = "pair_fixed_true",
  PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO, seed = NULL, verbose = FALSE)

round(af_res_12$BESTDESIGN, 3) # The resulting design
af_res_12$BESTVAL # The T-optimal criterion value
af_res_12$CPUTIME # CPU time

# Test optimality by equivalence theorem
af_eqv_12 <- equivalence(ngrid = 100, PSO_RESULT = af_res_12, MODEL_INFO = af_info_12,
  DISTANCE = sq_diff, dsLower = -1, dsUpper = 1, crit_type = "pair_fixed_true",
  PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

# Draw the directional derivative curve
plot(af_eqv_12$Grid_1, af_eqv_12$DirDeriv, type = "l", col = "blue",
  main = "af_res_12", xlab = "x", ylab = "Directional Derivative"); abline(h = 0)
points(af_res_12$BESTDESIGN[,1], rep(0, nrow(af_res_12$BESTDESIGN)), pch = 16)

# Following above, we add the 3rd model in Atkinson and Fedorov (1975b)
af3 <- function(x, p) p[1] + p[2]*sin(0.5*pi*x) + p[3]*cos(0.5*pi*x) + p[4]*sin(pi*x)
af_info_13 <- list(
  list(model = af1, para = AF_para_af1),
  list(model = af3, paraLower = rep(-10, 4), paraUpper = rep(10, 4))
)
# Find another T-optimal design for models af1 and af3
af_res_13 <- DiscrimOD(MODEL_INFO = af_info_13, DISTANCE = sq_diff,
  nSupp = 5, dsLower = -1.0, dsUpper = 1.0, crit_type = "pair_fixed_true",
  PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO, seed = NULL, verbose = FALSE)

# Re-organize model list for finding the max-min T-optimal design
af_info_maxmin <- list(af_info_12[[1]], af_info_12[[2]], af_info_13[[2]])
# Define the vector of optimal criterion values for efficiency computations
af_vals_pair <- c(af_res_12$BESTVAL, af_res_13$BESTVAL)

# Search for max-min T-optimal design for discriminating af1, af2 and af3
af_res_maxmin <- DiscrimOD(MODEL_INFO = af_info_maxmin, DISTANCE = sq_diff,
  nSupp = 5, dsLower = -1, dsUpper = 1, crit_type = "maxmin_fixed_true",
  MaxMinStdVals = af_vals_pair, PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO,
  seed = NULL, verbose = FALSE)

round(af_res_maxmin$BESTDESIGN, 3) # The resulting design
af_res_maxmin$BESTVAL # The T-optimal criterion value
af_res_maxmin$CPUTIME # CPU time

# Test optimality by equivalence theorem
af_eqv_maxmin <- equivalence(ngrid = 100, PSO_RESULT = af_res_maxmin,
  MODEL_INFO = af_info_maxmin, DISTANCE = sq_diff, dsLower = -1, dsUpper = 1,
  crit_type = "maxmin_fixed_true", MaxMinStdVals = af_vals_pair,
  PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

af_eqv_maxmin$alpha # The weight of efficiency values

# Draw the directional derivative curve
plot(af_eqv_maxmin$Grid_1, af_eqv_maxmin$DirDeriv, type = "l", col = "blue",
  main = "af_res_maxmin", xlab = "x", ylab = "Directional Derivative"); abline(h = 0)
points(af_res_maxmin$BESTDESIGN[,1], rep(0, nrow(af_res_maxmin$BESTDESIGN)), pch = 16)

# Not Run
# For other distance measure, here are a few

# Heteroscedastic Normal model
# heter_norm <- function(xt, xr) {
#   var_t <- xt^2
#   var_r <- xr^2
#   (var_t + (xt - xr)^2)/var_r - log(var_t/var_r)
# }

# Logistic regression model
# logit_diff <- function(xt, xr) {
#   exp_t <- exp(xt)
#   exp_r <- exp(xr)
#   mu_t <- exp_t/(1 + exp_t)
#   mu_r <- exp_r/(1 + exp_r)
#   mu_t*(log(mu_t) - log(mu_r)) + (1 - mu_t)*(log(1.0 - mu_t) - log(1.0 - mu_r))
# }

# Gamma regression model
# gamma_diff <- function(xt, xr) log(xr/xt) + (xt - xr)/xr
```
