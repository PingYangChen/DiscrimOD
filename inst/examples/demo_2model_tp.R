# File Name: demo_2model_algorithm.r
# Description:
#   The codes for the implementation of Section 3.4.2, the case of three
#   models with Gaussian errors.  We shown in this case how to use the PSO-QN
#   algorithm to find the T-optimal designs for each pair discrimination.
#   Then, we use the PSO-S-QN algorithm to find the max-min T-optimal design
#   for discriminating among all three models.
# Reference:
#   Atkinson, A. C. and Fedorov, V. V. (1975a). The design of experiments for
#     discriminating between two rival models. Biometrika, 62(1):57-70.
#   Atkinson, A. C. and Fedorov, V. V. (1975b). Optimal design: experiments for
#     discriminating between several models. Biometrika 62 (2), 289–303.
# ----------------------------------------------------------------------------

# Load packages
library(DiscrimOD); library(Rcpp)

# ----------------------------------------------------------------------------
# 1. Create competing models using C++ codes (faster)
# ----------------------------------------------------------------------------
# af1975_1 <- cppFunction('
#   Rcpp::NumericVector af1975_1(Rcpp::NumericMatrix x, Rcpp::NumericVector p) {
#     Rcpp::NumericVector eta(x.nrow());
#     for (int i = 0; i < x.nrow(); i++) {
#       eta(i) = p(0) + p(1)*std::exp(x(i,0)) + p(2)*std::exp(-1.0*x(i,0)); }
#     return eta;
# }')
# af1975_2 <- cppFunction('
#   Rcpp::NumericVector af1975_2(Rcpp::NumericMatrix x, Rcpp::NumericVector p) {
#     Rcpp::NumericVector eta(x.nrow());
#     for (int i = 0; i < x.nrow(); i++) {
#       eta(i) = p(0) + p(1)*x(i,0) + p(2)*x(i,0)*x(i,0); }
#     return eta;
# }')
#
# af1975_disp <- cppFunction('
#   Rcpp::NumericVector af1975_disp(Rcpp::NumericMatrix x, Rcpp::NumericVector p) {
#     Rcpp::NumericVector eta(x.nrow());
#     for (int i = 0; i < x.nrow(); i++) { eta(i) = p(0); }
#     return eta;
# }')

# Or, create competing models using R codes
af1975_1 <- function(x, p) p[1] + p[2]*exp(x) + p[3]*exp(-x)
af1975_2 <- function(x, p) p[1] + p[2]*x + p[3]*x^2
af1975_disp <- function(x, p) rep(p[1], length(x))

# ----------------------------------------------------------------------------
# 2. Set the discrimination design problem
# ----------------------------------------------------------------------------
# The true model is 'af1975_1' and we set the nominal values
mean_af1975_1 <- c(4.5, -1.5, -2)
disp_af1975_1 <- c(1.0)
# Create a model list as the input of our PSO-QN algorithm
model_af1975_12 <- list(
  # The fitst object should always be the true model
  # Input the model function to the label 'model'
  # Input the nominal values to the labal 'para'
  list(mean = af1975_1, disp = af1975_disp, meanPara = mean_af1975_1, dispPara = disp_af1975_1),
  # Starting from the second object, specify the rival models
  # Input the model function to the label 'model'
  # Input the lower and upper bounds for the rival parameters to the labal
  # 'paraLower' and 'paraUpper' (based on experiences)
  list(mean = af1975_2, disp = af1975_disp, meanParaLower = rep(-10, 3), meanParaUpper = rep(10, 3), dispParaLower = c(1), dispParaUpper = c(1))
)


# ----------------------------------------------------------------------------
# 3. Create distance function using C++ codes (faster)
# ----------------------------------------------------------------------------
# Squared difference between two models
# sq_diff <- cppFunction('
#   Rcpp::NumericVector gamma_diff(Rcpp::NumericVector xt, Rcpp::NumericVector xr, Rcpp::NumericVector st, Rcpp::NumericVector sr) {
#     Rcpp::NumericVector div(xt.size()); double diff;
#     for (int i = 0; i < xt.size(); i++) {
#       diff = xt(i) - xr(i); div(i) = diff*diff;
#     }
#     return div;
# }')

# Or, create distance function using R codes
sq_diff <- function(xt, xr, st, sr) (xt - xr)^2


# ----------------------------------------------------------------------------
# 4-1. Run PSO-QN Algorithm
# ----------------------------------------------------------------------------

# Set PSO options for pariwise discrimination design cases
PSO_INFO <- getPSOInfo(nSwarm = 32, maxIter = 100)
# Set L-BFGS algorithm options
LBFGS_INFO <- getLBFGSInfo(LBFGS_RETRY = 4)
# Run PSO-QN algorithm

#MODEL_INFO = model_af1975_12; DISTANCE = sq_diff; nSupp = 4; dsLower = -1; dsUpper = 1; MODEL_PAIR = MODEL_PAIR; WT_PAIR = NULL
#crit_type = "pair_fixed_true"; minWt = 0.0; MaxMinStdVals = NULL
#PSO_INFO = PSO_INFO; LBFGS_INFO = LBFGS_INFO; seed = 10; verbose = TRUE

#MODEL_PAIR <- rbind(c(1, 2), c(2, 1))
MODEL_PAIR <- rbind(c(1, 2))

PSOQN_af1975_12 <- DiscrimOD(MODEL_INFO = model_af1975_12, DISTANCE = sq_diff, nSupp = 4,
                             dsLower = -1, dsUpper = 1, crit_type = "pair_fixed_true",
                             MODEL_PAIR = MODEL_PAIR,
                             PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO,
                             seed = 10, verbose = TRUE)

PSOQN_af1975_12$BESTDESIGN
#           dim_1     weight
#obs_1 -1.0000000 0.25267570
#obs_2 -0.6692933 0.42772012
#obs_3  0.1438411 0.24732430
#obs_4  0.9569755 0.07227988

# ----------------------------------------------------------------------------
# 6. Check equivalecne theorem for each result
# ----------------------------------------------------------------------------

# Check PSO-QN result by the equivalence theorem
eqv_PSOQN <- equivalence(ngrid = 100, PSO_RESULT = PSOQN_af1975_12,
                         MODEL_INFO = model_af1975_12, DISTANCE = sq_diff,
                         dsLower = -1, dsUpper = 1, crit_type = "pair_fixed_true",
                         MODEL_PAIR = MODEL_PAIR,
                         PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

# Draw the curve of directional derivative function
plot(eqv_PSOQN$eqv$Grid_1, eqv_PSOQN$eqv$DirDeriv, type = "l", col = "black",
     xlab = "x", ylab = "Directional Derivative",
     main = "PSO-QN result for case af1975_2 vs. af1975_1")
abline(h = 0, col = "grey50", lty = 2)
points(PSOQN_af1975_12$BESTDESIGN[,1], rep(0, nrow(PSOQN_af1975_12$BESTDESIGN)), pch = 16)
