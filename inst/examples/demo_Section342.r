# File Name: demo_Section342.r
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
af1975_1 <- cppFunction('
  Rcpp::NumericVector af1975_1(Rcpp::NumericMatrix x, Rcpp::NumericVector p) {
    Rcpp::NumericVector eta(x.nrow());
    for (int i = 0; i < x.nrow(); i++) {
      eta(i) = p(0) + p(1)*std::exp(x(i,0)) + p(2)*std::exp(-1.0*x(i,0)); }
    return eta;
}')
af1975_2 <- cppFunction('
  Rcpp::NumericVector af1975_2(Rcpp::NumericMatrix x, Rcpp::NumericVector p) {
    Rcpp::NumericVector eta(x.nrow());
    for (int i = 0; i < x.nrow(); i++) {
      eta(i) = p(0) + p(1)*x(i,0) + p(2)*x(i,0)*x(i,0); }
    return eta;
}')
af1975_3 <- cppFunction('
  Rcpp::NumericVector af1975_3(Rcpp::NumericMatrix x, Rcpp::NumericVector p) {
    double pi_val = 3.141592653589793;
    Rcpp::NumericVector eta(x.nrow());
    for (int i = 0; i < x.nrow(); i++) {
      eta(i) = p(0) + p(1)*std::sin(0.5*pi_val*x(i,0)) + p(2)*std::cos(0.5*pi_val*x(i,0)) + p(3)*std::sin(pi_val*x(i,0)); }
    return eta;
}')

# Or, create competing models using R codes
#af1975_1 <- function(x, p) p[1] + p[2]*exp(x) + p[3]*exp(-x)
#af1975_2 <- function(x, p) p[1] + p[2]*x + p[3]*x^2
#af1975_3 <- function(x, p) p[1] + p[2]*sin(0.5*pi*x) + p[3]*cos(0.5*pi*x) + p[4]*sin(pi*x)


# ----------------------------------------------------------------------------
# 2. Set the discrimination design problem
# ----------------------------------------------------------------------------
# The true model is 'af1975_1' and we set the nominal values
para_af1975_1 <- c(4.5, -1.5, -2)
# Create a model list as the input of our PSO-QN algorithm
model_af1975 <- list(
  # The fitst object should always be the true model
  # Input the model function to the label 'model'
  # Input the nominal values to the labal 'para'
  list(model = af1975_1, para = para_af1975_1),
  # Starting from the second object, specify the rival models
  # Input the model function to the label 'model'
  # Input the lower and upper bounds for the rival parameters to the labal
  # 'paraLower' and 'paraUpper' (based on experiences)
  list(model = af1975_2, paraLower = rep(-10, 3), paraUpper = rep(10, 3)),
  # The third model
  list(model = af1975_3, paraLower = rep(-10, 4), paraUpper = rep(10, 4))
)


# ----------------------------------------------------------------------------
# 3. Create distance function using C++ codes (faster)
# ----------------------------------------------------------------------------
# Squared difference between two models
sq_diff <- cppFunction('
  Rcpp::NumericVector gamma_diff(Rcpp::NumericVector xt, Rcpp::NumericVector xr) {
    Rcpp::NumericVector div(xt.size()); double diff;
    for (int i = 0; i < xt.size(); i++) {
      diff = xt(i) - xr(i); div(i) = diff*diff;
    }
    return div;
}')

# Or, create distance function using R codes
#sq_diff <- function(xt, xr) (xt - xr)^2


# ----------------------------------------------------------------------------
# 4. Set PSO-QN Algorithm Parameters
# ----------------------------------------------------------------------------

# Set PSO options for pariwise discrimination design cases
PSO_INFO <- getPSOInfo(nSwarm = 32, maxIter = 200)
# Set L-BFGS algorithm options
LBFGS_INFO <- getLBFGSInfo(LBFGS_RETRY = 4)


# ----------------------------------------------------------------------------
# 5-1. Run PSO-QN for the case: af1975_2 vs. af1975_1
# ----------------------------------------------------------------------------

# Get list for two models
two_model_af1975_12 <- list(model_af1975[[1]], model_af1975[[2]])
# Run PSO-QN algorithm
res_af1975_12 <- DiscrimOD(MODEL_INFO = two_model_af1975_12, DISTANCE = sq_diff, nSupp = 4,
                           dsLower = -1, dsUpper = 1, crit_type = "pair_fixed_true",
                           PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO,
                           seed = 10, verbose = TRUE)
# Check the result
names(res_af1975_12)
# $BESTDESIGN: the best design found by the PSO-QN algorithm
# $BESTVAL:    the criterion value of the best design
# $GBESTHIST:  the updating path of the cirterion value in PSO loop
# $CPUTIME:    the computing time
res_af1975_12$BESTDESIGN
#           dim_1     weight
#obs_1 -1.0000000 0.25267570
#obs_2 -0.6692933 0.42772012
#obs_3  0.1438411 0.24732430
#obs_4  0.9569755 0.07227988
res_af1975_12$BESTVAL
# 0.001086725

# Check by the equivalence theorem
eqv_af1975_12 <- equivalence(ngrid = 100, PSO_RESULT = res_af1975_12,
                             MODEL_INFO = two_model_af1975_12, DISTANCE = sq_diff,
                             dsLower = -1, dsUpper = 1, crit_type = "pair_fixed_true",
                             PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

names(eqv_af1975_12)
# $Grid_1:   the grid on the design space
# $Grid_2:   (currently not used)
# $DirDeriv: the values of directional derivative function on the grid
# $MAX_DD:   the maximal value of the directional derivative function
eqv_af1975_12$MAX_DD
# -1.329804e-08

# Draw the curve of directional derivative function
plot(eqv_af1975_12$Grid_1, eqv_af1975_12$DirDeriv, type = "l", col = "black",
     xlab = "x", ylab = "Directional Derivative",
     main = "T-optimal design for case af1975_2 vs. af1975_1")
abline(h = 0, col = "grey50", lty = 2)
points(res_af1975_12$BESTDESIGN[,1], rep(0, nrow(res_af1975_12$BESTDESIGN)), pch = 16)


# ----------------------------------------------------------------------------
# 5-2. Run PSO-QN for the case: af1975_3 vs. af1975_1
# ----------------------------------------------------------------------------

# Get list for two models
two_model_af1975_13 <- list(model_af1975[[1]], model_af1975[[3]])
# Run PSO-QN algorithm
res_af1975_13 <- DiscrimOD(MODEL_INFO = two_model_af1975_13, DISTANCE = sq_diff, nSupp = 5,
                           dsLower = -1, dsUpper = 1, crit_type = "pair_fixed_true",
                           PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO,
                           seed = 10, verbose = TRUE)
# Check the result
names(res_af1975_13)
# $BESTDESIGN: the best design found by the PSO-QN algorithm
# $BESTVAL:    the criterion value of the best design
# $GBESTHIST:  the updating path of the cirterion value in PSO loop
# $CPUTIME:    the computing time
res_af1975_13$BESTDESIGN
#           dim_1     weight
#obs_1 -1.0000000 0.19159011
#obs_2 -0.7405088 0.32281553
#obs_3 -0.1044060 0.22737157
#obs_4  0.6339742 0.17718308
#obs_5  1.0000000 0.08103971
res_af1975_13$BESTVAL
# 0.005715191

# Check by the equivalence theorem
eqv_af1975_13 <- equivalence(ngrid = 100, PSO_RESULT = res_af1975_13,
                             MODEL_INFO = two_model_af1975_13, DISTANCE = sq_diff,
                             dsLower = -1, dsUpper = 1, crit_type = "pair_fixed_true",
                             PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

names(eqv_af1975_13)
# $Grid_1:   the grid on the design space
# $Grid_2:   (currently not used)
# $DirDeriv: the values of directional derivative function on the grid
# $MAX_DD:   the maximal value of the directional derivative function
eqv_af1975_13$MAX_DD
# -3.039192e-08

# Draw the curve of directional derivative function
plot(eqv_af1975_13$Grid_1, eqv_af1975_13$DirDeriv, type = "l", col = "black",
     xlab = "x", ylab = "Directional Derivative",
     main = "T-optimal design for case af1975_3 vs. af1975_1")
abline(h = 0, col = "grey50", lty = 2)
points(res_af1975_13$BESTDESIGN[,1], rep(0, nrow(res_af1975_13$BESTDESIGN)), pch = 16)


# ----------------------------------------------------------------------------
# 6. Run PSO-S-QN for max-min T-optimal design for all models
# ----------------------------------------------------------------------------

# Set PSO options for pariwise discrimination design cases
PSO_MM_INFO <- getPSOInfo(nSwarm = 32, maxIter = 400)

# Get the values of T-optimal criterion found previouly which will be
# the denominators of the T-efficiencies
eff_denom <- c(res_af1975_12$BESTVAL, res_af1975_13$BESTVAL)

res_af1975 <- DiscrimOD(MODEL_INFO = model_af1975, DISTANCE = sq_diff, nSupp = 5,
                        dsLower = -1, dsUpper = 1, crit_type = "maxmin_fixed_true",
                        MaxMinStdVals = eff_denom,
                        PSO_INFO = PSO_MM_INFO, LBFGS_INFO = LBFGS_INFO,
                        seed = 10, verbose = TRUE)
# Check the result
names(res_af1975)
# $BESTDESIGN: the best design found by the PSO-QN algorithm
# $BESTVAL:    the criterion value of the best design
# $GBESTHIST:  the updating path of the cirterion value in PSO loop
# $CPUTIME:    the computing time
res_af1975$BESTDESIGN
#            dim_1     weight
#obs_1 -1.00000000 0.22829853
#obs_2 -0.70494003 0.38283585
#obs_3 -0.02059612 0.21781345
#obs_4  0.56840642 0.11185337
#obs_5  1.00000000 0.05919881
res_af1975$BESTVAL
# 0.8061735

eqv_af1975 <- equivalence(ngrid = 100, PSO_RESULT = res_af1975,
                          MODEL_INFO = model_af1975, DISTANCE = sq_diff,
                          dsLower = -1, dsUpper = 1, crit_type = "maxmin_fixed_true",
                          MaxMinStdVals = eff_denom,
                          PSO_INFO = PSO_MM_INFO, LBFGS_INFO = LBFGS_INFO)

names(eqv_af1975)
# $Grid_1:   the grid on the design space
# $Grid_2:   (currently not used)
# $DirDeriv: the values of directional derivative function on the grid
# $MAX_DD:   the maximal value of the directional derivative function
# $alpha:    the best weight vector of the efficiency values
eqv_af1975$MAX_DD
# 0.005839902

# Draw the curve of directional derivative function
plot(eqv_af1975$Grid_1, eqv_af1975$DirDeriv, type = "l", col = "black",
     xlab = "x", ylab = "Directional Derivative",
     main = "max-min T-optimal design for af1975 example")
abline(h = 0, col = "grey50", lty = 2)
points(res_af1975$BESTDESIGN[,1], rep(0, nrow(res_af1975$BESTDESIGN)), pch = 16)

