# File Name: demo_Section4_T.r
# Description:
#   The codes for the implementation of Section 4, the case of five toxicology
#   models with lognromal errors.  We shown in this case how to use the PSO-QN
#   algorithm to find the KL-optimal designs for each pair discrimination.
#   Then, we use the PSO-S-QN algorithm to find the max-min KL-optimal design
#   for discriminating among all five models.
# Reference:

# ----------------------------------------------------------------------------

# Load packages
library(DiscrimOD); library(Rcpp)

# ----------------------------------------------------------------------------
# 1. Create competing models using C++ codes (faster)
# ----------------------------------------------------------------------------
tox5 <- cppFunction('
  Rcpp::NumericVector tox5(Rcpp::NumericMatrix x, Rcpp::NumericVector p) {
    Rcpp::NumericVector eta(x.nrow());
    for (int i = 0; i < x.nrow(); i++) {
      eta(i) = p(0)*(p(2) - (p(2) - 1.0)*std::exp((-1.0)*std::pow(x(i,0)/p(1), p(3)))); }
    return eta;
}')
tox4 <- cppFunction('
  Rcpp::NumericVector tox4(Rcpp::NumericMatrix x, Rcpp::NumericVector p) {
    Rcpp::NumericVector eta(x.nrow());
    for (int i = 0; i < x.nrow(); i++) {
      eta(i) = p(0)*(p(2) - (p(2) - 1.0)*std::exp((-1.0)*x(i,0)/p(1))); }
    return eta;
}')
tox3 <- cppFunction('
  Rcpp::NumericVector tox3(Rcpp::NumericMatrix x, Rcpp::NumericVector p) {
    Rcpp::NumericVector eta(x.nrow());
    for (int i = 0; i < x.nrow(); i++) {
      eta(i) = p(0)*std::exp((-1.0)*std::pow(x(i,0)/p(1), p(2))); }
    return eta;
}')
tox2 <- cppFunction('
  Rcpp::NumericVector tox2(Rcpp::NumericMatrix x, Rcpp::NumericVector p) {
    Rcpp::NumericVector eta(x.nrow());
    for (int i = 0; i < x.nrow(); i++) {
      eta(i) = p(0)*std::exp((-1.0)*x(i,0)/p(1)); }
    return eta;
}')
tox1 <- cppFunction('
  Rcpp::NumericVector tox1(Rcpp::NumericMatrix x, Rcpp::NumericVector p) {
    Rcpp::NumericVector eta(x.nrow());
    for (int i = 0; i < x.nrow(); i++) { eta(i) = p(0); }
    return eta;
}')


# Or, create competing models using R codes
#tox5 <- function(x, p) p[1]*(p[3] - (p[3] - 1)*exp(-(x/p[2])^p[4]))
#tox4 <- function(x, p) p[1]*(p[3] - (p[3] - 1)*exp(-(x/p[2])))
#tox3 <- function(x, p) p[1]*exp(-(x/p[2])^p[3])
#tox2 <- function(x, p) p[1]*exp(-(x/p[2]))
#tox1 <- function(x, p) rep(p[1], length(x))


# ----------------------------------------------------------------------------
# 2. Set the discrimination design problem
# ----------------------------------------------------------------------------
# The true model is 'linlogi4' and we set the nominal values
para_tox_5 <- c(4.282, 835.571, 0.739, 3.515)
# Create a model list as the input of our PSO-QN algorithm
model_tox <- list(
  # The fitst object should always be the true model
  # Input the model function to the label 'model'
  # Input the nominal values to the labal 'para'
  list(model = tox5, para = para_tox_5),
  # Starting from the second object, specify the rival models
  # Input the model function to the label 'model'
  # Input the lower and upper bounds for the rival parameters to the labal
  # 'paraLower' and 'paraUpper' (based on experiences)
  list(model = tox4, paraLower = c(0, 0, 0), paraUpper = c(20, 5000, 1)),
  # The third model
  list(model = tox3, paraLower = c(0, 0, 1), paraUpper = c(20, 5000, 15)),
  # The fourth model
  list(model = tox2, paraLower = c(0, 0), paraUpper = c(20, 5000)),
  # The fifth model
  list(model = tox1, paraLower = 0, paraUpper = 20)
)


# ----------------------------------------------------------------------------
# 3. Create distance function using C++ codes (faster)
# ----------------------------------------------------------------------------
# KL divergence between two models with lognorm errors
kl_div_lognorm <- cppFunction('
  Rcpp::NumericVector kl_div_lognorm(Rcpp::NumericVector xt, Rcpp::NumericVector xr) {
    double sigsq = 0.01; // Assumed constant
    Rcpp::NumericVector div(xt.size());
    double xt2, xr2, vt, vr, mt, mr;
    for (int i = 0; i < xt.size(); i++) {
      xt2 = xt(i)*xt(i);
      xr2 = xr(i)*xr(i);
      vt = (std::exp(sigsq) - 1.0)*xt2;
      vr = (std::exp(sigsq) - 1.0)*xr2;
      mt = std::log(xt(i)) - 0.5*std::log(1.0 + (vt/xt2));
      mr = std::log(xr(i)) - 0.5*std::log(1.0 + (vr/xr2));
      div(i) = 0.5*((mt - mr)*(mt - mr))/sigsq;
    }
    return div;
}')

# Or, create distance function using R codes
#kl_div_lognorm <- function(xt, xr) {
#    sigsq <- 0.01 # Assumed constant
#    var_t <- (exp(sigsq) - 1.0)*(xt^2)
#    var_r <- (exp(sigsq) - 1.0)*(xr^2)
#    mu_t <- log(xt) - 0.5*log(1.0 + (var_t/(xt^2)))
#    mu_r <- log(xr) - 0.5*log(1.0 + (var_r/(xr^2)))
#    ((mu_r - mu_t)^2)/(2*sigsq)
#  }


# ----------------------------------------------------------------------------
# 4. Set PSO-QN Algorithm Parameters
# ----------------------------------------------------------------------------

# Set PSO options for pariwise discrimination design cases
PSO_INFO <- getPSOInfo(nSwarm = 32, maxIter = 200)
# Set L-BFGS algorithm options
LBFGS_INFO <- getLBFGSInfo(LBFGS_RETRY = 4)


# ----------------------------------------------------------------------------
# 5-1. Run PSO-QN for the case: tox4 vs. tox5
# ----------------------------------------------------------------------------

# Get list for two models
two_model_tox_54 <- list(model_tox[[1]], model_tox[[2]])
# Run PSO-QN algorithm
res_tox_54 <- DiscrimOD(MODEL_INFO = two_model_tox_54, DISTANCE = kl_div_lognorm, nSupp = 3,
                        dsLower = 0, dsUpper = 1250, crit_type = "pair_fixed_true",
                        PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO,
                        seed = 10, verbose = TRUE)
# Check the result
names(res_tox_54)
# $BESTDESIGN: the best design found by the PSO-QN algorithm
# $BESTVAL:    the criterion value of the best design
# $GBESTHIST:  the updating path of the cirterion value in PSO loop
# $CPUTIME:    the computing time
res_tox_54$BESTDESIGN
#         dim_1    weight
#obs_1    0.000 0.2712305
#obs_2  487.448 0.5000000
#obs_3 1065.369 0.2287695
res_tox_54$BESTVAL
# 0.09263272

# Check by the equivalence theorem
eqv_tox_54 <- equivalence(ngrid = 100, PSO_RESULT = res_tox_54,
                          MODEL_INFO = two_model_tox_54, DISTANCE = kl_div_lognorm,
                          dsLower = 0, dsUpper = 1250, crit_type = "pair_fixed_true",
                          PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

names(eqv_tox_54)
# $Grid_1:   the grid on the design space
# $Grid_2:   (currently not used)
# $DirDeriv: the values of directional derivative function on the grid
# $MAX_DD:   the maximal value of the directional derivative function
eqv_tox_54$MAX_DD
# -1.610767e-05

# Draw the curve of directional derivative function
plot(eqv_tox_54$Grid_1, eqv_tox_54$DirDeriv, type = "l", col = "black",
     xlab = "x", ylab = "Directional Derivative",
     main = "KL-optimal design for case tox_4 vs. tox_5")
abline(h = 0, col = "grey50", lty = 2)
points(res_tox_54$BESTDESIGN[,1], rep(0, nrow(res_tox_54$BESTDESIGN)), pch = 16)


# ----------------------------------------------------------------------------
# 5-2. Run PSO-QN for the case: tox_3 vs. tox_5
# ----------------------------------------------------------------------------

# Get list for two models
two_model_tox_53 <- list(model_tox[[1]], model_tox[[3]])
# Run PSO-QN algorithm
res_tox_53 <- DiscrimOD(MODEL_INFO = two_model_tox_53, DISTANCE = kl_div_lognorm, nSupp = 4,
                        dsLower = 0, dsUpper = 1250, crit_type = "pair_fixed_true",
                        PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO,
                        seed = 10, verbose = TRUE)

res_tox_53$BESTDESIGN
#          dim_1     weight
#obs_1    0.0000 0.09250735
#obs_2  498.8916 0.29044494
#obs_3  979.7211 0.40749728
#obs_4 1250.0000 0.20955043
res_tox_53$BESTVAL
# 0.0316106

# Check by the equivalence theorem
eqv_tox_53 <- equivalence(ngrid = 100, PSO_RESULT = res_tox_53,
                          MODEL_INFO = two_model_tox_53, DISTANCE = kl_div_lognorm,
                          dsLower = 0, dsUpper = 1250, crit_type = "pair_fixed_true",
                          PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

eqv_tox_53$MAX_DD
# -1.044966e-06

# Draw the curve of directional derivative function
plot(eqv_tox_53$Grid_1, eqv_tox_53$DirDeriv, type = "l", col = "black",
     xlab = "x", ylab = "Directional Derivative",
     main = "KL-optimal design for case tox_3 vs. tox_5")
abline(h = 0, col = "grey50", lty = 2)
points(res_tox_53$BESTDESIGN[,1], rep(0, nrow(res_tox_53$BESTDESIGN)), pch = 16)


# ----------------------------------------------------------------------------
# 5-3. Run PSO-QN for the case: tox_2 vs. tox_5
# ----------------------------------------------------------------------------

# Get list for two models
two_model_tox_52 <- list(model_tox[[1]], model_tox[[4]])
# Run PSO-QN algorithm
res_tox_52 <- DiscrimOD(MODEL_INFO = two_model_tox_52, DISTANCE = kl_div_lognorm, nSupp = 3,
                        dsLower = 0, dsUpper = 1250, crit_type = "pair_fixed_true",
                        PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO,
                        seed = 1, verbose = TRUE)

res_tox_52$BESTDESIGN
#          dim_1    weight
#obs_1    0.0000 0.2712298
#obs_2  487.4482 0.5000010
#obs_3 1065.3698 0.2287692
res_tox_52$BESTVAL
# 0.09263269

# Check by the equivalence theorem
eqv_tox_52 <- equivalence(ngrid = 100, PSO_RESULT = res_tox_52,
                          MODEL_INFO = two_model_tox_52, DISTANCE = kl_div_lognorm,
                          dsLower = 0, dsUpper = 1250, crit_type = "pair_fixed_true",
                          PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

eqv_tox_52$MAX_DD
# -4.949537e-07

# Draw the curve of directional derivative function
plot(eqv_tox_52$Grid_1, eqv_tox_52$DirDeriv, type = "l", col = "black",
     xlab = "x", ylab = "Directional Derivative",
     main = "KL-optimal design for case tox_2 vs. tox_5")
abline(h = 0, col = "grey50", lty = 2)
points(res_tox_52$BESTDESIGN[,1], rep(0, nrow(res_tox_52$BESTDESIGN)), pch = 16)


# ----------------------------------------------------------------------------
# 5-4. Run PSO-QN for the case: tox_1 vs. tox_5
# ----------------------------------------------------------------------------

# Get list for two models
two_model_tox_51 <- list(model_tox[[1]], model_tox[[5]])
# Run PSO-QN algorithm
res_tox_51 <- DiscrimOD(MODEL_INFO = two_model_tox_51, DISTANCE = kl_div_lognorm, nSupp = 2,
                        dsLower = 0, dsUpper = 1250, crit_type = "pair_fixed_true",
                        PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO,
                        seed = 10, verbose = TRUE)

res_tox_51$BESTDESIGN
#      dim_1 weight
#obs_1     0    0.5
#obs_2  1250    0.5
res_tox_51$BESTVAL
# 1.100645

# Check by the equivalence theorem
eqv_tox_51 <- equivalence(ngrid = 100, PSO_RESULT = res_tox_51,
                          MODEL_INFO = two_model_tox_51, DISTANCE = kl_div_lognorm,
                          dsLower = 0, dsUpper = 1250, crit_type = "pair_fixed_true",
                          PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

eqv_tox_51$MAX_DD
# 1.593241e-07

# Draw the curve of directional derivative function
plot(eqv_tox_51$Grid_1, eqv_tox_51$DirDeriv, type = "l", col = "black",
     xlab = "x", ylab = "Directional Derivative",
     main = "KL-optimal design for case tox_1 vs. tox_5")
abline(h = 0, col = "grey50", lty = 2)
points(res_tox_51$BESTDESIGN[,1], rep(0, nrow(res_tox_51$BESTDESIGN)), pch = 16)

# ----------------------------------------------------------------------------
# 6. Run PSO-S-QN for max-min KL-optimal design for all models
# ----------------------------------------------------------------------------

# Set PSO options for pariwise discrimination design cases
PSO_MM_INFO <- getPSOInfo(nSwarm = 32, maxIter = 400)

# Get the values of KL-optimal criterion found previouly which will be
# the denominators of the T-efficiencies
eff_denom <- c(res_tox_54$BESTVAL, res_tox_53$BESTVAL, res_tox_52$BESTVAL, res_tox_51$BESTVAL)

res_tox <- DiscrimOD(MODEL_INFO = model_tox, DISTANCE = kl_div_lognorm, nSupp = 4,
                     dsLower = 0, dsUpper = 1250, crit_type = "maxmin_fixed_true",
                     MaxMinStdVals = eff_denom,
                     PSO_INFO = PSO_MM_INFO, LBFGS_INFO = LBFGS_INFO,
                     seed = 10, verbose = TRUE)
# Check the result
names(res_tox)
# $BESTDESIGN: the best design found by the PSO-QN algorithm
# $BESTVAL:    the criterion value of the best design
# $GBESTHIST:  the updating path of the cirterion value in PSO loop
# $CPUTIME:    the computing time
res_tox$BESTDESIGN
#          dim_1    weight
#obs_1    0.0000 0.2225314
#obs_2  453.3533 0.3395145
#obs_3 1044.5112 0.2506704
#obs_4 1250.0000 0.1872837
res_tox$BESTVAL
# 0.7678414

eqv_tox <- equivalence(ngrid = 100, PSO_RESULT = res_tox,
                       MODEL_INFO = model_tox, DISTANCE = kl_div_lognorm,
                       dsLower = 0, dsUpper = 1250, crit_type = "maxmin_fixed_true",
                       MaxMinStdVals = eff_denom,
                       PSO_INFO = PSO_MM_INFO, LBFGS_INFO = LBFGS_INFO)

names(eqv_tox)
# $Grid_1:   the grid on the design space
# $Grid_2:   (currently not used)
# $DirDeriv: the values of directional derivative function on the grid
# $MAX_DD:   the maximal value of the directional derivative function
# $alpha:    the best weight vector of the efficiency values
eqv_tox$MAX_DD
# 0.002781319

# Draw the curve of directional derivative function
plot(eqv_tox$Grid_1, eqv_tox$DirDeriv, type = "l", col = "black",
     xlab = "x", ylab = "Directional Derivative",
     main = "max-min KL-optimal design for toxicology example")
abline(h = 0, col = "grey50", lty = 2)
points(res_tox$BESTDESIGN[,1], rep(0, nrow(res_tox$BESTDESIGN)), pch = 16)

