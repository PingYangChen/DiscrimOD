# File Name: demo_Section341.r
# Description:
#   The codes for the implementation of Section 3.4.1, the case of two
#   pharmacokinetic models.  We shown in this case how to use the PSO-QN
#   algorithm to find the KL-optimal designs for different model assumptions.
# Reference:
#   Lopez-Fidalgo, J., Tommasi, C., and Trandafir, P. C. (2007). An optimal
#     experimental design criterion for discriminating between non-normal models.
#     Journal of the Royal Statistical Society: Series B (Statistical Methodology),
#     69(2):231-242.
# ----------------------------------------------------------------------------

# Load packages
library(DiscrimOD); library(Rcpp)

# ----------------------------------------------------------------------------
# 1. Create competing models using C++ codes (faster)
# ----------------------------------------------------------------------------
# Modified Michaelis-Menten
enzyme2 <- cppFunction('
  Rcpp::NumericVector enzyme2(Rcpp::NumericMatrix x, Rcpp::NumericVector p) {
    Rcpp::NumericVector eta(x.nrow());
    for (int i = 0; i < x.nrow(); i++) {
      eta(i) = p(0)*x(i,0)/(p(1) + x(i,0)) + p(2)*x(i,0); }
    return eta;
}')
# Michaelise-Menten
enzyme1 <- cppFunction('
  Rcpp::NumericVector enzyme1(Rcpp::NumericMatrix x, Rcpp::NumericVector p) {
    Rcpp::NumericVector eta(x.nrow());
    for (int i = 0; i < x.nrow(); i++) {
      eta(i) = p(0)*x(i,0)/(p(1) + x(i,0)); }
    return eta;
}')

# Or, create competing models using R codes
#enzyme2 <- function(x, p) p[1]*x/(p[2] + x) + p[3]*x  # Modified Michaelis-Menten
#enzyme1 <- function(x, p) p[1]*x/(p[2] + x)           # Michaelise-Menten


# ----------------------------------------------------------------------------
# 2. Set the discrimination design problem
# ----------------------------------------------------------------------------
# The true model is 'enzyme2' and we set the nominal values (para)
para_enzyme_2 <- c(1, 1, 1)
# Create a model list as the input of our PSO-QN algorithm
model_enzyme <- list(
  # The fitst object should always be the true model
  # Input the model function to the label 'model'
  # Input the nominal values to the labal 'para'
  list(model = enzyme2, para = para_enzyme_2),
  # Starting from the second object, specify the rival models
  # Input the model function to the label 'model'
  # Input the lower and upper bounds for the rival parameters to the labal
  # 'paraLower' and 'paraUpper' (based on experiences)
  list(model = enzyme1, paraLower = c(1e-4, 1e-4), paraUpper = c(20, 20))
)


# ----------------------------------------------------------------------------
# 3. Create distance function using C++ codes (faster)
# ----------------------------------------------------------------------------
# KL divergence for two lognormal distributions
kl_div_lognorm <- cppFunction('
  Rcpp::NumericVector kl_div_lognorm(Rcpp::NumericVector xt, Rcpp::NumericVector xr) {
    double sigsq = 1.0; // Assumed constant
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
# KL divergence for two gamma distributions
kl_div_gamma <- cppFunction('
  Rcpp::NumericVector kl_div_gamma(Rcpp::NumericVector xt, Rcpp::NumericVector xr) {
    Rcpp::NumericVector div(xt.size());
    for (int i = 0; i < xt.size(); i++) {
      div(i) = std::log(xr(i)/xt(i)) + (xt(i) - xr(i))/xr(i);
    }
    return div;
}')

# Or, create distance function using R codes
#kl_div_lognorm <- function(xt, xr) {
#    sigsq <- 1.0 # Assumed constant
#    var_t <- (exp(sigsq) - 1.0)*(xt^2)
#    var_r <- (exp(sigsq) - 1.0)*(xr^2)
#    mu_t <- log(xt) - 0.5*log(1.0 + (var_t/(xt^2)))
#    mu_r <- log(xr) - 0.5*log(1.0 + (var_r/(xr^2)))
#    ((mu_r - mu_t)^2)/(2*sigsq)
#  }
#kl_div_gamma <- function(xt, xr) { log(xr/xt) + (xt - xr)/xr }


# ----------------------------------------------------------------------------
# 4. Set PSO-QN Algorithm Parameters
# ----------------------------------------------------------------------------

# Set PSO options for pariwise discrimination design cases
PSO_INFO <- getPSOInfo(nSwarm = 32, maxIter = 200)
# Set L-BFGS algorithm options
LBFGS_INFO <- getLBFGSInfo(LBFGS_RETRY = 4)


# ----------------------------------------------------------------------------
# 5-1. Run PSO-QN for the lognormal case
# ----------------------------------------------------------------------------

# Run PSO-QN algorithm
res_lognorm <- DiscrimOD(MODEL_INFO = model_enzyme, DISTANCE = kl_div_lognorm, nSupp = 3,
                         dsLower = 0.1, dsUpper = 5, crit_type = "pair_fixed_true",
                         PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO,
                         seed = 100, verbose = TRUE)
# Check the result
names(res_lognorm)
# $BESTDESIGN: the best design found by the PSO-QN algorithm
# $BESTVAL:    the criterion value of the best design
# $GBESTHIST:  the updating path of the cirterion value in PSO loop
# $CPUTIME:    the computing time
res_lognorm$BESTDESIGN
#          dim_1    weight
# obs_1 0.100000 0.2940052
# obs_2 1.569043 0.5000002
# obs_3 5.000000 0.2059946
res_lognorm$BESTVAL
# 0.00256509

# Check by the equivalence theorem
eqv_lognorm <- equivalence(ngrid = 100, PSO_RESULT = res_lognorm,
                           MODEL_INFO = model_enzyme, DISTANCE = kl_div_lognorm,
                           dsLower = 0.1, dsUpper = 5, crit_type = "pair_fixed_true",
                           PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

names(eqv_lognorm)
# $Grid_1:   the grid on the design space
# $Grid_2:   (currently not used)
# $DirDeriv: the values of directional derivative function on the grid
# $MAX_DD:   the maximal value of the directional derivative function
eqv_lognorm$MAX_DD
# 3.011799e-08

# Draw the curve of directional derivative function
plot(eqv_lognorm$Grid_1, eqv_lognorm$DirDeriv, type = "l", col = "black",
     xlab = "x", ylab = "Directional Derivative",
     main = "KL-optimal design for lognormal error")
abline(h = 0, col = "grey50", lty = 2)
points(res_lognorm$BESTDESIGN[,1], rep(0, nrow(res_lognorm$BESTDESIGN)), pch = 16)



# ----------------------------------------------------------------------------
# 5-2. Run PSO-QN for the gamma case
# ----------------------------------------------------------------------------

# Run PSO-QN algorithm
res_gamma <- DiscrimOD(MODEL_INFO = model_enzyme, DISTANCE = kl_div_gamma, nSupp = 3,
                       dsLower = 0.1, dsUpper = 5, crit_type = "pair_fixed_true",
                       PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO,
                       seed = 100, verbose = TRUE)
# Check the result
names(res_gamma)
# $BESTDESIGN: the best design found by the PSO-QN algorithm
# $BESTVAL:    the criterion value of the best design
# $GBESTHIST:  the updating path of the cirterion value in PSO loop
# $CPUTIME:    the computing time
res_gamma$BESTDESIGN
#         dim_1    weight
# obs_1 0.10000 0.2869955
# obs_2 1.56901 0.5119481
# obs_3 5.00000 0.2010564
res_gamma$BESTVAL
# 0.002564359

# Check by the equivalence theorem
eqv_gamma <- equivalence(ngrid = 100, PSO_RESULT = res_gamma,
                         MODEL_INFO = model_enzyme, DISTANCE = kl_div_gamma,
                         dsLower = 0.1, dsUpper = 5, crit_type = "pair_fixed_true",
                         PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

names(eqv_gamma)
# $Grid_1:   the grid on the design space
# $Grid_2:   (currently not used)
# $DirDeriv: the values of directional derivative function on the grid
# $MAX_DD:   the maximal value of the directional derivative function
eqv_gamma$MAX_DD
# 3.011799e-08

# Draw the curve of directional derivative function
plot(eqv_gamma$Grid_1, eqv_gamma$DirDeriv, type = "l", col = "black",
     xlab = "x", ylab = "Directional Derivative",
     main = "KL-optimal design for gamma error")
abline(h = 0, col = "grey50", lty = 2)
points(res_gamma$BESTDESIGN[,1], rep(0, nrow(res_gamma$BESTDESIGN)), pch = 16)








