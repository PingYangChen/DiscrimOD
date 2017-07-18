# File Name: demo_Section343.r
# Description:
#   The codes for the implementation of Section 3.4.3, the case of four
#   logistic regression models.  We shown in this case how to use the PSO-QN
#   algorithm to find the KL-optimal designs for each pair discrimination.
#   Then, we use the PSO-S-QN algorithm to find the max-min KL-optimal design
#   for discriminating among all four models.
# Reference:
#   Tommasi, C., Martin-Martin, R., and Lopez-Fidalgo, J. (2016). Max-min
#     optimal discriminating designs for several statistical models.
#     Statistics and Computing, 26(6):1163-1172.
# ----------------------------------------------------------------------------

# Load packages
library(DiscrimOD); library(Rcpp)

# ----------------------------------------------------------------------------
# 1. Create competing models using C++ codes (faster)
# ----------------------------------------------------------------------------
logit_4 <- cppFunction('
  Rcpp::NumericVector logit_4(Rcpp::NumericMatrix x, Rcpp::NumericVector p) {
    Rcpp::NumericVector eta(x.nrow());
    for (int i = 0; i < x.nrow(); i++) {
      eta(i) = p(0) + p(1)*x(i,0) + p(2)*x(i,0)*x(i,0); }
    return eta;
}')
logit_3 <- cppFunction('
  Rcpp::NumericVector logit_3(Rcpp::NumericMatrix x, Rcpp::NumericVector p) {
    Rcpp::NumericVector eta(x.nrow());
    for (int i = 0; i < x.nrow(); i++) { eta(i) = x(i,0)*(p(0) + p(1)*x(i,0)); }
    return eta;
}')
logit_2 <- cppFunction('
  Rcpp::NumericVector logit_2(Rcpp::NumericMatrix x, Rcpp::NumericVector p) {
    Rcpp::NumericVector eta(x.nrow());
    for (int i = 0; i < x.nrow(); i++) { eta(i) = p(0) + p(1)*x(i,0); }
    return eta;
}')
logit_1 <- cppFunction('
  Rcpp::NumericVector logit_1(Rcpp::NumericMatrix x, Rcpp::NumericVector p) {
    Rcpp::NumericVector eta(x.nrow());
    for (int i = 0; i < x.nrow(); i++) { eta(i) = p(0)*x(i,0); }
    return eta;
}')

# Or, create competing models using R codes
#logit_4 <- function(x, p) p[1] + p[2]*x + p[3]*(x^2)  # Tommasi et al. (2016) Logistic 4
#logit_3 <- function(x, p) x*(p[1] + p[2]*x)           # Tommasi et al. (2016) Logistic 3
#logit_2 <- function(x, p) p[1] + p[2]*x               # Tommasi et al. (2016) Logistic 2
#logit_1 <- function(x, p) p[1]*x                      # Tommasi et al. (2016) Logistic 1


# ----------------------------------------------------------------------------
# 2. Set the discrimination design problem
# ----------------------------------------------------------------------------
# The true model is 'linlogi4' and we set the nominal values
para_logit_4 <- c(1, 1, 1)
# Create a model list as the input of our PSO-QN algorithm
model_logit <- list(
  # The fitst object should always be the true model
  # Input the model function to the label 'model'
  # Input the nominal values to the labal 'para'
  list(model = logit_4, para = para_logit_4),
  # Starting from the second object, specify the rival models
  # Input the model function to the label 'model'
  # Input the lower and upper bounds for the rival parameters to the labal
  # 'paraLower' and 'paraUpper' (based on experiences)
  list(model = logit_3, paraLower = c(-10, -10), paraUpper = c(10, 10)),
  # The third model
  list(model = logit_2, paraLower = c(-10, -10), paraUpper = c(10, 10)),
  # The fourth model
  list(model = logit_1, paraLower = -10, paraUpper = 10)
)


# ----------------------------------------------------------------------------
# 3. Create distance function using C++ codes (faster)
# ----------------------------------------------------------------------------
# KL divergence between two logistic regression models
kl_div_logit <- cppFunction('
  Rcpp::NumericVector logit_diff(Rcpp::NumericVector xt, Rcpp::NumericVector xr) {
    Rcpp::NumericVector div(xt.size());
    double et, er, mt, mr;
    for (int i = 0; i < xt.size(); i++) {
      et = std::exp(xt(i)); mt = et/(1.0 + et);
      er = std::exp(xr(i)); mr = er/(1.0 + er);
      div(i) = mt*std::log(mt/mr) + (1.0 - mt)*std::log((1.0 - mt)/(1.0 - mr));
    }
    return div;
}')

# Or, create distance function using R codes
#kl_div_logit <- function(xt, xr) {
#  exp_t <- exp(xt)
#  exp_r <- exp(xr)
#  mu_t <- exp_t/(1 + exp_t)
#  mu_r <- exp_r/(1 + exp_r)
#  mu_t*(log(mu_t) - log(mu_r)) + (1 - mu_t)*(log(1.0 - mu_t) - log(1.0 - mu_r))
#}


# ----------------------------------------------------------------------------
# 4. Set PSO-QN Algorithm Parameters
# ----------------------------------------------------------------------------

# Set PSO options for pariwise discrimination design cases
PSO_INFO <- getPSOInfo(nSwarm = 32, maxIter = 200)
# Set L-BFGS algorithm options
LBFGS_INFO <- getLBFGSInfo(LBFGS_RETRY = 4)


# ----------------------------------------------------------------------------
# 5-1. Run PSO-QN for the case: logit_3 vs. logit_4
# ----------------------------------------------------------------------------

# Get list for two models
two_model_logit_43 <- list(model_logit[[1]], model_logit[[2]])
# Run PSO-QN algorithm
res_logit_43 <- DiscrimOD(MODEL_INFO = two_model_logit_43, DISTANCE = kl_div_logit, nSupp = 3,
                          dsLower = 0, dsUpper = 1, crit_type = "pair_fixed_true",
                          PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO,
                          seed = 10, verbose = TRUE)
# Check the result
names(res_logit_43)
# $BESTDESIGN: the best design found by the PSO-QN algorithm
# $BESTVAL:    the criterion value of the best design
# $GBESTHIST:  the updating path of the cirterion value in PSO loop
# $CPUTIME:    the computing time
res_logit_43$BESTDESIGN
#          dim_1       weight
#obs_1 0.0000000 2.413067e-01
#obs_2 0.0000000 7.586933e-01
#obs_3 0.6676636 3.749152e-33
res_logit_43$BESTVAL
# 0.1109441

# Check by the equivalence theorem
eqv_logit_43 <- equivalence(ngrid = 100, PSO_RESULT = res_logit_43,
                            MODEL_INFO = two_model_logit_43, DISTANCE = kl_div_logit,
                            dsLower = 0, dsUpper = 1, crit_type = "pair_fixed_true",
                            PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

names(eqv_logit_43)
# $Grid_1:   the grid on the design space
# $Grid_2:   (currently not used)
# $DirDeriv: the values of directional derivative function on the grid
# $MAX_DD:   the maximal value of the directional derivative function
eqv_logit_43$MAX_DD
# -2.775558e-17

# Draw the curve of directional derivative function
plot(eqv_logit_43$Grid_1, eqv_logit_43$DirDeriv, type = "l", col = "black",
     xlab = "x", ylab = "Directional Derivative",
     main = "KL-optimal design for case logit_3 vs. logit_4")
abline(h = 0, col = "grey50", lty = 2)
points(res_logit_43$BESTDESIGN[,1], rep(0, nrow(res_logit_43$BESTDESIGN)), pch = 16)


# ----------------------------------------------------------------------------
# 5-2. Run PSO-QN for the case: logit_2 vs. logit_4
# ----------------------------------------------------------------------------

# Get list for two models
two_model_logit_42 <- list(model_logit[[1]], model_logit[[3]])
# Run PSO-QN algorithm
res_logit_42 <- DiscrimOD(MODEL_INFO = two_model_logit_42, DISTANCE = kl_div_logit, nSupp = 3,
                          dsLower = 0, dsUpper = 1, crit_type = "pair_fixed_true",
                          PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO,
                          seed = 10, verbose = TRUE)
# Check the result
res_logit_42$BESTDESIGN
#          dim_1    weight
#obs_1 0.0000000 0.2202030
#obs_2 0.4173956 0.4646236
#obs_3 1.0000000 0.3151734
res_logit_42$BESTVAL
# 0.0008454829

# Check by the equivalence theorem
eqv_logit_42 <- equivalence(ngrid = 100, PSO_RESULT = res_logit_42,
                            MODEL_INFO = two_model_logit_42, DISTANCE = kl_div_logit,
                            dsLower = 0, dsUpper = 1, crit_type = "pair_fixed_true",
                            PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

eqv_logit_42$MAX_DD
# 1.745368e-08

# Draw the curve of directional derivative function
plot(eqv_logit_42$Grid_1, eqv_logit_42$DirDeriv, type = "l", col = "black",
     xlab = "x", ylab = "Directional Derivative",
     main = "KL-optimal design for case logit_2 vs. logit_4")
abline(h = 0, col = "grey50", lty = 2)
points(res_logit_42$BESTDESIGN[,1], rep(0, nrow(res_logit_42$BESTDESIGN)), pch = 16)


# ----------------------------------------------------------------------------
# 5-3. Run PSO-QN for the case: logit_1 vs. logit_4
# ----------------------------------------------------------------------------

# Get list for two models
two_model_logit_41 <- list(model_logit[[1]], model_logit[[4]])
# Run PSO-QN algorithm
res_logit_41 <- DiscrimOD(MODEL_INFO = two_model_logit_41, DISTANCE = kl_div_logit, nSupp = 3,
                          dsLower = 0, dsUpper = 1, crit_type = "pair_fixed_true",
                          PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO,
                          seed = 10, verbose = TRUE)
# Check the result
res_logit_41$BESTDESIGN
#          dim_1       weight
#obs_1 0.0000000 5.422216e-03
#obs_2 0.0000000 9.945778e-01
#obs_3 0.6402028 3.749152e-33
res_logit_41$BESTVAL
# 0.1109441

# Check by the equivalence theorem
eqv_logit_41 <- equivalence(ngrid = 100, PSO_RESULT = res_logit_41,
                            MODEL_INFO = two_model_logit_41, DISTANCE = kl_div_logit,
                            dsLower = 0, dsUpper = 1, crit_type = "pair_fixed_true",
                            PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

eqv_logit_41$MAX_DD
# -2.775558e-17

# Draw the curve of directional derivative function
plot(eqv_logit_41$Grid_1, eqv_logit_41$DirDeriv, type = "l", col = "black",
     xlab = "x", ylab = "Directional Derivative",
     main = "KL-optimal design for case logit_2 vs. logit_4")
abline(h = 0, col = "grey50", lty = 2)
points(res_logit_41$BESTDESIGN[,1], rep(0, nrow(res_logit_41$BESTDESIGN)), pch = 16)


# ----------------------------------------------------------------------------
# 6. Run PSO-S-QN for max-min KL-optimal design for all models
# ----------------------------------------------------------------------------

# Set PSO options for pariwise discrimination design cases
PSO_MM_INFO <- getPSOInfo(nSwarm = 32, maxIter = 400)

# Get the values of T-optimal criterion found previouly which will be
# the denominators of the T-efficiencies
eff_denom <- c(res_logit_43$BESTVAL, res_logit_42$BESTVAL, res_logit_41$BESTVAL)

res_logit <- DiscrimOD(MODEL_INFO = model_logit, DISTANCE = kl_div_logit, nSupp = 3,
                       dsLower = 0, dsUpper = 1, crit_type = "maxmin_fixed_true",
                       MaxMinStdVals = eff_denom,
                       PSO_INFO = PSO_MM_INFO, LBFGS_INFO = LBFGS_INFO,
                       seed = 10, verbose = TRUE)
# Check the result
names(res_logit)
# $BESTDESIGN: the best design found by the PSO-QN algorithm
# $BESTVAL:    the criterion value of the best design
# $GBESTHIST:  the updating path of the cirterion value in PSO loop
# $CPUTIME:    the computing time
res_logit$BESTDESIGN
#         dim_1    weight
#obs_1 0.000000 0.6184616
#obs_2 0.359819 0.2392604
#obs_3 1.000000 0.1422780
res_logit$BESTVAL
# 0.6184616

eqv_logit <- equivalence(ngrid = 100, PSO_RESULT = res_logit,
                         MODEL_INFO = model_logit, DISTANCE = kl_div_logit,
                         dsLower = 0, dsUpper = 1, crit_type = "maxmin_fixed_true",
                         MaxMinStdVals = eff_denom,
                         PSO_INFO = PSO_MM_INFO, LBFGS_INFO = LBFGS_INFO)

names(eqv_logit)
# $Grid_1:   the grid on the design space
# $Grid_2:   (currently not used)
# $DirDeriv: the values of directional derivative function on the grid
# $MAX_DD:   the maximal value of the directional derivative function
# $alpha:    the best weight vector of the efficiency values
eqv_logit$MAX_DD
# 1.465643e-06

# Draw the curve of directional derivative function
plot(eqv_logit$Grid_1, eqv_logit$DirDeriv, type = "l", col = "black",
     xlab = "x", ylab = "Directional Derivative",
     main = "max-min KL-optimal design for logit example")
abline(h = 0, col = "grey50", lty = 2)
points(res_logit$BESTDESIGN[,1], rep(0, nrow(res_logit$BESTDESIGN)), pch = 16)

