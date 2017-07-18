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

# Or, create competing models using R codes
#af1975_1 <- function(x, p) p[1] + p[2]*exp(x) + p[3]*exp(-x)
#af1975_2 <- function(x, p) p[1] + p[2]*x + p[3]*x^2


# ----------------------------------------------------------------------------
# 2. Set the discrimination design problem
# ----------------------------------------------------------------------------
# The true model is 'af1975_1' and we set the nominal values
para_af1975_1 <- c(4.5, -1.5, -2)
# Create a model list as the input of our PSO-QN algorithm
model_af1975_12 <- list(
  # The fitst object should always be the true model
  # Input the model function to the label 'model'
  # Input the nominal values to the labal 'para'
  list(model = af1975_1, para = para_af1975_1),
  # Starting from the second object, specify the rival models
  # Input the model function to the label 'model'
  # Input the lower and upper bounds for the rival parameters to the labal
  # 'paraLower' and 'paraUpper' (based on experiences)
  list(model = af1975_2, paraLower = rep(-10, 3), paraUpper = rep(10, 3))
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
# 4-1. Run PSO-QN Algorithm
# ----------------------------------------------------------------------------

# Set PSO options for pariwise discrimination design cases
PSO_INFO <- getPSOInfo(nSwarm = 32, maxIter = 200)
# Set L-BFGS algorithm options
LBFGS_INFO <- getLBFGSInfo(LBFGS_RETRY = 4)
# Run PSO-QN algorithm
PSOQN_af1975_12 <- DiscrimOD(MODEL_INFO = model_af1975_12, DISTANCE = sq_diff, nSupp = 4,
                             dsLower = -1, dsUpper = 1, crit_type = "pair_fixed_true",
                             PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO,
                             seed = 10, verbose = TRUE)

PSOQN_af1975_12$BESTDESIGN
#           dim_1     weight
#obs_1 -1.0000000 0.25267570
#obs_2 -0.6692933 0.42772012
#obs_3  0.1438411 0.24732430
#obs_4  0.9569755 0.07227988

# ----------------------------------------------------------------------------
# 4-2. Run NestedPSO Algorithm
# ----------------------------------------------------------------------------

# Specify options for outer and inner PSO loops
NESTEDPSO_INFO <- getPSOInfo(nSwarm = c(32, 32), maxIter = c(200, 100))
# Set 'IF_INNER_LBFGS = FALSE' to disable the L-BFGS algorithm
LBFGS_NOTRUN <- getLBFGSInfo(IF_INNER_LBFGS = FALSE)
# Run PSO-QN algorithm
NestedPSO_af1975_12 <- DiscrimOD(MODEL_INFO = model_af1975_12, DISTANCE = sq_diff, nSupp = 4,
                                 dsLower = -1, dsUpper = 1, crit_type = "pair_fixed_true",
                                 PSO_INFO = NESTEDPSO_INFO, LBFGS_INFO = LBFGS_NOTRUN,
                                 seed = 10, verbose = TRUE)

NestedPSO_af1975_12$BESTDESIGN
#           dim_1    weight
#obs_1 0.08800934 0.1490641
#obs_2 0.71710883 0.2694418
#obs_3 1.00000000 0.2742851
#obs_4 1.00000000 0.3072091

# ----------------------------------------------------------------------------
# 4-3. Run Fedorov-Wynn Algorithm
# ----------------------------------------------------------------------------

# Set Fedorov-Wynn options
FED_INFO <- getFEDInfo(FED_MAXIT = 200, FED_TRIM = 3, FED_TRIM_EPS = 1e-2,
                       freeRun = 1.0, FED_EPS = 1e-6, FED_ALPHA_GRID = 20)

FedWynn_af1975_12 <- DiscrimFedWynn(MODEL_INFO = model_af1975_12, DISTANCE = sq_diff,
                                    dsLower = -1, dsUpper = 1,
                                    FED_INFO = FED_INFO, LBFGS_INFO = LBFGS_INFO,
                                    seed = 10, verbose = TRUE)

FedWynn_af1975_12$BESTDESIGN
#             dim_1      weight
#obs_1  -0.99936418 0.248895086
#obs_2  -0.67426521 0.052631579
#obs_3  -0.66344320 0.354205848
#obs_4  -0.60458921 0.008473679
#obs_5  -0.38646299 0.022562752
#obs_6  -0.14618467 0.022562752
#obs_7   0.01495641 0.022562752
#obs_8   0.14054055 0.128702036
#obs_9   0.15984290 0.062394582
#obs_10  0.38620416 0.022562752
#obs_11  0.96532781 0.043238541
#obs_12  0.99620636 0.011207642

# ----------------------------------------------------------------------------
# 4-4. Run Remes Algorithm
# ----------------------------------------------------------------------------

Remes_af1975_12 <- DiscrimUnifApproxT(MODEL_INFO = model_af1975_12, nSupp = 4,
                                      dsLower = -1, dsUpper = 1,
                                      REMES_MAXIT = 200, REMES_FreeRun = 1.0, REMES_EPS = 1e-2,
                                      LBFGS_INFO = LBFGS_INFO, seed = 10, verbose = TRUE)

Remes_af1975_12$BESTDESIGN # may not be a design
#           dim_1      weight
#obs_1 -1.0000000  1.31989632
#obs_2 -0.6673253 -0.18147928
#obs_3  0.1497666 -0.09985344
#obs_4  1.0000000 -0.03856361


# ----------------------------------------------------------------------------
# 5. Use L-BFGS algorithm to calculate the criterion values for each result
# ----------------------------------------------------------------------------

PSOQN_ACCVALUE <- designCriterion(PSOQN_af1975_12$BESTDESIGN,
                                  MODEL_INFO = model_af1975_12, DISTANCE = sq_diff,
                                  dsLower = -1, dsUpper = 1, crit_type = "pair_fixed_true",
                                  PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

PSOQN_af1975_12$BESTVAL
# 0.001086725
PSOQN_ACCVALUE$cri_val
# 0.001086725

NestedPSO_ACCVALUE <- designCriterion(NestedPSO_af1975_12$BESTDESIGN,
                                      MODEL_INFO = model_af1975_12, DISTANCE = sq_diff,
                                      dsLower = -1, dsUpper = 1, crit_type = "pair_fixed_true",
                                      PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

NestedPSO_af1975_12$BESTVAL
# 0.3964421
NestedPSO_ACCVALUE$cri_val
# 1.224024e-12 (The correct value is much smaller)

FedWynn_ACCVALUE <- designCriterion(FedWynn_af1975_12$BESTDESIGN,
                                    MODEL_INFO = model_af1975_12, DISTANCE = sq_diff,
                                    dsLower = -1, dsUpper = 1, crit_type = "pair_fixed_true",
                                    PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)
FedWynn_af1975_12$BESTVAL
# 0.001017412
FedWynn_ACCVALUE$cri_val
# 0.001017297

if (is.character(Remes_af1975_12$BESTVAL)) {
  print("the result is not a experimental design")
} else {
  Remes_ACCVALUE <- designCriterion(Remes_af1975_12$BESTDESIGN,
                                    MODEL_INFO = model_af1975_12, DISTANCE = sq_diff,
                                    dsLower = -1, dsUpper = 1, crit_type = "pair_fixed_true",
                                    PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)
}

Remes_af1975_12$BESTVAL
# "NOT A DESIGN"
# "the result is not a experimental design"

# ----------------------------------------------------------------------------
# 6. Check equivalecne theorem for each result
# ----------------------------------------------------------------------------

# Check PSO-QN result by the equivalence theorem
eqv_PSOQN <- equivalence(ngrid = 100, PSO_RESULT = PSOQN_af1975_12,
                         MODEL_INFO = model_af1975_12, DISTANCE = sq_diff,
                         dsLower = -1, dsUpper = 1, crit_type = "pair_fixed_true",
                         PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

# Draw the curve of directional derivative function
plot(eqv_PSOQN$Grid_1, eqv_PSOQN$DirDeriv, type = "l", col = "black",
     xlab = "x", ylab = "Directional Derivative",
     main = "PSO-QN result for case af1975_2 vs. af1975_1")
abline(h = 0, col = "grey50", lty = 2)
points(PSOQN_af1975_12$BESTDESIGN[,1], rep(0, nrow(PSOQN_af1975_12$BESTDESIGN)), pch = 16)

# Check NestedPSO result by the equivalence theorem
eqv_NestedPSO <- equivalence(ngrid = 100, PSO_RESULT = NestedPSO_af1975_12,
                             MODEL_INFO = model_af1975_12, DISTANCE = sq_diff,
                             dsLower = -1, dsUpper = 1, crit_type = "pair_fixed_true",
                             PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

# Draw the curve of directional derivative function
plot(eqv_NestedPSO$Grid_1, eqv_NestedPSO$DirDeriv, type = "l", col = "black",
     xlab = "x", ylab = "Directional Derivative",
     main = "NestedPSO result for case af1975_2 vs. af1975_1")
abline(h = 0, col = "grey50", lty = 2)
points(NestedPSO_af1975_12$BESTDESIGN[,1], rep(0, nrow(NestedPSO_af1975_12$BESTDESIGN)), pch = 16)

# Check Fedorov-Wynn result by the equivalence theorem
eqv_FedWynn <- equivalence(ngrid = 100, DESIGN = FedWynn_af1975_12$BESTDESIGN,
                           MODEL_INFO = model_af1975_12, DISTANCE = sq_diff,
                           dsLower = -1, dsUpper = 1, crit_type = "pair_fixed_true",
                           PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

# Draw the curve of directional derivative function
plot(eqv_FedWynn$Grid_1, eqv_FedWynn$DirDeriv, type = "l", col = "black",
     xlab = "x", ylab = "Directional Derivative",
     main = "Fedorov-Wynn result for case af1975_2 vs. af1975_1")
abline(h = 0, col = "grey50", lty = 2)
points(FedWynn_af1975_12$BESTDESIGN[,1], rep(0, nrow(FedWynn_af1975_12$BESTDESIGN)), pch = 16)


if (is.character(Remes_af1975_12$BESTVAL)) {
  print("the result is not a experimental design")
} else {
  # Check Remes result by the equivalence theorem
  eqv_Remes <- equivalence(ngrid = 100, PSO_RESULT = Remes_af1975_12,
                           MODEL_INFO = model_af1975_12, DISTANCE = sq_diff,
                           dsLower = -1, dsUpper = 1, crit_type = "pair_fixed_true",
                           PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

  # Draw the curve of directional derivative function
  plot(eqv_Remes$Grid_1, eqv_Remes$DirDeriv, type = "l", col = "black",
       xlab = "x", ylab = "Directional Derivative",
       main = "Remes result for case af1975_2 vs. af1975_1")
  abline(h = 0, col = "grey50", lty = 2)
  points(Remes_af1975_12$BESTDESIGN[,1], rep(0, nrow(Remes_af1975_12$BESTDESIGN)), pch = 16)
}


# ----------------------------------------------------------------------------
# 7. Compare the inner loops of PSO-QN and NestedPSO
# ----------------------------------------------------------------------------

# Set stopping critria for PSO and L-BFGS algorithms
# Focus on the setting of the second component of each parameter (inner PSO loop)
NPSO_STOP_EARLY <- getPSOInfo(nSwarm = c(32, 32), maxIter = c(200, 100),
                              freeRun = c(1, 0.25), tol = c(1e-6, 1e-06))
# Set L-BFGS algorithm options with stopping criterion
LBFGS_STOP_EARLY <- getLBFGSInfo(LBFGS_RETRY = 4, FVAL_EPS = 1e-6)

optiamlDesign <- PSOQN_af1975_12$BESTDESIGN

valMat <- matrix(0, 100, 2); colnames(valMat) <- c("L-BFGS", "PSO")
for (i in 1:nrow(valMat)) {

  LFBSG_VAL <- designCriterion(optiamlDesign, MODEL_INFO = model_af1975_12, DISTANCE = sq_diff,
                  dsLower = -1, dsUpper = 1, crit_type = "pair_fixed_true",
                  PSO_INFO = NPSO_STOP_EARLY, LBFGS_INFO = LBFGS_STOP_EARLY)

  PSO_VAL <- designCriterion(optiamlDesign, MODEL_INFO = model_af1975_12, DISTANCE = sq_diff,
                               dsLower = -1, dsUpper = 1, crit_type = "pair_fixed_true",
                               PSO_INFO = NPSO_STOP_EARLY, LBFGS_INFO = LBFGS_NOTRUN)

  valMat[i,] <- c(LFBSG_VAL$cri_val, PSO_VAL$cri_val)
}

sumTable <- t(sapply(1:2, function(i) {
  c(range(valMat[,i]), sd(valMat[,i]))
}))
dimnames(sumTable) <- list(colnames(valMat), c("Min.", "Max.", "Std.Err."))
sumTable
#               Min.        Max.     Std.Err.
#L-BFGS 0.001086725 0.001086725 1.347162e-14
#PSO    0.004124766 0.180649305 4.103751e-02

