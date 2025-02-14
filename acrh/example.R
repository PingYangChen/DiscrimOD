library(DiscrimOD)

# Atkinson and Fedorov (1975a): T-optimal
# Two R functions of competing models are given by
af1 <- function(x, p) p[1] + p[2]*exp(x) + p[3]*exp(-x)
af2 <- function(x, p) p[1] + p[2]*x + p[3]*x^2

af1_v <- function(x, p) p[1]
af2_v <- function(x, p) p[1]

# Set the model information
# The nominla value in 'm1' is 4.5, -1.5, -2.0
# For 'af2', we set the parameter space to be [-10, 10]^3 and
# the initial guess (for LBFGS) of the rival model parameter is zero vector
AF_para_af1 <- c(4.5, -1.5, -2)
af_info_12 <- list(
  # The first list should be the true model and the specified nominal values
  list(mean = af1, disp = af1_v, meanPara = AF_para_af1, dispPara = 1),
  # Then the rival models are listed accordingly. We also need to specify the model space.
  list(mean = af2, disp = af2_v,
       meanParaLower = rep(-10, 3), meanParaUpper = rep(10, 3),
       dispParaLower = c(1), dispParaUpper = c(5) )
)
# Define the R function for the distance measure in T-optimal criterion
# xt is the mean values of the true model
# xr is the mean values of the rival model
sq_diff <- function(xt, xr, vt, vr) (vt/vr)*(xt - xr)^2

# Initialize PSO and BFGS options
PSO_INFO <- getPSOInfo(nSwarm = 32, maxIter = 100)
LBFGS_INFO <- getLBFGSInfo(LBFGS_RETRY = 2)


# getDesignInfo(D_TYPE = "approx", MODEL_INFO = af_info_12, dist_func = sq_diff,
#                         crit_type = "pair_fixed_true", MaxMinStdVals = 0, minWt = .0,
#                         dSupp = 1, nSupp = 4, dsLower = -1, dsUpper = 1)


# Find T-optimal design for models af1 and af2
af_res_12 <- DiscrimOD(MODEL_INFO = af_info_12, DISTANCE = sq_diff,
                       nSupp = 4, dsLower = -1, dsUpper = 1, crit_type = "pair_fixed_true",
                       PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO, seed = NULL, verbose = TRUE)

round(af_res_12$BESTDESIGN, 3) # The resulting design
af_res_12$BESTVAL # The T-optimal criterion value
af_res_12$CPUTIME # CPU time


designCriterion(af_res_12$BESTDESIGN, MODEL_INFO = af_info_12, DISTANCE = sq_diff,
                dsLower = -1, dsUpper = 1, crit_type = "pair_fixed_true", MaxMinStdVals = NULL,
                PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)


# Test optimality by equivalence theorem
af_eqv_12 <- equivalence(ngrid = 100, PSO_RESULT = af_res_12, MODEL_INFO = af_info_12,
                         DISTANCE = sq_diff, dsLower = -1, dsUpper = 1, crit_type = "pair_fixed_true",
                         PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

# Draw the directional derivative curve
plot(af_eqv_12$Grid_1, af_eqv_12$DirDeriv, type = "l", col = "blue",
     main = "af_res_12", xlab = "x", ylab = "Directional Derivative"); abline(h = 0)
points(af_res_12$BESTDESIGN[,1], rep(0, nrow(af_res_12$BESTDESIGN)), pch = 16)
