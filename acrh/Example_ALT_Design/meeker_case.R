library(DiscrimOD)

#source("./divFunctions.R")

# Case: Pascual, F. G., & Meeker, W. Q. (1997). Analysis of fatigue data with runouts
#       based on a model with nonconstant standard deviation and a fatigue limit parameter.
#       Journal of testing and evaluation, 25(3), 292-301.

# Given LN is TRUE
# Use the objective function 'kldiv_censored_lnIsTrue_tc_1000' in divFunctions.R
ln_mean <- function(x, p) p[1] + p[2]*log(x - 75.71 + 1e-12)
we_mean <- function(x, p) p[1] + p[2]*log(x - 75.71 + 1e-12)
ln_disp <- function(x, p) exp(p[1] + p[2]*log(x + 1e-12))
we_disp <- function(x, p) exp(p[1] + p[2]*log(x + 1e-12))
#
ln_mean_para <- c(14.75, -1.39)
ln_disp_para <- c(10.97, -2.50)
#
model_info <- list(
  # The first list should be the true model and the specified nominal values
  list(mean = ln_mean, disp = ln_disp, meanPara = ln_mean_para, dispPara = ln_disp_para),
  # Then the rival models are listed accordingly. We also need to specify the model space.
  list(mean = we_mean, disp = we_disp,
       meanParaLower = c(13, -2.5), meanParaUpper = c(20, -0.01),
       dispParaLower = c(10, -3.0), dispParaUpper = c(20, -0.01) )
)
#
nSupp <- 4
dsRange <- c(76, 150)
# Initialize PSO and BFGS options
PSO_INFO <- getPSOInfo(nSwarm = 32, maxIter = 100)
LBFGS_INFO <- getLBFGSInfo(LBFGS_RETRY = 5)

# Find KL-optimal design for models
res <- DiscrimOD(MODEL_INFO = model_info, DISTANCE = kldiv_censored_lnIsTrue_tc_1000,
                 nSupp = nSupp, dsLower = dsRange[1], dsUpper = dsRange[2], crit_type = "pair_fixed_true",
                 PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO, seed = 1, verbose = TRUE)


round(res$BESTDESIGN, 3) # The resulting design
res$BESTVAL # The CKL-optimal criterion value
res$CPUTIME # CPU time

#testDesign <- cbind(c(86.29, 89.456, 90, 142), c(0,0,1,0))

designCriterion(res$BESTDESIGN, MODEL_INFO = model_info, DISTANCE = kldiv_censored_lnIsTrue_tc_1000,
                dsLower = dsRange[1], dsUpper = dsRange[2], crit_type = "pair_fixed_true", MaxMinStdVals = NULL,
                PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

# Test optimality by equivalence theorem
eqv <- equivalence(ngrid = 100, PSO_RESULT = res, MODEL_INFO = model_info,
                   DISTANCE = kldiv_censored_lnIsTrue_tc_1000,
                   dsLower = dsRange[1], dsUpper = dsRange[2], crit_type = "pair_fixed_true",
                   PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

# Draw the directional derivative curve
plot(eqv$Grid_1, eqv$DirDeriv, type = "l", col = "blue",
     main = "res", xlab = "x", ylab = "Directional Derivative"); abline(h = 0)
points(res$BESTDESIGN[,1], rep(0, nrow(res$BESTDESIGN)), pch = 16)

