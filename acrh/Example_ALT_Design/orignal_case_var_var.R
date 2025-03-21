library(DiscrimOD)

source("./divFunctions.R")

# Case: Orginal case with fixed variance

# Given LN is TRUE
# Use the objective function 'kldiv_censored_lnIsTrue_tc_5000' in divFunctions.R
# OR
# Given WB is TRUE
# Use the objective function 'kldiv_censored_wbIsTrue_tc_5000' in divFunctions.R
ln_mean <- function(x, p) p[1] + p[2]*(11605/(x + 273.15))
we_mean <- function(x, p) p[1] + p[2]*(11605/(x + 273.15))
ln_disp <- function(x, p) rep(p[1], length(x)) # make sure output length is equal to length of x (nSupp)
we_disp <- function(x, p) rep(p[1], length(x)) # make sure output length is equal to length of x (nSupp)
#
lm_para <- c(-13.46341, 0.6277041)
#
model_info <- list(
  # The first list should be the true model and the specified nominal values
  list(mean = ln_mean, disp = ln_disp, meanPara = lm_para, dispPara = 0.9780103),
  # Then the rival models are listed accordingly. We also need to specify the model space.
  list(mean = we_mean, disp = we_disp,
       meanParaLower = c(-15, 0.4), meanParaUpper = c(-10, 0.8),
       dispParaLower = c(0.9780103), dispParaUpper = c(4.9780103) )
)
#
nSupp <- 3
dsRange <- c(10, 40)
# Initialize PSO and BFGS options
PSO_INFO <- getPSOInfo(nSwarm = 64, maxIter = 100)
LBFGS_INFO <- getLBFGSInfo(LBFGS_RETRY = 5)

# Find KL-optimal design for models
res <- DiscrimOD(MODEL_INFO = model_info, DISTANCE = kldiv_censored_lnIsTrue_tc_5000,
                 nSupp = nSupp, dsLower = dsRange[1], dsUpper = dsRange[2], crit_type = "pair_fixed_true",
                 PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO, seed = 100, verbose = TRUE)


round(res$BESTDESIGN, 3) # The resulting design
res$BESTVAL # The CKL-optimal criterion value
res$CPUTIME # CPU time

designCriterion(res$BESTDESIGN, MODEL_INFO = model_info, DISTANCE = kldiv_censored_lnIsTrue_tc_5000,
                dsLower = dsRange[1], dsUpper = dsRange[2], crit_type = "pair_fixed_true", MaxMinStdVals = NULL,
                PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)


# Test optimality by equivalence theorem
eqv <- equivalence(ngrid = 100, PSO_RESULT = res, MODEL_INFO = model_info,
                   DISTANCE = kldiv_censored_lnIsTrue_tc_5000,
                   dsLower = dsRange[1], dsUpper = dsRange[2], crit_type = "pair_fixed_true",
                   PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

# Draw the directional derivative curve
plot(eqv$Grid_1, eqv$DirDeriv, type = "l", col = "blue",
     main = "res", xlab = "x", ylab = "Directional Derivative"); abline(h = 0)
points(res$BESTDESIGN[,1], rep(0, nrow(res$BESTDESIGN)), pch = 16)

