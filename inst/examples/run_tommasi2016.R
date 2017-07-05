# Example: A and Fedorov (1975a, b)
library(DiscrimOD)

# Set PSO options for pariwise discrimination design cases
PSO_INFO <- getPSOInfo(nSwarm = c(16, 32), maxIter = c(100, 100))
# Set PSO options for max-min discrimination design cases
PSO_MAXMIN <- getPSOInfo(nSwarm = c(32, 32), maxIter = c(200, 100))
# Set L-BFGS algorithm options
LBFGS_INFO <- getLBFGSInfo(LBFGS_RETRY = 2)
# Set a NOT-RUN L-BFGS algorithm for trying NestedPSO (for fun)
LBFGS_NOTRUN <- getLBFGSInfo(IF_INNER_LBFGS = FALSE)

# Create competing models
linlogi4 <- function(x, p) p[1] + p[2]*x + p[3]*(x^2)  # Tommasi et al. (2016) Logistic 4
linlogi3 <- function(x, p) x*(p[1] + p[2]*x)           # Tommasi et al. (2016) Logistic 3
linlogi2 <- function(x, p) p[1] + p[2]*x               # Tommasi et al. (2016) Logistic 2
linlogi1 <- function(x, p) p[1]*x                      # Tommasi et al. (2016) Logistic 1

# Set the nominal values for the first model (null model)
para_linlogi_4 <- c(1, 1, 1)
# Create the model list
model_linearLogistic <- list(
  list(model = linlogi4, para = para_linlogi_4),
  list(model = linlogi3, paraLower = c(-10, -10), paraUpper = c(10, 10)),
  list(model = linlogi2, paraLower = c(-10, -10), paraUpper = c(10, 10)),
  list(model = linlogi1, paraLower = c(-10), paraUpper = c(10))
)

# Create the lists for pairwise discrimination designs
two_model_linearLogistic <- list(
  list(model_linearLogistic[[1]], model_linearLogistic[[2]]),
  list(model_linearLogistic[[1]], model_linearLogistic[[3]]),
  list(model_linearLogistic[[1]], model_linearLogistic[[4]])
)
# Specify the number fo support points for pairwise discrimination designs
two_nSupp <- c(3, 3, 3)

# Set distance function for KL-optimal design
logit_diff <- function(xt, xr) {
  exp_t <- exp(xt)
  exp_r <- exp(xr)
  mu_t <- exp_t/(1 + exp_t)
  mu_r <- exp_r/(1 + exp_r)
  mu_t*(log(mu_t) - log(mu_r)) + (1 - mu_t)*(log(1.0 - mu_t) - log(1.0 - mu_r))
}

# Set the number of regeneration when the resulting design is not optimal
reGen <- 4; tolMaxDD <- 1e-6;

# Get optimal designs with two different approaches
pairRes_linearLogistic_q <- pairRes_linearLogistic_n <- vector("list", length(two_model_linearLogistic))

# Start for each pairwise discrimination design
for (iC in 1:length(two_model_linearLogistic)) {

  MODEL_INFO <- two_model_linearLogistic[[iC]]
  DISTANCE <- logit_diff

  # PSO-QN
  genCount <- 0; testMaxDD <- 1e20
  while ((genCount < reGen) & (testMaxDD > tolMaxDD)) {
    genCount <- genCount + 1
    # Use PSO-QN algorithm to find the optimal deisgn
    out_q <- DiscrimOD(MODEL_INFO, DISTANCE, two_nSupp[iC], dsLower = 0, dsUpper = 1,
                       crit_type = "pair_fixed_true",
                       PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)
    # equivalence theorem
    eqv_q <- equivalence(ngrid = 100, PSO_RESULT = out_q, MODEL_INFO = MODEL_INFO, DISTANCE = DISTANCE,
                         dsLower = 0, dsUpper = 1,
                         crit_type = "pair_fixed_true",
                         PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

    testMaxDD <- eqv_q$MAX_DD
    # If not optimal, regenerate
    if (genCount > 1) {
      if (testMaxDD < pairRes_linearLogistic_q[[iC]]$eqv$MAX_DD) {
        pairRes_linearLogistic_q[[iC]] <- list(out = out_q, eqv = eqv_q)
      }
    } else {
      pairRes_linearLogistic_q[[iC]] <- list(out = out_q, eqv = eqv_q)
    }
  }

  # NestedPSO
  genCount <- 0; testMaxDD <- 1e20
  while ((genCount < reGen) & (testMaxDD > tolMaxDD)) {
    genCount <- genCount + 1
    # Use NestedPSO algorithm to find the optimal deisgn
    out_n <- DiscrimOD(MODEL_INFO, DISTANCE, two_nSupp[iC], dsLower = 0, dsUpper = 1,
                       crit_type = "pair_fixed_true",
                       PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_NOTRUN)
    # equivalence theorem
    eqv_n <- equivalence(ngrid = 100, PSO_RESULT = out_n, MODEL_INFO = MODEL_INFO, DISTANCE = DISTANCE,
                         dsLower = 0, dsUpper = 1,
                         crit_type = "pair_fixed_true",
                         PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

    testMaxDD <- eqv_n$MAX_DD
    # If not optimal, regenerate
    if (genCount > 1) {
      if (testMaxDD < pairRes_linearLogistic_n[[iC]]$eqv$MAX_DD) {
        pairRes_linearLogistic_n[[iC]] <- list(out = out_n, eqv = eqv_n)
      }
    } else {
      pairRes_linearLogistic_n[[iC]] <- list(out = out_n, eqv = eqv_n)
    }
  }
}

# Get results from the PSO-QN algorithm
for (iC in 1:length(two_model_linearLogistic)) {
  tmp <- pairRes_linearLogistic_q[[iC]]
  # The max-min optimal discrimination design
  print(round(tmp$out$BESTDESIGN, 3))
  # The best efficiency value
  print(tmp$out$BESTVAL)
  # The computing time
  print(tmp$out$CPUTIME)
}

# Draw the plot for checking the equivalence theorem
par(mfrow = c(2, 2))
for (iC in 1:length(two_model_linearLogistic)) {
  tmp <- pairRes_linearLogistic_q[[iC]]
  plot(tmp$eqv$Grid_1, tmp$eqv$DirDeriv, type = "l", col = "black",
       xlab = "x", ylab = "Directional Derivative", main = paste0("Result from PSO-QN: ", iC))
  abline(h = 0, col = "grey50", lty = 2)
  points(tmp$out$BESTDESIGN[,1], rep(0, nrow(tmp$out$BESTDESIGN)), pch = 16)
}

par(mfrow = c(2, 2))
for (iC in 1:length(two_model_linearLogistic)) {
  tmp <- pairRes_linearLogistic_n[[iC]]
  plot(tmp$eqv$Grid_1, tmp$eqv$DirDeriv, type = "l", col = "black",
       xlab = "x", ylab = "Directional Derivative", main = paste0("Result from NestedPSO: ", iC))
  abline(h = 0, col = "grey50", lty = 2)
  points(tmp$out$BESTDESIGN[,1], rep(0, nrow(tmp$out$BESTDESIGN)), pch = 16)
}

par(mfrow = c(1, 1))

# Get max-min T-optimal optimal design
MODEL_INFO <- model_linearLogistic
DISTANCE <- logit_diff
eff_denom_linearLogistic <- sapply(1:length(pairRes_linearLogistic_q), function(k) pairRes_linearLogistic_q[[k]]$out$BESTVAL)
mm_nSupp <- 3

mmRes_linearLogistic <- NULL

genCount <- 0; testMaxDD <- 1e20
while ((genCount < reGen) & (testMaxDD > tolMaxDD)) {
  genCount <- genCount + 1
  # PSO-S-QN algorithm
  out <- DiscrimOD(MODEL_INFO, DISTANCE, mm_nSupp, dsLower = 0, dsUpper = 1,
                   crit_type = "maxmin_fixed_true", MaxMinStdVals = eff_denom_linearLogistic,
                   PSO_INFO = PSO_MAXMIN, LBFGS_INFO = LBFGS_INFO, seed = NULL, verbose = TRUE)

  eqv <- equivalence(ngrid = 100, PSO_RESULT = out, MODEL_INFO = MODEL_INFO, DISTANCE = DISTANCE,
                     dsLower = 0, dsUpper = 1,
                     crit_type = "maxmin_fixed_true", MaxMinStdVals = eff_denom_linearLogistic,
                     PSO_INFO = PSO_MAXMIN, LBFGS_INFO = LBFGS_INFO)

  testMaxDD <- eqv$MAX_DD

  if (genCount > 1) {
    if (testMaxDD < mmRes_linearLogistic$eqv$MAX_DD) { mmRes_linearLogistic <- list(out = out, eqv = eqv) }
  } else {
    mmRes_linearLogistic <- list(out = out, eqv = eqv)
  }
}

# The max-min optimal discrimination design
round(mmRes_linearLogistic$out$BESTDESIGN, 3)
# The best efficiency value
mmRes_linearLogistic$out$BESTVAL
# The computing time
mmRes_linearLogistic$out$CPUTIME

# Draw the plot for checking the equivalence theorem for max-min discrimination design
tmp <- mmRes_linearLogistic
plot(tmp$eqv$Grid_1, tmp$eqv$DirDeriv, type = "l", col = "black",
     xlab = "x", ylab = "Directional Derivative", main = "Max-min Discrimination Design")
abline(h = 0, col = "grey50", lty = 2)
points(tmp$out$BESTDESIGN[,1], rep(0, nrow(tmp$out$BESTDESIGN)), pch = 16)
