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
af1975_1 <- function(x, p) p[1] + p[2]*exp(x) + p[3]*exp(-x)
af1975_2 <- function(x, p) p[1] + p[2]*x + p[3]*x^2
af1975_3 <- function(x, p) p[1] + p[2]*sin(0.5*pi*x) + p[3]*cos(0.5*pi*x) + p[4]*sin(pi*x)

# Set the nominal values for the first model (null model)
para_af1975_1 <- c(4.5, -1.5, -2)
# Create the model list
model_af1975 <- list(
  list(model = af1975_1, para = para_af1975_1),
  list(model = af1975_2, paraLower = rep(-10, 3), paraUpper = rep(10, 3)),
  list(model = af1975_3, paraLower = rep(-10, 4), paraUpper = rep(10, 4))
)

# Create the lists for pairwise discrimination designs
two_model_af1975 <- list(
  list(model_af1975[[1]], model_af1975[[2]]),
  list(model_af1975[[1]], model_af1975[[3]])
)
# Specify the number fo support points for pairwise discrimination designs
two_nSupp <- c(4, 5)

# Set distance function
sq_diff <- function(xt, xr) (xt - xr)^2

# Set the number of regeneration when the resulting design is not optimal
reGen <- 4; tolMaxDD <- 1e-6;

# Get optimal designs with two different approaches
pairRes_af1975_q <- pairRes_af1975_n <- vector("list", length(two_model_af1975))

# Start for each pairwise discrimination design
for (iC in 1:length(two_model_af1975)) {

  MODEL_INFO <- two_model_af1975[[iC]]
  DISTANCE <- sq_diff

  # PSO-QN
  genCount <- 0; testMaxDD <- 1e20
  while ((genCount < reGen) & (testMaxDD > tolMaxDD)) {
    genCount <- genCount + 1
    # Use PSO-QN algorithm to find the optimal deisgn
    out_q <- DiscrimOD(MODEL_INFO, DISTANCE, two_nSupp[iC], dsLower = -1, dsUpper = 1,
                       crit_type = "pair_fixed_true",
                       PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)
    # equivalence theorem
    eqv_q <- equivalence(ngrid = 100, PSO_RESULT = out_q, MODEL_INFO = MODEL_INFO, DISTANCE = DISTANCE,
                         dsLower = -1, dsUpper = 1,
                         crit_type = "pair_fixed_true",
                         PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

    testMaxDD <- eqv_q$MAX_DD
    # If not optimal, regenerate
    if (genCount > 1) {
      if (testMaxDD < pairRes_af1975_q[[iC]]$eqv$MAX_DD) { pairRes_af1975_q[[iC]] <- list(out = out_q, eqv = eqv_q) }
    } else {
      pairRes_af1975_q[[iC]] <- list(out = out_q, eqv = eqv_q)
    }
  }

  # NestedPSO
  genCount <- 0; testMaxDD <- 1e20
  while ((genCount < reGen) & (testMaxDD > tolMaxDD)) {
    genCount <- genCount + 1
    # Use NestedPSO algorithm to find the optimal deisgn
    out_n <- DiscrimOD(MODEL_INFO, DISTANCE, two_nSupp[iC], dsLower = -1, dsUpper = 1,
                       crit_type = "pair_fixed_true",
                       PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_NOTRUN)
    # equivalence theorem
    eqv_n <- equivalence(ngrid = 100, PSO_RESULT = out_n, MODEL_INFO = MODEL_INFO, DISTANCE = DISTANCE,
                         dsLower = -1, dsUpper = 1,
                         crit_type = "pair_fixed_true",
                         PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

    testMaxDD <- eqv_n$MAX_DD
    # If not optimal, regenerate
    if (genCount > 1) {
      if (testMaxDD < pairRes_af1975_n[[iC]]$eqv$MAX_DD) { pairRes_af1975_n[[iC]] <- list(out = out_n, eqv = eqv_n) }
    } else {
      pairRes_af1975_n[[iC]] <- list(out = out_n, eqv = eqv_n)
    }
  }
}

# Get results from the PSO-QN algorithm
for (iC in 1:length(two_model_af1975)) {
  tmp <- pairRes_af1975_q[[iC]]
  # The max-min optimal discrimination design
  print(round(tmp$out$BESTDESIGN, 3))
  # The best efficiency value
  print(tmp$out$BESTVAL)
  # The computing time
  print(tmp$out$CPUTIME)
}

# Draw the plot for checking the equivalence theorem
par(mfrow = c(1, 2))
for (iC in 1:length(two_model_af1975)) {
  tmp <- pairRes_af1975_q[[iC]]
  plot(tmp$eqv$Grid_1, tmp$eqv$DirDeriv, type = "l", col = "black",
       xlab = "x", ylab = "Directional Derivative", main = paste0("Result from PSO-QN: ", iC))
  abline(h = 0, col = "grey50", lty = 2)
  points(tmp$out$BESTDESIGN[,1], rep(0, nrow(tmp$out$BESTDESIGN)), pch = 16)
}

par(mfrow = c(1, 2))
for (iC in 1:length(two_model_af1975)) {
  tmp <- pairRes_af1975_n[[iC]]
  plot(tmp$eqv$Grid_1, tmp$eqv$DirDeriv, type = "l", col = "black",
       xlab = "x", ylab = "Directional Derivative", main = paste0("Result from NestedPSO: ", iC))
  abline(h = 0, col = "grey50", lty = 2)
  points(tmp$out$BESTDESIGN[,1], rep(0, nrow(tmp$out$BESTDESIGN)), pch = 16)
}

par(mfrow = c(1, 1))

# Get max-min T-optimal optimal design
MODEL_INFO <- model_af1975
DISTANCE <- sq_diff
eff_denom_af1975 <- sapply(1:length(pairRes_af1975_q), function(k) pairRes_af1975_q[[k]]$out$BESTVAL)
mm_nSupp <- 5

mmRes_af1975 <- NULL

genCount <- 0; testMaxDD <- 1e20
while ((genCount < reGen) & (testMaxDD > tolMaxDD)) {
  genCount <- genCount + 1
  out <- DiscrimOD(MODEL_INFO, DISTANCE, mm_nSupp, dsLower = -1, dsUpper = 1,
                   crit_type = "maxmin_fixed_true", MaxMinStdVals = eff_denom_af1975,
                   PSO_INFO = PSO_MAXMIN, LBFGS_INFO = LBFGS_INFO, seed = NULL, verbose = TRUE)

  eqv <- equivalence(ngrid = 100, PSO_RESULT = out, MODEL_INFO = MODEL_INFO, DISTANCE = DISTANCE,
                     dsLower = -1, dsUpper = 1,
                     crit_type = "maxmin_fixed_true", MaxMinStdVals = eff_denom_af1975,
                     PSO_INFO = PSO_MAXMIN, LBFGS_INFO = LBFGS_INFO)

  testMaxDD <- eqv$MAX_DD

  if (genCount > 1) {
    if (testMaxDD < mmRes_af1975$eqv$MAX_DD) { mmRes_af1975 <- list(out = out, eqv = eqv) }
  } else {
    mmRes_af1975 <- list(out = out, eqv = eqv)
  }
}

# The max-min optimal discrimination design
round(mmRes_af1975$out$BESTDESIGN, 3)
# The best efficiency value
mmRes_af1975$out$BESTVAL
# The computing time
mmRes_af1975$out$CPUTIME

# Draw the plot for checking the equivalence theorem for max-min discrimination design
tmp <- mmRes_af1975
plot(tmp$eqv$Grid_1, tmp$eqv$DirDeriv, type = "l", col = "black",
     xlab = "x", ylab = "Directional Derivative", main = "Max-min Discrimination Design")
abline(h = 0, col = "grey50", lty = 2)
points(tmp$out$BESTDESIGN[,1], rep(0, nrow(tmp$out$BESTDESIGN)), pch = 16)
