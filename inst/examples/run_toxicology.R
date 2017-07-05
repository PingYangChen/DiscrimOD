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
tox5 <- function(x, p) p[1]*(p[3] - (p[3] - 1)*exp(-(x/p[2])^p[4]))
tox4 <- function(x, p) p[1]*(p[3] - (p[3] - 1)*exp(-(x/p[2])))
tox3 <- function(x, p) p[1]*exp(-(x/p[2])^p[3])
tox2 <- function(x, p) p[1]*exp(-(x/p[2]))
tox1 <- function(x, p) rep(p[1], length(x))

# Set the nominal values for the first model (null model)
para_tox_5 <- c(4.282, 835.571, 0.739, 3.515)
# Create the model list
model_tox <- list(
  list(model = tox5, para = para_tox_5),
  list(model = tox4, paraLower = c(0, 0, 0), paraUpper = c(20, 5000, 1)),
  list(model = tox3, paraLower = c(0, 0, 1), paraUpper = c(20, 5000, 15)),
  list(model = tox2, paraLower = c(0, 0), paraUpper = c(20, 5000)),
  list(model = tox1, paraLower = c(0), paraUpper = c(20))
)

# Create the lists for pairwise discrimination designs
two_model_tox <- list(
  list(model_tox[[1]], model_tox[[2]]),
  list(model_tox[[1]], model_tox[[3]]),
  list(model_tox[[1]], model_tox[[4]]),
  list(model_tox[[1]], model_tox[[5]])
)
# Specify the number fo support points for pairwise discrimination designs
two_nSupp <- c(3, 4, 3, 2)

# Set distance function
# T-optimal design
#sq_diff <- function(xt, xr) (xt - xr)^2
# KL-optimal design
log_norm_B <- function(xt, xr) {
  sigsq <- 0.01
  var_t <- (exp(sigsq) - 1.0)*(xt^2)
  var_r <- (exp(sigsq) - 1.0)*(xr^2)
  mu_t <- log(xt) - 0.5*log(1.0 + (var_t/(xt^2)))
  mu_r <- log(xr) - 0.5*log(1.0 + (var_r/(xr^2)))
  ((mu_r - mu_t)^2)/(2*sigsq)
}


# Set the number of regeneration when the resulting design is not optimal
reGen <- 4; tolMaxDD <- 1e-6;

# Get optimal designs with two different approaches
pairRes_tox_q <- pairRes_tox_n <- vector("list", length(two_model_tox))

# Start for each pairwise discrimination design
for (iC in 1:length(two_model_tox)) {

  MODEL_INFO <- two_model_tox[[iC]]
  #DISTANCE <- sq_diff
  DISTANCE <- log_norm_B

  # PSO-QN
  genCount <- 0; testMaxDD <- 1e20
  while ((genCount < reGen) & (testMaxDD > tolMaxDD)) {
    genCount <- genCount + 1
    # Use PSO-QN algorithm to find the optimal deisgn
    out_q <- DiscrimOD(MODEL_INFO, DISTANCE, two_nSupp[iC], dsLower = 0, dsUpper = 1250,
                       crit_type = "pair_fixed_true",
                       PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)
    # equivalence theorem
    eqv_q <- equivalence(ngrid = 100, PSO_RESULT = out_q, MODEL_INFO = MODEL_INFO, DISTANCE = DISTANCE,
                         dsLower = 0, dsUpper = 1250,
                         crit_type = "pair_fixed_true",
                         PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

    testMaxDD <- eqv_q$MAX_DD
    # If not optimal, regenerate
    if (genCount > 1) {
      if (testMaxDD < pairRes_tox_q[[iC]]$eqv$MAX_DD) { pairRes_tox_q[[iC]] <- list(out = out_q, eqv = eqv_q) }
    } else {
      pairRes_tox_q[[iC]] <- list(out = out_q, eqv = eqv_q)
    }
  }

  # NestedPSO
  genCount <- 0; testMaxDD <- 1e20
  while ((genCount < reGen) & (testMaxDD > tolMaxDD)) {
    genCount <- genCount + 1
    # Use NestedPSO algorithm to find the optimal deisgn
    out_n <- DiscrimOD(MODEL_INFO, DISTANCE, two_nSupp[iC], dsLower = 0, dsUpper = 1250,
                       crit_type = "pair_fixed_true",
                       PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_NOTRUN)
    # equivalence theorem
    eqv_n <- equivalence(ngrid = 100, PSO_RESULT = out_n, MODEL_INFO = MODEL_INFO, DISTANCE = DISTANCE,
                         dsLower = 0, dsUpper = 1250,
                         crit_type = "pair_fixed_true",
                         PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

    testMaxDD <- eqv_n$MAX_DD
    # If not optimal, regenerate
    if (genCount > 1) {
      if (testMaxDD < pairRes_tox_n[[iC]]$eqv$MAX_DD) { pairRes_tox_n[[iC]] <- list(out = out_n, eqv = eqv_n) }
    } else {
      pairRes_tox_n[[iC]] <- list(out = out_n, eqv = eqv_n)
    }
  }
}

# Get results from the PSO-QN algorithm
for (iC in 1:length(two_model_tox)) {
  tmp <- pairRes_tox_q[[iC]]
  # The max-min optimal discrimination design
  print(round(tmp$out$BESTDESIGN, 3))
  # The best efficiency value
  print(tmp$out$BESTVAL)
  # The computing time
  print(tmp$out$CPUTIME)
}

# Draw the plot for checking the equivalence theorem
par(mfrow = c(2, 2))
for (iC in 1:length(two_model_tox)) {
  tmp <- pairRes_tox_q[[iC]]
  plot(tmp$eqv$Grid_1, tmp$eqv$DirDeriv, type = "l", col = "black",
       xlab = "x", ylab = "Directional Derivative", main = paste0("Result from PSO-QN: ", iC))
  abline(h = 0, col = "grey50", lty = 2)
  points(tmp$out$BESTDESIGN[,1], rep(0, nrow(tmp$out$BESTDESIGN)), pch = 16)
}

par(mfrow = c(2, 2))
for (iC in 1:length(two_model_tox)) {
  tmp <- pairRes_tox_n[[iC]]
  plot(tmp$eqv$Grid_1, tmp$eqv$DirDeriv, type = "l", col = "black",
       xlab = "x", ylab = "Directional Derivative", main = paste0("Result from NestedPSO: ", iC))
  abline(h = 0, col = "grey50", lty = 2)
  points(tmp$out$BESTDESIGN[,1], rep(0, nrow(tmp$out$BESTDESIGN)), pch = 16)
}

par(mfrow = c(1, 1))

# Get max-min T-optimal optimal design
MODEL_INFO <- model_tox
#DISTANCE <- sq_diff
DISTANCE <- log_norm_B
eff_denom_tox <- sapply(1:length(pairRes_tox_q), function(k) pairRes_tox_q[[k]]$out$BESTVAL)
mm_nSupp <- 4

mmRes_tox <- NULL

genCount <- 0; testMaxDD <- 1e20
while ((genCount < reGen) & (testMaxDD > tolMaxDD)) {
  genCount <- genCount + 1
  # PSO-S-QN algorithm
  out <- DiscrimOD(MODEL_INFO, DISTANCE, mm_nSupp, dsLower = 0, dsUpper = 1250,
                   crit_type = "maxmin_fixed_true", MaxMinStdVals = eff_denom_tox,
                   PSO_INFO = PSO_MAXMIN, LBFGS_INFO = LBFGS_INFO, seed = NULL, verbose = TRUE)

  eqv <- equivalence(ngrid = 100, PSO_RESULT = out, MODEL_INFO = MODEL_INFO, DISTANCE = DISTANCE,
                     dsLower = 0, dsUpper = 1250,
                     crit_type = "maxmin_fixed_true", MaxMinStdVals = eff_denom_tox,
                     PSO_INFO = PSO_MAXMIN, LBFGS_INFO = LBFGS_INFO)

  testMaxDD <- eqv$MAX_DD

  if (genCount > 1) {
    if (testMaxDD < mmRes_tox$eqv$MAX_DD) { mmRes_tox <- list(out = out, eqv = eqv) }
  } else {
    mmRes_tox <- list(out = out, eqv = eqv)
  }
}

# The max-min optimal discrimination design
round(mmRes_tox$out$BESTDESIGN, 3)
# The best efficiency value
mmRes_tox$out$BESTVAL
# The computing time
mmRes_tox$out$CPUTIME

# Draw the plot for checking the equivalence theorem for max-min discrimination design
tmp <- mmRes_tox
plot(tmp$eqv$Grid_1, tmp$eqv$DirDeriv, type = "l", col = "black",
     xlab = "x", ylab = "Directional Derivative", main = "Max-min Discrimination Design")
abline(h = 0, col = "grey50", lty = 2)
points(tmp$out$BESTDESIGN[,1], rep(0, nrow(tmp$out$BESTDESIGN)), pch = 16)
