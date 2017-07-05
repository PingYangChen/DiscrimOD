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
enzyme2 <- function(x, p) p[1]*x/(p[2] + x) + p[3]*x  # Modified Michaelis-Menten
enzyme1 <- function(x, p) p[1]*x/(p[2] + x)           # Michaelise-Menten

# Set the nominal values for the first model (null model)
para_mmm_2 <- c(1, 1, 1)
# Create the model list
model_mmm <- list(
  list(model = enzyme2, para = para_mmm_2),
  list(model = enzyme1, paraLower = c(-20, -20), paraUpper = c(20, 20))
)

# Define the distance functions
distFuncList <- list(
  sq_diff = { function(xt, xr) (xt - xr)^2 }
  ,
  heter_norm = { function(xt, xr) {
    var_t <- xt^2
    var_r <- xr^2
    (var_t + (xt - xr)^2)/var_r - log(var_t/var_r)
  }},
  log_norm_A = { function(xt, xr) {
    vSQ <- 1.0
    s_t <- log(1.0 + (vSQ/(xt^2)))
    s_r <- log(1.0 + (vSQ/(xr^2)))
    mu_t <- log(xt) - s_t
    mu_r <- log(xr) - s_r
    0.5*log(s_r/s_t) - (s_r - s_t - (mu_r - mu_t)^2)/(2*s_r)
  }},
  log_norm_B = { function(xt, xr) {
    sigsq <- 1.0
    var_t <- (exp(sigsq) - 1.0)*(xt^2)
    var_r <- (exp(sigsq) - 1.0)*(xr^2)
    mu_t <- log(xt) - 0.5*log(1.0 + (var_t/(xt^2)))
    mu_r <- log(xr) - 0.5*log(1.0 + (var_r/(xr^2)))
    ((mu_r - mu_t)^2)/(2*sigsq)
  }},
  log_norm_C = { function(xt, xr) {
    c <- 1; d <- 1
    var_t <- d*exp(c*xt)
    var_r <- d*exp(c*xr)
    s_t <- log(1.0 + (var_t/(xt^2)))
    s_r <- log(1.0 + (var_r/(xr^2)))
    mu_t <- log(xt) - s_t
    mu_r <- log(xr) - s_r
    0.5*log(s_r/s_t) - (s_r - s_t - (mu_r - mu_t)*(mu_r - mu_t))/(2*s_r)
  }},
  logit_diff = { function(xt, xr) {
    exp_t <- exp(xt)
    exp_r <- exp(xr)
    mu_t <- exp_t/(1 + exp_t)
    mu_r <- exp_r/(1 + exp_r)
    mu_t*(log(mu_t) - log(mu_r)) + (1 - mu_t)*(log(1.0 - mu_t) - log(1.0 - mu_r))
  }},
  gamma_diff = { function(xt, xr) log(xr/xt) + (xt - xr)/xr }
)

# Set the number of regeneration when the resulting design is not optimal
reGen <- 4; tolMaxDD <- 1e-6;

# Get T-optimal designs with two different approaches
pairRes_mmm_q <- pairRes_mmm_n <- vector("list", length(distFuncList))
names(pairRes_mmm_q) <- names(distFuncList)
names(pairRes_mmm_n) <- names(distFuncList)

# Start for each distance function
for (iC in 1:length(distFuncList)) {

  MODEL_INFO <- model_mmm
  DISTANCE <- distFuncList[[iC]]

  genCount <- 0; testMaxDD <- 1e20
  while ((genCount < reGen) & (testMaxDD > tolMaxDD)) {
    genCount <- genCount + 1
    # Use PSO-QN algorithm to find the optimal deisgn
    out_q <- DiscrimOD(MODEL_INFO, DISTANCE, 3, dsLower = 0.1, dsUpper = 5,
                       crit_type = "pair_fixed_true",
                       PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)
    # equivalence theorem
    eqv_q <- equivalence(ngrid = 100, PSO_RESULT = out_q, MODEL_INFO = MODEL_INFO, DISTANCE = DISTANCE,
                         dsLower = 0.1, dsUpper = 5,
                         crit_type = "pair_fixed_true",
                         PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

    testMaxDD <- eqv_q$MAX_DD
    # If not optimal, regenerate
    if (genCount > 1) {
      if (testMaxDD < pairRes_mmm_q[[iC]]$eqv$MAX_DD) { pairRes_mmm_q[[iC]] <- list(out = out_q, eqv = eqv_q) }
    } else {
      pairRes_mmm_q[[iC]] <- list(out = out_q, eqv = eqv_q)
    }
  }

  # NestedPSO
  genCount <- 0; testMaxDD <- 1e20
  while ((genCount < reGen) & (testMaxDD > tolMaxDD)) {
    genCount <- genCount + 1
    # Use NestedPSO algorithm to find the optimal deisgn
    out_n <- DiscrimOD(MODEL_INFO, DISTANCE, 3, dsLower = 0.1, dsUpper = 5,
                       crit_type = "pair_fixed_true",
                       PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_NOTRUN)
    # equivalence theorem
    eqv_n <- equivalence(ngrid = 100, PSO_RESULT = out_n, MODEL_INFO = MODEL_INFO, DISTANCE = DISTANCE,
                         dsLower = 0.1, dsUpper = 5,
                         crit_type = "pair_fixed_true",
                         PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)
    pairRes_mmm_n[[iC]] <- list(out = out_n, eqv = eqv_n)

    testMaxDD <- eqv_n$MAX_DD
    # If not optimal, regenerate
    if (genCount > 1) {
      if (testMaxDD < pairRes_mmm_n[[iC]]$eqv$MAX_DD) { pairRes_mmm_n[[iC]] <- list(out = out_n, eqv = eqv_n) }
    } else {
      pairRes_mmm_n[[iC]] <- list(out = out_n, eqv = eqv_n)
    }
  }
}


# Get results from the PSO-QN algorithm
for (iC in 1:length(distFuncList)) {
  tmp <- pairRes_mmm_q[[iC]]
  # The max-min optimal discrimination design
  print(round(tmp$out$BESTDESIGN, 3))
  # The best efficiency value
  print(tmp$out$BESTVAL)
  # The computing time
  print(tmp$out$CPUTIME)
}

# Draw the plot for checking the equivalence theorem
par(mfrow = c(4, 2))
for (iC in 1:length(distFuncList)) {
  tmp <- pairRes_mmm_q[[iC]]
  plot(tmp$eqv$Grid_1, tmp$eqv$DirDeriv, type = "l", col = "black",
       xlab = "x", ylab = "Directional Derivative", main = paste0("Result from PSO-QN: ", iC))
  abline(h = 0, col = "grey50", lty = 2)
  points(tmp$out$BESTDESIGN[,1], rep(0, nrow(tmp$out$BESTDESIGN)), pch = 16)
}




