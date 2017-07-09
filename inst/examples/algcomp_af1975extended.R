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
# Set Fedorov-Wynn options
FED_INFO <- getFEDInfo(FED_MAXIT = 200, FED_TRIM = 3, FED_TRIM_EPS = 1e-3,
                       freeRun = 0.5, FED_EPS = 1e-6, FED_ALPHA_GRID = 20)


# Create competing models
af1975ex_1 <- function(x, p) p[1] + p[2]*exp(p[3]*x) + p[4]*exp(-p[5]*x)
af1975ex_2 <- function(x, p) p[1] + p[2]*x + p[3]*x^2 + p[4]*x^3
af1975ex_3 <- function(x, p) p[1] + p[2]*sin(0.5*pi*x) + p[3]*cos(0.5*pi*x) + p[4]*sin(pi*x)

# Set the nominal values for the first model (null model)
para_af1975ex_1 <- c(4.5, -1.5, 0.5, -2, 0.5)
# Create the model list
model_af1975ex <- list(
  list(model = af1975ex_1, para = para_af1975ex_1),
  list(model = af1975ex_2, paraLower = rep(-10, 4), paraUpper = rep(10, 4)),
  list(model = af1975ex_3, paraLower = rep(-10, 4), paraUpper = rep(10, 4))
)

# Create the lists for pairwise discrimination designs
two_model_af1975ex <- list(
  list(model_af1975ex[[1]], model_af1975ex[[2]]),
  list(model_af1975ex[[1]], model_af1975ex[[3]])
)
# Specify the number fo support points for pairwise discrimination designs
two_nSupp <- c(5, 5)

# Set distance function
sq_diff <- function(xt, xr) (xt - xr)^2

# Set the number of regeneration when the resulting design is not optimal
nRep <- 1

# Get optimal designs with two different approaches
pairRes_af1975ex_q <- pairRes_af1975ex_n <- pairRes_af1975ex_f <- pairRes_af1975ex_r <- vector("list", nRep)

# Start for each pairwise discrimination design
for (iC in 1:nRep) {

  MODEL_INFO <- two_model_af1975ex[[1]]
  DISTANCE <- sq_diff
  #DISTANCE <- logit_diff

  # PSO-QN
  out_q <- DiscrimOD(MODEL_INFO, DISTANCE, two_nSupp[1], dsLower = -1, dsUpper = 1,
                     crit_type = "pair_fixed_true",
                     PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

  pairRes_af1975ex_q[[iC]] <- out_q

  # NestedPSO
  out_n <- DiscrimOD(MODEL_INFO, DISTANCE, two_nSupp[1], dsLower = -1, dsUpper = 1,
                     crit_type = "pair_fixed_true",
                     PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_NOTRUN)

  pairRes_af1975ex_n[[iC]] <- out_n

  # Fedorov-Wynn
  out_f <- DiscrimFedWynn(MODEL_INFO, DISTANCE, dsLower = -1, dsUpper = 1,
                          FED_INFO = FED_INFO, LBFGS_INFO = LBFGS_INFO)

  pairRes_af1975ex_f[[iC]] <- out_f

  # Remes
  out_r <- DiscrimUnifApprox(MODEL_INFO, DISTANCE, two_nSupp[1], dsLower = -1, dsUpper = 1,
                             REMES_MAXIT = 50, REMES_FreeRun = 1.0, REMES_EPS = 1e-6,
                             LBFGS_INFO = LBFGS_INFO)

  pairRes_af1975ex_r[[iC]] <- out_r

}





