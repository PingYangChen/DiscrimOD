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
FED_INFO <- getFEDInfo(FED_MAXIT = 100, FED_TRIM = 5, FED_TRIM_EPS = 1e-2,
                       freeRun = 1.0, FED_EPS = 1e-6, FED_ALPHA_GRID = 20)


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
nRep <- 1

# Get optimal designs with two different approaches
algCompRes <- lapply(1:nRep, function(k) {
  a <- vector("list", 4); names(a) <- c("PSOQN", "NESTEDPSO", "FEDWYNN", "REMES")
  a
})

# Start for each pairwise discrimination design
for (iC in 1:nRep) {

  MODEL_INFO <- two_model_af1975[[1]]
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

  out_f$BESTVAL/out_q$BESTVAL

  pairRes_af1975ex_f[[iC]] <- out_f

  # Remes
  out_r <- DiscrimUnifApprox(MODEL_INFO, DISTANCE, two_nSupp[1], dsLower = -1, dsUpper = 1,
                             REMES_MAXIT = 100, REMES_FreeRun = 1.0, REMES_EPS = 1e-2,
                             LBFGS_INFO = LBFGS_INFO)
  out_r
  out_r$BESTVAL/out_q$BESTVAL

  pairRes_af1975ex_r[[iC]] <- out_r

}





