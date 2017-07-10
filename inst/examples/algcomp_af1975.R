# Example: A and Fedorov (1975a, b)
library(DiscrimOD)
gdpath <- "D:/Ping_Yang/Google Drive/PYChen_Statistics_NCKU"
projPath <- file.path(gdpath, "Researches/2015 min-max optimal discriminating designs")
#projPath <- "./2017_PSOQN"
outputPath <- file.path(projPath, "pkgOutput_algComp_af1975")
if (!dir.exists(outputPath)) { dir.create(outputPath) }

caseName <- "af1975"

nIter <- 200; nRep <- 50
# Set PSO options for pariwise discrimination design cases
PSO_INFO <- getPSOInfo(nSwarm = c(32, 32), maxIter = c(nIter, 100))
# Set L-BFGS algorithm options
LBFGS_INFO <- getLBFGSInfo(LBFGS_RETRY = 2)
# Set a NOT-RUN L-BFGS algorithm for trying NestedPSO (for fun)
LBFGS_NOTRUN <- getLBFGSInfo(IF_INNER_LBFGS = FALSE)
# Set Fedorov-Wynn options
FED_INFO <- getFEDInfo(FED_MAXIT = nIter, FED_TRIM = 3, FED_TRIM_EPS = 1e-2,
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

DL <- -1; DU <- 1
# Create the lists for pairwise discrimination designs
two_model <- list(
  list(model_af1975[[1]], model_af1975[[2]]),
  list(model_af1975[[1]], model_af1975[[3]])
)
# Specify the number fo support points for pairwise discrimination designs
two_nSupp <- c(4, 5)

two_optimal <- list(
  cbind(c(-1.0000, -0.6690, 0.1440, 0.9570), c(0.2530,  0.4280, 0.2470, 0.0720)),
  cbind(c(-1.0000, -0.7405, -0.1044, 0.6340, 1.0000), c(0.1916, 0.3228, 0.2274, 0.1772, 0.0810))
)

# Set distance function
sq_diff <- function(xt, xr) (xt - xr)^2

# Get optimal designs with two different approaches
algCompRes <- lapply(1:length(two_model), function(i) {
  lapply(1:nRep, function(k) {
    a <- vector("list", 4); names(a) <- c("PSOQN", "NESTEDPSO", "FEDWYNN", "REMES")
    a
  })
})

# Start for each pairwise discrimination design
DISTANCE <- sq_diff
for (iC in 1:length(two_model)) {

  cat(paste0("Case: ", caseName, " ;Sub: ", iC, "\n"))

  for (iR in 1:nRep) {

    MODEL_INFO <- two_model[[iC]]
    OPTIMAL <- two_optimal[[iC]]
    OPT_VAL <- designCriterion(OPTIMAL, MODEL_INFO, DISTANCE, dsLower = DL, dsUpper = DU,
                               crit_type = "pair_fixed_true",
                               PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

    # PSO-QN
    out_q <- DiscrimOD(MODEL_INFO, DISTANCE, two_nSupp[iC], dsLower = DL, dsUpper = DU,
                       crit_type = "pair_fixed_true",
                       PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO, verbose = TRUE)
    eff_q <- out_q$BESTVAL/OPT_VAL$cri_val

    algCompRes[[iC]][[iR]][[1]] <- list(RES = out_q, EFF = eff_q)

    # NestedPSO
    out_n <- DiscrimOD(MODEL_INFO, DISTANCE, two_nSupp[iC], dsLower = DL, dsUpper = DU,
                       crit_type = "pair_fixed_true",
                       PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_NOTRUN, verbose = TRUE)

    eff_n <- out_n$BESTVAL/OPT_VAL$cri_val

    algCompRes[[iC]][[iR]][[2]] <- list(RES = out_n, EFF = eff_n)

    # Fedorov-Wynn
    out_f <- DiscrimFedWynn(MODEL_INFO, DISTANCE, dsLower = DL, dsUpper = DU,
                            FED_INFO = FED_INFO, LBFGS_INFO = LBFGS_INFO, verbose = TRUE)

    eff_f <- out_f$BESTVAL/OPT_VAL$cri_val


    algCompRes[[iC]][[iR]][[3]] <- list(RES = out_f, EFF = eff_f)

    # Remes
    out_r <- DiscrimUnifApproxT(MODEL_INFO, two_nSupp[iC], dsLower = DL, dsUpper = DU,
                                REMES_MAXIT = nIter, REMES_FreeRun = 1.0, REMES_EPS = 1e-2,
                                LBFGS_INFO = LBFGS_INFO, verbose = TRUE)
    eff_r <- -99
    if (is.numeric(out_r$BESTVAL)) {
      eff_r <- out_r$BESTVAL/OPT_VAL$cri_val
    }
    algCompRes[[iC]][[iR]][[4]] <- list(RES = out_r, EFF = eff_r)
  }
  # SAVE RESULT
  tmp <- algCompRes[[iC]]
  effvals <- matrix(0, nRep, 4*2)
  colnames(effvals) <- paste0(rep(c("PSOQN", "NESTEDPSO", "FEDWYNN", "REMES"), 2), rep(c("EFF", "CPU"), each = 4))
  for (iR in 1:nRep) {
    effvals[iR,] <- c(tmp[[iR]][[1]]$EFF, tmp[[iR]][[2]]$EFF, tmp[[iR]][[3]]$EFF, tmp[[iR]][[4]]$EFF,
                      tmp[[iR]][[1]]$RES$CPUTIME, tmp[[iR]][[2]]$RES$CPUTIME, tmp[[iR]][[3]]$RES$CPUTIME, tmp[[iR]][[4]]$RES$CPUTIME)
  }
  write.csv(effvals, file.path(outputPath, paste0("algComp_Summary_", caseName, "_", iC, ".csv")))
}

save.image(file.path(outputPath, paste0("algComp_", caseName, ".Rdata")))
