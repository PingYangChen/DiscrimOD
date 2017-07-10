# Example: A and Fedorov (1975a, b)
library(DiscrimOD)
gdpath <- "D:/Ping_Yang/Google Drive/PYChen_Statistics_NCKU"
projPath <- file.path(gdpath, "Researches/2015 min-max optimal discriminating designs")
#projPath <- "./2017_PSOQN"
outputPath <- file.path(projPath, "pkgOutput_algComp_michaelismenten")
if (!dir.exists(outputPath)) { dir.create(outputPath) }

caseName <- "michaelismenten"

nIter <- 200; nRep <- 50
# Set PSO options for pariwise discrimination design cases
PSO_INFO <- getPSOInfo(nSwarm = c(32, 32), maxIter = c(nIter, 100))
# Set L-BFGS algorithm options
LBFGS_INFO <- getLBFGSInfo(LBFGS_RETRY = 2); LBFGS_CRIT <- getLBFGSInfo(LBFGS_RETRY = 8)
# Set a NOT-RUN L-BFGS algorithm for trying NestedPSO (for fun)
LBFGS_NOTRUN <- getLBFGSInfo(IF_INNER_LBFGS = FALSE)
# Set Fedorov-Wynn options
FED_INFO <- getFEDInfo(FED_MAXIT = nIter, FED_TRIM = 3, FED_TRIM_EPS = 1e-2,
                       freeRun = 1.0, FED_EPS = 1e-6, FED_ALPHA_GRID = 20)

# Create competing models
enzyme2 <- function(x, p) p[1]*x/(p[2] + x) + p[3]*x  # Modified Michaelis-Menten 
enzyme1 <- function(x, p) p[1]*x/(p[2] + x)           # Michaelise-Menten

# Set the nominal values for the first model (null model)
para_mmm_2 <- c(1, 1, 1)
model_mmm <- list(
  list(model = enzyme2, para = para_mmm_2),
  list(model = enzyme1, paraLower = c(-20, -20), paraUpper = c(20, 20))
)

DL <- 0.1; DU <- 5.0
# Create the lists for pairwise discrimination designs
two_model <- model_mmm
# Specify the number fo support points for pairwise discrimination designs
two_nSupp <- 3

two_optimal <- list(
  cbind(c(0.1000, 1.5730, 5.0000), c(0.2935, 0.4996, 0.2069)),
  cbind(c(0.1000, 1.5730, 5.0000), c(0.2868, 0.5116, 0.2016))
)

# Set distance function
distFunSet <- list(
  log_norm_B = { function(xt, xr) {
    sigsq <- 1
    var_t <- (exp(sigsq) - 1.0)*(xt^2)
    var_r <- (exp(sigsq) - 1.0)*(xr^2)
    mu_t <- log(xt) - 0.5*log(1.0 + (var_t/(xt^2)))
    mu_r <- log(xr) - 0.5*log(1.0 + (var_r/(xr^2)))
    ((mu_r - mu_t)^2)/(2*sigsq)
  }},
 gamma_diff = { function(xt, xr) log(xr/xt) + (xt - xr)/xr }
)

# Set the number of regeneration when the resulting design is not optimal
nRep <- 50

# Get optimal designs with two different approaches
algCompRes <- lapply(1:length(distFunSet), function(i) {
  lapply(1:nRep, function(k) {
    a <- vector("list", 3); names(a) <- c("PSOQN", "NESTEDPSO", "FEDWYNN")
    a
  })
})

# Start for each pairwise discrimination design
MODEL_INFO <- two_model

for (iC in 1:length(distFunSet)) {

  effvals <- matrix(0, nRep, 3*2)
  colnames(effvals) <- paste0(rep(c("PSOQN", "NESTEDPSO", "FEDWYNN"), 2), rep(c("EFF", "CPU"), each = 3))

  for (iR in 1:nRep) {
    cat(paste0("Case: ", caseName, "; Sub: ", iC, "; Rep: ", iR, "\n"))
    DISTANCE <- distFunSet[[iC]]
    OPTIMAL <- two_optimal[[iC]]

    OPT_VAL <- designCriterion(OPTIMAL, MODEL_INFO, DISTANCE, dsLower = DL, dsUpper = DU, 
                               crit_type = "pair_fixed_true", 
                               PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_CRIT)

    eachRep <- vector("list", 3)
    # PSO-QN
    out_q <- DiscrimOD(MODEL_INFO, DISTANCE, two_nSupp, dsLower = DL, dsUpper = DU,
                       crit_type = "pair_fixed_true",
                       PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO, verbose = TRUE)

    cri_q <- designCriterion(out_q$BESTDESIGN, MODEL_INFO, DISTANCE, dsLower = DL, dsUpper = DU, 
                             crit_type = "pair_fixed_true", 
                             PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_CRIT)

    eff_q <- cri_q$cri_val/OPT_VAL$cri_val

    #algCompRes[[iC]][[iR]][[1]] <- list(RES = out_q, EFF = eff_q)
    eachRep[[1]] <- list(RES = out_q, EFF = eff_q)

    algCompRes[[iC]][[iR]][[1]] <- list(RES = out_q, EFF = eff_q)

    # NestedPSO
    out_n <- DiscrimOD(MODEL_INFO, DISTANCE, two_nSupp, dsLower = DL, dsUpper = DU,
                       crit_type = "pair_fixed_true",
                       PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_NOTRUN, verbose = TRUE)

    cri_n <- designCriterion(out_n$BESTDESIGN, MODEL_INFO, DISTANCE, dsLower = DL, dsUpper = DU, 
                             crit_type = "pair_fixed_true", 
                             PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_CRIT)

    eff_n <- cri_n$cri_val/OPT_VAL$cri_val

    #algCompRes[[iC]][[iR]][[2]] <- list(RES = out_n, EFF = eff_n)
    eachRep[[2]] <- list(RES = out_n, EFF = eff_n)


    # Fedorov-Wynn
    out_f <- DiscrimFedWynn(MODEL_INFO, DISTANCE, dsLower = DL, dsUpper = DU,
                            FED_INFO = FED_INFO, LBFGS_INFO = LBFGS_INFO, verbose = TRUE)

    cri_f <- designCriterion(out_f$BESTDESIGN, MODEL_INFO, DISTANCE, dsLower = DL, dsUpper = DU, 
                             crit_type = "pair_fixed_true", 
                             PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_CRIT)

    eff_f <- cri_f$cri_val/OPT_VAL$cri_val


    #algCompRes[[iC]][[iR]][[3]] <- list(RES = out_f, EFF = eff_f)
    eachRep[[3]] <- list(RES = out_f, EFF = eff_f)

    # No Remes
    effvals[iR,] <- c(eachRep[[1]]$EFF, eachRep[[2]]$EFF, eachRep[[3]]$EFF, 
                      eachRep[[1]]$RES$CPUTIME, eachRep[[2]]$RES$CPUTIME, eachRep[[3]]$RES$CPUTIME)
  }
  # SAVE RESULT
  write.csv(effvals, file.path(outputPath, paste0("algComp_Summary_", caseName, "_", iC, ".csv")))
}

save.image(file.path(outputPath, paste0("algComp_", caseName, ".Rdata")))
