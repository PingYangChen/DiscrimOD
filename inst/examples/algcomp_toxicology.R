# Example: A and Fedorov (1975a, b)
library(DiscrimOD)
#gdpath <- "D:/Ping_Yang/Google Drive/PYChen_Statistics_NCKU"
#projPath <- file.path(gdpath, "Researches/2015 min-max optimal discriminating designs")
projPath <- "./2017_PSOQN"
outputPath <- file.path(projPath, "pkgOutput_algComp_toxicology")
if (!dir.exists(outputPath)) { dir.create(outputPath) }

caseName <- "toxicology"

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
tox5 <- function(x, p) p[1]*(p[3] - (p[3] - 1)*exp(-(x/p[2])^p[4]))
tox4 <- function(x, p) p[1]*(p[3] - (p[3] - 1)*exp(-(x/p[2])))
tox3 <- function(x, p) p[1]*exp(-(x/p[2])^p[3])
tox2 <- function(x, p) p[1]*exp(-(x/p[2]))
tox1 <- function(x, p) rep(p[1], length(x))
# Set the nominal values for the first model (null model)
para_tox_5 <- c(4.282, 835.571, 0.739, 3.515)
model_tox <- list(
  list(model = tox5, para = para_tox_5),
  list(model = tox4, paraLower = c(0, 0, 0), paraUpper = c(20, 5000, 1)),
  list(model = tox3, paraLower = c(0, 0, 1), paraUpper = c(20, 5000, 15)),
  list(model = tox2, paraLower = c(0, 0), paraUpper = c(20, 5000)),
  list(model = tox1, paraLower = c(0), paraUpper = c(20))
)
DL <- 0; DU <- 1250
# Create the lists for pairwise discrimination designs
two_model <- list(
  list(model_tox[[1]], model_tox[[2]]),
  list(model_tox[[1]], model_tox[[3]]),
  list(model_tox[[1]], model_tox[[4]]),
  list(model_tox[[1]], model_tox[[5]])
)
# Specify the number fo support points for pairwise discrimination designs
two_nSupp <- c(3, 4, 3, 2)

two_optimal <- list(
  cbind(c(0.000, 468.186, 1064.179), c(0.249, 0.498, 0.253)),
  cbind(c(0.000, 484.213, 963.141, 1250.000), c(0.092, 0.280, 0.407, 0.221)),
  cbind(c(0.000, 468.156, 1064.178), c(0.249, 0.498, 0.253)),
  cbind(c(0.000, 1250.000), c(0.500, 0.500))
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

  effvals <- matrix(0, nRep, 4*2)
  colnames(effvals) <- paste0(rep(c("PSOQN", "NESTEDPSO", "FEDWYNN", "REMES"), 2), rep(c("EFF", "CPU"), each = 4))

  for (iR in 1:nRep) {
    cat(paste0("Case: ", caseName, "; Sub: ", iC, "; Rep: ", iR, "\n"))
    MODEL_INFO <- two_model[[iC]]
    OPTIMAL <- two_optimal[[iC]]
    OPT_VAL <- designCriterion(OPTIMAL, MODEL_INFO, DISTANCE, dsLower = DL, dsUpper = DU, 
                               crit_type = "pair_fixed_true", 
                               PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

    eachRep <- vector("list", 4)
     # PSO-QN
    out_q <- DiscrimOD(MODEL_INFO, DISTANCE, two_nSupp[iC], dsLower = DL, dsUpper = DU,
                       crit_type = "pair_fixed_true",
                       PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO, verbose = TRUE)

    cri_q <- designCriterion(out_q$BESTDESIGN, MODEL_INFO, DISTANCE, dsLower = DL, dsUpper = DU, 
                             crit_type = "pair_fixed_true", 
                             PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_CRIT)

    eff_q <- cri_q$cri_val/OPT_VAL$cri_val

    #algCompRes[[iC]][[iR]][[1]] <- list(RES = out_q, EFF = eff_q)
    eachRep[[1]] <- list(RES = out_q, EFF = eff_q)

    # NestedPSO
    out_n <- DiscrimOD(MODEL_INFO, DISTANCE, two_nSupp[iC], dsLower = DL, dsUpper = DU,
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

    # Remes
    out_r <- DiscrimUnifApproxT(MODEL_INFO, two_nSupp[iC], dsLower = DL, dsUpper = DU,
                                REMES_MAXIT = nIter, REMES_FreeRun = 1.0, REMES_EPS = 1e-2,
                                LBFGS_INFO = LBFGS_INFO, verbose = TRUE)
    eff_r <- -99
    if (is.numeric(out_r$BESTVAL)) {
      cri_r <- designCriterion(out_r$BESTDESIGN, MODEL_INFO, DISTANCE, dsLower = DL, dsUpper = DU, 
                             crit_type = "pair_fixed_true", 
                             PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_CRIT)

      eff_r <- cri_r$cri_val/OPT_VAL$cri_val
    }
    #algCompRes[[iC]][[iR]][[4]] <- list(RES = out_r, EFF = eff_r)
    eachRep[[4]] <- list(RES = out_r, EFF = eff_r)

    effvals[iR,] <- c(eachRep[[1]]$EFF, eachRep[[2]]$EFF, eachRep[[3]]$EFF, eachRep[[4]]$EFF,
                      eachRep[[1]]$RES$CPUTIME, eachRep[[2]]$RES$CPUTIME, eachRep[[3]]$RES$CPUTIME, eachRep[[4]]$RES$CPUTIME)
  }
  # SAVE RESULT
  #tmp <- algCompRes[[iC]]
  #effvals <- matrix(0, nRep, 4*2)
  #colnames(effvals) <- paste0(rep(c("PSOQN", "NESTEDPSO", "FEDWYNN", "REMES"), 2), rep(c("EFF", "CPU"), each = 4))
  #for (iR in 1:nRep) {
  #  effvals[iR,] <- c(tmp[[iR]][[1]]$EFF, tmp[[iR]][[2]]$EFF, tmp[[iR]][[3]]$EFF, tmp[[iR]][[4]]$EFF,
  #                    tmp[[iR]][[1]]$RES$CPUTIME, tmp[[iR]][[2]]$RES$CPUTIME, tmp[[iR]][[3]]$RES$CPUTIME, tmp[[iR]][[4]]$RES$CPUTIME)    
  #}
  write.csv(effvals, file.path(outputPath, paste0("algComp_Summary_", caseName, "_", iC, ".csv")))
}

save.image(file.path(outputPath, paste0("algComp_", caseName, ".Rdata")))
