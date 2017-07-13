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
LBFGS_INFO <- getLBFGSInfo(LBFGS_RETRY = 4); LBFGS_CRIT <- getLBFGSInfo(LBFGS_RETRY = 8)
# Set a NOT-RUN L-BFGS algorithm for trying NestedPSO (for fun)
LBFGS_NOTRUN <- getLBFGSInfo(IF_INNER_LBFGS = FALSE)
# Set Fedorov-Wynn options
FED_INFO <- getFEDInfo(FED_MAXIT = nIter, FED_TRIM = 3, FED_TRIM_EPS = 1e-2,
                       freeRun = 1.0, FED_EPS = 1e-6, FED_ALPHA_GRID = 20)

# Create competing models
af1975_1 <- cppFunction('
  Rcpp::NumericVector enzyme2(Rcpp::NumericMatrix x, Rcpp::NumericVector p) {
    Rcpp::NumericVector eta(x.nrow());
    for (int i = 0; i < x.nrow(); i++) { 
      eta(i) = p(0) + p(1)*std::exp(x(i,0)) + p(2)*std::exp(-1.0*x(i,0)); }
    return eta; 
}')
af1975_2 <- cppFunction('
  Rcpp::NumericVector enzyme2(Rcpp::NumericMatrix x, Rcpp::NumericVector p) {
    Rcpp::NumericVector eta(x.nrow());
    for (int i = 0; i < x.nrow(); i++) { 
      eta(i) = p(0) + p(1)*x(i,0) + p(2)*x(i,0)*x(i,0); }
    return eta; 
}')
af1975_3 <- cppFunction('
  Rcpp::NumericVector enzyme2(Rcpp::NumericMatrix x, Rcpp::NumericVector p) {
    double pi_val = 3.141592653589793;
    Rcpp::NumericVector eta(x.nrow());
    for (int i = 0; i < x.nrow(); i++) { 
      eta(i) = p(0) + p(1)*std::sin(0.5*pi_val*x(i,0)) + p(2)*std::cos(0.5*pi_val*x(i,0)) + p(3)*std::sin(pi_val*x(i,0)); }
    return eta; 
}')

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
sq_diff <- cppFunction('
  Rcpp::NumericVector gamma_diff(Rcpp::NumericVector xt, Rcpp::NumericVector xr) {
    Rcpp::NumericVector div(xt.size()); double diff;
    for (int i = 0; i < xt.size(); i++) {
      diff = xt(i) - xr(i); div(i) = diff*diff;
    }
    return div;
}')

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
                               PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_CRIT)

    eachRep <- vector("list", 4)
    # PSO-QN
    cat("PSO-QN\n")
    out_q <- DiscrimOD(MODEL_INFO, DISTANCE, two_nSupp[iC], dsLower = DL, dsUpper = DU,
                       crit_type = "pair_fixed_true",
                       PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO, verbose = TRUE)
    
    cri_q <- designCriterion(out_q$BESTDESIGN, MODEL_INFO, DISTANCE, dsLower = DL, dsUpper = DU, 
                             crit_type = "pair_fixed_true", 
                             PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_CRIT)

    eff_q <- cri_q$cri_val/OPT_VAL$cri_val
    eachRep[[1]] <- list(RES = out_q, EFF = eff_q)

    # NestedPSO
    cat("NestedPSO\n")
    out_n <- DiscrimOD(MODEL_INFO, DISTANCE, two_nSupp[iC], dsLower = DL, dsUpper = DU,
                       crit_type = "pair_fixed_true",
                       PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_NOTRUN, verbose = TRUE)

    cri_n <- designCriterion(out_n$BESTDESIGN, MODEL_INFO, DISTANCE, dsLower = DL, dsUpper = DU, 
                             crit_type = "pair_fixed_true", 
                             PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_CRIT)

    eff_n <- cri_n$cri_val/OPT_VAL$cri_val
    eachRep[[2]] <- list(RES = out_n, EFF = eff_n)

    # Fedorov-Wynn
    cat("Fedorov\n")
    out_f <- DiscrimFedWynn(MODEL_INFO, DISTANCE, dsLower = DL, dsUpper = DU,
                            FED_INFO = FED_INFO, LBFGS_INFO = LBFGS_INFO, verbose = TRUE)

    cri_f <- designCriterion(out_f$BESTDESIGN, MODEL_INFO, DISTANCE, dsLower = DL, dsUpper = DU, 
                             crit_type = "pair_fixed_true", 
                             PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_CRIT)

    eff_f <- cri_f$cri_val/OPT_VAL$cri_val
    eachRep[[3]] <- list(RES = out_f, EFF = eff_f)

    # Remes
    cat("Remes\n")
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
    eachRep[[4]] <- list(RES = out_r, EFF = eff_r)

    effvals[iR,] <- c(eachRep[[1]]$EFF, eachRep[[2]]$EFF, eachRep[[3]]$EFF, eachRep[[4]]$EFF,
                      eachRep[[1]]$RES$CPUTIME, eachRep[[2]]$RES$CPUTIME, eachRep[[3]]$RES$CPUTIME, eachRep[[4]]$RES$CPUTIME)
  }
  # SAVE RESULT
  write.csv(effvals, file.path(outputPath, paste0("algComp_Summary_", caseName, "_", iC, ".csv")), row.names = FALSE, quote = FALSE)
}

#save.image(file.path(outputPath, paste0("algComp_", caseName, ".Rdata")))
