# Example: A and Fedorov (1975a, b)
library(DiscrimOD)
gdpath <- "D:/Ping_Yang/Google Drive/PYChen_Statistics_NCKU"
projPath <- file.path(gdpath, "Researches/2015 min-max optimal discriminating designs")
#projPath <- "./2017_PSOQN"
outputPath <- file.path(projPath, "pkgOutput_algComp_logistic")
if (!dir.exists(outputPath)) { dir.create(outputPath) }

caseName <- "logistic"

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
linlogi4 <- cppFunction('
  Rcpp::NumericVector enzyme2(Rcpp::NumericMatrix x, Rcpp::NumericVector p) {
    Rcpp::NumericVector eta(x.nrow());
    for (int i = 0; i < x.nrow(); i++) { 
      eta(i) = p(0) + p(1)*x(i,0) + p(2)*x(i,0)*x(i,0); }
    return eta; 
}')
linlogi3 <- cppFunction('
  Rcpp::NumericVector enzyme2(Rcpp::NumericMatrix x, Rcpp::NumericVector p) {
    Rcpp::NumericVector eta(x.nrow());
    for (int i = 0; i < x.nrow(); i++) { eta(i) = x(i,0)*(p(0) + p(1)*x(i,0)); }
    return eta; 
}')
linlogi2 <- cppFunction('
  Rcpp::NumericVector enzyme2(Rcpp::NumericMatrix x, Rcpp::NumericVector p) {
    Rcpp::NumericVector eta(x.nrow());
    for (int i = 0; i < x.nrow(); i++) { eta(i) = p(0) + p(1)*x(i,0); }
    return eta; 
}')
linlogi1 <- cppFunction('
  Rcpp::NumericVector enzyme2(Rcpp::NumericMatrix x, Rcpp::NumericVector p) {
    Rcpp::NumericVector eta(x.nrow());
    for (int i = 0; i < x.nrow(); i++) { eta(i) = p(0)*x(i,0); }
    return eta; 
}')

# Set the nominal values for the first model (null model)
para_linlogi_4 <- c(1, 1, 1)
# Create the model list
model_linearLogistic <- list(
  list(model = linlogi4, para = para_linlogi_4),
  list(model = linlogi3, paraLower = c(-10, -10), paraUpper = c(10, 10)),
  list(model = linlogi2, paraLower = c(-10, -10), paraUpper = c(10, 10)),
  list(model = linlogi1, paraLower = c(-10), paraUpper = c(10))
)

DL <- 0; DU <- 1
# Create the lists for pairwise discrimination designs
two_model <- list(
  list(model_linearLogistic[[1]], model_linearLogistic[[2]]),
  list(model_linearLogistic[[1]], model_linearLogistic[[3]]),
  list(model_linearLogistic[[1]], model_linearLogistic[[4]])
)

# Specify the number fo support points for pairwise discrimination designs
two_nSupp <- c(3, 3, 3)

#
two_optimal <- list(

)

# Set distance function
logit_diff_cpp <- cppFunction('
  Rcpp::NumericVector logit_diff(Rcpp::NumericVector xt, Rcpp::NumericVector xr) {
    Rcpp::NumericVector div(xt.size());
    double et, er, mt, mr;
    for (int i = 0; i < xt.size(); i++) {
      et = std::exp(xt(i)); mt = et/(1.0 + et);
      er = std::exp(xr(i)); mr = er/(1.0 + er);
      div(i) = mt*std::log(mt/mr) + (1.0 - mt)*std::log((1.0 - mt)/(1.0 - mr));
    }
    return div;
}')


# Start for each pairwise discrimination design
DISTANCE <- logit_diff_cpp
for (iC in 1:length(two_model)) {
  iC <- 3

  effvals <- matrix(0, nRep, 3*2)
  colnames(effvals) <- paste0(rep(c("PSOQN", "NESTEDPSO", "FEDWYNN"), 2), rep(c("EFF", "CPU"), each = 3))

  for (iR in 1:nRep) {
    cat(paste0("Case: ", caseName, "; Sub: ", iC, "; Rep: ", iR, "\n"))
    MODEL_INFO <- two_model[[iC]]
    OPTIMAL <- two_optimal[[iC]]
    OPT_VAL <- designCriterion(OPTIMAL, MODEL_INFO, DISTANCE, dsLower = DL, dsUpper = DU, 
                               crit_type = "pair_fixed_true", 
                               PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_CRIT)

    eachRep <- vector("list", 3)
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

    # No Remes
    effvals[iR,] <- c(eachRep[[1]]$EFF, eachRep[[2]]$EFF, eachRep[[3]]$EFF, 
                      eachRep[[1]]$RES$CPUTIME, eachRep[[2]]$RES$CPUTIME, eachRep[[3]]$RES$CPUTIME)
  }
  # SAVE RESULT
  write.csv(effvals, file.path(outputPath, paste0("algComp_Summary_", caseName, "_", iC, ".csv")), row.names = FALSE, quote = FALSE)
}

#save.image(file.path(outputPath, paste0("algComp_", caseName, ".Rdata")))
