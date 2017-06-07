library(DiscrimOD)
library(Rcpp); library(RcppArmadillo); #library(assertthat)

#path <- "D:/Ping_Yang/Google Drive/PYChen_Statistics_NCKU/Researches/2015 min-max optimal discriminating designs/DiscrimOD"
#cppPath <- file.path(path, "src")
#rPath <- file.path(path, "R")
#sourceCpp(file.path(cppPath, "cppPSO.cpp"))
#source(file.path(rPath, "DiscrimOD_assist.R"))
#source(file.path(rPath, "DiscrimOD_tools.R"))
#source(file.path(rPath, "DiscrimOD.R"))

## Define R functions
afex1 <- function(x, p) p[1] + p[2]*exp(p[3]*x) + p[4]*exp(-p[5]*x)
afex2 <- function(x, p) p[1] + p[2]*x + p[3]*(x^2) + p[4]*(x^3)
sq_diff <- function(xt, xr) (xt - xr)^2

## Define C++ functions
# AF extension 1
afex1_cpp <- cppFunction(
  'Rcpp::NumericVector afex1(SEXP xx, SEXP pp) {
    arma::rowvec x = Rcpp::as<arma::rowvec>(xx);
    arma::rowvec p = Rcpp::as<arma::rowvec>(pp);
    arma::rowvec eta = p(0) + p(1)*arma::exp(p(2)*x) + p(3)*arma::exp((-1.0)*p(4)*x);
    return Rcpp::wrap(eta);
  }', depends = "RcppArmadillo")
# AF extension 2
afex2_cpp <- cppFunction(
  'Rcpp::NumericVector afex2(SEXP xx, SEXP pp) {
    arma::rowvec x = Rcpp::as<arma::rowvec>(xx);
    arma::rowvec p = Rcpp::as<arma::rowvec>(pp);
    arma::rowvec x2 = x % x;
    arma::rowvec eta = p(0) + p(1)*x + p(2)*x2 + p(3)*(x2 % x);
    return Rcpp::wrap(eta);
  }', depends = "RcppArmadillo")

afex1_para <- c(4.5, -1.5, 0.5, -2.0, 0.5)

afex1_list <- list(
  list(model = afex1, para = afex1_para),
  list(model = afex2, paraLower = rep(-10, 4), paraUpper = rep(10, 4))
)

afex1_list_cpp <- list(
  list(model = afex1_cpp, para = afex1_para),
  list(model = afex2_cpp, paraLower = rep(-10, 4), paraUpper = rep(10, 4))
)


sq_diff_cpp <- cppFunction('Rcpp::NumericVector l2_diff(SEXP xt, SEXP xr) {
  arma::rowvec eta_t = Rcpp::as<arma::rowvec>(xt);
  arma::rowvec eta_r = Rcpp::as<arma::rowvec>(xr);
  arma::rowvec div = (eta_t - eta_r) % (eta_t - eta_r);
  return Rcpp::wrap(div);
}', depends = "RcppArmadillo")

#
PSO_INFO <- getPSOInfo(nSwarm = c(16, 128),
                       maxIter = c(100, 200))
LBFGS_INFO <- getLBFGSInfo(LBFGS_RETRY = 2)

dsLower <- -1
dsUpper <- 1
nSupp <- 5

afex_res <- DiscrimOD(afex1_list, sq_diff, nSupp, dsLower, dsUpper,
                      crit_type = "pair_fixed_true", MaxMinStdVals = NULL,
                      PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO,
                      seed = 100, verbose = TRUE)

afex_res$BESTVAL
round(afex_res$BESTDESIGN, 3)

afex_res_cpp <- DiscrimOD(afex1_list_cpp, sq_diff_cpp, nSupp, dsLower, dsUpper,
                          crit_type = "pair_fixed_true", MaxMinStdVals = NULL,
                          PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO,
                          seed = 100, verbose = TRUE)

afex_res_cpp$BESTVAL
round(afex_res_cpp$BESTDESIGN, 3)

DESIGN1 <- out$BESTDESIGN
dc <- designCriterion(DESIGN1, afex1_list_cpp, sq_diff_cpp, dsLower, dsUpper,
                      crit_type = "pair_fixed_true", MaxMinStdVals = NULL,
                      PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)
dc

curve(afex1_cpp(x, dc$theta2[1,]))
curve(afex2_cpp(x, dc$theta2[2,]), add = T, col = 2)

eqv <- equivalence(PSO_RESULT = out, MODEL_INFO = afex1_list_cpp, DISTANCE = sq_diff_cpp,
                   ngrid = 100, dsLower = dsLower, dsUpper = dsUpper,
                   crit_type = "pair_fixed_true", MaxMinStdVals = NULL,
                   PSO_INFO = PSO_INFO, LBFGS_INFO = LBFGS_INFO)

plot(eqv$Grid_1, eqv$DirDeriv, type = "l", col = "blue"); abline(h = 0);
points(out$BESTDESIGN[,1], rep(0, nrow(out$BESTDESIGN)), pch = 19)



