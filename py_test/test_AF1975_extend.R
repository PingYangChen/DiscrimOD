library(DiscrimOD)
library(Rcpp); library(RcppArmadillo); #library(assertthat)

#path <- "D:/Ping_Yang/Google Drive/PYChen_Statistics_NCKU/Researches/2015 min-max optimal discriminating designs/DiscrimOD"
#cppPath <- file.path(path, "src")
#rPath <- file.path(path, "R")
#sourceCpp(file.path(cppPath, "cppPSO.cpp"))
#source(file.path(rPath, "DiscrimOD_assist.R"))
#source(file.path(rPath, "DiscrimOD_tools.R"))
#source(file.path(rPath, "DiscrimOD.R"))

## Define C++ functions
# Competitive 
m1_cpp <- cppFunction(
  'Rcpp::NumericVector m1(SEXP xx, SEXP pp) {
    Rcpp::NumericVector x = Rcpp::as<Rcpp::NumericVector>(xx);
    Rcpp::NumericVector p = Rcpp::as<Rcpp::NumericVector>(pp);
    Rcpp::NumericVector eta(x.size());
    for (int i = 0; i < x.size(); i++) {
      eta[i] = p[0] + p[1]*std::exp(p[2]*x[i]) + p[3]*std::exp(-p[4]*x[i]);
    }
    return eta;
  }')
# Uncompetitive 
m2_cpp <- cppFunction(
  'Rcpp::NumericVector m2(SEXP xx, SEXP pp) {
    Rcpp::NumericVector x = Rcpp::as<Rcpp::NumericVector>(xx);
    Rcpp::NumericVector p = Rcpp::as<Rcpp::NumericVector>(pp);
    Rcpp::NumericVector eta(x.size());
    for (int i = 0; i < x.size(); i++) {
      eta[i] = p[0] + p[1]*x[i] + p(2)*x[i]*x[i] + p(3)*x[i]*x[i]*x[i];
    }
    return eta;
  }')

m1_para <- c(4.5, -1.5, 0.5, -2.0, 0.5)
modelList <- list(
  list(model = m1_cpp, para = m1_para),
  list(model = m2_cpp,
       paraLower = rep(-10, 4),
       paraUpper = rep(10, 4),
       paraInit = rep(0, 4))
)


dist_cpp <- cppFunction('Rcpp::NumericVector dist(SEXP xt, SEXP xr) {
  // Re-define variable types
  Rcpp::NumericVector val_t = Rcpp::as<Rcpp::NumericVector>(xt);
  Rcpp::NumericVector val_r = Rcpp::as<Rcpp::NumericVector>(xr);
  Rcpp::NumericVector div(val_t.size());
  for (int i = 0; i < val_t.size(); i++) {
    div[i] = (val_t[i] - val_r[i])*(val_t[i] - val_r[i]);
  }
  return div;
}')

ALG_INFO <- getAlgInfo(nSwarm = 16, maxIter = 100,
                       LBFGS_RETRY = 6,
                       FVAL_EPS = 0, GRAD_EPS = 1e-10,
                       LINESEARCH_MAX = 1, LINESEARCH_ARMIJO = 0.1)

dsLower <- -1
dsUpper <- 1
nSupp <- 5

out <- DiscrimOD(modelList, dist_cpp,
                 nSupp, dsLower, dsUpper,
                 crit_type = "pair_fixed_true",
                 MaxMinStdVals = NULL, ALG_INFO, seed = NULL, verbose = TRUE)
out$BESTVAL
round(out$BESTDESIGN, 3)

eqv <- equivalence(PSO_RESULT = out, MODEL_INFO = modelList,
                   DISTANCE = dist_cpp, ngrid = 100, crit_type = "pair_fixed_true",
                   dsLower = dsLower, dsUpper = dsUpper, MaxMinStdVals = NULL, 
                   ALG_INFO = ALG_INFO)

plot(eqv$Grid_1, eqv$DirDeriv, type = "l", col = "blue"); abline(h = 0);
points(out$BESTDESIGN[,1], rep(0, nrow(out$BESTDESIGN)), pch = 19)



