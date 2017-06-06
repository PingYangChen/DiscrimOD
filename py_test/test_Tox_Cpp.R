library(DiscrimOD)
library(Rcpp); library(RcppArmadillo)

# Define C++ functions
m5_cpp <- cppFunction(
  'Rcpp::NumericVector m5(SEXP xx, SEXP pp) {
    Rcpp::NumericVector x = Rcpp::as<Rcpp::NumericVector>(xx);
    Rcpp::NumericVector p = Rcpp::as<Rcpp::NumericVector>(pp);
    Rcpp::NumericVector eta(x.size());
    for (int i = 0; i < x.size(); i++) {
      eta[i] = p[0]*(p[2] - (p[2] - 1)*std::exp((-1.0)*std::pow(x[i]/p[1], p[3])));
    }
    return eta;
  }')
m4_cpp <- cppFunction(
  'Rcpp::NumericVector m4(SEXP xx, SEXP pp) {
    Rcpp::NumericVector x = Rcpp::as<Rcpp::NumericVector>(xx);
    Rcpp::NumericVector p = Rcpp::as<Rcpp::NumericVector>(pp);
    Rcpp::NumericVector eta(x.size());
    for (int i = 0; i < x.size(); i++) {
      eta[i] = p[0]*(p[2] - (p[2] - 1)*std::exp((-1.0)*(x[i]/p[1])));
    }
    return eta;
  }')
m3_cpp <- cppFunction(
  'Rcpp::NumericVector m3(SEXP xx, SEXP pp) {
    Rcpp::NumericVector x = Rcpp::as<Rcpp::NumericVector>(xx);
    Rcpp::NumericVector p = Rcpp::as<Rcpp::NumericVector>(pp);
    Rcpp::NumericVector eta(x.size());
    for (int i = 0; i < x.size(); i++) {
      eta[i] = p[0]*std::exp((-1.0)*std::pow(x[i]/p[1], p[2]));
    }
    return eta;
  }')
m2_cpp <- cppFunction(
  'Rcpp::NumericVector m2(SEXP xx, SEXP pp) {
    Rcpp::NumericVector x = Rcpp::as<Rcpp::NumericVector>(xx);
    Rcpp::NumericVector p = Rcpp::as<Rcpp::NumericVector>(pp);
    Rcpp::NumericVector eta(x.size());
    for (int i = 0; i < x.size(); i++) {
      eta[i] = p[0]*std::exp((-1.0)*(x[i]/p[1]));
    }
    return eta;
  }')
m1_cpp <- cppFunction(
  'Rcpp::NumericVector m1(SEXP xx, SEXP pp) {
    Rcpp::NumericVector x = Rcpp::as<Rcpp::NumericVector>(xx);
    Rcpp::NumericVector p = Rcpp::as<Rcpp::NumericVector>(pp);
    Rcpp::NumericVector eta(x.size(), p[0]);
    return eta;
  }')


# R Function
m5_par <- c(4.282, 835.571, 0.739, 3.515)
MODEL_INFO_cpp <- list(
  # m5
  list(model = m5_cpp,
       para = m5_par),
  # m4
  list(model = m4_cpp,
       paraLower = c(0, 0, 0),
       paraUpper = c(20, 5000, 1)),
  # m3
  list(model = m3_cpp,
       paraLower = c(0, 0, 1),
       paraUpper = c(20, 5000, 15)),
  # m2
  list(model = m2_cpp,
       paraLower = c(0, 0),
       paraUpper = c(20, 5000)),
  # m1
  list(model = m1_cpp,
       paraLower = c(0),
       paraUpper = c(20))
)

DISTANCE_cpp <- cppFunction('Rcpp::NumericVector dist(SEXP xt, SEXP xr) {
  // Re-define variable types
  Rcpp::NumericVector val_t = Rcpp::as<Rcpp::NumericVector>(xt);
  Rcpp::NumericVector val_r = Rcpp::as<Rcpp::NumericVector>(xr);
  Rcpp::NumericVector div(val_t.size());
  for (int i = 0; i < val_t.size(); i++) {
    div[i] = (val_t[i] - val_r[i])*(val_t[i] - val_r[i]);
  }
  return div;
}')

#
ALG_INFO <- getAlgInfo(nSwarm = c(16, 128),
                       maxIter = c(100, 200),
                       LBFGS_RETRY = rep(1, 2))

dsLower <- 0
dsUpper <- 1250

TWO_M_INFO_cpp <- list(
  list(TWO_M = list(MODEL_INFO_cpp[[1]], MODEL_INFO_cpp[[2]]), nSupp = 3),
  list(TWO_M = list(MODEL_INFO_cpp[[1]], MODEL_INFO_cpp[[3]]), nSupp = 4),
  list(TWO_M = list(MODEL_INFO_cpp[[1]], MODEL_INFO_cpp[[4]]), nSupp = 3),
  list(TWO_M = list(MODEL_INFO_cpp[[1]], MODEL_INFO_cpp[[5]]), nSupp = 2)
)


outALL <- NULL
outAll <- lapply(1:length(TWO_M_INFO_cpp), function(CaseID) {
  #CaseID <- 1
  out <- DiscrimOD(TWO_M_INFO_cpp[[CaseID]]$TWO_M, DISTANCE_cpp,
                   TWO_M_INFO_cpp[[CaseID]]$nSupp, dsLower, dsUpper,
                   crit_type = "pair_fixed_true",
                   MaxMinStdVals = NULL, IF_INNER_LBFGS = TRUE,
                   ALG_INFO = ALG_INFO, seed = NULL, verbose = TRUE)
  
  round(out$BESTDESIGN, 3)
  out$BESTVAL
  
  eqv <- equivalence(PSO_RESULT = out, MODEL_INFO = TWO_M_INFO_cpp[[CaseID]]$TWO_M,
                     DISTANCE = DISTANCE_cpp, ngrid = 100,
                     dsLower = dsLower, dsUpper = dsUpper, MaxMinStdVals = NULL, 
                     ALG_INFO = ALG_INFO)
  plot(eqv$Grid_1, eqv$DirDeriv, type = "l", col = "blue"); abline(h = 0);
  points(out$BESTDESIGN[,1], rep(0, nrow(out$BESTDESIGN)), pch = 16)
  DESIGN1 <- out$BESTDESIGN
  #DESIGN1 <- cbind(c(0, 468.186, 1064.179), c(.249, .498, .253))
  designCriterion(DESIGN1, TWO_M_INFO_cpp[[CaseID]]$TWO_M, DISTANCE_cpp, dsLower, dsUpper,
                  IF_INNER_LBFGS = TRUE,
                  crit_type = "pair_fixed_true", MaxMinStdVals = NULL, ALG_INFO = ALG_INFO)
  list(res = out, eqv = eqv)
})

CaseID <- 3
plot(outAll[[CaseID]]$eqv$Grid_1, outAll[[CaseID]]$eqv$DirDeriv, type = "l", col = "blue"); abline(h = 0);
points(outAll[[CaseID]]$res$BESTDESIGN[,1], rep(0, nrow(outAll[[CaseID]]$res$BESTDESIGN)), pch = 19)

vals_MM <- sapply(1:length(outAll), function(i) outAll[[i]]$res$BESTVAL)
nSupp_MM <- 4
out <- DiscrimOD(MODEL_INFO_cpp, DISTANCE_cpp,
                 nSupp_MM, dsLower, dsUpper, crit_type = "maxmin_fixed_true",
                 IF_INNER_LBFGS = TRUE,
                 MaxMinStdVals = vals_MM, ALG_INFO, seed = NULL, verbose = TRUE)
round(out$BESTDESIGN, 3)

eqv <- equivalence(PSO_RESULT = out, MODEL_INFO = MODEL_INFO_cpp,
                   DISTANCE = DISTANCE_cpp, ngrid = 100, crit_type = "maxmin_fixed_true",
                   dsLower = dsLower, dsUpper = dsUpper, MaxMinStdVals = vals_MM, 
                   ALG_INFO = ALG_INFO)

plot(eqv$Grid_1, eqv$DirDeriv, type = "l", col = "blue"); abline(h = 0);
points(out$BESTDESIGN[,1], rep(0, nrow(out$BESTDESIGN)), pch = 16)

DESIGN1 <- out$BESTDESIGN
designCriterion(DESIGN1, MODEL_INFO_cpp, DISTANCE_cpp, dsLower, dsUpper,
                IF_INNER_LBFGS = TRUE, crit_type = "maxmin_fixed_true", 
                MaxMinStdVals = vals_MM, ALG_INFO = ALG_INFO)



