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
m1_cpp <- cppFunction('Rcpp::NumericVector m1(SEXP xx, SEXP pp) {
  Rcpp::NumericMatrix x = Rcpp::as<Rcpp::NumericMatrix>(xx);
  Rcpp::NumericVector p = Rcpp::as<Rcpp::NumericVector>(pp);
  Rcpp::NumericVector eta(x.nrow());
  for (int i = 0; i < x.nrow(); i++) {
    eta[i] = p[0]*x[i,0]/(p[1]*(1.0 + x[i,1]/p[2]) + x[i,0]);
  }
  return eta;
}')
# Uncompetitive 
m2_cpp <- cppFunction('Rcpp::NumericVector m2(SEXP xx, SEXP pp) {
  Rcpp::NumericMatrix x = Rcpp::as<Rcpp::NumericMatrix>(xx);
  Rcpp::NumericVector p = Rcpp::as<Rcpp::NumericVector>(pp);
  Rcpp::NumericVector eta(x.nrow());
  for (int i = 0; i < x.nrow(); i++) {
    eta[i] = p[0]*x[i,0]/(p[1] + x[i,0]*(1.0 + x[i,1]/p[2]));
  }
  return eta;
}')
# Noncompetitive 
m3_cpp <- cppFunction('Rcpp::NumericVector m3(SEXP xx, SEXP pp) {
  Rcpp::NumericMatrix x = Rcpp::as<Rcpp::NumericMatrix>(xx);
  Rcpp::NumericVector p = Rcpp::as<Rcpp::NumericVector>(pp);
  Rcpp::NumericVector eta(x.nrow());
  for (int i = 0; i < x.nrow(); i++) {
    eta[i] = p[0]*x[i,0]/((p[1] + x[i,0])*(1.0 + x[i,1]/p[2]));
  }
  return eta;
}')
# Mixed
m4_cpp <- cppFunction('Rcpp::NumericVector m4(SEXP xx, SEXP pp) {
  Rcpp::NumericMatrix x = Rcpp::as<Rcpp::NumericMatrix>(xx);
  Rcpp::NumericVector p = Rcpp::as<Rcpp::NumericVector>(pp);
  Rcpp::NumericVector eta(x.nrow());
  for (int i = 0; i < x.nrow(); i++) {
    eta[i] = p[0]*x[i,0]/(p[1]*(1.0 + x[i,1]/p[2]) + x[i,0]*(1.0 + x[i,1]/p[3]));
  }
  return eta;
}')

mixed_para <- c(1513, 6.59, 1.063, 1.849) # Vmax, Km, Kc, Ku
modelList <- list(
  list(model = m4_cpp, para = mixed_para),
  list(model = m3_cpp,
       paraLower = rep(1e-4, 3),
       paraUpper = c(2500, rep(60, 2)),
       paraInit = mixed_para[1:3]),
  list(model = m2_cpp,
       paraLower = rep(1e-4, 3),
       paraUpper = c(2500, rep(60, 2)),
       paraInit = mixed_para[c(1,2,4)]),
	list(model = m1_cpp,
       paraLower = rep(1e-4, 3),
       paraUpper = c(2500, rep(60, 2)),
       paraInit = mixed_para[1:3])
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

ALG_INFO <- getAlgInfo(nSwarm = 16, maxIter = 200,
                       LBFGS_RETRY = 6,
                       FVAL_EPS = 0, GRAD_EPS = 1e-10,
                       LINESEARCH_MAX = 1, LINESEARCH_ARMIJO = 0.1)

dsLower <- c(0, 0)
dsUpper <- c(100, 100)

modelList_2 <- list(
  list(TWO_M = list(modelList[[1]], modelList[[2]]), nSupp = 4),
  list(TWO_M = list(modelList[[1]], modelList[[3]]), nSupp = 5),
	list(TWO_M = list(modelList[[1]], modelList[[4]]), nSupp = 5)
)


outALL <- vector("list", 2)
for (CaseID in 1:2) {

  #sourceCpp(file.path(cppPath, "cppPSO.cpp"))
  CaseID <- 1
  out <- DiscrimOD(modelList_2[[CaseID]]$TWO_M, dist_cpp,
                   modelList_2[[CaseID]]$nSupp, dsLower, dsUpper,
                   crit_type = "pair_fixed_true",
                   MaxMinStdVals = NULL, ALG_INFO, seed = NULL, verbose = TRUE)
  out$BESTVAL
  round(out$BESTDESIGN, 3)

  eqv <- equivalence(PSO_RESULT = out, MODEL_INFO = modelList_2[[CaseID]]$TWO_M,
                     DISTANCE = dist_cpp, ngrid = 100, crit_type = "pair_fixed_true",
                     dsLower = dsLower, dsUpper = dsUpper, MaxMinStdVals = NULL, 
                     ALG_INFO = ALG_INFO)
	
  #eqv						
  
  image(eqv$Grid_1, eqv$Grid_2, t(eqv$DirDeriv))
  points(out$BESTDESIGN[,1:2], pch = 19)
  
  range(eqv$DirDeriv)
  eqv$DirDeriv[89:100, 1:5]
  
  #plot(eqv$Grid_1, eqv$DirDeriv, type = "l", col = "blue"); abline(h = 0);
  #points(out$BESTDESIGN[,1], rep(0, nrow(out$BESTDESIGN)), pch = 19)
  DESIGN1 <- out$BESTDESIGN
  #DESIGN1 <- cbind(c(0, 468.186, 1064.179), c(.249, .498, .253))
  designCriterion(DESIGN1, modelList_2[[CaseID]]$TWO_M, dist_cpp, dsLower, dsUpper,
                  crit_type = "pair_fixed_true", MaxMinStdVals = NULL, ALG_INFO)

  outALL[[CaseID]] <- list(res = out, eqv = eqv)
}


CaseID <- 1
outALL[[CaseID]]$res
plot(outALL[[CaseID]]$eqv$Grid_1, outALL[[CaseID]]$eqv$DirDeriv, type = "l", col = "blue"); abline(h = 0);
points(outALL[[CaseID]]$res$BESTDESIGN[,1], rep(0, nrow(outALL[[CaseID]]$res$BESTDESIGN)), pch = 19)


vals_MM <- sapply(1:length(outALL), function(i) outALL[[i]]$res$BESTVAL)
#vals_MM <- c(0.1,0.1)
nSupp_MM <- 5
out_MM <- DiscrimOD(MODEL_INFO_Cpp, DISTANCE_Cpp(),
                    nSupp_MM, dsLower, dsUpper, crit_type = "maxmin_fixed_true",
                    MaxMinStdVals = vals_MM, ALG_INFO, seed = NULL, verbose = TRUE)
round(out_MM$BESTDESIGN, 3)

PARA_SET <- designCriterion(out_MM$BESTDESIGN, MODEL_INFO_Cpp, DISTANCE_Cpp(), dsLower, dsUpper,
                crit_type = "maxmin_fixed_true", MaxMinStdVals = vals_MM, ALG_INFO)$theta2

eqv <- equivalence(PSO_RESULT = out_MM, MODEL_INFO = MODEL_INFO_Cpp, DISTANCE = DISTANCE_Cpp(),
                   ngrid = 100, crit_type = "maxmin_fixed_true",
                   dsLower = dsLower, dsUpper = dsUpper,
                   MaxMinStdVals = vals_MM, ALG_INFO = ALG_INFO)

plot(eqv$Grid_1, eqv$DirDeriv, type = "l", col = "blue"); abline(h = 0);
points(out_MM$BESTDESIGN[,1], rep(0, nrow(out_MM$BESTDESIGN)), pch = 19)


