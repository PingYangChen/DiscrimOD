library(DiscrimOD)
library(Rcpp)

AF_para_m1 <-c (4.5, -1.5, -2)
# R Function
MODEL_INFO <- list(
  list(model = function(x, p) p[1] + p[2]*exp(x) + p[3]*exp(-x), para = AF_para_m1),
  list(model = function(x, p) p[1] + p[2]*x + p[3]*x^2,
       paraLower = rep(-10, 3),
       paraUpper = rep(10, 3),
       paraInit = rep(0, 3))
)
DISTANCE <- function(xt, xr) (xt - xr)^2


# C++ Funciton
TM_CppCode <- 'double TM_Cpp(SEXP xx, SEXP pp)
{
  // Re-define variable types
  Rcpp::NumericVector x(xx); Rcpp::NumericVector p(pp);
  double val = p(0) + p(1)*std::exp(x(0)) + p(2)*std::exp(-x(0));
  return val;
}'
RM_CppCode <- 'double RM_Cpp(SEXP xx, SEXP pp)
{
  // Re-define variable types
  Rcpp::NumericVector x(xx); Rcpp::NumericVector p(pp);
  double val = p(0) + p(1)*x(0) + p(2)*x(0)*x(0);
  return val;
}'
MODEL_INFO_Cpp <- list(
  list(model = cppFunction(TM_CppCode), para = AF_para_m1),
  list(model = cppFunction(RM_CppCode),
       paraLower = rep(-10, 3),
       paraUpper = rep(10, 3),
       paraInit = AF_para_m1)
)

DISTANCE_CppCode <- '
double DISTANCE_Cpp(SEXP xt, SEXP xr)
{
  // Re-define variable types
  double val_t = Rcpp::as<double>(xt);
  double val_r = Rcpp::as<double>(xr);
  double val = (val_t - val_r)*(val_t - val_r);
  return val;
}
'
DISTANCE_Cpp <- cppFunction(DISTANCE_CppCode)

#
ALG_INFO <- getAlgInfo(nSwarm = 32, maxIter = 100, typePSO = 2)

nSupp <- 4
dsLower <- -1
dsUpper <- 1
MaxMinStdVals = NULL; seed = NULL; verbose = TRUE

out <- DiscrimOD(MODEL_INFO_Cpp, DISTANCE_Cpp, nSupp, dsLower, dsUpper,
                 MaxMinStdVals, ALG_INFO, seed, verbose)
round(out$BESTDESIGN, 3)

nn <- 10
system.time(
for (i in 1:nn) {
DESIGN1 <- cbind(c(-1,-.669,.144,.957), c(.253,.428,.247,.072))
designCriterion(DESIGN1, MODEL_INFO_Cpp, DISTANCE_Cpp, dsLower, dsUpper,
                MaxMinStdVals, ALG_INFO)
})[3]/nn


