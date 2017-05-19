library(DiscrimOD)
library(Rcpp)

# R Function
MODEL_INFO <- list(
  list(model = function(x, p) p[1] + p[2]*exp(x) + p[3]*exp(-x), para = c(4.5, -1.5, -2)),
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
  return p(0) + p(1)*std::exp(x(0)) + p(2)*std::exp(-x(0));
}
'
RM_CppCode <- 'double RM_Cpp(SEXP xx, SEXP pp)
{
  // Re-define variable types
  Rcpp::NumericVector x(xx); Rcpp::NumericVector p(pp);
  return p(0) + p(1)*x(0) + p(2)*x(0)*x(0);
}
'
MODEL_INFO_Cpp <- list(
  list(model = cppFunction(TM_CppCode), para = c(4.5, -1.5, -2)),
  list(model = cppFunction(RM_CppCode),
       paraLower = rep(-10, 3),
       paraUpper = rep(10, 3),
       paraInit = rep(0, 3))
)

DISTANCE_CppCode <- '
double DISTANCE_Cpp(SEXP xt, SEXP xr)
{
  // Re-define variable types
  double val_t = Rcpp::as<double>(xt);
  double val_r = Rcpp::as<double>(xr);
  return (val_t - val_r)*(val_t - val_r);
}
'
DISTANCE_Cpp <- cppFunction(DISTANCE_CppCode)

#
ALG_INFO <- getAlgInfo(nSwarm = 32, maxIter = 200, typePSO = 2)

nSupp <- 4
dsLower <- -1
dsUpper <- 1
MaxMinStdVals = NULL; seed = NULL; verbose = TRUE

out <- DiscrimOD(MODEL_INFO_Cpp, DISTANCE_Cpp, nSupp, dsLower, dsUpper,
                 MaxMinStdVals, ALG_INFO, seed, verbose)
round(out$BESTDESIGN, 3)

DESIGN1 <- cbind(c(-1, 0, 1), rep(1/3, 3))
designCriterion(DESIGN1, MODEL_INFO, DISTANCE, 3, dsLower, dsUpper,
                MaxMinStdVals, ALG_INFO)



