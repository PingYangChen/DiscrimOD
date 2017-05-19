library(DiscrimOD)
library(Rcpp)

# R Function
m5_par <- c(4.282, 835.571, 0.739, 3.515)
MODEL_IINFO <- list(
  # m5
  list(model = function(x, p) p[1]*(p[3] - (p[3] - 1)*exp(-(x/p[2])^p[4])),
       para = m5_par),
  # m4
  list(model = function(x, p) p[1]*(p[3] - (p[3] - 1)*exp(-(x/p[2]))),
       paraLower = c(0, 0, 0),
       paraUpper = c(10, 10, 1),
       paraInit = m5_par[1:3]),
  # m3
  list(model = function(x, p) p[1]*exp(-(x/p[2])^p[3]),
       paraLower = c(0, 0, 1),
       paraUpper = c(10, 10, 10),
       paraInit = m5_par[c(1, 2, 4)]),
  # m2
  list(model = function(x, p) p[1]*exp(-(x/p[2])),
       paraLower = c(0, 0),
       paraUpper = c(10, 10),
       paraInit = m5_par[1:2]),
  # m1
  list(model = function(x, p) p[1],
       paraLower = c(0),
       paraUpper = c(10),
       paraInit = m5_par[1])
)
DISTANCE <- function(xt, xr) (xt - xr)^2


# C++ Funciton
tox5_CppCode <- 'double tox5_Cpp(SEXP xx, SEXP pp)
{
  // Re-define variable types
  Rcpp::NumericVector x(xx); Rcpp::NumericVector p(pp);
  return p(0)*(p(2) - (p(2) - 1.0)*std::exp(-std::pow(x(0)/p(1), p(3))));
}'
tox4_CppCode <- 'double tox4_Cpp(SEXP xx, SEXP pp)
{
  // Re-define variable types
  Rcpp::NumericVector x(xx); Rcpp::NumericVector p(pp);
  return p(0)*(p(2) - (p(2) - 1.0)*std::exp(-(x(0)/p(1))));
}'
tox3_CppCode <- 'double tox3_Cpp(SEXP xx, SEXP pp)
{
  // Re-define variable types
  Rcpp::NumericVector x(xx); Rcpp::NumericVector p(pp);
  return p(0)*std::exp(-std::pow(x(0)/p(1), p(2)));
}'
tox2_CppCode <- 'double tox2_Cpp(SEXP xx, SEXP pp)
{
  // Re-define variable types
  Rcpp::NumericVector x(xx); Rcpp::NumericVector p(pp);
  return p(0)*std::exp(-(x(0)/p(1)));
}'
tox1_CppCode <- 'double tox1_Cpp(SEXP xx, SEXP pp)
{
  // Re-define variable types
  Rcpp::NumericVector x(xx); Rcpp::NumericVector p(pp);
  return p(0);
}'


MODEL_INFO_Cpp <- list(
  # m5
  list(model = cppFunction(tox5_CppCode),
       para = m5_par),
  # m4
  list(model = cppFunction(tox4_CppCode),
       paraLower = c(0, 0, 0),
       paraUpper = c(10, 10, 1),
       paraInit = m5_par[1:3]),
  # m3
  list(model = cppFunction(tox3_CppCode),
       paraLower = c(0, 0, 1),
       paraUpper = c(10, 10, 10),
       paraInit = m5_par[c(1, 2, 4)]),
  # m2
  list(model = cppFunction(tox2_CppCode),
       paraLower = c(0, 0),
       paraUpper = c(10, 10),
       paraInit = m5_par[1:2]),
  # m1
  list(model = cppFunction(tox1_CppCode),
       paraLower = c(0),
       paraUpper = c(10),
       paraInit = m5_par[1])
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

nSupp <- 3
dsLower <- 0
dsUpper <- 1250

SUB_MODEL_INFO_Cpp <- list(MODEL_INFO_Cpp[[1]], MODEL_INFO_Cpp[[2]])

out <- DiscrimOD(SUB_MODEL_INFO_Cpp, DISTANCE_Cpp, nSupp, dsLower, dsUpper,
                 MaxMinStdVals = NULL, ALG_INFO, seed = NULL, verbose = TRUE)
round(out$BESTDESIGN, 3)

DESIGN1 <- cbind(c(-1, 0, 1), rep(1/3, 3))
designCriterion(out$BESTDESIGN, MODEL_INFO, DISTANCE, dsLower, dsUpper,
                MaxMinStdVals, ALG_INFO)








