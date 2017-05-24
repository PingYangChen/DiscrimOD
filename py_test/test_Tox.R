library(DiscrimOD)
library(Rcpp)

# R Function
m5_par <- c(4.282, 835.571, 0.739, 3.515)
MODEL_INFO <- list(
  # m5
  list(model = function(x, p) p[1]*(p[3] - (p[3] - 1)*exp(-(x/p[2])^p[4])),
       para = m5_par),
  # m4
  list(model = function(x, p) p[1]*(p[3] - (p[3] - 1)*exp(-(x/p[2]))),
       paraLower = c(0, 800, 0),
       paraUpper = c(10, 5000, 1),
       paraInit = m5_par[1:3]),
  # m3
  list(model = function(x, p) p[1]*exp(-(x/p[2])^p[3]),
       paraLower = c(0, 800, 1),
       paraUpper = c(10, 5000, 15),
       paraInit = m5_par[c(1, 2, 4)]),
  # m2
  list(model = function(x, p) p[1]*exp(-(x/p[2])),
       paraLower = c(0, 800),
       paraUpper = c(10, 5000),
       paraInit = m5_par[1:2]),
  # m1
  list(model = function(x, p) p[1],
       paraLower = c(0),
       paraUpper = c(10),
       paraInit = m5_par[1])
)
DISTANCE <- function(xt, xr) (xt - xr)^2


# C++ Funciton
tox5_CppCode <- 'Rcpp::NumericVector tox5_Cpp(SEXP xx, SEXP pp)
{
  // Re-define variable types
  Rcpp::NumericVector x(xx); Rcpp::NumericVector p(pp);
  Rcpp::NumericVector val(x.size());
  for (int i = 0; i < x.size(); i++) {
    val(i) = p(0)*(p(2) - (p(2) - 1.0)*std::exp(-std::pow(x(i)/p(1), p(3))));
  }
  return val;
}'
tox4_CppCode <- 'Rcpp::NumericVector tox4_Cpp(SEXP xx, SEXP pp)
{
  // Re-define variable types
  Rcpp::NumericVector x(xx); Rcpp::NumericVector p(pp);
  Rcpp::NumericVector val(x.size());
  for (int i = 0; i < x.size(); i++) {
    val(i) = p(0)*(p(2) - (p(2) - 1.0)*std::exp(-(x(i)/p(1))));
  }
  return val;
}'
tox3_CppCode <- 'Rcpp::NumericVector tox3_Cpp(SEXP xx, SEXP pp)
{
  // Re-define variable types
  Rcpp::NumericVector x(xx); Rcpp::NumericVector p(pp);
  Rcpp::NumericVector val(x.size());
  for (int i = 0; i < x.size(); i++) {
    val(i) = p(0)*std::exp(-std::pow(x(i)/p(1), p(2)));
  }
  return val;
}'
tox2_CppCode <- 'Rcpp::NumericVector tox2_Cpp(SEXP xx, SEXP pp)
{
  // Re-define variable types
  Rcpp::NumericVector x(xx); Rcpp::NumericVector p(pp);
  Rcpp::NumericVector val(x.size());
  for (int i = 0; i < x.size(); i++) {
    val(i) = p(0)*std::exp(-(x(i)/p(1)));
  }
  return val;
}'
tox1_CppCode <- 'Rcpp::NumericVector tox1_Cpp(SEXP xx, SEXP pp)
{
  // Re-define variable types
  Rcpp::NumericVector x(xx); Rcpp::NumericVector p(pp);
  Rcpp::NumericVector val(x.size(), p(0));
  return val;
}'


MODEL_INFO_Cpp <- list(
  # m5
  list(model = cppFunction(tox5_CppCode),
       para = m5_par),
  # m4
  list(model = cppFunction(tox4_CppCode),
       paraLower = c(0, 800, 0),
       paraUpper = c(10, 5000, 1),
       paraInit = m5_par[1:3]),
  # m3
  list(model = cppFunction(tox3_CppCode),
       paraLower = c(0, 800, 1),
       paraUpper = c(10, 5000, 15),
       paraInit = m5_par[c(1, 2, 4)]),
  # m2
  list(model = cppFunction(tox2_CppCode),
       paraLower = c(0, 800),
       paraUpper = c(10, 5000),
       paraInit = m5_par[1:2]),
  # m1
  list(model = cppFunction(tox1_CppCode),
       paraLower = c(0),
       paraUpper = c(10),
       paraInit = m5_par[1])
)

DISTANCE_CppCode <- 'Rcpp::NumericVector DISTANCE_Cpp(SEXP xt, SEXP xr)
{
  // Re-define variable types
  Rcpp::NumericVector val_t(xt);
  Rcpp::NumericVector val_r(xr);
  Rcpp::NumericVector val(val_t.size());
  for (int i = 0; i < val.size(); i++) val(i) = (val_t(i) - val_r(i))*(val_t(i) - val_r(i));
  return val;
}
'
DISTANCE_Cpp <- cppFunction(DISTANCE_CppCode)

#
ALG_INFO <- getAlgInfo(nSwarm = 32, maxIter = 200, typePSO = 2,
                       LBFGS_RETRY = 3,
                       FVAL_EPS = 1e-6, GRAD_EPS = 1e-8, LINESEARCH_MAX = 1e20)

nSupp <- 3
dsLower <- 0
dsUpper <- 1250

SUB_MODEL_INFO_Cpp <- list(MODEL_INFO_Cpp[[1]], MODEL_INFO_Cpp[[2]])
SUB_MODEL_INFO <- list(MODEL_INFO[[1]], MODEL_INFO[[2]])

out <- DiscrimOD(SUB_MODEL_INFO_Cpp, DISTANCE_Cpp, nSupp, dsLower, dsUpper,
                 MaxMinStdVals = NULL, ALG_INFO, seed = NULL, verbose = TRUE)
round(out$BESTDESIGN, 3)

eqv <- equivalence(PSO_RESULT = out, MODEL_INFO = SUB_MODEL_INFO, DISTANCE = DISTANCE, ngrid = 100,
                   dsLower = dsLower, dsUpper = dsUpper, MaxMinStdVals = NULL, ALG_INFO = ALG_INFO)

plot(eqv$Grid_1, eqv$DirDeriv, type = "l", col = "blue"); abline(h = 0);
points(out$BESTDESIGN[,1], rep(0, nrow(out$BESTDESIGN)), pch = 19)


DESIGN1 <- out$BESTDESIGN
DESIGN1 <- cbind(c(0, 468.186, 1064.179), c(.249, .498, .253))
designCriterion(DESIGN1, SUB_MODEL_INFO, DISTANCE, dsLower, dsUpper,
                MaxMinStdVals = NULL, ALG_INFO)



