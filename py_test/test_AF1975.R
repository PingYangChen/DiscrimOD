library(DiscrimOD)
library(Rcpp)

AF_para_m1 <-c (4.5, -1.5, -2)
# R Function
MODEL_INFO <- list(
  list(model = function(x, p) p[1] + p[2]*exp(x) + p[3]*exp(-x), para = AF_para_m1),
  list(model = function(x, p) p[1] + p[2]*x + p[3]*x^2,
       paraLower = rep(-10, 3),
       paraUpper = rep(10, 3),
       paraInit = runif(3, -10,10))
)
DISTANCE <- function(xt, xr) (xt - xr)^2


# C++ Funciton
TM_CppCode <- 'Rcpp::NumericVector TM_Cpp(SEXP xx, SEXP pp)
{
  // Re-define variable types
  Rcpp::NumericVector x(xx); Rcpp::NumericVector p(pp);
  Rcpp::NumericVector val(x.size());
  for (int i = 0; i < x.size(); i++) val(i) = p(0) + p(1)*std::exp(x(i)) + p(2)*std::exp(-x(i));
  return val;
}'
RM_CppCode <- 'Rcpp::NumericVector RM_Cpp(SEXP xx, SEXP pp)
{
  // Re-define variable types
  Rcpp::NumericVector x(xx); Rcpp::NumericVector p(pp);
  Rcpp::NumericVector val(x.size());
  for (int i = 0; i < x.size(); i++) val(i) = p(0) + p(1)*x(i) + p(2)*(x(i)*x(i));
  return val;
}'
MODEL_INFO_Cpp <- list(
  list(model = cppFunction(TM_CppCode), para = AF_para_m1),
  list(model = cppFunction(RM_CppCode),
       paraLower = rep(-10, 3),
       paraUpper = rep(10, 3),
       paraInit = runif(3,-10,10))
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
ALG_INFO <- getAlgInfo(nSwarm = 32, maxIter = 100, typePSO = 0,
                       LBFGS_RETRY = 3,
                       FVAL_EPS = 1e-6, GRAD_EPS = 1e-9, LINESEARCH_C = 1e-4)

nSupp <- 4
dsLower <- -1
dsUpper <- 1
MaxMinStdVals = NULL; seed = NULL; verbose = TRUE

out <- DiscrimOD(MODEL_INFO_Cpp, DISTANCE_Cpp, nSupp, dsLower, dsUpper,
                 MaxMinStdVals, ALG_INFO, seed, verbose)
round(out$BESTDESIGN, 3)



nn <- 100
aa <- numeric(nn)
system.time(
for (i in 1:nn) {
DESIGN1 <- out$BESTDESIGN
#uu <- runif(4); DESIGN1 <- cbind(runif(4, -1, 1), uu/sum(uu))
#DESIGN1 <- cbind(c(-1,-.669,.144,.957), c(.253,.428,.247,.072))
aa[i] <- designCriterion(DESIGN1, MODEL_INFO_Cpp, DISTANCE_Cpp, dsLower, dsUpper,
                         MaxMinStdVals, ALG_INFO)$cri_val
})[3]/nn

aa
sum(aa < -1)

