
#include <math.h>
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

//' Distance measurement function for the T-optiaml criterion
//'
//' @param xt vector of design points.
//' @param xr vector of design weight. Its length must be equal to the length of \code{x} and \code{sum(w)} should be 1.
//' @return Fisher information matrix.
//' @references Dette, H., Kiss, C., Bevanda, M., & Bretz, F. (2010). Optimal designs for the EMAX, log-linear and exponential models. Biometrika, 97(2), 513-518.
//' @details
//' 
//' @export
//[[Rcpp::export]]
arma::rowvec T_optimal(SEXP xt, SEXP xr)
{
  // Re-define variable type
  arma::rowvec eta_t = Rcpp::as<arma::rowvec>(xt);
  arma::rowvec eta_r = Rcpp::as<arma::rowvec>(xr);
  if (eta_t.n_elem != eta_r.n_elem) { Rcpp::stop("two mean vectors must have the same length."); }
  return (eta_t - eta_r) % (eta_t - eta_r);
}

//' @export
//[[Rcpp::export]]
arma::rowvec KL_normal_hetero(SEXP xt, SEXP xr)
{
  // Re-define variable type
  arma::rowvec eta_t = Rcpp::as<arma::rowvec>(xt);
  arma::rowvec eta_r = Rcpp::as<arma::rowvec>(xr);
  if (eta_t.n_elem != eta_r.n_elem) { Rcpp::stop("two mean vectors must have the same length."); }
  arma::rowvec var_t = eta_t % eta_t;
  arma::rowvec var_r = eta_r % eta_r;
  return (var_t + arma::pow(eta_t - eta_r, 2))/var_r - arma::log(var_t/var_r);
}

//' @export
//[[Rcpp::export]]
arma::rowvec KL_logistic(SEXP xt, SEXP xr)
{
  // Re-define variable type
  arma::rowvec eta_t = Rcpp::as<arma::rowvec>(xt);
  arma::rowvec eta_r = Rcpp::as<arma::rowvec>(xr);
  if (eta_t.n_elem != eta_r.n_elem) { Rcpp::stop("two mean vectors must have the same length."); }
  arma::rowvec exp_t = arma::exp(eta_t);
  arma::rowvec exp_r = arma::exp(eta_r);
  arma::rowvec mu_t = exp_t/(1.0 + exp_t); 
  arma::rowvec mu_r = exp_r/(1.0 + exp_r);
  return mu_t % (arma::log(mu_t) - arma::log(mu_r)) + (1.0 - mu_t) % (arma::log(1.0 - mu_t) - arma::log(1.0 - mu_r));
}

//' @export
//[[Rcpp::export]]
arma::rowvec KL_lognormal_a(SEXP xt, SEXP xr)
{
  double vSQ = 1.0;
  // Re-define variable type
  arma::rowvec eta_t = Rcpp::as<arma::rowvec>(xt);
  arma::rowvec eta_r = Rcpp::as<arma::rowvec>(xr);
  if (eta_t.n_elem != eta_r.n_elem) { Rcpp::stop("two mean vectors must have the same length."); }
  arma::rowvec s_t = arma::log(1.0 + (vSQ/(eta_t % eta_t)));
  arma::rowvec s_r = arma::log(1.0 + (vSQ/(eta_r % eta_r)));
  arma::rowvec mu_t = arma::log(eta_t) - s_t;
  arma::rowvec mu_r = arma::log(eta_r) - s_r;
  return 0.5*arma::log(s_r/s_t) - (s_r - s_t - (mu_r - mu_t) % (mu_r - mu_t))/(2*s_r);
}

//' @export
//[[Rcpp::export]]
arma::rowvec KL_lognormal_b(SEXP xt, SEXP xr)
{
  double sigsq = 1.0;
  // Re-define variable type
  arma::rowvec eta_t = Rcpp::as<arma::rowvec>(xt);
  arma::rowvec eta_r = Rcpp::as<arma::rowvec>(xr);
  if (eta_t.n_elem != eta_r.n_elem) { Rcpp::stop("two mean vectors must have the same length."); }
  arma::rowvec var_t = (std::exp(sigsq) - 1.0)*(eta_t % eta_t);
  arma::rowvec var_r = (std::exp(sigsq) - 1.0)*(eta_r % eta_r);
  arma::rowvec mu_t = arma::log(eta_t) - 0.5*arma::log(1.0 + (var_t/(eta_t % eta_t)));
  arma::rowvec mu_r = arma::log(eta_r) - 0.5*arma::log(1.0 + (var_r/(eta_r % eta_r)));
  return ((mu_r - mu_t) % (mu_r - mu_t))/(2*sigsq);
}

//' @export
//[[Rcpp::export]]
arma::rowvec KL_lognormal_c(SEXP xt, SEXP xr)
{
  double c = 1.0; double d = 1.0;
  // Re-define variable type
  arma::rowvec eta_t = Rcpp::as<arma::rowvec>(xt);
  arma::rowvec eta_r = Rcpp::as<arma::rowvec>(xr);
  if (eta_t.n_elem != eta_r.n_elem) { Rcpp::stop("two mean vectors must have the same length."); }
  arma::rowvec var_t = d * arma::exp(c*eta_t);
  arma::rowvec var_r = d * arma::exp(c*eta_r);
  arma::rowvec s_t = arma::log(1.0 + (var_t/(eta_t % eta_t)));
  arma::rowvec s_r = arma::log(1.0 + (var_r/(eta_r % eta_r)));
  arma::rowvec mu_t = arma::log(eta_t) - s_t;
  arma::rowvec mu_r = arma::log(eta_r) - s_r;
  return 0.5*arma::log(s_r/s_t) - (s_r - s_t - (mu_r - mu_t) % (mu_r - mu_t))/(2*s_r);
}

//' @export
//[[Rcpp::export]]
arma::rowvec KL_gamma_b(SEXP xt, SEXP xr)
{
  // Re-define variable type
  arma::rowvec eta_t = Rcpp::as<arma::rowvec>(xt);
  arma::rowvec eta_r = Rcpp::as<arma::rowvec>(xr);
  if (eta_t.n_elem != eta_r.n_elem) { Rcpp::stop("two mean vectors must have the same length."); }
  return arma::log(eta_r/eta_t) + (eta_t - eta_r)/eta_r;
}
