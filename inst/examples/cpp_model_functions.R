library(Rcpp)

# Michaelis-Menten
enzyme2 <- cppFunction('
  Rcpp::NumericVector enzyme2(Rcpp::NumericMatrix x, Rcpp::NumericVector p) {
    Rcpp::NumericVector eta(x.nrow());
    for (int i = 0; i < x.nrow(); i++) { 
      eta(i) = p(0)*x(i,0)/(p(1) + x(i,0)) + p(2)*x(i,0); }
    return eta; 
}')

enzyme1 <- cppFunction('
  Rcpp::NumericVector enzyme2(Rcpp::NumericMatrix x, Rcpp::NumericVector p) {
    Rcpp::NumericVector eta(x.nrow());
    for (int i = 0; i < x.nrow(); i++) { 
      eta(i) = p(0)*x(i,0)/(p(1) + x(i,0)); }
    return eta; 
}')

# AF1975
af1975_1 <- cppFunction('
  Rcpp::NumericVector enzyme2(Rcpp::NumericMatrix x, Rcpp::NumericVector p) {
    Rcpp::NumericVector eta(x.nrow());
    for (int i = 0; i < x.nrow(); i++) { 
      eta(i) = p(0) + p(1)*std::exp(x(i,0)) + p(2)*std::exp(-1.0*x(i,0)); }
    return eta; 
}')

af1975_2 <- cppFunction('
  Rcpp::NumericVector enzyme2(Rcpp::NumericMatrix x, Rcpp::NumericVector p) {
    Rcpp::NumericVector eta(x.nrow());
    for (int i = 0; i < x.nrow(); i++) { 
      eta(i) = p(0) + p(1)*x(i,0) + p(2)*x(i,0)*x(i,0); }
    return eta; 
}')

af1975_3 <- cppFunction('
  Rcpp::NumericVector enzyme2(Rcpp::NumericMatrix x, Rcpp::NumericVector p) {
  	double pi_val = 3.141592653589793;
    Rcpp::NumericVector eta(x.nrow());
    for (int i = 0; i < x.nrow(); i++) { 
      eta(i) = p(0) + p(1)*std::sin(0.5*pi_val*x(i,0)) + p(2)*std::cos(0.5*pi_val*x(i,0)) + p(3)*std::sin(pi_val*x(i,0)); }
    return eta; 
}')

# tox
tox5 <- cppFunction('
  Rcpp::NumericVector enzyme2(Rcpp::NumericMatrix x, Rcpp::NumericVector p) {
    Rcpp::NumericVector eta(x.nrow());
    for (int i = 0; i < x.nrow(); i++) { 
      eta(i) = p(0)*(p(2) - (p(2) - 1.0)*std::exp((-1.0)*std::pow(x(i,0)/p(1), p(3)))); }
    return eta; 
}')

tox4 <- cppFunction('
  Rcpp::NumericVector enzyme2(Rcpp::NumericMatrix x, Rcpp::NumericVector p) {
    Rcpp::NumericVector eta(x.nrow());
    for (int i = 0; i < x.nrow(); i++) { 
      eta(i) = p(0)*(p(2) - (p(2) - 1.0)*std::exp((-1.0)*x(i,0)/p(1))); }
    return eta; 
}')

tox3 <- cppFunction('
  Rcpp::NumericVector enzyme2(Rcpp::NumericMatrix x, Rcpp::NumericVector p) {
    Rcpp::NumericVector eta(x.nrow());
    for (int i = 0; i < x.nrow(); i++) { 
      eta(i) = p(0)*std::exp((-1.0)*std::pow(x(i,0)/p(1), p(2))); }
    return eta; 
}')

tox2 <- cppFunction('
  Rcpp::NumericVector enzyme2(Rcpp::NumericMatrix x, Rcpp::NumericVector p) {
    Rcpp::NumericVector eta(x.nrow());
    for (int i = 0; i < x.nrow(); i++) { 
      eta(i) = p(0)*std::exp((-1.0)*x(i,0)/p(1)); }
    return eta; 
}')

tox1 <- cppFunction('
  Rcpp::NumericVector enzyme2(Rcpp::NumericMatrix x, Rcpp::NumericVector p) {
    Rcpp::NumericVector eta(x.nrow());
    for (int i = 0; i < x.nrow(); i++) { eta(i) = p(0); }
    return eta; 
}')

# logit
linlogi4 <- cppFunction('
  Rcpp::NumericVector enzyme2(Rcpp::NumericMatrix x, Rcpp::NumericVector p) {
    Rcpp::NumericVector eta(x.nrow());
    for (int i = 0; i < x.nrow(); i++) { 
      eta(i) = p(0) + p(1)*x(i,0) + p(2)*x(i,0)*x(i,0); }
    return eta; 
}')

linlogi3 <- cppFunction('
  Rcpp::NumericVector enzyme2(Rcpp::NumericMatrix x, Rcpp::NumericVector p) {
    Rcpp::NumericVector eta(x.nrow());
    for (int i = 0; i < x.nrow(); i++) { eta(i) = x(i,0)*(p(0) + p(1)*x(i,0)); }
    return eta; 
}')

linlogi2 <- cppFunction('
  Rcpp::NumericVector enzyme2(Rcpp::NumericMatrix x, Rcpp::NumericVector p) {
    Rcpp::NumericVector eta(x.nrow());
    for (int i = 0; i < x.nrow(); i++) { eta(i) = p(0) + p(1)*x(i,0); }
    return eta; 
}')

linlogi1 <- cppFunction('
  Rcpp::NumericVector enzyme2(Rcpp::NumericMatrix x, Rcpp::NumericVector p) {
    Rcpp::NumericVector eta(x.nrow());
    for (int i = 0; i < x.nrow(); i++) { eta(i) = p(0)*x(i,0); }
    return eta; 
}')