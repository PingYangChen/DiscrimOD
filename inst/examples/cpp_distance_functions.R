# DISTANCE FUNCTION LIST
library(Rcpp)

l2_diff_cpp <- cppFunction('
	Rcpp::NumericVector gamma_diff(Rcpp::NumericVector xt, Rcpp::NumericVector xr) {
		Rcpp::NumericVector div(xt.size()); double diff;
		for (int i = 0; i < xt.size(); i++) {
			diff = xt(i) - xr(i); div(i) = diff*diff;
		}
  	return div;
}')

lognorm_for_mm_cpp <- cppFunction('
	Rcpp::NumericVector log_norm_B(Rcpp::NumericVector xt, Rcpp::NumericVector xr) {
	  double sigsq = 1.0;
	  Rcpp::NumericVector div(xt.size());
	  double xt2, xr2, vt, vr, mt, mr;
	  for (int i = 0; i < xt.size(); i++) {
	  	xt2 = xt(i)*xt(i); 
	  	xr2 = xr(i)*xr(i);
	  	vt = (std::exp(sigsq) - 1.0)*xt2;
	  	vr = (std::exp(sigsq) - 1.0)*xr2;
		  mt = std::log(xt(i)) - 0.5*std::log(1.0 + (vt/xt2));
		  mr = std::log(xr(i)) - 0.5*std::log(1.0 + (vr/xr2));
		  div(i) = 0.5*((mt - mr)*(mt - mr))/sigsq;
	  }  
	  return div;
}')

lognorm_for_tox_cpp <- cppFunction('
	Rcpp::NumericVector log_norm_B(Rcpp::NumericVector xt, Rcpp::NumericVector xr) {
	  double sigsq = 0.01;
	  Rcpp::NumericVector div(xt.size());
	  double xt2, xr2, vt, vr, mt, mr;
	  for (int i = 0; i < xt.size(); i++) {
	  	xt2 = xt(i)*xt(i); 
	  	xr2 = xr(i)*xr(i);
	  	vt = (std::exp(sigsq) - 1.0)*xt2;
	  	vr = (std::exp(sigsq) - 1.0)*xr2;
		  mt = std::log(xt(i)) - 0.5*std::log(1.0 + (vt/xt2));
		  mr = std::log(xr(i)) - 0.5*std::log(1.0 + (vr/xr2));
		  div(i) = 0.5*((mt - mr)*(mt - mr))/sigsq;
	  }  
	  return div;
}')

gamma_diff_cpp <- cppFunction('
	Rcpp::NumericVector gamma_diff(Rcpp::NumericVector xt, Rcpp::NumericVector xr) {
		Rcpp::NumericVector div(xt.size());
		for (int i = 0; i < xt.size(); i++) {
	  	div(i) = std::log(xr(i)/xt(i)) + (xt(i) - xr(i))/xr(i);
		}
  	return div;
}')

logit_diff_cpp <- cppFunction('
	Rcpp::NumericVector logit_diff(Rcpp::NumericVector xt, Rcpp::NumericVector xr) {
		Rcpp::NumericVector div(xt.size());
		double et, er, mt, mr;
		for (int i = 0; i < xt.size(); i++) {
			et = std::exp(xt(i)); mt = et/(1.0 + et);
			er = std::exp(xr(i)); mr = er/(1.0 + er);
	  	div(i) = mt*std::log(mt/mr) + (1.0 - mt)*std::log((1.0 - mt)/(1.0 - mr));
		}
  	return div;
}')
