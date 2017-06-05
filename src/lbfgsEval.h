//
arma::irowvec boxCheck(arma::rowvec &R_PARA, const arma::rowvec &R_UPPER, const arma::rowvec &R_LOWER, const arma::irowvec &R_NBD)
{
  int len = R_PARA.n_elem;
  int nbd;
  arma::irowvec BDD(len, fill::zeros); // 0: central; 1: forward; 2: backward
  for (int i = 0; i < len; i++) {
    nbd = R_NBD(i);
    switch(nbd) {
      case 1: { // Lower
        if ((R_PARA(i) - 5e-4) < R_LOWER(i)) { BDD(i) = 1; }
        if (R_PARA(i) < R_LOWER(i)) { R_PARA(i) = R_LOWER(i); }
        break; 
      } 
      case 2: { // Both
        if ((R_PARA(i) - 5e-4) < R_LOWER(i)) { BDD(i) = 1; }
        if ((R_PARA(i) + 5e-4) > R_UPPER(i)) { BDD(i) = 2; }
        if (R_PARA(i) < R_LOWER(i)) { R_PARA(i) = R_LOWER(i); }
        if (R_PARA(i) > R_UPPER(i)) { R_PARA(i) = R_UPPER(i); }
        break; 
      } 
      case 3: { // Upper
        if ((R_PARA(i) + 5e-4) > R_UPPER(i)) { BDD(i) = 2; }
        if (R_PARA(i) > R_UPPER(i)) R_PARA(i) = R_UPPER(i); 
        break; 
      } 
    }
  }
  return BDD;
}
//
double domainMapping(const int INV, const double par, const int nbd, const double upper, const double lower)
{
  double tmp;
  double out = par;
  if (INV == 0) {
    switch(nbd) {
      case 1: // Lower
      {
        tmp = par - lower; 
        if (tmp < 1e-20) { out = -1e12; } else { out = std::log(tmp); } 
        break;
      }
      case 2: // Both
      {
        tmp = (par - lower)/(upper - lower);
        if (tmp < 1e-20) { out = -1e12; } else if (tmp > (1 - 1e-20)) { out = 1e12; } else { out = std::log((tmp/(1-tmp))); } break;
      }
      case 3: // Upper
      {
        tmp = upper - par; 
        if (tmp < 1e-20) { out = -1e12; } else { out = std::log(tmp); } break;
      }
    }  
  } else {
    switch(nbd) {
      case 1:
      {
        tmp = std::exp(par) + lower; out = tmp; break;
      }
      case 2:
      {
        tmp = std::exp(par); tmp = tmp/(1+tmp);
        tmp = tmp*(upper - lower) + lower; out = tmp; break;
      }
      case 3:
      {
        tmp = upper - std::exp(par); out = tmp; break;
      }
    }  
  }
  return out;
}

//
double f_fn(const arma::rowvec &R_PARA, const arma::rowvec &T_PARA, const arma::mat &DESIGN, const arma::rowvec &WT, model_diff_func *func_input,
  const arma::rowvec &R_UPPER, const arma::rowvec &R_LOWER, const arma::irowvec &R_NBD)
{
  
  arma::rowvec R_PARA_OS(R_PARA.n_elem, fill::zeros);
  for (uword i = 0; i < R_PARA.n_elem; i++) { R_PARA_OS(i) = domainMapping(1, R_PARA(i), R_NBD(i), R_UPPER(i), R_LOWER(i)); }

  Rcpp::EvalBase *m1_func = (Rcpp::EvalBase *) func_input->M1_FUNC;
  Rcpp::EvalBase *m2_func = (Rcpp::EvalBase *) func_input->M2_FUNC;
  Rcpp::EvalBase *distFunc = (Rcpp::EvalBase *) func_input->DISTFUNC;

  double fvalTmp = 1e20;
 
  Rcpp::NumericVector eta_T_Rform, eta_R_Rform, DIV_Rform;
  
  eta_T_Rform = (Rcpp::NumericVector) m1_func->eval(Rcpp::wrap(DESIGN), Rcpp::wrap(T_PARA)); 
  eta_R_Rform = (Rcpp::NumericVector) m2_func->eval(Rcpp::wrap(DESIGN), Rcpp::wrap(R_PARA_OS)); 
  DIV_Rform = (Rcpp::NumericVector) distFunc->eval(Rcpp::wrap(eta_T_Rform), Rcpp::wrap(eta_R_Rform));  

  arma::rowvec DIV(DIV_Rform.begin(), DIV_Rform.size(), false);

  fvalTmp = arma::accu(WT % DIV);

  if (std::isnan(fvalTmp)) { fvalTmp = 1e20; }
  if (!(arma::is_finite(fvalTmp))) { fvalTmp = 1e20; }
  //Rprintf("%4.4f\n", fvalTmp);
  return fvalTmp;
}

arma::rowvec f_gr(const double &FVAL, const arma::rowvec &R_PARA, const arma::rowvec &T_PARA, const arma::mat &DESIGN, const arma::rowvec &WT, 
                  model_diff_func *func_input, const arma::irowvec &BDD, const arma::rowvec &R_UPPER, const arma::rowvec &R_LOWER, const arma::irowvec &R_NBD)
{
  //forward difference method
  arma::rowvec gr(R_PARA.n_elem, fill::zeros);
  arma::rowvec fr_para(R_PARA.n_elem, fill::zeros); 
  arma::rowvec bk_para(R_PARA.n_elem, fill::zeros);  
  double fr_val = 0.0; 
  double bk_val = 0.0;
  for (uword i = 0; i < R_PARA.n_elem; i++) {
    fr_para = R_PARA; bk_para = R_PARA;
    if (BDD(i) == 1) { // forward
      fr_para(i) += 1e-3; 
      fr_val = f_fn(fr_para, T_PARA, DESIGN, WT, func_input, R_UPPER, R_LOWER, R_NBD); 
      bk_val = FVAL;
    } else if (BDD(i) == 2) {
      bk_para(i) -= 1e-3; 
      bk_val = f_fn(bk_para, T_PARA, DESIGN, WT, func_input, R_UPPER, R_LOWER, R_NBD); 
      fr_val = FVAL;
    } else {
      fr_para(i) += 5e-4; bk_para(i) -= 5e-4;  
      fr_val = f_fn(fr_para, T_PARA, DESIGN, WT, func_input, R_UPPER, R_LOWER, R_NBD); 
      bk_val = f_fn(bk_para, T_PARA, DESIGN, WT, func_input, R_UPPER, R_LOWER, R_NBD); 
    }
    gr(i) = fr_val*1e3 - bk_val*1e3;
  }
  //rvecPrintf(gr);
  return gr;
}

//
lbfgsfloatval_t evaluate(void *ex, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step)
{
  
  lbfgs_eval EXT = *(lbfgs_eval*)(ex);

  model_diff_func* func_input = EXT.func_input;
  arma::mat DESIGN = EXT.DESIGN;
  arma::rowvec R_UPPER = EXT.R_UPPER; 
  arma::rowvec R_LOWER = EXT.R_LOWER;
  arma::rowvec T_PARA = EXT.T_PARA;
  arma::rowvec WT = EXT.WT;
  arma::irowvec R_NBD = EXT.R_NBD;

  arma::rowvec R_PARA(n);
  for (int i = 0; i < n; i++) { R_PARA(i) = x[i]; }

  //arma::irowvec BDD = boxCheck(R_PARA, R_UPPER, R_LOWER, R_NBD);  
  arma::irowvec BDD(n, fill::zeros);

  double fx = f_fn(R_PARA, T_PARA, DESIGN, WT, func_input, R_UPPER, R_LOWER, R_NBD); 

  arma::rowvec gr = f_gr(fx, R_PARA, T_PARA, DESIGN, WT, func_input, BDD, R_UPPER, R_LOWER, R_NBD); 
  for (int i = 0; i < n; i++) { g[i] = gr(i); }

  return (lbfgsfloatval_t)fx;
}