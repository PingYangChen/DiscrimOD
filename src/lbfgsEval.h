
double domainMapping(const int INV, const double par, const int nbd, const double upper, const double lower)
{
  double tmp;
  double out = par;
  if (INV == 0) {
    switch(nbd) {
      case 1: // Lower
      {
        tmp = par - lower;
        if (tmp < 1e-12) { out = -1e8; } else { out = std::log(tmp); }
        break;
      }
      case 2: // Both
      { 
        if (std::abs(upper - lower) < 1e-12) {
          out = upper;
        } else {
          tmp = (par - lower)/(upper - lower);
          if (tmp < 1e-12) { out = -1e8; } else if (tmp > (1 - 1e-12)) { out = 1e8; } else { out = std::log((tmp/(1.0 - tmp))); }  
        }
        break; 
      }
      case 3: // Upper
      {
        tmp = upper - par;
        if (tmp < 1e-12) { out = -1e8; } else { out = std::log(tmp); } break;
      }
    }
  } else {
    switch(nbd) {
      case 1:
      {
        tmp = std::exp(par);
        if (!std::isfinite(tmp)) { out = 1e8; } else { out = tmp + lower; } break;
      }
      case 2:
      {
        if (std::abs(upper - lower) < 1e-12) {
          out = upper;
        } else {
          tmp = std::exp(par);
          if (!std::isfinite(tmp)) { out = upper; } else { out = (tmp/(1.0 + tmp))*(upper - lower) + lower; } 
        }
        break;
      }
      case 3:
      {
        tmp = std::exp(par); if (!std::isfinite(tmp)) { out = -1e8; } else { out = upper - tmp; } break;
      }
    }
  }
  return out;
}

//
double f_fn(const arma::rowvec R_PARA, const arma::rowvec T_PARA, const arma::mat DESIGN, const arma::rowvec WT, model_diff_func *func_input,
  const arma::rowvec R_UPPER, const arma::rowvec R_LOWER, const arma::irowvec R_NBD,
  const int T_VAR_PARA_L0, const int R_VAR_PARA_L0)
{

  arma::rowvec R_PARA_OS(R_PARA.n_elem, arma::fill::zeros);
  for (arma::uword i = 0; i < R_PARA.n_elem; i++) { R_PARA_OS(i) = domainMapping(1, R_PARA(i), R_NBD(i), R_UPPER(i), R_LOWER(i)); }

  arma::rowvec T_PARA_M = T_PARA.subvec(0, T_VAR_PARA_L0 - 1);
  arma::rowvec T_PARA_V = T_PARA.subvec(T_VAR_PARA_L0, T_PARA.n_elem - 1);
  arma::rowvec R_PARA_M = R_PARA_OS.subvec(0, R_VAR_PARA_L0 - 1);
  arma::rowvec R_PARA_V = R_PARA_OS.subvec(R_VAR_PARA_L0, R_PARA.n_elem - 1);

  double fvalTmp = 1e10;

  if (R_PARA_OS.is_finite()) {
    Rcpp::EvalBase2 *m1_func = (Rcpp::EvalBase2 *) func_input->M1_FUNC;
    Rcpp::EvalBase2 *m2_func = (Rcpp::EvalBase2 *) func_input->M2_FUNC;
    Rcpp::EvalBase2 *v1_func = (Rcpp::EvalBase2 *) func_input->V1_FUNC;
    Rcpp::EvalBase2 *v2_func = (Rcpp::EvalBase2 *) func_input->V2_FUNC;
    Rcpp::EvalBase4 *distFunc = (Rcpp::EvalBase4 *) func_input->DISTFUNC;

    Rcpp::Shield<SEXP> DESIGN_SEXP(Rcpp::wrap(DESIGN));
    Rcpp::Shield<SEXP> T_PARA_M_SEXP(Rcpp::wrap(T_PARA_M));
    Rcpp::Shield<SEXP> T_PARA_V_SEXP(Rcpp::wrap(T_PARA_V));
    Rcpp::Shield<SEXP> R_PARA_M_SEXP(Rcpp::wrap(R_PARA_M));
    Rcpp::Shield<SEXP> R_PARA_V_SEXP(Rcpp::wrap(R_PARA_V));

    Rcpp::NumericMatrix DESIGN_Rform = Rcpp::as<Rcpp::NumericMatrix>(DESIGN_SEXP);
    Rcpp::NumericVector T_PARA_M_Rform = Rcpp::as<Rcpp::NumericVector>(T_PARA_M_SEXP);
    Rcpp::NumericVector T_PARA_V_Rform = Rcpp::as<Rcpp::NumericVector>(T_PARA_V_SEXP);
    Rcpp::NumericVector R_PARA_M_Rform = Rcpp::as<Rcpp::NumericVector>(R_PARA_M_SEXP);
    Rcpp::NumericVector R_PARA_V_Rform = Rcpp::as<Rcpp::NumericVector>(R_PARA_V_SEXP);

    Rcpp::NumericVector eta_T_Rform((int)WT.n_elem), eta_R_Rform((int)WT.n_elem);
    Rcpp::NumericVector var_T_Rform((int)WT.n_elem), var_R_Rform((int)WT.n_elem);
    Rcpp::NumericVector DIV_Rform((int)WT.n_elem);

    eta_T_Rform = (Rcpp::NumericVector) m1_func->eval(DESIGN_Rform, T_PARA_M_Rform);
    eta_R_Rform = (Rcpp::NumericVector) m2_func->eval(DESIGN_Rform, R_PARA_M_Rform);
    var_T_Rform = (Rcpp::NumericVector) v1_func->eval(DESIGN_Rform, T_PARA_V_Rform);
    var_R_Rform = (Rcpp::NumericVector) v2_func->eval(DESIGN_Rform, R_PARA_V_Rform);
    DIV_Rform = (Rcpp::NumericVector) distFunc->eval(eta_T_Rform, eta_R_Rform, var_T_Rform, var_R_Rform);

    bool eta_T_finite = Rcpp::all(Rcpp::is_finite(eta_T_Rform));
    bool eta_R_finite = Rcpp::all(Rcpp::is_finite(eta_R_Rform));
    bool var_T_finite = Rcpp::all(Rcpp::is_finite(var_T_Rform));
    bool var_R_finite = Rcpp::all(Rcpp::is_finite(var_R_Rform));
    bool DIV_finite = Rcpp::all(Rcpp::is_finite(DIV_Rform));

    bool checkFinite = eta_T_finite && eta_R_finite && var_T_finite && var_R_finite && DIV_finite;

    if (checkFinite) {
      arma::rowvec DIV(DIV_Rform.begin(), DIV_Rform.size());
      fvalTmp = 0;
      for (arma::uword i = 0; i < WT.n_elem; i++) { fvalTmp += WT(i)*DIV(i); }
      //fvalTmp = arma::accu(WT % DIV);
    }
  }
  return fvalTmp;
}

arma::rowvec f_gr(const double FVAL, const arma::rowvec R_PARA, const arma::rowvec T_PARA, const arma::mat DESIGN, const arma::rowvec WT,
                  model_diff_func *func_input, const arma::irowvec BDD, const arma::rowvec R_UPPER, const arma::rowvec R_LOWER,
                  const arma::irowvec R_NBD,
                  const int T_VAR_PARA_L0, const int R_VAR_PARA_L0, const double FD_DELTA)
{
  // finite difference method
  double inv_delta = 1.0/FD_DELTA;
  arma::rowvec gr(R_PARA.n_elem, arma::fill::zeros);
  arma::rowvec fr_para(R_PARA.n_elem, arma::fill::zeros);
  arma::rowvec bk_para(R_PARA.n_elem, arma::fill::zeros);
  double fr_val = FVAL;
  double bk_val = FVAL;
  for (arma::uword i = 0; i < R_PARA.n_elem; i++) {
    fr_para = R_PARA; bk_para = R_PARA;
    if (BDD(i) == 1) { // forward
      fr_para(i) += FD_DELTA;
      fr_val = f_fn(fr_para, T_PARA, DESIGN, WT, func_input, R_UPPER, R_LOWER, R_NBD, T_VAR_PARA_L0, R_VAR_PARA_L0);
    } else if (BDD(i) == 2) {
      bk_para(i) -= FD_DELTA;
      bk_val = f_fn(bk_para, T_PARA, DESIGN, WT, func_input, R_UPPER, R_LOWER, R_NBD, T_VAR_PARA_L0, R_VAR_PARA_L0);
    } else {
      fr_para(i) += 0.5*FD_DELTA; bk_para(i) -= 0.5*FD_DELTA;
      fr_val = f_fn(fr_para, T_PARA, DESIGN, WT, func_input, R_UPPER, R_LOWER, R_NBD, T_VAR_PARA_L0, R_VAR_PARA_L0);
      bk_val = f_fn(bk_para, T_PARA, DESIGN, WT, func_input, R_UPPER, R_LOWER, R_NBD, T_VAR_PARA_L0, R_VAR_PARA_L0);
    }
    gr(i) = fr_val*inv_delta - bk_val*inv_delta;
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
  arma::rowvec R_UPPER = EXT.UPPER;
  arma::rowvec R_LOWER = EXT.LOWER;
  arma::rowvec T_PARA = EXT.T_PARA;
  int T_VAR_PARA_L0 = EXT.T_VAR_PARA_L0;
  int R_VAR_PARA_L0 = EXT.R_VAR_PARA_L0;
  arma::rowvec WT = EXT.WT;
  arma::irowvec R_NBD = EXT.NBD;
  double FD_DELTA = EXT.FD_DELTA;

  arma::rowvec R_PARA(n);
  for (int i = 0; i < n; i++) { R_PARA(i) = (double)(x[i]); }

  double fx = 1e10;
  arma::rowvec gr(n, fill::zeros);

  if (R_PARA.is_finite()) {
    arma::irowvec BDD(n, fill::zeros);
    fx = f_fn(R_PARA, T_PARA, DESIGN, WT, func_input, R_UPPER, R_LOWER, R_NBD, T_VAR_PARA_L0, R_VAR_PARA_L0);
    gr = f_gr(fx, R_PARA, T_PARA, DESIGN, WT, func_input, BDD, R_UPPER, R_LOWER, R_NBD, T_VAR_PARA_L0, R_VAR_PARA_L0, FD_DELTA);
  }
  for (int i = 0; i < n; i++) { g[i] = (lbfgsfloatval_t)(gr(i)); }
  return (lbfgsfloatval_t)fx;
}


// FEDOROV-WYNN NEXT POINT
double dd_fn(const arma::mat X_VEC, const arma::rowvec R_PARA, const arma::rowvec T_PARA,
             const int T_VAR_PARA_L0, const int R_VAR_PARA_L0, model_diff_func *func_input,
             const arma::rowvec X_UPPER, const arma::rowvec X_LOWER, const arma::irowvec X_NBD)
{

  arma::mat X_OS(1, X_VEC.n_cols, fill::zeros);
  for (uword i = 0; i < X_VEC.n_cols; i++) { X_OS(0, i) = domainMapping(1, X_VEC(0, i), X_NBD(i), X_UPPER(i), X_LOWER(i)); }

  arma::rowvec T_PARA_M = T_PARA.subvec(0, T_VAR_PARA_L0 - 1);
  arma::rowvec T_PARA_V = T_PARA.subvec(T_VAR_PARA_L0, T_PARA.n_elem - 1);
  arma::rowvec R_PARA_M = R_PARA.subvec(0, R_VAR_PARA_L0 - 1);
  arma::rowvec R_PARA_V = R_PARA.subvec(R_VAR_PARA_L0, R_PARA.n_elem - 1);

  double fvalTmp = 1e10;

  if (X_OS.is_finite()) {
    Rcpp::EvalBase2 *m1_func = (Rcpp::EvalBase2 *) func_input->M1_FUNC;
    Rcpp::EvalBase2 *m2_func = (Rcpp::EvalBase2 *) func_input->M2_FUNC;
    Rcpp::EvalBase2 *v1_func = (Rcpp::EvalBase2 *) func_input->V1_FUNC;
    Rcpp::EvalBase2 *v2_func = (Rcpp::EvalBase2 *) func_input->V2_FUNC;
    Rcpp::EvalBase4 *distFunc = (Rcpp::EvalBase4 *) func_input->DISTFUNC;

    Rcpp::Shield<SEXP> DESIGN_SEXP(Rcpp::wrap(X_OS));
    Rcpp::Shield<SEXP> T_PARA_M_SEXP(Rcpp::wrap(T_PARA_M));
    Rcpp::Shield<SEXP> T_PARA_V_SEXP(Rcpp::wrap(T_PARA_V));
    Rcpp::Shield<SEXP> R_PARA_M_SEXP(Rcpp::wrap(R_PARA_M));
    Rcpp::Shield<SEXP> R_PARA_V_SEXP(Rcpp::wrap(R_PARA_V));

    Rcpp::NumericMatrix DESIGN_Rform = Rcpp::as<Rcpp::NumericMatrix>(DESIGN_SEXP);
    Rcpp::NumericVector T_PARA_M_Rform = Rcpp::as<Rcpp::NumericVector>(T_PARA_M_SEXP);
    Rcpp::NumericVector T_PARA_V_Rform = Rcpp::as<Rcpp::NumericVector>(T_PARA_V_SEXP);
    Rcpp::NumericVector R_PARA_M_Rform = Rcpp::as<Rcpp::NumericVector>(R_PARA_M_SEXP);
    Rcpp::NumericVector R_PARA_V_Rform = Rcpp::as<Rcpp::NumericVector>(R_PARA_V_SEXP);

    Rcpp::NumericVector eta_T_Rform((int)X_VEC.n_elem), eta_R_Rform((int)X_VEC.n_elem);
    Rcpp::NumericVector var_T_Rform((int)X_VEC.n_elem), var_R_Rform((int)X_VEC.n_elem);
    Rcpp::NumericVector DIV_Rform((int)X_VEC.n_elem);

    eta_T_Rform = (Rcpp::NumericVector) m1_func->eval(DESIGN_Rform, T_PARA_M_Rform);
    eta_R_Rform = (Rcpp::NumericVector) m2_func->eval(DESIGN_Rform, R_PARA_M_Rform);
    var_T_Rform = (Rcpp::NumericVector) v1_func->eval(DESIGN_Rform, T_PARA_V_Rform);
    var_R_Rform = (Rcpp::NumericVector) v2_func->eval(DESIGN_Rform, R_PARA_V_Rform);
    DIV_Rform = (Rcpp::NumericVector) distFunc->eval(eta_T_Rform, eta_R_Rform, var_T_Rform, var_R_Rform);

    bool eta_T_finite = Rcpp::all(Rcpp::is_finite(eta_T_Rform));
    bool eta_R_finite = Rcpp::all(Rcpp::is_finite(eta_R_Rform));
    bool var_T_finite = Rcpp::all(Rcpp::is_finite(var_T_Rform));
    bool var_R_finite = Rcpp::all(Rcpp::is_finite(var_R_Rform));
    bool DIV_finite = Rcpp::all(Rcpp::is_finite(DIV_Rform));

    bool checkFinite = eta_T_finite && eta_R_finite && var_T_finite && var_R_finite && DIV_finite;

    if (checkFinite) {
      if (Rcpp::all(Rcpp::is_finite(DIV_Rform))) {
        fvalTmp = (-1.0)*(double)DIV_Rform[0];
      }
    }
  }
  return fvalTmp;
}

arma::rowvec dd_gr(const double FVAL, const arma::mat X_VEC, const arma::rowvec R_PARA, const arma::rowvec T_PARA,
                   const int T_VAR_PARA_L0, const int R_VAR_PARA_L0,
                   model_diff_func *func_input, const arma::irowvec BDD, const arma::rowvec X_UPPER, const arma::rowvec X_LOWER,
                   const arma::irowvec X_NBD, const double FD_DELTA)
{
  // finite difference method
  double inv_delta = 1.0/FD_DELTA;
  arma::rowvec gr(X_VEC.n_cols, fill::zeros);
  arma::mat fr_x(1, X_VEC.n_cols, fill::zeros);
  arma::mat bk_x(1, X_VEC.n_cols, fill::zeros);
  double fr_val = FVAL;
  double bk_val = FVAL;
  for (uword i = 0; i < X_VEC.n_cols; i++) {
    fr_x = X_VEC; bk_x = X_VEC;
    if (BDD(i) == 1) { // forward
      fr_x(0,i) += FD_DELTA;
      fr_val = dd_fn(fr_x, R_PARA, T_PARA, T_VAR_PARA_L0, R_VAR_PARA_L0, func_input, X_UPPER, X_LOWER, X_NBD);
    } else if (BDD(i) == 2) {
      bk_x(0,i) -= FD_DELTA;
      bk_val = dd_fn(bk_x, R_PARA, T_PARA, T_VAR_PARA_L0, R_VAR_PARA_L0, func_input, X_UPPER, X_LOWER, X_NBD);
    } else {
      fr_x(0,i) += 0.5*FD_DELTA; bk_x(0,i) -= 0.5*FD_DELTA;
      fr_val = dd_fn(fr_x, R_PARA, T_PARA, T_VAR_PARA_L0, R_VAR_PARA_L0, func_input, X_UPPER, X_LOWER, X_NBD);
      bk_val = dd_fn(bk_x, R_PARA, T_PARA, T_VAR_PARA_L0, R_VAR_PARA_L0, func_input, X_UPPER, X_LOWER, X_NBD);
    }
    gr(i) = fr_val*inv_delta - bk_val*inv_delta;
  }
  //rvecPrintf(gr);
  return gr;
}
//
lbfgsfloatval_t evaluate_dirdev(void *ex, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step)
{

  lbfgs_eval EXT = *(lbfgs_eval*)(ex);

  model_diff_func* func_input = EXT.func_input;

  arma::rowvec X_UPPER = EXT.UPPER;
  arma::rowvec X_LOWER = EXT.LOWER;
  arma::rowvec T_PARA = EXT.T_PARA;
  arma::rowvec R_PARA = EXT.R_PARA;
  int T_VAR_PARA_L0 = EXT.T_VAR_PARA_L0;
  int R_VAR_PARA_L0 = EXT.R_VAR_PARA_L0;
  arma::irowvec X_NBD = EXT.NBD;
  double FD_DELTA = EXT.FD_DELTA;

  arma::mat X_VEC(1, n, fill::zeros);
  for (int i = 0; i < n; i++) { X_VEC(0,i) = (double)(x[i]); }

  double fx = 1e10;
  arma::rowvec gr(n, fill::zeros);

  if (X_VEC.is_finite()) {
    arma::irowvec BDD(n, fill::zeros);
    fx = dd_fn(X_VEC, R_PARA, T_PARA, T_VAR_PARA_L0, R_VAR_PARA_L0, func_input, X_UPPER, X_LOWER, X_NBD);
    gr = dd_gr(fx, X_VEC, R_PARA, T_PARA, T_VAR_PARA_L0, R_VAR_PARA_L0, func_input, BDD, X_UPPER, X_LOWER, X_NBD, FD_DELTA);
  }
  for (int i = 0; i < n; i++) { g[i] = (lbfgsfloatval_t)(gr(i)); }
  return (lbfgsfloatval_t)fx;
}

// Remes'a algorithm
double cpoly(const arma::rowvec R_PARA, const arma::rowvec T_PARA, const int T_VAR_PARA_L0, const int R_VAR_PARA_L0,
             const arma::mat DESIGN_PT, model_diff_func *func_input)
{
  arma::rowvec T_PARA_M = T_PARA.subvec(0, T_VAR_PARA_L0 - 1);
  arma::rowvec T_PARA_V = T_PARA.subvec(T_VAR_PARA_L0, T_PARA.n_elem - 1);
  arma::rowvec R_PARA_M = R_PARA.subvec(0, R_VAR_PARA_L0 - 1);
  arma::rowvec R_PARA_V = R_PARA.subvec(R_VAR_PARA_L0, R_PARA.n_elem - 1);

  double cployVal = 0.0;

  Rcpp::EvalBase2 *m1_func = (Rcpp::EvalBase2 *) func_input->M1_FUNC;
  Rcpp::EvalBase2 *m2_func = (Rcpp::EvalBase2 *) func_input->M2_FUNC;
  Rcpp::EvalBase2 *v1_func = (Rcpp::EvalBase2 *) func_input->V1_FUNC;
  Rcpp::EvalBase2 *v2_func = (Rcpp::EvalBase2 *) func_input->V2_FUNC;

  Rcpp::Shield<SEXP> DESIGN_SEXP(Rcpp::wrap(DESIGN_PT));
  Rcpp::Shield<SEXP> T_PARA_M_SEXP(Rcpp::wrap(T_PARA_M));
  Rcpp::Shield<SEXP> T_PARA_V_SEXP(Rcpp::wrap(T_PARA_V));
  Rcpp::Shield<SEXP> R_PARA_M_SEXP(Rcpp::wrap(R_PARA_M));
  Rcpp::Shield<SEXP> R_PARA_V_SEXP(Rcpp::wrap(R_PARA_V));

  Rcpp::NumericMatrix DESIGN_Rform = Rcpp::as<Rcpp::NumericMatrix>(DESIGN_SEXP);
  Rcpp::NumericVector T_PARA_M_Rform = Rcpp::as<Rcpp::NumericVector>(T_PARA_M_SEXP);
  Rcpp::NumericVector T_PARA_V_Rform = Rcpp::as<Rcpp::NumericVector>(T_PARA_V_SEXP);
  Rcpp::NumericVector R_PARA_M_Rform = Rcpp::as<Rcpp::NumericVector>(R_PARA_M_SEXP);
  Rcpp::NumericVector R_PARA_V_Rform = Rcpp::as<Rcpp::NumericVector>(R_PARA_V_SEXP);

  Rcpp::NumericVector eta_T_Rform((int)DESIGN_PT.n_elem), eta_R_Rform((int)DESIGN_PT.n_elem);
  Rcpp::NumericVector var_T_Rform((int)DESIGN_PT.n_elem), var_R_Rform((int)DESIGN_PT.n_elem);
  Rcpp::NumericVector DIV_Rform((int)DESIGN_PT.n_elem);

  eta_T_Rform = (Rcpp::NumericVector) m1_func->eval(DESIGN_Rform, T_PARA_M_Rform);
  eta_R_Rform = (Rcpp::NumericVector) m2_func->eval(DESIGN_Rform, R_PARA_M_Rform);
  var_T_Rform = (Rcpp::NumericVector) v1_func->eval(DESIGN_Rform, T_PARA_V_Rform);
  var_R_Rform = (Rcpp::NumericVector) v2_func->eval(DESIGN_Rform, R_PARA_V_Rform);

  bool eta_T_finite = Rcpp::all(Rcpp::is_finite(eta_T_Rform));
  bool eta_R_finite = Rcpp::all(Rcpp::is_finite(eta_R_Rform));
  bool var_T_finite = Rcpp::all(Rcpp::is_finite(var_T_Rform));
  bool var_R_finite = Rcpp::all(Rcpp::is_finite(var_R_Rform));

  bool checkFinite = eta_T_finite && eta_R_finite && var_T_finite && var_R_finite;

  if (checkFinite) {
    DIV_Rform = eta_T_Rform - eta_R_Rform;
    cployVal = (double)(DIV_Rform[0]);
  }
  return cployVal;
}

arma::rowvec rival_gr(const arma::rowvec R_PARA, const int R_VAR_PARA_L0, const arma::mat DESIGN_PT,
                      model_diff_func *func_input, const arma::irowvec BDD, const double FD_DELTA)
{
  arma::rowvec R_PARA_M = R_PARA.subvec(0, R_VAR_PARA_L0 - 1);
  arma::rowvec R_PARA_V = R_PARA.subvec(R_VAR_PARA_L0, R_PARA.n_elem - 1);

  Rcpp::EvalBase2 *m2_func = (Rcpp::EvalBase2 *) func_input->M2_FUNC;
  // finite difference method
  double inv_delta = 1.0/FD_DELTA;
  arma::rowvec gr(R_PARA_M.n_elem, fill::zeros);
  arma::rowvec fr_para(R_PARA_M.n_elem, fill::zeros);
  arma::rowvec bk_para(R_PARA_M.n_elem, fill::zeros);
  double fr_val, bk_val;

  Shield<SEXP> DESIGN_SEXP(Rcpp::wrap(DESIGN_PT));
  Rcpp::NumericMatrix DESIGN_Rform = Rcpp::as<Rcpp::NumericMatrix>(DESIGN_SEXP);

  Rcpp::NumericVector fr_eta_R_Rform((int)DESIGN_PT.n_rows), bk_eta_R_Rform((int)DESIGN_PT.n_rows);

  for (uword i = 0; i < R_PARA_M.n_elem; i++) {
    fr_para = R_PARA_M; bk_para = R_PARA_M;
    if (BDD(i) == 1) { // forward
      fr_para(i) += FD_DELTA;
    } else if (BDD(i) == 2) {
      bk_para(i) -= FD_DELTA;
    } else {
      fr_para(i) += 0.5*FD_DELTA; bk_para(i) -= 0.5*FD_DELTA;
    }

    Shield<SEXP> fr_R_PARA_SEXP(Rcpp::wrap(fr_para));
    Shield<SEXP> bk_R_PARA_SEXP(Rcpp::wrap(bk_para));

    Rcpp::NumericVector fr_R_PARA_Rform = Rcpp::as<Rcpp::NumericVector>(fr_R_PARA_SEXP);
    Rcpp::NumericVector bk_R_PARA_Rform = Rcpp::as<Rcpp::NumericVector>(bk_R_PARA_SEXP);

    fr_eta_R_Rform = (Rcpp::NumericVector) m2_func->eval(DESIGN_Rform, fr_R_PARA_Rform);
    bk_eta_R_Rform = (Rcpp::NumericVector) m2_func->eval(DESIGN_Rform, bk_R_PARA_Rform);
    fr_val = (double)(fr_eta_R_Rform[0]); bk_val = (double)(bk_eta_R_Rform[0]);
    gr(i) = fr_val*inv_delta - bk_val*inv_delta;
  }
  return gr;
}


//
double remes_r_fn(const arma::rowvec R_PARA, const arma::rowvec T_PARA, const int T_VAR_PARA_L0, const int R_VAR_PARA_L0,
                  const arma::mat DESIGN, model_diff_func *func_input,
                  const arma::rowvec R_UPPER, const arma::rowvec R_LOWER, const arma::irowvec R_NBD)
{
  arma::rowvec R_PARA_OS(R_PARA.n_elem, fill::zeros);
  for (uword i = 0; i < R_PARA.n_elem; i++) { R_PARA_OS(i) = domainMapping(1, R_PARA(i), R_NBD(i), R_UPPER(i), R_LOWER(i)); }

  double fvalTmp = 1e10;
  if (R_PARA_OS.is_finite()) {
    arma::vec cpolyVals(DESIGN.n_rows);
    for (uword i = 0; i < DESIGN.n_rows; i++) {
      arma::mat DESIGN_PT = DESIGN.row(i);
      cpolyVals(i) = std::abs(cpoly(R_PARA_OS, T_PARA, T_VAR_PARA_L0, R_VAR_PARA_L0, DESIGN_PT, func_input));
    }
    fvalTmp = cpolyVals.max();
  }
  return fvalTmp;
}

arma::rowvec remes_r_gr(const double FVAL, const arma::rowvec R_PARA, const arma::rowvec T_PARA, const int T_VAR_PARA_L0, const int R_VAR_PARA_L0,
                        const arma::mat DESIGN, model_diff_func *func_input,
                        const arma::irowvec BDD, const arma::rowvec R_UPPER, const arma::rowvec R_LOWER,
                        const arma::irowvec R_NBD, const double FD_DELTA)
{
  // finite difference method
  double inv_delta = 1.0/FD_DELTA;
  arma::rowvec gr(R_PARA.n_elem, fill::zeros);
  arma::rowvec fr_para(R_PARA.n_elem, fill::zeros);
  arma::rowvec bk_para(R_PARA.n_elem, fill::zeros);
  double fr_val = FVAL;
  double bk_val = FVAL;
  for (uword i = 0; i < R_PARA.n_elem; i++) {
    fr_para = R_PARA; bk_para = R_PARA;
    if (BDD(i) == 1) { // forward
      fr_para(i) += FD_DELTA;
      fr_val = remes_r_fn(fr_para, T_PARA, T_VAR_PARA_L0, R_VAR_PARA_L0, DESIGN, func_input, R_UPPER, R_LOWER, R_NBD);
    } else if (BDD(i) == 2) {
      bk_para(i) -= FD_DELTA;
      bk_val = remes_r_fn(bk_para, T_PARA, T_VAR_PARA_L0, R_VAR_PARA_L0, DESIGN, func_input, R_UPPER, R_LOWER, R_NBD);
    } else {
      fr_para(i) += 0.5*FD_DELTA; bk_para(i) -= 0.5*FD_DELTA;
      fr_val = remes_r_fn(fr_para, T_PARA, T_VAR_PARA_L0, R_VAR_PARA_L0, DESIGN, func_input, R_UPPER, R_LOWER, R_NBD);
      bk_val = remes_r_fn(bk_para, T_PARA, T_VAR_PARA_L0, R_VAR_PARA_L0, DESIGN, func_input, R_UPPER, R_LOWER, R_NBD);
    }
    gr(i) = fr_val*inv_delta - bk_val*inv_delta;
  }
  return gr;
}

//
double remes_x_fn(const arma::mat X_VEC, const arma::rowvec R_PARA, const arma::rowvec T_PARA, const int T_VAR_PARA_L0, const int R_VAR_PARA_L0,
                  model_diff_func *func_input, const arma::rowvec X_UPPER, const arma::rowvec X_LOWER, const arma::irowvec X_NBD)
{

  arma::mat X_OS(1, X_VEC.n_cols, fill::zeros);
  for (uword i = 0; i < X_VEC.n_cols; i++) { X_OS(0, i) = domainMapping(1, X_VEC(0, i), X_NBD(i), X_UPPER(i), X_LOWER(i)); }

  double fvalTmp = 1e10;
  if (X_OS.is_finite()) {
    fvalTmp = (-1.0)*std::abs(cpoly(R_PARA, T_PARA, T_VAR_PARA_L0, R_VAR_PARA_L0, X_OS, func_input));
  }
  return fvalTmp;
}

arma::rowvec remes_x_gr(const double FVAL, const arma::mat X_VEC, const arma::rowvec R_PARA, const arma::rowvec T_PARA,
                        const int T_VAR_PARA_L0, const int R_VAR_PARA_L0, model_diff_func *func_input,
                        const arma::irowvec BDD, const arma::rowvec X_UPPER, const arma::rowvec X_LOWER,
                        const arma::irowvec X_NBD, const double FD_DELTA)
{
  // finite difference method
  double inv_delta = 1.0/FD_DELTA;
  arma::rowvec gr(X_VEC.n_cols, fill::zeros);
  arma::mat fr_x(1, X_VEC.n_cols, fill::zeros);
  arma::mat bk_x(1, X_VEC.n_cols, fill::zeros);
  double fr_val = FVAL;
  double bk_val = FVAL;
  for (uword i = 0; i < X_VEC.n_cols; i++) {
    fr_x = X_VEC; bk_x = X_VEC;
    if (BDD(i) == 1) { // forward
      fr_x(0,i) += FD_DELTA;
      fr_val = remes_x_fn(fr_x, R_PARA, T_PARA, T_VAR_PARA_L0, R_VAR_PARA_L0, func_input, X_UPPER, X_LOWER, X_NBD);
    } else if (BDD(i) == 2) {
      bk_x(0,i) -= FD_DELTA;
      bk_val = remes_x_fn(bk_x, R_PARA, T_PARA, T_VAR_PARA_L0, R_VAR_PARA_L0, func_input, X_UPPER, X_LOWER, X_NBD);
    } else {
      fr_x(0,i) += 0.5*FD_DELTA; bk_x(0,i) -= 0.5*FD_DELTA;
      fr_val = remes_x_fn(fr_x, R_PARA, T_PARA, T_VAR_PARA_L0, R_VAR_PARA_L0, func_input, X_UPPER, X_LOWER, X_NBD);
      bk_val = remes_x_fn(bk_x, R_PARA, T_PARA, T_VAR_PARA_L0, R_VAR_PARA_L0, func_input, X_UPPER, X_LOWER, X_NBD);
    }
    gr(i) = fr_val*inv_delta - bk_val*inv_delta;
  }
  //rvecPrintf(gr);
  return gr;
}

//
lbfgsfloatval_t evaluate_remes_r(void *ex, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step)
{
  lbfgs_eval EXT = *(lbfgs_eval*)(ex);

  model_diff_func* func_input = EXT.func_input;
  arma::mat DESIGN = EXT.DESIGN;
  arma::rowvec R_UPPER = EXT.UPPER;
  arma::rowvec R_LOWER = EXT.LOWER;
  arma::rowvec T_PARA = EXT.T_PARA;
  int T_VAR_PARA_L0 = EXT.T_VAR_PARA_L0;
  int R_VAR_PARA_L0 = EXT.R_VAR_PARA_L0;
  arma::irowvec R_NBD = EXT.NBD;
  double FD_DELTA = EXT.FD_DELTA;

  arma::rowvec R_PARA(n);
  for (int i = 0; i < n; i++) { R_PARA(i) = (double)(x[i]); }

  double fx = 1e10;
  arma::rowvec gr(n, fill::zeros);

  if (R_PARA.is_finite()) {
    arma::irowvec BDD(n, fill::zeros);
    fx = remes_r_fn(R_PARA, T_PARA, T_VAR_PARA_L0, R_VAR_PARA_L0, DESIGN, func_input, R_UPPER, R_LOWER, R_NBD);
    gr = remes_r_gr(fx, R_PARA, T_PARA, T_VAR_PARA_L0, R_VAR_PARA_L0, DESIGN, func_input, BDD, R_UPPER, R_LOWER, R_NBD, FD_DELTA);
  }
  for (int i = 0; i < n; i++) { g[i] = (lbfgsfloatval_t)(gr(i)); }
  return (lbfgsfloatval_t)fx;
}

lbfgsfloatval_t evaluate_remes_x(void *ex, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step)
{

  lbfgs_eval EXT = *(lbfgs_eval*)(ex);

  model_diff_func* func_input = EXT.func_input;

  arma::rowvec X_UPPER = EXT.UPPER;
  arma::rowvec X_LOWER = EXT.LOWER;
  arma::rowvec T_PARA = EXT.T_PARA;
  arma::rowvec R_PARA = EXT.R_PARA;
  int T_VAR_PARA_L0 = EXT.T_VAR_PARA_L0;
  int R_VAR_PARA_L0 = EXT.R_VAR_PARA_L0;
  arma::irowvec X_NBD = EXT.NBD;
  double FD_DELTA = EXT.FD_DELTA;

  arma::mat X_VEC(1, n, fill::zeros);
  for (int i = 0; i < n; i++) { X_VEC(0,i) = (double)(x[i]); }

  double fx = 1e10;
  arma::rowvec gr(n, fill::zeros);

  if (X_VEC.is_finite()) {
    arma::irowvec BDD(n, fill::zeros);
    fx = remes_x_fn(X_VEC, R_PARA, T_PARA, T_VAR_PARA_L0, R_VAR_PARA_L0, func_input, X_UPPER, X_LOWER, X_NBD);
    gr = remes_x_gr(fx, X_VEC, R_PARA, T_PARA, T_VAR_PARA_L0, R_VAR_PARA_L0, func_input, BDD, X_UPPER, X_LOWER, X_NBD, FD_DELTA);
  }
  for (int i = 0; i < n; i++) { g[i] = (lbfgsfloatval_t)(gr(i)); }
  return (lbfgsfloatval_t)fx;
}
