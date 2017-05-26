
// INCLUDE HEADER FILES
#include "lbfgsKernel.h"

// DECLARE FUNCTIONS
double criterionList(const int &LOOPID, const PSO_OPTIONS PSO_OPTS[],
                     const OBJ_INFO &OBJ, const Rcpp::List MODEL_INFO_LIST, const SEXP DIST_FUNC_SEXP, const SEXP env, 
										 const arma::mat &DESIGN, const arma::rowvec &WT, arma::mat &R_PARA);
double minDistCalc(const int &LOOPID, const PSO_OPTIONS PSO_OPTS[],
                   const int &tid, const int &rid, const OBJ_INFO &OBJ, const Rcpp::List MODEL_INFO_LIST, const SEXP DIST_FUNC_SEXP, const SEXP env, 
									 const arma::mat &DESIGN, const arma::rowvec &WT, arma::rowvec &R_PARA_OUT);
arma::rowvec directionalDerivative(const OBJ_INFO &OBJ, const arma::mat &dsGrid, const arma::mat &PARA_SET, const arma::rowvec &alpha, 
                                   const Rcpp::List MODEL_INFO_LIST, const SEXP DIST_FUNC_SEXP, const SEXP env);
arma::rowvec distCalc(const OBJ_INFO &OBJ, const arma::mat &x, const int &tid, const int &rid, const arma::mat &PARA_SET, 
                      const Rcpp::List MODEL_INFO_LIST, const SEXP DIST_FUNC_SEXP, const SEXP env);

// BODY
double DesignCriterion(const int &LOOPID, const PSO_OPTIONS PSO_OPTS[], const OBJ_INFO &OBJ, const Rcpp::List MODEL_INFO_LIST, const SEXP DIST_FUNC_SEXP, 
											 const SEXP env, const arma::rowvec &FIXEDVALUE, const rowvec &x, arma::mat &R_PARA)
{
	int nSupp = OBJ.nSupp;
	int dSupp = OBJ.dSupp;
	int d_type = OBJ.d_type;
	//Rprintf("get design\n");

  double val = 0.0;
	switch (d_type) {
		// Exact Design Module
		case 0:
		{
      arma::mat DESIGN(nSupp, dSupp, fill::zeros);
      for (int i = 0; i < nSupp; i++) { DESIGN.row(i) = x.subvec(i*dSupp, (i+1)*dSupp - 1); }
      arma::rowvec WT(nSupp);
			WT.fill(1.0/(double)nSupp);
      val = criterionList(LOOPID, PSO_OPTS, OBJ, MODEL_INFO_LIST, DIST_FUNC_SEXP, env, DESIGN, WT, R_PARA);  
      val *= -1.0;  
			break;
		}
		// Approximation Design Module
		case 1:
		{
      arma::mat DESIGN(nSupp, dSupp, fill::zeros);
      for (int i = 0; i < nSupp; i++) { DESIGN.row(i) = x.subvec(i*dSupp, (i+1)*dSupp - 1); }
      arma::rowvec WT(nSupp);
			// Get Design Weight
			arma::rowvec ang = x.subvec(nSupp*dSupp, nSupp*dSupp + nSupp - 2);
			arma::rowvec wcumsin(nSupp), wcos(nSupp);
			arma::rowvec wsin = arma::sin(ang);
			wcumsin(0) = 1.0; 
			for (int i = 1; i < nSupp; i++) { wcumsin(i) = wcumsin(i-1)*wsin(i-1); }
			wcos(nSupp - 1) = 1.0; wcos.subvec(0, nSupp - 2) = arma::cos(ang);
			WT = wcumsin % wcos;
			WT = WT % WT;
      val = criterionList(LOOPID, PSO_OPTS, OBJ, MODEL_INFO_LIST, DIST_FUNC_SEXP, env, DESIGN, WT, R_PARA);  
      val *= -1.0;  
			break;
		}
    // Find Optimal Weight for Equivalence Theorem on Max-min Optimal Design
    case 1001:
    {
      arma::mat DESIGN(nSupp, dSupp, fill::zeros);
      for (int i = 0; i < nSupp; i++) { DESIGN.row(i) = FIXEDVALUE.subvec(i*dSupp, (i+1)*dSupp - 1); }
      double CRIT_VAL = FIXEDVALUE(FIXEDVALUE.n_elem - 1);
      int n_model = (int)x.n_elem;
      arma::rowvec alpha(n_model + 1);
      arma::rowvec ang = x;
      arma::rowvec wcumsin(n_model + 1), wcos(n_model + 1);
      arma::rowvec wsin = arma::sin(ang);
      wcumsin(0) = 1.0; 
      for (int i = 1; i < (n_model + 1); i++) { wcumsin(i) = wcumsin(i-1)*wsin(i-1); }
      wcos(n_model) = 1.0; wcos.subvec(0, n_model - 1) = arma::cos(ang);
      alpha = wcumsin % wcos;
      alpha = alpha % alpha;
      arma::rowvec DIV = directionalDerivative(OBJ, DESIGN, OBJ.paras, alpha, MODEL_INFO_LIST, DIST_FUNC_SEXP, env);
      DIV -= CRIT_VAL;
      val = arma::accu(DIV % DIV);
      break;
    }
	}
	return val;
}

double criterionList(const int &LOOPID, const PSO_OPTIONS PSO_OPTS[],
                     const OBJ_INFO &OBJ, const Rcpp::List MODEL_INFO_LIST, const SEXP DIST_FUNC_SEXP, const SEXP env, 
										 const arma::mat &DESIGN, const arma::rowvec &WT, arma::mat &R_PARA)
{
	int crit_type = OBJ.crit_type;
	int N_model = OBJ.N_model;

  R_PARA.reset(); R_PARA.set_size(N_model, OBJ.dParas.max());
	double val = 1e20;
	switch (crit_type) {
		case 0: // Fixed True
		{
      R_PARA.submat(0, 0, 0, OBJ.dParas(0) - 1) = OBJ.paras.submat(0, 0, 0, OBJ.dParas(0) - 1);
			arma::rowvec R_PARA_tmp(OBJ.dParas(1));
			val = minDistCalc(LOOPID, PSO_OPTS, 0, 1, OBJ, MODEL_INFO_LIST, DIST_FUNC_SEXP, env, DESIGN, WT, R_PARA_tmp);
      R_PARA.submat(1, 0, 1, OBJ.dParas(1) - 1) = R_PARA_tmp;
			break;
		}
		case 1: // Max-min, Fixed True
		{
      R_PARA.submat(0, 0, 0, OBJ.dParas(0) - 1) = OBJ.paras.submat(0, 0, 0, OBJ.dParas(0) - 1);
			arma::rowvec std_vals = OBJ.std_vals;
			arma::rowvec eff_vals(N_model - 1);
			for (int i = 1; i < N_model; i++) {
				arma::rowvec R_PARA_tmp(OBJ.dParas(i));
				eff_vals(i-1) = minDistCalc(LOOPID, PSO_OPTS, 0, i, OBJ, MODEL_INFO_LIST, DIST_FUNC_SEXP, env, DESIGN, WT, R_PARA_tmp);	
        R_PARA.submat(i, 0, i, OBJ.dParas(i) - 1) = R_PARA_tmp;
			}
			eff_vals = eff_vals/std_vals;
			val = eff_vals.min();
			break;
		}
	}
	return val;
}

// SUBFUNCTIONS

// Minimal Distance Between Two Models
double minDistCalc(const int &LOOPID, const PSO_OPTIONS PSO_OPTS[], const int &tid, const int &rid, const OBJ_INFO &OBJ, 
                   const Rcpp::List MODEL_INFO_LIST, const SEXP DIST_FUNC_SEXP, const SEXP env, 
									 const arma::mat &DESIGN, const arma::rowvec &WT, arma::rowvec &R_PARA_OUT)
{
	arma::rowvec T_PARA = OBJ.paras.submat(tid, 0, tid, OBJ.dParas(tid) - 1);

	int dParas = OBJ.dParas(rid);
  arma::rowvec R_PARA_INI = OBJ.parasInit.submat(rid, 0, rid, dParas - 1);
  arma::rowvec R_UPPER = OBJ.parasUpper.submat(rid, 0, rid, dParas - 1);
  arma::rowvec R_LOWER = OBJ.parasLower.submat(rid, 0, rid, dParas - 1);

  arma::irowvec R_NBD(dParas); 
  arma::rowvec R_PARA_EX = R_PARA_INI;
  arma::rowvec R_PARA_EX_1 = R_PARA_INI;
  for (int d = 0; d < dParas; d++) { R_NBD(d) = (int)OBJ.parasBdd(rid, d); }

  double fx, fx1;
  int CONV_STATUS;
  CONV_STATUS = lbfgsKernel(LOOPID, PSO_OPTS, fx, R_PARA_EX, T_PARA, DESIGN, WT, tid, rid, MODEL_INFO_LIST, DIST_FUNC_SEXP, env, R_UPPER, R_LOWER, R_NBD);

  int LBFGS_RETRY = PSO_OPTS[LOOPID].LBFGS_RETRY + (1 - CONV_STATUS);
  int count = 0;
  while ((count < LBFGS_RETRY)) {
    R_PARA_EX_1 = randu(1, dParas) % (R_UPPER - R_LOWER) + R_LOWER;
    CONV_STATUS = lbfgsKernel(LOOPID, PSO_OPTS, fx1, R_PARA_EX_1, T_PARA, DESIGN, WT, tid, rid, MODEL_INFO_LIST, DIST_FUNC_SEXP, env, R_UPPER, R_LOWER, R_NBD);
    if (std::abs(fx1 - fx) > 0.1) { count--; }
    if (fx1 < fx) { count--; fx = fx1; R_PARA_EX = R_PARA_EX_1; } 
    count++;
  }
  R_PARA_OUT = R_PARA_EX;
	return fx;
}

// Equivalence Theorem
arma::rowvec directionalDerivative(const OBJ_INFO &OBJ, const arma::mat &dsGrid, const arma::mat &PARA_SET, const arma::rowvec &alpha, 
                                   const Rcpp::List MODEL_INFO_LIST, const SEXP DIST_FUNC_SEXP, const SEXP env)
{
  int crit_type = OBJ.crit_type;
  arma::rowvec dirDer(dsGrid.n_rows, fill::zeros);
  switch (crit_type) {
    case 0: // Fixed True
    {
      dirDer = distCalc(OBJ, dsGrid, 0, 1, PARA_SET, MODEL_INFO_LIST, DIST_FUNC_SEXP, env);
      break;
    }
    case 1: // Max-min, Fixed True
    {
      arma::rowvec std_vals = OBJ.std_vals;
      for (int i = 1; i < OBJ.N_model; i++) {
        arma::rowvec DIV = distCalc(OBJ, dsGrid, 0, i, PARA_SET, MODEL_INFO_LIST, DIST_FUNC_SEXP, env);
        dirDer += (alpha(i-1)/std_vals(i-1))*DIV;
      }
      break;
    }
  }
  return dirDer;
}

arma::rowvec distCalc(const OBJ_INFO &OBJ, const arma::mat &x, const int &tid, const int &rid, const arma::mat &PARA_SET, 
                      const Rcpp::List MODEL_INFO_LIST, const SEXP DIST_FUNC_SEXP, const SEXP env)
{
  Rcpp::EvalBase *distFunc = NULL;
  if (TYPEOF(DIST_FUNC_SEXP) == EXTPTRSXP) {   
    distFunc = new Rcpp::EvalCompiled(DIST_FUNC_SEXP, env); 
  } else {
    distFunc = new Rcpp::EvalStandard(DIST_FUNC_SEXP, env);  
  } 

  Rcpp::EvalBase *m1_func = NULL;
  SEXP tmp1 = as<SEXP>(MODEL_INFO_LIST[tid]);
  if (TYPEOF(tmp1) == EXTPTRSXP) {   
    m1_func = new Rcpp::EvalCompiled(tmp1, env);
  } else {                                                
    m1_func = new Rcpp::EvalStandard(tmp1, env);
  }  
  
  Rcpp::EvalBase *m2_func = NULL;
  SEXP tmp2 = as<SEXP>(MODEL_INFO_LIST[rid]);
  if (TYPEOF(tmp2) == EXTPTRSXP) {   
    m2_func = new Rcpp::EvalCompiled(tmp2, env);
  } else {                                                
    m2_func = new Rcpp::EvalStandard(tmp2, env);
  } 

  arma::rowvec eta_T(x.n_rows), eta_R(x.n_rows), DIV(x.n_rows);
  eta_T = (arma::rowvec) m1_func->eval(Rcpp::wrap(x), Rcpp::wrap(PARA_SET.submat(tid, 0, tid, OBJ.dParas(tid) - 1))); 
  eta_R = (arma::rowvec) m2_func->eval(Rcpp::wrap(x), Rcpp::wrap(PARA_SET.submat(rid, 0, rid, OBJ.dParas(rid) - 1))); 
  DIV = (arma::rowvec) distFunc->eval(Rcpp::wrap(eta_T), Rcpp::wrap(eta_R)); 
  delete m1_func; delete m2_func; delete distFunc; m1_func = NULL; m2_func = NULL; distFunc = NULL;
  
  return DIV;
}
