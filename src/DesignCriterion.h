
// INCLUDE HEADER FILES
#include "lbfgsKernel.h"

// DECLARE FUNCTIONS
double criterionList(const int &LOOPID, const PSO_OPTIONS PSO_OPTS[], const OBJ_INFO &OBJ, model_diff_func *MODEL_COLLECTOR[],  
                     const arma::mat &DESIGN, const arma::rowvec &WT, arma::mat &R_PARA);
double minDistCalc(const int &LOOPID, const PSO_OPTIONS PSO_OPTS[], const OBJ_INFO &OBJ, model_diff_func *MODEL_COLLECTOR[], const int &PAIRID, 
                   const arma::mat &DESIGN, const arma::rowvec &WT, arma::rowvec &R_PARA_OUT);
arma::rowvec directionalDerivative(const OBJ_INFO &OBJ, const arma::mat &dsGrid, const arma::mat &PARA_SET, const arma::rowvec &alpha, 
                                   model_diff_func *MODEL_COLLECTOR[]);
arma::rowvec distCalc(const OBJ_INFO &OBJ, const arma::mat &x, const arma::mat &PARA_SET, 
                      model_diff_func *MODEL_COLLECTOR[], const int &PAIRID);

// BODY
double DesignCriterion(const int &LOOPID, const PSO_OPTIONS PSO_OPTS[], const OBJ_INFO &OBJ, model_diff_func *MODEL_COLLECTOR[], 
                       const arma::rowvec &FIXEDVALUE, const rowvec &x, arma::mat &R_PARA)
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
      arma::rowvec WT(nSupp, fill::zeros);
			WT.fill(1.0/(double)nSupp);
      val = criterionList(LOOPID, PSO_OPTS, OBJ, MODEL_COLLECTOR, DESIGN, WT, R_PARA);  
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
			arma::rowvec wcumsin(nSupp, fill::zeros), wcos(nSupp, fill::zeros);
			arma::rowvec wsin = arma::sin(ang);
			wcumsin(0) = 1.0; 
			for (int i = 1; i < nSupp; i++) { wcumsin(i) = wcumsin(i-1)*wsin(i-1); }
			wcos(nSupp - 1) = 1.0; wcos.subvec(0, nSupp - 2) = arma::cos(ang);
			WT = wcumsin % wcos;
			WT = WT % WT;
      val = criterionList(LOOPID, PSO_OPTS, OBJ, MODEL_COLLECTOR, DESIGN, WT, R_PARA);  
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
      arma::rowvec alpha(n_model + 1, fill::zeros);
      arma::rowvec ang = x;
      arma::rowvec wcumsin(n_model + 1, fill::zeros), wcos(n_model + 1, fill::zeros);
      arma::rowvec wsin = arma::sin(ang);
      wcumsin(0) = 1.0; 
      for (int i = 1; i < (n_model + 1); i++) { wcumsin(i) = wcumsin(i-1)*wsin(i-1); }
      wcos(n_model) = 1.0; wcos.subvec(0, n_model - 1) = arma::cos(ang);
      alpha = wcumsin % wcos;
      alpha = alpha % alpha;
      arma::rowvec DIV = directionalDerivative(OBJ, DESIGN, OBJ.paras, alpha, MODEL_COLLECTOR);
      DIV -= CRIT_VAL;
      val = arma::accu(DIV % DIV);
      break;
    }
	}
	return val;
}

double criterionList(const int &LOOPID, const PSO_OPTIONS PSO_OPTS[], const OBJ_INFO &OBJ, model_diff_func *MODEL_COLLECTOR[],  
										 const arma::mat &DESIGN, const arma::rowvec &WT, arma::mat &R_PARA)
{
	int crit_type = OBJ.crit_type;
	int N_PAIR = OBJ.N_PAIR;

  R_PARA.reset(); R_PARA.set_size(OBJ.dParas.n_elem, OBJ.dParas.max());
	double val = 1e20;
	switch (crit_type) {
		case 0: // Fixed True
		{
      R_PARA.submat(0, 0, 0, OBJ.dParas(0) - 1) = OBJ.paras.submat(0, 0, 0, OBJ.dParas(0) - 1);
			arma::rowvec R_PARA_tmp(OBJ.dParas(1));
			val = minDistCalc(LOOPID, PSO_OPTS, OBJ, MODEL_COLLECTOR, 0, DESIGN, WT, R_PARA_tmp);
      R_PARA.submat(1, 0, 1, OBJ.dParas(1) - 1) = R_PARA_tmp;
			break;
		}
		case 1: // Max-min, Fixed True
		{
      R_PARA.submat(0, 0, 0, OBJ.dParas(0) - 1) = OBJ.paras.submat(0, 0, 0, OBJ.dParas(0) - 1);
			arma::rowvec std_vals = OBJ.std_vals;
			arma::rowvec eff_vals(N_PAIR, fill::zeros);
			for (int i = 0; i < N_PAIR; i++) {
				arma::rowvec R_PARA_tmp(OBJ.dParas(i+1));
				eff_vals(i) = minDistCalc(LOOPID, PSO_OPTS, OBJ, MODEL_COLLECTOR, i, DESIGN, WT, R_PARA_tmp);	
        R_PARA.submat(i+1, 0, i+1, OBJ.dParas(i+1) - 1) = R_PARA_tmp;
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
double minDistCalc(const int &LOOPID, const PSO_OPTIONS PSO_OPTS[], const OBJ_INFO &OBJ, model_diff_func *MODEL_COLLECTOR[], const int &PAIRID, 
									 const arma::mat &DESIGN, const arma::rowvec &WT, arma::rowvec &R_PARA_OUT)
{
  model_diff_func* func_input = MODEL_COLLECTOR[PAIRID];

  arma::imat MODEL_PAIR = OBJ.MODEL_PAIR;
  int tmID = MODEL_PAIR(PAIRID, 0);
  int rmID = MODEL_PAIR(PAIRID, 1);

	arma::rowvec T_PARA = OBJ.paras.submat(tmID, 0, tmID, OBJ.dParas(tmID) - 1);

	int dParas = OBJ.dParas(rmID);
  arma::rowvec R_PARA_INI = OBJ.parasInit.submat(rmID, 0, rmID, dParas - 1);
  arma::rowvec R_UPPER = OBJ.parasUpper.submat(rmID, 0, rmID, dParas - 1);
  arma::rowvec R_LOWER = OBJ.parasLower.submat(rmID, 0, rmID, dParas - 1);

  arma::irowvec R_NBD(dParas); 
  arma::rowvec R_PARA_EX = R_PARA_INI;
  arma::rowvec R_PARA_EX_1 = R_PARA_INI;
  for (int d = 0; d < dParas; d++) { R_NBD(d) = (int)OBJ.parasBdd(rmID, d); }

  double fx = 1e20, fx1 = 1e20;
  int CONV_STATUS;
  CONV_STATUS = lbfgsKernel(LOOPID, PSO_OPTS, fx, R_PARA_EX, T_PARA, DESIGN, WT, func_input, R_UPPER, R_LOWER, R_NBD);

  int LBFGS_RETRY = PSO_OPTS[LOOPID].LBFGS_RETRY + (1 - CONV_STATUS);
  if ((!std::isfinite(fx)) | std::isnan(fx)) { LBFGS_RETRY++; }
  int count = 0;
  while ((count < LBFGS_RETRY)) {
    R_PARA_EX_1 = randu(1, dParas) % (R_UPPER - R_LOWER) + R_LOWER;
    CONV_STATUS = lbfgsKernel(LOOPID, PSO_OPTS, fx1, R_PARA_EX_1, T_PARA, DESIGN, WT, func_input, R_UPPER, R_LOWER, R_NBD);
    //if ((std::abs(fx1 - fx) > 0.1) | (CONV_STATUS == 0)) { count--; }
    if (std::isfinite(fx1) & (!std::isnan(fx1))) {
      if ((!std::isfinite(fx)) | std::isnan(fx)) { 
        fx = fx1; R_PARA_EX = R_PARA_EX_1; 
      }  else {
        if (fx1 < fx) { fx = fx1; R_PARA_EX = R_PARA_EX_1; }     
      }
    } 
    count++;
  }
  R_PARA_OUT = R_PARA_EX;
	return fx;
}

// Equivalence Theorem
arma::rowvec directionalDerivative(const OBJ_INFO &OBJ, const arma::mat &dsGrid, const arma::mat &PARA_SET, const arma::rowvec &alpha, 
                                   model_diff_func *MODEL_COLLECTOR[])
{
  int crit_type = OBJ.crit_type;
  arma::rowvec dirDer(dsGrid.n_rows, fill::zeros);
  switch (crit_type) {
    case 0: // Fixed True
    {
      dirDer = distCalc(OBJ, dsGrid, PARA_SET, MODEL_COLLECTOR, 0);
      break;
    }
    case 1: // Max-min, Fixed True
    {
      arma::rowvec std_vals = OBJ.std_vals;
      for (int i = 0; i < OBJ.N_PAIR; i++) {
        arma::rowvec DIV = distCalc(OBJ, dsGrid, PARA_SET, MODEL_COLLECTOR, i);
        dirDer += (alpha(i)/std_vals(i))*DIV;
      }
      break;
    }
  }
  return dirDer;
}

arma::rowvec distCalc(const OBJ_INFO &OBJ, const arma::mat &x, const arma::mat &PARA_SET, 
                      model_diff_func *MODEL_COLLECTOR[], const int &PAIRID)
{
  model_diff_func* func_input = MODEL_COLLECTOR[PAIRID];
  Rcpp::EvalBase *m1_func = (Rcpp::EvalBase *) func_input->M1_FUNC;
  Rcpp::EvalBase *m2_func = (Rcpp::EvalBase *) func_input->M2_FUNC;
  Rcpp::EvalBase *distFunc = (Rcpp::EvalBase *) func_input->DISTFUNC;

  arma::imat MODEL_PAIR = OBJ.MODEL_PAIR;
  int tmID = MODEL_PAIR(PAIRID, 0);
  int rmID = MODEL_PAIR(PAIRID, 1);

  arma::rowvec eta_T(x.n_rows, fill::zeros), eta_R(x.n_rows, fill::zeros), DIV(x.n_rows, fill::zeros);
  eta_T = (arma::rowvec) m1_func->eval(Rcpp::wrap(x), Rcpp::wrap(PARA_SET.submat(tmID, 0, tmID, OBJ.dParas(tmID) - 1))); 
  eta_R = (arma::rowvec) m2_func->eval(Rcpp::wrap(x), Rcpp::wrap(PARA_SET.submat(rmID, 0, rmID, OBJ.dParas(rmID) - 1))); 
  DIV = (arma::rowvec) distFunc->eval(Rcpp::wrap(eta_T), Rcpp::wrap(eta_R)); 
    
  return DIV;
}
