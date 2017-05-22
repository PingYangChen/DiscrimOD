
// INCLUDE HEADER FILES
#include "lbfgsKernel.h"

// DECLARE FUNCTIONS
double criterionList(const OBJ_INFO &OBJ, const MODEL_SET MODELS[], Rcpp::EvalBase* distFunc, 
										 const arma::mat &DESIGN, const arma::rowvec &WT, arma::rowvec &R_PARA);
double minDistCalc(const int rid, const OBJ_INFO &OBJ, const MODEL_SET MODELS[], Rcpp::EvalBase* distFunc, 
									 const arma::mat &DESIGN, const arma::rowvec &WT, arma::rowvec &R_PARA_OUT);
double paraTransform(const int INV, const double par, const int nbd, const double upper, const double lower);

// BODY
double DesignCriterion(const int &LOOPID, const PSO_OPTIONS PSO_OPTS[], const OBJ_INFO &OBJ, 
											 const MODEL_SET MODELS[], Rcpp::EvalBase* distFunc, const rowvec &p, const rowvec &x, arma::rowvec &R_PARA)
{
	int nSupp = OBJ.nSupp;
	int dSupp = OBJ.dSupp;
	//int d_type = OBJ.d_type;
	//Rprintf("get design\n");
	arma::mat DESIGN(nSupp, dSupp, fill::zeros);
	for (int i = 0; i < nSupp; i++) { DESIGN.row(i) = x.subvec(i*dSupp, (i+1)*dSupp - 1); }
	arma::rowvec WT(nSupp);

	//switch (d_type) {
		// Exact Design Module
		/*case 0:
		{
			WT.fill(1.0/(double)nSupp);
			break;
		}*/
		// Approximation Design Module
		//case 1:
		//{
			// Get Design Weight
			arma::rowvec ang = x.subvec(nSupp*dSupp, nSupp*dSupp + nSupp - 2);
			arma::rowvec wcumsin(nSupp), wcos(nSupp);
			arma::rowvec wsin = arma::sin(ang);
			wcumsin(0) = 1.0; 
			for (int i = 1; i < nSupp; i++) { wcumsin(i) = wcumsin(i-1)*wsin(i-1); }
			wcos(nSupp - 1) = 1.0; wcos.subvec(0, nSupp - 2) = arma::cos(ang);
			WT = wcumsin % wcos;
			WT = WT % WT;
			//break;
		//}
	//}
	double val = criterionList(OBJ, MODELS, distFunc, DESIGN, WT, R_PARA);		
	//R_PARA.set_size(OBJ.dParas(1));
	//val = minDistCalc(1, OBJ, MODELS, DESIGN, WT, R_PARA);		
	return (-1.0)*val;
}

double criterionList(const OBJ_INFO &OBJ, const MODEL_SET MODELS[], Rcpp::EvalBase* distFunc, 
										 const arma::mat &DESIGN, const arma::rowvec &WT, arma::rowvec &R_PARA)
{
	int crit_type = OBJ.crit_type;
	int N_model = OBJ.N_model;

	double val = 1e20;
	switch (crit_type) {
		case 0: // Fixed True
		{
			R_PARA.reset(); R_PARA.set_size(OBJ.dParas(1));
			val = minDistCalc(1, OBJ, MODELS, distFunc, DESIGN, WT, R_PARA);
			break;
		}
		case 1: // Max-min, Fixed True
		{
			arma::rowvec std_vals = OBJ.std_vals;
			arma::rowvec eff_vals(N_model - 1);
			for (int i = 1; i < N_model; i++) {
				R_PARA.reset();	R_PARA.set_size(OBJ.dParas(i));
				eff_vals(i) = minDistCalc(i, OBJ, MODELS, distFunc, DESIGN, WT, R_PARA);	
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
double minDistCalc(const int rid, const OBJ_INFO &OBJ, const MODEL_SET MODELS[], Rcpp::EvalBase* distFunc, 
									 const arma::mat &DESIGN, const arma::rowvec &WT, arma::rowvec &R_PARA_OUT)
{
	arma::rowvec T_PARA = OBJ.paras.submat(0, 0, 0, OBJ.dParas(0) - 1);

	int dParas = OBJ.dParas(rid);
  arma::rowvec R_PARA_INI = OBJ.parasInit.submat(rid, 0, rid, dParas - 1);
  arma::rowvec R_UPPER = OBJ.parasUpper.submat(rid, 0, rid, dParas - 1);
  arma::rowvec R_LOWER = OBJ.parasLower.submat(rid, 0, rid, dParas - 1);

  Rcpp::EvalBase *m1_func = MODELS[0].modelFunc;
  Rcpp::EvalBase *m2_func = MODELS[rid].modelFunc;
  
  arma::irowvec R_NBD(dParas); arma::rowvec R_PARA_EX(dParas), R_PARA_EX_1(dParas);
  for (int d = 0; d < dParas; d++) { 
    R_NBD(d) = (int)OBJ.parasBdd(rid, d); 
    R_PARA_EX(d) = paraTransform(0, R_PARA_INI(d), R_NBD(d), R_UPPER(d), R_LOWER(d));
  }

  double fx, fx1;
  fx = lbfgsKernel(R_PARA_EX, T_PARA, DESIGN, WT, m1_func, m2_func, distFunc, R_UPPER, R_LOWER, R_NBD);

  int count = 0;
  while ((count < 2)) {
    R_PARA_EX_1 = randu(1, dParas) % (R_UPPER - R_LOWER) + R_LOWER;
    for (int d = 0; d < dParas; d++) { R_PARA_EX_1(d) = paraTransform(0, R_PARA_EX_1(d), R_NBD(d), R_UPPER(d), R_LOWER(d)); }
    fx1 = lbfgsKernel(R_PARA_EX_1, T_PARA, DESIGN, WT, m1_func, m2_func, distFunc, R_UPPER, R_LOWER, R_NBD);
    if (fx1 < fx) { count--; fx = fx1; R_PARA_EX = R_PARA_EX_1; }
    count++;
  }
  for (int d = 0; d < dParas; d++) R_PARA_OUT(d) = paraTransform(1, R_PARA_EX(d), R_NBD(d), R_UPPER(d), R_LOWER(d));
	return fx;
}

