
// INCLUDE HEADER FILES
#include <malloc.h>
#include "lbfgsKernel.h"

// DEFINE STUCTURES
typedef struct {
  Rcpp::EvalBase* dist_func;
  Rcpp::EvalBase* t_func;
  Rcpp::EvalBase* r_func;
  arma::rowvec T_PARA;
  arma::mat DESIGN;
  arma::rowvec WT;
  double* rUpper;
  double* rLower;
  int* rNBD;
} LBFGS_EX, *Ptr_LBFGS_EX;


// DECLARE FUNCTIONS
double criterionList(const OBJ_INFO &OBJ, const MODEL_SET MODELS[], const arma::mat &DESIGN, const arma::rowvec &WT, arma::rowvec &R_PARA);
lbfgsfloatval_t evaluate(void *ex, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step);
double fminfn(int n, double *par, void *ex);
void fmingr(int n, double *par, lbfgsfloatval_t *gr, void *ex);
double minDistCalc(const int rid, const OBJ_INFO &OBJ,  const MODEL_SET MODELS[], const arma::mat &DESIGN, const arma::rowvec &WT, 
									 arma::rowvec &R_PARA_OUT);
double paraTransform(const int INV, const double par, const int nbd, const double upper, const double lower);

// BODY
double DesignCriterion(const int &LOOPID, const PSO_OPTIONS PSO_OPTS[], const OBJ_INFO &OBJ,
											 const MODEL_SET MODELS[], const rowvec &p, const rowvec &x, arma::rowvec &R_PARA)
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
	double val = criterionList(OBJ, MODELS, DESIGN, WT, R_PARA);		
	//R_PARA.set_size(OBJ.dParas(1));
	//val = minDistCalc(1, OBJ, MODELS, DESIGN, WT, R_PARA);		
	return (-1.0)*val;
}

double criterionList(const OBJ_INFO &OBJ, const MODEL_SET MODELS[], const arma::mat &DESIGN, const arma::rowvec &WT, arma::rowvec &R_PARA) 
{
	int crit_type = OBJ.crit_type;
	int N_model = OBJ.N_model;

	double val = 1e20;
	switch (crit_type) {
		case 0: // Fixed True
		{
			R_PARA.reset(); R_PARA.set_size(OBJ.dParas(1));
			val = minDistCalc(1, OBJ, MODELS, DESIGN, WT, R_PARA);
			break;
		}
		case 1: // Max-min, Fixed True
		{
			arma::rowvec std_vals = OBJ.std_vals;
			arma::rowvec eff_vals(N_model - 1);
			for (int i = 1; i < N_model; i++) {
				R_PARA.reset();	R_PARA.set_size(OBJ.dParas(i));
				eff_vals(i) = minDistCalc(i, OBJ, MODELS, DESIGN, WT, R_PARA);	
			}
			eff_vals = eff_vals/std_vals;
			val = eff_vals.min();
			break;
		}
	}
	return val;
}

// SUBFUNCTIONS
lbfgsfloatval_t evaluate(void *ex, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step)
{
	lbfgsfloatval_t fx = fminfn(n, (double*)x, ex);
	fmingr(n, (double*)x, g, ex);
  return fx;
}

// ====== The Objective function used in Quasi-Newton Algorithm ====== //
double fminfn(int n, double *par, void *ex)
{
  LBFGS_EX exList = *(LBFGS_EX*)(ex);
  arma::rowvec R_PARA(n);
  for (int q = 0; q < n; q++) { R_PARA(q) = paraTransform(1, par[q], exList.rNBD[q], exList.rUpper[q], exList.rLower[q]); }

  arma::rowvec T_PARA = exList.T_PARA;
  arma::mat DESIGN = exList.DESIGN;
	int nSupp = DESIGN.n_rows; //int dSupp = DESIGN.n_cols;
	arma::rowvec WT = exList.WT;

  Rcpp::EvalBase* t_func = (Rcpp::EvalBase *) exList.t_func;
	Rcpp::EvalBase* r_func = (Rcpp::EvalBase *) exList.r_func;
	Rcpp::EvalBase* dist_func = (Rcpp::EvalBase *) exList.dist_func;

  double fvalTmp = 0;
  double eta_T, eta_R, DIV;
  for (int i = 0; i < nSupp; i++) {
  	arma::rowvec x = DESIGN.row(i);
		eta_T = (double) t_func->eval(Rcpp::wrap(x), Rcpp::wrap(T_PARA)); 
  	eta_R = (double) r_func->eval(Rcpp::wrap(x), Rcpp::wrap(R_PARA)); 
  	DIV = (double) dist_func->eval(Rcpp::wrap(eta_T), Rcpp::wrap(eta_R));  
    fvalTmp += WT(i)*DIV;
  }
  if (std::isnan(fvalTmp)) { fvalTmp = 1e20; }
  if (!(arma::is_finite(fvalTmp))) { fvalTmp = 1e20; }
  return fvalTmp;
}

// ====== The Gradient function used in Quasi-Newton Algorithm ====== //
void fmingr(int n, double *par, lbfgsfloatval_t *gr, void *ex)
{
  double par1[n], par0[n];
  double fv, bv;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) { par1[j] = par[j]; par0[j] = par[j]; }
    par1[i] = par1[i] + 5e-4; par0[i] = par0[i] - 5e-4;
    //central difference
    fv = fminfn(n, par1, ex); bv = fminfn(n, par0, ex);
    gr[i] = fv*1e3 - bv*1e3;
  }
}

// Minimal Distance Between Two Models
double minDistCalc(const int rid, const OBJ_INFO &OBJ, const MODEL_SET MODELS[], 
									 const arma::mat &DESIGN, const arma::rowvec &WT, arma::rowvec &R_PARA_OUT)
{
	arma::rowvec T_PARA = OBJ.paras.submat(0, 0, 0, OBJ.dParas(0) - 1);

	int len_rPara = OBJ.dParas(rid);
  rowvec R_PARA_INI(len_rPara);
  double lower[len_rPara], upper[len_rPara]; int nbd[len_rPara];

  for (int d = 0; d < len_rPara; d++) {
    upper[d] = OBJ.parasUpper(rid, d);	lower[d] = OBJ.parasLower(rid, d);
    nbd[d] = (int)OBJ.parasBdd(rid, d);
    R_PARA_INI(d) = OBJ.parasInit(rid, d);
  }

	lbfgsfloatval_t *rpar = lbfgs_malloc(len_rPara);
	lbfgsfloatval_t *rpar1 = lbfgs_malloc(len_rPara);

	lbfgs_parameter_t LBFGS_PAR;
	lbfgs_parameter_init(&LBFGS_PAR);
	LBFGS_PAR.epsilon = 1e-10; // default 1e-5
	LBFGS_PAR.gtol = 1e-1; // default 0.9
	LBFGS_PAR.ftol = 1e-6; // default 1e-4
	LBFGS_PAR.delta = 1e-6; // default 0

  LBFGS_EX ex = {};
  ex.t_func = MODELS[0].model_func; ex.r_func = MODELS[rid].model_func; ex.dist_func = OBJ.dist_func;
	ex.T_PARA = T_PARA; ex.DESIGN = DESIGN; ex.WT = WT;
  ex.rUpper = upper; ex.rLower = lower; ex.rNBD = nbd;

  for (int d = 0; d < len_rPara; d++) rpar[d] = paraTransform(0, R_PARA_INI(d), nbd[d], upper[d], lower[d]);

  //int ret;
  lbfgsfloatval_t fx, fx1;
  //ret = lbfgs(D_INFO.drPara(rid), rpar, &fx, evaluate, NULL, &ex, &LBFGS_PAR);
  lbfgs(len_rPara, rpar, &fx, evaluate, NULL, &ex, &LBFGS_PAR);

  int count = 0;
  while ((count < 1)) {
    for (int d = 0; d < len_rPara; d++) {
      rpar1[d] = as_scalar(randu(1))*(upper[d]-lower[d])+lower[d];
      rpar1[d] = paraTransform(0, rpar1[d], nbd[d], upper[d], lower[d]);
    }
    //ret = lbfgs(D_INFO.drPara(rid), rpar1, &fx1, evaluate, NULL, &ex, &LBFGS_PAR);
    lbfgs(len_rPara, rpar1, &fx1, evaluate, NULL, &ex, &LBFGS_PAR);
    if (fx1 < fx) {
      count--;
      fx = fx1; for (int d = 0; d < len_rPara; d++) { rpar[d] = rpar1[d]; }
    }
    count++;
  }
	for (int d = 0; d < len_rPara; d++) R_PARA_OUT(d) = paraTransform(1, rpar[d], nbd[d], upper[d], lower[d]);

	//Fmin[0] = (double)fx;
	lbfgs_free(rpar); lbfgs_free(rpar1);
	return (double)fx;
}

//
double paraTransform(const int INV, const double par, const int nbd, const double upper, const double lower)
{
  double tmp; double out = par;
  if (INV == 0) {
    switch(nbd) {
      case 1: { tmp = std::log((par - lower)); out = tmp; break; } // Lower
      case 2: { tmp = (par - lower)/(upper - lower); tmp = std::log((tmp/(1-tmp))); out = tmp; break; } // Both
      case 3: { tmp = std::log((upper - par)); out = tmp; break; } // Upper
    }
  } else {
    switch(nbd) {
      case 1: { tmp = std::exp(par) + lower; out = tmp; break; }
      case 2: { tmp = std::exp(par); tmp = tmp/(1+tmp); tmp = tmp*(upper - lower) + lower; out = tmp; break; }
      case 3: { tmp = upper - std::exp(par); out = tmp; break; }
    }
  }
  return out;
}
