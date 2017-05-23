
// DECLARE FUNCTIONS

double f_fn(const int dPara, const arma::rowvec &R_PARA_EX, const arma::rowvec &T_PARA, const arma::mat &DESIGN, const arma::rowvec &WT,
						Rcpp::EvalBase* m1_func, Rcpp::EvalBase* m2_func, Rcpp::EvalBase* dist_func,
						const arma::rowvec &R_UPPER, const arma::rowvec &R_LOWER, const arma::irowvec &R_NBD);
arma::rowvec f_gr(const int dPara, const arma::rowvec &R_PARA_EX, const arma::rowvec &T_PARA, const arma::mat &DESIGN, const arma::rowvec &WT,
									Rcpp::EvalBase* m1_func, Rcpp::EvalBase* m2_func, Rcpp::EvalBase* dist_func,
									const arma::rowvec &R_UPPER, const arma::rowvec &R_LOWER, const arma::irowvec &R_NBD);
double paraTransform(const int INV, const double par, const int nbd, const double upper, const double lower);

// BODY
int lbfgsKernel(const int &LOOPID, const PSO_OPTIONS PSO_OPTS[], double &FVAL, arma::rowvec &R_PARA_EX, 
	 						  const arma::rowvec &T_PARA, const arma::mat &DESIGN, const arma::rowvec &WT,
								Rcpp::EvalBase* m1_func, Rcpp::EvalBase* m2_func, Rcpp::EvalBase* dist_func,
								const arma::rowvec &R_UPPER, const arma::rowvec &R_LOWER, const arma::irowvec &R_NBD)
{
	// LBFGS SETTINGS (Temporary)
	int LBFGS_MAXIT = PSO_OPTS[LOOPID].LBFGS_MAXIT;
	// CORRECTION SETTINGS
	int LBFGS_LM = PSO_OPTS[LOOPID].LBFGS_LM;
	arma::rowvec ai(LBFGS_LM), bi(LBFGS_LM);
	// STOPPING CRITERION SETTING
	double FVAL_EPS = PSO_OPTS[LOOPID].FVAL_EPS; // DEFAULT 0.0
	double GRAD_EPS = PSO_OPTS[LOOPID].GRAD_EPS; // DEFAULT 1e-5

	// LINE SEARCH PARAMETERS
	int LINESEARCH_MAXTRIAL = PSO_OPTS[LOOPID].LINESEARCH_MAXTRIAL;
	double LINESEARCH_MAX = PSO_OPTS[LOOPID].LINESEARCH_MAX;
	double LINESEARCH_MIN = PSO_OPTS[LOOPID].LINESEARCH_MIN;
	double LINESEARCH_C = PSO_OPTS[LOOPID].LINESEARCH_C;
	double LINESEARCH_TAU = PSO_OPTS[LOOPID].LINESEARCH_TAU;
	double alpha = LINESEARCH_MAX;

	// Initialization
	int dPara = R_PARA_EX.n_elem;
	// New Position
	arma::rowvec R_PARA_EX_NEW(dPara);
	// Objective Function Value
	double FVAL_NEW;
	FVAL = f_fn(dPara, R_PARA_EX, T_PARA, DESIGN, WT, m1_func, m2_func, dist_func, R_UPPER, R_LOWER, R_NBD); 
	// Gradient
	arma::rowvec GRAD(dPara), GRAD_NEW(dPara);	
	GRAD = f_gr(dPara, R_PARA_EX, T_PARA, DESIGN, WT, m1_func, m2_func, dist_func, R_UPPER, R_LOWER, R_NBD); 
	// Hessian Matrix
	arma::mat EYE(dPara, dPara); EYE.eye();
	arma::mat H_0 = EYE;
	// Direction Correction
	arma::rowvec DIR(dPara);
	arma::rowvec S_ONE(dPara), Y_ONE(dPara);
	arma::mat S_MAT(LBFGS_LM, dPara, fill::zeros), Y_MAT(LBFGS_LM, dPara);
	double RHO = 1.0; arma::rowvec RHO_VEC(LBFGS_LM);
	// Stopping Criterion
	double FVAL_TSET = 1e20, GRAD_TEST = 1e20;
	bool LBFGS_STOP = FALSE; int CONV = 0;
	// START L-BFGS LOOP
	for (int t = 0; t < LBFGS_MAXIT; t++) {
		// Update Directions
		if (t == 0) {
			DIR = (-1.0)*(GRAD * H_0);
		} else if (t < LBFGS_LM) {
			DIR = GRAD; ai.zeros(); bi.zeros();
			for (int i = 0; i < t; i++) {
				ai(i) = RHO_VEC(t - 1 - i)*arma::as_scalar(DIR * S_MAT.row(t - 1 - i).t()); // newest to oldest
				DIR -= ai(i)*Y_MAT.row(t - 1 - i);
			}
			for (int i = 0; i < t; i++) {
				bi(i) = RHO_VEC(i)*arma::as_scalar(DIR * Y_MAT.row(i).t()); // oldest to newest
				DIR += (ai(t - 1 - i) - bi(i))*S_MAT.row(i);
			}
			DIR *= -1.0;
		} else {
			DIR = GRAD; ai.zeros(); bi.zeros();
			for (int i = 0; i < LBFGS_LM; i++) {
				ai(i) = RHO_VEC(LBFGS_LM - 1 - i)*arma::as_scalar(DIR * S_MAT.row(LBFGS_LM - 1 - i).t()); // newest to oldest
				DIR -= ai(i)*Y_MAT.row(LBFGS_LM - 1 - i);
			}
			for (int i = 0; i < LBFGS_LM; i++) {
				bi(i) = RHO_VEC(i)*arma::as_scalar(DIR * Y_MAT.row(i).t()); // oldest to newest
				DIR += (ai(LBFGS_LM - 1 - i) - bi(i))*S_MAT.row(i);
			}
			DIR *= -1.0;
		}
		// Perform the Backtracking Line Search for Suitable Step Size 
		double DIR_VAL = arma::as_scalar(DIR * GRAD.t());
		if (DIR_VAL < 0) {
			alpha = LINESEARCH_MAX;
			double LS_VAL; int LS_COUNTER = 0; bool FLAG = TRUE;
			while ((LS_COUNTER < LINESEARCH_MAXTRIAL) & FLAG) {
				arma::rowvec R_PARA_EX_FOR = R_PARA_EX + alpha*DIR;
				LS_VAL = f_fn(dPara, R_PARA_EX_FOR, T_PARA, DESIGN, WT, m1_func, m2_func, dist_func, R_UPPER, R_LOWER, R_NBD); 
				if (LS_VAL > (FVAL + alpha*LINESEARCH_C*DIR_VAL)) {
					alpha *= LINESEARCH_TAU;
					if (alpha < LINESEARCH_MIN) { alpha = LINESEARCH_MIN; FLAG = FALSE; }
				} else {
					FLAG = FALSE;
				}
				LS_COUNTER++;
			}
		} else {
			LBFGS_STOP = TRUE;
		} 

		if (LBFGS_STOP) { t = LBFGS_MAXIT + 1; } else {
			// Update New Position
			R_PARA_EX_NEW = R_PARA_EX + alpha*DIR;
			FVAL_NEW = f_fn(dPara, R_PARA_EX_NEW, T_PARA, DESIGN, WT, m1_func, m2_func, dist_func, R_UPPER, R_LOWER, R_NBD); 
			GRAD_NEW = f_gr(dPara, R_PARA_EX_NEW, T_PARA, DESIGN, WT, m1_func, m2_func, dist_func, R_UPPER, R_LOWER, R_NBD); 
			// Update S and Y
			S_ONE = R_PARA_EX_NEW - R_PARA_EX;
			Y_ONE = GRAD_NEW - GRAD;
			RHO = 1.0/as_scalar(S_ONE*Y_ONE.t());
			// Update Memory
			if (t < LBFGS_LM) {
			 	S_MAT.row(t) = S_ONE; Y_MAT.row(t) = Y_ONE; RHO_VEC(t) = RHO;
			} else {
			 	arma::mat S_TMP = S_MAT; arma::mat Y_TMP = Y_MAT; arma::rowvec RHO_TMP = RHO_VEC;
			 	S_MAT.rows(0, LBFGS_LM - 2) = S_TMP.rows(1, LBFGS_LM - 1); S_MAT.row(LBFGS_LM - 1) = S_ONE;
			 	Y_MAT.rows(0, LBFGS_LM - 2) = Y_TMP.rows(1, LBFGS_LM - 1); Y_MAT.row(LBFGS_LM - 1) = Y_ONE;
			 	RHO_VEC.subvec(0, LBFGS_LM - 2) = RHO_TMP.subvec(1, LBFGS_LM - 1); RHO_VEC(LBFGS_LM - 1) = RHO;
			}
			// Examine Stopping Criterion
			if (FVAL_EPS > 0) {
				FVAL_TSET = std::abs(FVAL - FVAL_NEW)/FVAL_NEW;
				if (FVAL_TSET < FVAL_EPS) { LBFGS_STOP = TRUE; CONV = 1; }	
			}
			if (GRAD_EPS > 0) {
				GRAD_TEST = arma::as_scalar(GRAD_NEW * GRAD_NEW.t());
				double XNORM = arma::as_scalar(R_PARA_EX_NEW * R_PARA_EX_NEW.t());
				if (XNORM < 1.0) { XNORM = 1.0; }
				GRAD_TEST = GRAD_TEST/XNORM;
				if (GRAD_TEST < GRAD_EPS) { LBFGS_STOP = TRUE; CONV = 1; }	
			}
			// Update Approximated Hessian Matrix
		  H_0 = (EYE - RHO*(Y_ONE.t() * S_ONE))	* H_0 * (EYE - RHO*(S_ONE.t() * Y_ONE)) + RHO*(S_ONE.t() * S_ONE);
			// Replace Current Values
			R_PARA_EX = R_PARA_EX_NEW;
			FVAL = FVAL_NEW;
			GRAD = GRAD_NEW;
		}
		// DEBUGGING 
		/*
		rowvec CHECK_PARA(dPara);
		for (int q = 0; q < dPara; q++) { CHECK_PARA(q) = paraTransform(1, R_PARA_EX(q), R_NBD(q), R_UPPER(q), R_LOWER(q)); } 
		Rprintf("Iteration %d :\n", t);
		Rprintf("F = %4.8f\n", FVAL);
		Rprintf("G = "); for (int q = 0; q < dPara; q++) { Rprintf("%4.4f", GRAD(q)); if (q < (dPara - 1)) Rprintf(", "); else Rprintf("\n"); } 
		Rprintf("alpha = %4.4f\n", alpha);
		Rprintf("DIR = "); for (int q = 0; q < dPara; q++) { Rprintf("%4.4f", DIR(q)); if (q < (dPara - 1)) Rprintf(", "); else Rprintf("\n"); } 
		Rprintf("P = "); for (int q = 0; q < dPara; q++) { Rprintf("%4.4f", CHECK_PARA(q)); if (q < (dPara - 1)) Rprintf(", "); else Rprintf("\n"); } 
		Rprintf("RHO = %4.4f\n", RHO);
		Rprintf("S = "); for (int q = 0; q < dPara; q++) { Rprintf("%4.4f", S_ONE(q)); if (q < (dPara - 1)) Rprintf(", "); else Rprintf("\n"); } 
		Rprintf("Y = "); for (int q = 0; q < dPara; q++) { Rprintf("%4.4f", Y_ONE(q)); if (q < (dPara - 1)) Rprintf(", "); else Rprintf("\n"); } 
		Rprintf("\n");
		*/
	}
	// END L-BFGS LOOP
	//Rprintf("\n");
	return CONV;
}

// SUBFUNCTIONS
double f_fn(const int dPara, const arma::rowvec &R_PARA_EX, const arma::rowvec &T_PARA, const arma::mat &DESIGN, const arma::rowvec &WT,
						Rcpp::EvalBase* m1_func, Rcpp::EvalBase* m2_func, Rcpp::EvalBase* dist_func,
						const arma::rowvec &R_UPPER, const arma::rowvec &R_LOWER, const arma::irowvec &R_NBD)
{
	// Transform to original scale
	arma::rowvec R_PARA(dPara);
	for (int q = 0; q < dPara; q++) { R_PARA(q) = paraTransform(1, R_PARA_EX(q), R_NBD(q), R_UPPER(q), R_LOWER(q)); }  // (-Inf, Inf) -> Original
	//int nSupp = DESIGN.n_rows;
  double fvalTmp = 0;
  arma::rowvec eta_T, eta_R, DIV;

  eta_T = (arma::rowvec) m1_func->eval(Rcpp::wrap(DESIGN), Rcpp::wrap(T_PARA)); 
  eta_R = (arma::rowvec) m2_func->eval(Rcpp::wrap(DESIGN), Rcpp::wrap(R_PARA)); 
  DIV = (arma::rowvec) dist_func->eval(Rcpp::wrap(eta_T), Rcpp::wrap(eta_R));  
 
  fvalTmp = arma::accu(WT % DIV);
  /*for (int i = 0; i < nSupp; i++) {
  	arma::rowvec x = DESIGN.row(i);
		eta_T = (double) m1_func->eval(Rcpp::wrap(x), Rcpp::wrap(T_PARA)); 
  	eta_R = (double) m2_func->eval(Rcpp::wrap(x), Rcpp::wrap(R_PARA)); 
  	DIV = (double) dist_func->eval(Rcpp::wrap(eta_T), Rcpp::wrap(eta_R));  
    fvalTmp += WT(i)*DIV;
  }*/
  if (std::isnan(fvalTmp)) { fvalTmp = 1e20; }
  if (!(arma::is_finite(fvalTmp))) { fvalTmp = 1e20; }
  return fvalTmp;
}

arma::rowvec f_gr(const int dPara, const arma::rowvec &R_PARA_EX, const arma::rowvec &T_PARA, const arma::mat &DESIGN, const arma::rowvec &WT,
									Rcpp::EvalBase* m1_func, Rcpp::EvalBase* m2_func, Rcpp::EvalBase* dist_func,
									const arma::rowvec &R_UPPER, const arma::rowvec &R_LOWER, const arma::irowvec &R_NBD)
{
	//central difference method
	arma::rowvec f_para(dPara), gr(dPara); arma::rowvec b_para(dPara);
  double f_val; double b_val;
  for (int i = 0; i < dPara; i++) {
  	f_para = R_PARA_EX; b_para = R_PARA_EX;
    f_para(i) += 5e-4; b_para(i) -= 5e-4;
    f_val = f_fn(dPara, f_para, T_PARA, DESIGN, WT, m1_func, m2_func, dist_func, R_UPPER, R_LOWER, R_NBD); 
    b_val = f_fn(dPara, b_para, T_PARA, DESIGN, WT, m1_func, m2_func, dist_func, R_UPPER, R_LOWER, R_NBD);
    gr(i) = f_val*1e3 - b_val*1e3;
  }
  return gr;
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
