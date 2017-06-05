

// DECLARE FUNCTIONS
double f_fn(const arma::rowvec &R_PARA_EX, const arma::rowvec &T_PARA, const arma::mat &DESIGN, const arma::rowvec &WT, model_diff_func *func_input);
arma::rowvec f_gr(const double &FVAL, const arma::rowvec &R_PARA_EX, const arma::rowvec &T_PARA, const arma::mat &DESIGN, const arma::rowvec &WT, 
									model_diff_func *func_input);
void boxCheck(arma::rowvec &R_PARA_EX, const arma::rowvec &R_UPPER, const arma::rowvec &R_LOWER, const arma::irowvec &R_NBD);

// BODY
int lbfgsKernel(const int &LOOPID, PSO_OPTIONS PSO_OPTS[], double &FVAL, arma::rowvec &R_PARA_EX, 
	 						  const arma::rowvec &T_PARA, const arma::mat &DESIGN, const arma::rowvec &WT,
								model_diff_func *func_input, 
								const arma::rowvec &R_UPPER, const arma::rowvec &R_LOWER, const arma::irowvec &R_NBD)
{
	// LBFGS SETTINGS (Temporary)
	int LBFGS_MAXIT = PSO_OPTS[LOOPID].LBFGS_MAXIT;
	// CORRECTION SETTINGS
	int LBFGS_LM = PSO_OPTS[LOOPID].LBFGS_LM;
	arma::rowvec ai(LBFGS_LM, fill::zeros), bi(LBFGS_LM, fill::zeros);
	// STOPPING CRITERION SETTING 
	double FVAL_EPS = PSO_OPTS[LOOPID].FVAL_EPS; // DEFAULT 0.0
	double GRAD_EPS = PSO_OPTS[LOOPID].GRAD_EPS; // DEFAULT 1e-5

	// LINE SEARCH PARAMETERS
	int LINESEARCH_MAXTRIAL = PSO_OPTS[LOOPID].LINESEARCH_MAXTRIAL;
	double LINESEARCH_MAX = PSO_OPTS[LOOPID].LINESEARCH_MAX;
	double LINESEARCH_MIN = PSO_OPTS[LOOPID].LINESEARCH_MIN;
	double LINESEARCH_RHO = PSO_OPTS[LOOPID].LINESEARCH_RHO;
	double LINESEARCH_ARMIJO = PSO_OPTS[LOOPID].LINESEARCH_ARMIJO;
	//double LINESEARCH_WOLFE = PSO_OPTS[LOOPID].LINESEARCH_WOLFE;
	double alpha = LINESEARCH_MAX;
	int LS_COUNTER = 0; bool FLAG = TRUE; 
	//double ZOOM_DIFF = 1.0; bool ZOOM_FLAG = TRUE;

	// Initialization
	boxCheck(R_PARA_EX, R_UPPER, R_LOWER, R_NBD);
	int dPara = R_PARA_EX.n_elem;
	// New Position
	arma::rowvec R_PARA_EX_NEW(dPara, fill::zeros);
	// Objective Function Value
	double FVAL_NEW = 1e20;
	FVAL = f_fn(R_PARA_EX, T_PARA, DESIGN, WT, func_input); 
	// Gradient
	arma::rowvec GRAD(dPara, fill::zeros), GRAD_NEW(dPara, fill::zeros);	
	GRAD = f_gr(FVAL, R_PARA_EX, T_PARA, DESIGN, WT, func_input); 
	// Hessian Matrix
	arma::mat EYE(dPara, dPara); EYE.eye();
	arma::mat H_0 = EYE;
	// Direction Correction
	arma::rowvec DIR(dPara, fill::zeros);
	arma::rowvec S_ONE(dPara, fill::zeros), Y_ONE(dPara, fill::zeros);
	arma::mat S_MAT(LBFGS_LM, dPara, fill::zeros), Y_MAT(LBFGS_LM, dPara, fill::zeros);
	double RHO = 1.0; arma::rowvec RHO_VEC(LBFGS_LM, fill::zeros);
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
			// backtracking
			alpha = LINESEARCH_MAX;
			LS_COUNTER = 0; FLAG = TRUE; 
			while ((LS_COUNTER < LINESEARCH_MAXTRIAL) & FLAG) {
				if (alpha < LINESEARCH_MIN) { alpha = LINESEARCH_MIN; FLAG = FALSE; }
				R_PARA_EX_NEW = R_PARA_EX + alpha*DIR;
				boxCheck(R_PARA_EX_NEW, R_UPPER, R_LOWER, R_NBD);
				FVAL_NEW = f_fn(R_PARA_EX_NEW, T_PARA, DESIGN, WT, func_input); 
				if (FVAL_NEW > (FVAL + alpha*LINESEARCH_ARMIJO*DIR_VAL)) {
					alpha = LINESEARCH_RHO*alpha;
				} else {
					FLAG = FALSE;
				} 
				LS_COUNTER++;				
			}
			GRAD_NEW = f_gr(FVAL_NEW, R_PARA_EX_NEW, T_PARA, DESIGN, WT, func_input); 
		} else {
			LBFGS_STOP = TRUE;
		} 

		if (LBFGS_STOP) { t += LBFGS_MAXIT; } else {
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
			
			// DEBUGGING   
			/*
			Rprintf("Iteration %d :\n", t);
			Rprintf("F = %4.8f\n", FVAL_NEW);
			if (GRAD_NEW.has_nan() | GRAD_NEW.has_inf()) { Rprintf("G has Inf or NaN\n"); } else {
				Rprintf("G = "); for (int q = 0; q < dPara; q++) { Rprintf("%4.4f", GRAD_NEW(q)); if (q < (dPara - 1)) Rprintf(", "); else Rprintf("\n"); } 	
			}
			Rprintf("alpha = %4.4f by %d iterations of LS\n", alpha, LS_COUNTER-1);
			if (DIR.has_nan() | DIR.has_inf()) { Rprintf("DIR has Inf or NaN\n"); } else {
				Rprintf("DIR = "); for (int q = 0; q < dPara; q++) { Rprintf("%4.4f", DIR(q)); if (q < (dPara - 1)) Rprintf(", "); else Rprintf("\n"); } 
			}
			if (R_PARA_EX_NEW.has_nan() | R_PARA_EX_NEW.has_inf()) { Rprintf("P has Inf or NaN\n"); } else {
				Rprintf("P = "); for (int q = 0; q < dPara; q++) { Rprintf("%4.4f", R_PARA_EX_NEW(q)); if (q < (dPara - 1)) Rprintf(", "); else Rprintf("\n"); } 
			}
			if (std::isfinite(RHO)) { Rprintf("RHO = %4.4f\n", RHO); } else { Rprintf("RHO = Inf or NaN\n"); }
			if (S_ONE.has_nan() | S_ONE.has_inf()) { Rprintf("S has Inf or NaN\n"); } else {
				Rprintf("S = "); for (int q = 0; q < dPara; q++) { Rprintf("%4.4f", S_ONE(q)); if (q < (dPara - 1)) Rprintf(", "); else Rprintf("\n"); } 
			}
			if (Y_ONE.has_nan() | Y_ONE.has_inf()) { Rprintf("Y has Inf or NaN\n"); } else {
				Rprintf("Y = "); for (int q = 0; q < dPara; q++) { Rprintf("%4.4f", Y_ONE(q)); if (q < (dPara - 1)) Rprintf(", "); else Rprintf("\n"); } 
			}
			Rprintf("\n");   
			*/
			// Examine Stopping Criterion
			arma::ivec TWOSTOP(2, fill::zeros);
			if (FVAL_EPS > 0) {
				FVAL_TSET = std::abs(FVAL - FVAL_NEW)/FVAL_NEW;
				if (FVAL_TSET < FVAL_EPS) { TWOSTOP(0) = 1; }	  
			} else {
				TWOSTOP(0) = 1;
			}
			if (GRAD_EPS > 0) {
				GRAD_TEST = arma::as_scalar(GRAD_NEW * GRAD_NEW.t());
				double XNORM = arma::as_scalar(R_PARA_EX_NEW * R_PARA_EX_NEW.t());
				if (XNORM < 1.0) { XNORM = 1.0; }
				GRAD_TEST = GRAD_TEST/XNORM;
				if (GRAD_TEST < GRAD_EPS) { TWOSTOP(1) = 1; }	
			} else {
				TWOSTOP(1) = 1;
			}
			if (all(TWOSTOP > 0)) { LBFGS_STOP = TRUE; CONV = 1; }
			// Update Approximated Hessian Matrix
			if (!LBFGS_STOP) H_0 = (EYE - RHO*(Y_ONE.t() * S_ONE))	* H_0 * (EYE - RHO*(S_ONE.t() * Y_ONE)) + RHO*(S_ONE.t() * S_ONE);	
			// Replace Current Values
			R_PARA_EX = R_PARA_EX_NEW;
			FVAL = FVAL_NEW;
			GRAD = GRAD_NEW;
		}
	}
	// END L-BFGS LOOP
	//Rprintf("\n");
	return CONV;
}

// SUBFUNCTIONS
double f_fn(const arma::rowvec &R_PARA_EX, const arma::rowvec &T_PARA, const arma::mat &DESIGN, const arma::rowvec &WT, model_diff_func *func_input)
{
	
  Rcpp::EvalBase *m1_func = (Rcpp::EvalBase *) func_input->M1_FUNC;
  Rcpp::EvalBase *m2_func = (Rcpp::EvalBase *) func_input->M2_FUNC;
  Rcpp::EvalBase *distFunc = (Rcpp::EvalBase *) func_input->DISTFUNC;

	arma::rowvec R_PARA = R_PARA_EX;
	//int nSupp = DESIGN.n_rows;
  double fvalTmp = 1e20;
  //arma::rowvec eta_T(DESIGN.n_rows, fill::zeros), eta_R(DESIGN.n_rows, fill::zeros);
  //arma::rowvec DIV(DESIGN.n_rows, fill::zeros);

  Rcpp::NumericVector eta_T_Rform, eta_R_Rform, DIV_Rform;
  
  eta_T_Rform = (Rcpp::NumericVector) m1_func->eval(Rcpp::wrap(DESIGN), Rcpp::wrap(T_PARA)); 
  eta_R_Rform = (Rcpp::NumericVector) m2_func->eval(Rcpp::wrap(DESIGN), Rcpp::wrap(R_PARA)); 
  DIV_Rform = (Rcpp::NumericVector) distFunc->eval(Rcpp::wrap(eta_T_Rform), Rcpp::wrap(eta_R_Rform));  

  arma::rowvec DIV(DIV_Rform.begin(), DIV_Rform.size(), false);

  fvalTmp = arma::accu(WT % DIV);

  if (std::isnan(fvalTmp)) { fvalTmp = 1e20; }
  if (!(arma::is_finite(fvalTmp))) { fvalTmp = 1e20; }
  
  return fvalTmp;
}

arma::rowvec f_gr(const double &FVAL, const arma::rowvec &R_PARA_EX, const arma::rowvec &T_PARA, const arma::mat &DESIGN, const arma::rowvec &WT, 
									model_diff_func *func_input)
{
	//forward difference method
	arma::rowvec gr(R_PARA_EX.n_elem, fill::zeros);
	arma::rowvec fr_para(R_PARA_EX.n_elem, fill::zeros); 
	arma::rowvec bk_para(R_PARA_EX.n_elem, fill::zeros);
  double fr_val = 0.0; 
  double bk_val = 0.0;
  for (uword i = 0; i < R_PARA_EX.n_elem; i++) {
  	fr_para = R_PARA_EX; 
  	bk_para = R_PARA_EX;
    fr_para(i) += 1e-3; bk_para(i) -= 1e-3;
    //fr_para(i) += 1e-4;
    fr_val = f_fn(fr_para, T_PARA, DESIGN, WT, func_input); 
    bk_val = f_fn(bk_para, T_PARA, DESIGN, WT, func_input); 
    //gr(i) = fr_val*1e4 - FVAL*1e4;
    gr(i) = fr_val*2e3 - bk_val*2e3;
  }
  return gr;
}

//
void boxCheck(arma::rowvec &R_PARA_EX, const arma::rowvec &R_UPPER, const arma::rowvec &R_LOWER, const arma::irowvec &R_NBD)
{
	int len = R_PARA_EX.n_elem;
	int nbd = -1;
	for (int i = 0; i < len; i++) {
		nbd = R_NBD(i);
		switch(nbd) {
      case 1: { if (R_PARA_EX(i) < R_LOWER(i)) R_PARA_EX(i) = R_LOWER(i); break; } // Lower
      case 2: { 
      	if (R_PARA_EX(i) < R_LOWER(i)) { R_PARA_EX(i) = R_LOWER(i); }
      	if (R_PARA_EX(i) > R_UPPER(i)) { R_PARA_EX(i) = R_UPPER(i); }
      	break; 
      } // Both
      case 3: { if (R_PARA_EX(i) > R_UPPER(i)) R_PARA_EX(i) = R_UPPER(i); break; } // Upper
    }
	}
}


void line_search()
{
	// strong Wolfe conditions.
	double LS_FVAL_0, LS_FVAL; LS_FVAL_0 = FVAL;
	arma::rowvec LS_GRAD_0, LS_GRAD; LS_GRAD_0 = GRAD;
	arma::rowvec LS_X;
	double alpha_0 = LINESEARCH_MIN; 
			alpha = alpha_0 + (LINESEARCH_MAX - alpha_0)*as_scalar(arma::randu(1)); //
			// alpha is alpha_i
			// alpha_0 is alpha_{i-1}
			LS_COUNTER = 0; FLAG = TRUE; 
			while ((LS_COUNTER < LINESEARCH_MAXTRIAL) & FLAG) {

				LS_X = R_PARA_EX + alpha*DIR;
				LS_FVAL = f_fn(LS_X, T_PARA, DESIGN, WT, func_input); 
				LS_GRAD = f_gr(LS_FVAL, LS_X, T_PARA, DESIGN, WT, func_input); 
				
				if ((LS_FVAL > (FVAL + alpha*LINESEARCH_ARMIJO*DIR_VAL)) | ((LS_COUNTER > 0) & (LS_FVAL > LS_FVAL_0))) {
					double aL = alpha_0; double aH = alpha;
					ZOOM_DIFF = 1.0; ZOOM_FLAG = TRUE;
					while ((ZOOM_DIFF > 1e-5) & ZOOM_FLAG) {
						alpha = 0.5*(aL + aH);
						LS_X = R_PARA_EX + alpha*DIR;
						LS_FVAL = f_fn(LS_X, T_PARA, DESIGN, WT, func_input); 
						LS_GRAD = f_gr(LS_FVAL, LS_X, T_PARA, DESIGN, WT, func_input); 
						rowvec xL = R_PARA_EX + aL*DIR;
						double fL = f_fn(xL, T_PARA, DESIGN, WT, func_input); 
						if ((LS_FVAL > (FVAL + alpha*LINESEARCH_ARMIJO*DIR_VAL)) | (LS_FVAL > fL)) {
							aH = alpha;
						} else {
							if (std::abs(arma::as_scalar(LS_GRAD * DIR.t())) <= ((-1.0)*LINESEARCH_WOLFE*DIR_VAL)) { ZOOM_FLAG = FALSE; } 
							if (ZOOM_FLAG & ((arma::as_scalar(LS_GRAD * DIR.t())*(aH - aL)) >= 0)) { aH = aL; }
							if (ZOOM_FLAG) { aL = alpha; }
						}
						ZOOM_DIFF = std::abs(aH - aL);
					}
					FLAG = FALSE;
				} 
				if (FLAG & (std::abs(arma::as_scalar(LS_GRAD * DIR.t())) <= ((-1.0)*LINESEARCH_WOLFE*DIR_VAL))) { FLAG = FALSE; }
				if (FLAG & (arma::as_scalar(LS_GRAD * DIR.t()) >= 0)) {
					double aL = alpha; double aH = alpha_0;
					ZOOM_DIFF = 1.0; ZOOM_FLAG = TRUE;
					while ((ZOOM_DIFF > 1e-5) & ZOOM_FLAG) {
						alpha = 0.5*(aL + aH);
						LS_X = R_PARA_EX + alpha*DIR;
						LS_FVAL = f_fn(dPara, LS_X, T_PARA, DESIGN, WT, m1_func, m2_func, dist_func, R_UPPER, R_LOWER, R_NBD); 
						LS_GRAD = f_gr(dPara, LS_X, T_PARA, DESIGN, WT, m1_func, m2_func, dist_func, R_UPPER, R_LOWER, R_NBD); 
						rowvec xL = R_PARA_EX + aL*DIR;
						double fL = f_fn(dPara, xL, T_PARA, DESIGN, WT, m1_func, m2_func, dist_func, R_UPPER, R_LOWER, R_NBD); 
						if ((LS_FVAL > (FVAL + alpha*LINESEARCH_ARMIJO*DIR_VAL)) | (LS_FVAL > fL)) {
							aH = alpha;
						} else {
							if (std::abs(arma::as_scalar(LS_GRAD * DIR.t())) <= ((-1.0)*LINESEARCH_WOLFE*DIR_VAL)) { ZOOM_FLAG = FALSE; } 
							if (ZOOM_FLAG & ((arma::as_scalar(LS_GRAD * DIR.t())*(aH - aL)) >= 0)) { aH = aL; }
							if (ZOOM_FLAG) { aL = alpha; }
						}
						ZOOM_DIFF = std::abs(aH - aL);
					}
					FLAG = FALSE;
				}
				if (FLAG) {
					alpha_0 = alpha;
					LS_FVAL_0 = LS_FVAL;
					LS_GRAD_0 = LS_GRAD;
					alpha = alpha_0 + (LINESEARCH_MAX - alpha_0)*as_scalar(arma::randu(1));	
				}
				LS_COUNTER++;				
			}
}
