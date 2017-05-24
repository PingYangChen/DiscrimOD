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
				LS_FVAL = f_fn(dPara, LS_X, T_PARA, DESIGN, WT, m1_func, m2_func, dist_func, R_UPPER, R_LOWER, R_NBD); 
				LS_GRAD = f_gr(dPara, LS_X, T_PARA, DESIGN, WT, m1_func, m2_func, dist_func, R_UPPER, R_LOWER, R_NBD); 
				
				if ((LS_FVAL > (FVAL + alpha*LINESEARCH_ARMIJO*DIR_VAL)) | ((LS_COUNTER > 0) & (LS_FVAL > LS_FVAL_0))) {
					double aL = alpha_0; double aH = alpha;
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