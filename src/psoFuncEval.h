// DECLARE FUNCTIONS
#include "DesignCriterion.h"

// BODY
void psoFuncEval(const bool &IF_PARALLEL, const int LOOPID, const PSO_OPTIONS PSO_OPTS[], const OBJ_INFO &OBJ, const PSO_DYN &PSO_DYN, 
								 model_diff_func *MODEL_COLLECTOR[], const rowvec &FIXEDVALUE, const mat &swarm, vec &fSwarm)
{	
  int nSwarm = (int)swarm.n_rows;
 /* if (IF_PARALLEL) { 
		// PARALLEL LOOP (DOES NOT WORK FOR Rcpp::EvalBase*. MAY REPLACE OpenMP BY RcppParallel IN THE FUTURE.)
		int iParallel;
		#pragma omp parallel private(iParallel) 
		{
		#pragma omp for
			for (iParallel = 0; iParallel < nSwarm; iParallel++) {
				rowvec PARTICLE = conv_to<rowvec>::from(swarm.row(iParallel));
				// Optimal Design Criteria
				arma::mat R_PARA;
				fSwarm(iParallel) = DesignCriterion(LOOPID, PSO_OPTS, OBJ, MODEL_INFO_LIST, DIST_FUNC_SEXP, env, FIXEDVALUE, PARTICLE, R_PARA);
			}
			#pragma omp barrier
		}
  } else {*/
		// NON-PARALLEL LOOP
		for (int iSwarm = 0; iSwarm < nSwarm; iSwarm++) {
			rowvec PARTICLE = conv_to<rowvec>::from(swarm.row(iSwarm));
			// Optimal Design Criteria
			arma::mat R_PARA;
			fSwarm(iSwarm) = DesignCriterion(LOOPID, PSO_OPTS, OBJ, MODEL_COLLECTOR, FIXEDVALUE, PARTICLE, R_PARA);
		}
  //}
}
