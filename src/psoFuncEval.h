// DECLARE FUNCTIONS
#include "DesignCriterion.h"

// BODY
void psoFuncEval(const bool &IF_PARALLEL, const int &LOOPID, const PSO_OPTIONS PSO_OPTS[], const OBJ_INFO &OBJ, const PSO_DYN &PSO_DYN, 
								 const MODEL_SET MODELS[], Rcpp::EvalBase* distFunc, const rowvec &FIXEDVALUE, const mat &swarm, vec &fSwarm)
{	
  // SET PSO PARAMETERS
  int nSwarm = (int)swarm.n_rows;
  /*if (IF_PARALLEL) { 
		// PARALLEL LOOP (DOES NOT WORK FOR Rcpp::EvalBase*. MAY REPLACE OpenMP BY RcppParallel IN THE FUTURE.)
		int iParallel;
		#pragma omp parallel private(iParallel) 
		{
		#pragma omp for
			for (iParallel = 0; iParallel < nSwarm; iParallel++) {
				rowvec PARTICLE = conv_to<rowvec>::from(swarm.row(iParallel));
				// Optimal Design Criteria
				fSwarm(iParallel) = DesignCriterion(LOOPID, PSO_OPTS, OBJ, meanFunc_core, corrFunc_core, FIXEDVALUE, PARTICLE);
			}
		}
  } else {*/
		// NON-PARALLEL LOOP
		for (int iSwarm = 0; iSwarm < nSwarm; iSwarm++) {
			rowvec PARTICLE = conv_to<rowvec>::from(swarm.row(iSwarm));
			// Optimal Design Criteria
			arma::mat R_PARA;
			fSwarm(iSwarm) = DesignCriterion(LOOPID, PSO_OPTS, OBJ, MODELS, distFunc, FIXEDVALUE, PARTICLE, R_PARA);
		}
  //}
}
