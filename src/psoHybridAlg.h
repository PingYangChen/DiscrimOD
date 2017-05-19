// DECLARE FUNCTIONS


// BODY
// hybridAlg(IF_PARALLEL, LOOPID, PSO_OPTS, OBJ, FIXEDVALUE, swarm, PBest, GBest, fSwarm, fPBest, fGBest);
void psoHybridAlg(const bool &IF_PARALLEL, const int &LOOPID, const PSO_OPTIONS PSO_OPTS[], const OBJ_INFO &OBJ, const PSO_DYN &PSO_DYN, const rowvec &FIXEDVALUE, 
					 		 		arma::mat &swarm, arma::mat &PBest, arma::rowvec &GBest, arma::vec &fSwarm, arma::vec &fPBest, double &fGBest)
{	
  // SET PSO PARAMETERS
  PSO_OPTS[LOOPID].hybrid_type;
  
  int nSwarm = (int)swarm.n_rows;
  if (IF_PARALLEL) {
		// PARALLEL LOOP
		int iParallel;
		#pragma omp parallel private(iParallel) 
		{
		#pragma omp for
			for (iParallel = 0; iParallel < nSwarm; iParallel++) {
				rowvec PARTICLE = conv_to<rowvec>::from(swarm.row(iParallel));
				// Optimal Design Criteria
				fSwarm(iParallel) = DesignCriterion(LOOPID, PSO_OPTS, OBJ, FIXEDVALUE, PARTICLE);
				//Benchmark Objective Functions
				//fSwarm(iParallel) = ObjectiveFunction(OBJ, FIXEDVALUE, PARTICLE);
			}
		}
  } else {
		// NON-PARALLEL LOOP
		for (int iSwarm = 0; iSwarm < nSwarm; iSwarm++) {
			rowvec PARTICLE = conv_to<rowvec>::from(swarm.row(iSwarm));
			// Optimal Design Criteria
			fSwarm(iSwarm) = DesignCriterion(LOOPID, PSO_OPTS, OBJ, FIXEDVALUE, PARTICLE);
			//Benchmark Objective Functions
			//fSwarm(iSwarm) = ObjectiveFunction(OBJ, FIXEDVALUE, PARTICLE);
		}
  }
}
