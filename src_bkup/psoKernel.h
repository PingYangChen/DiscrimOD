
// BODY
// PSO MAIN FUNCTIONS
void PSO_MAIN(const int &LOOPID, const PSO_OPTIONS PSO_OPTS[], const OBJ_INFO &OBJ, Rcpp::List MODEL_INFO_LIST, SEXP DIST_FUNC_SEXP, SEXP env,
              const rowvec &FIXEDVALUE, const bool &IF_PARALLEL, const bool COUNTER_ON, Ptr_PSO_Result Ptr_PSO_Result)
{
	/* -- BEGIN -- */
  // GET PSO PARAMETERS
	int nSwarm    = PSO_OPTS[LOOPID].nSwarm; 
	int dSwarm    = PSO_OPTS[LOOPID].dSwarm; 
  int nGroup    = PSO_OPTS[LOOPID].nGroup;
	int maxIter   = PSO_OPTS[LOOPID].maxIter; 
	int checkConv = PSO_OPTS[LOOPID].checkConv; 
	double freeRun   = PSO_OPTS[LOOPID].freeRun; 
	double tol       = PSO_OPTS[LOOPID].tol; 
  rowvec varUpper  = PSO_OPTS[LOOPID].varUpper;
  rowvec varLower  = PSO_OPTS[LOOPID].varLower;

	// DECLARE VARIABLES
  arma::mat swarm(nSwarm, dSwarm), vStep(nSwarm, dSwarm), PBest(nSwarm, dSwarm), GrBest(nGroup, dSwarm);
  arma::rowvec velMax(dSwarm);
  arma::rowvec GBest(dSwarm);
  arma::vec fSwarm(nSwarm), fPBest(nSwarm), fGrBest(nGroup);
  double fGBest;
  arma::uword GBestIdx;
  arma::rowvec fGBestHist(maxIter + 1, fill::zeros);
	PSO_DYN PSO_DYN;

  /* -- START INITIALIZATION -- */
  if (COUNTER_ON) { Rprintf("PSO Loop: Initializing .. "); }
  // GENERATE THE VMAX MATRIX
  double vk = PSO_OPTS[LOOPID].vk;
  velMax = (varUpper - varLower)/vk;
  // INITIALIZE RANDOM SWARM
  swarm = randu(nSwarm, dSwarm) % repmat(varUpper - varLower, nSwarm, 1) + repmat(varLower, nSwarm, 1);
  // INITIALIZE VELOCITY
  vStep.fill(0);
  // INITIALIZE OBJECTIVE FUNCTION VALUES
  psoFuncEval(IF_PARALLEL, LOOPID, PSO_OPTS, OBJ, PSO_DYN, MODEL_INFO_LIST, DIST_FUNC_SEXP, env, FIXEDVALUE, swarm, fSwarm); 
  // INITIALIZE LOCAL BEST
  fPBest = fSwarm;	PBest = swarm;
  // INITIALIZE GLOBAL BEST
  fGBest = fPBest.min(GBestIdx); 
	GBest = PBest.row(GBestIdx);	
  // INITIALIZE PSO DYNAMIC PARAMETERS
  psoUpdateDynPara(LOOPID, PSO_OPTS, -1, PSO_DYN, swarm, PBest, GBest, fSwarm, fPBest, fGBest);
	// SAVE INITIAL GLOBAL BEST VALUE
	fGBestHist(0) = fGBest;
	  // SET ITERATION COUNTER
  int t; 
  if (COUNTER_ON) Rprintf("OK \n"); 
  /* -- FINISH INITIALIZATION -- */

  /* -- START PSO LOOP -- */
  for (t = 0; t < maxIter; t++) {
  	// PRINT OUT PROGRESS
    if (COUNTER_ON) {
      if (t == 0)  Rprintf("PSO Loop: Updating ..    "); 
      Rprintf("\b\b\b%2.0f%%", (double)((t+1)*100/maxIter)); 
      if (t == (maxIter - 1)) Rprintf("\n"); 
    }
    // UPDATE VELOCITY
		psoUpdateParticle(LOOPID, PSO_OPTS, PSO_DYN, PBest, GBest, velMax, varUpper, varLower, vStep, swarm);
    // UPDATE SWARM POSITION
    psoCheckParticle(LOOPID, PSO_OPTS, PSO_DYN, varUpper, varLower, swarm);	
    // UPDATE OBJECTIVE FUNCTION VALUES
    psoFuncEval(IF_PARALLEL, LOOPID, PSO_OPTS, OBJ, PSO_DYN, MODEL_INFO_LIST, DIST_FUNC_SEXP, env, FIXEDVALUE, swarm, fSwarm); 
    // UPDATE THE LOCAL AND GLOBAL BEST
    if (any(fSwarm < fPBest)) {
      uvec RowChange = find(fSwarm < fPBest);
      fPBest.elem(RowChange) = fSwarm.elem(RowChange);
      PBest.rows(RowChange) = swarm.rows(RowChange);
    }
    if (min(fPBest) < fGBest) {
      fGBest = fPBest.min(GBestIdx); GBest = PBest.row(GBestIdx); PSO_DYN.succ_GB = 1;
    }		
    // UPDATE PSO DYNAMIC PARAMETERS
    psoUpdateDynPara(LOOPID, PSO_OPTS, t, PSO_DYN, swarm, PBest, GBest, fSwarm, fPBest, fGBest);
    // SAVE CURRENT GLOBAL BEST VALUE
    fGBestHist(t+1) = fGBest; 
    // CHECK STOPPING CRITERION
    if ((t > (int)(freeRun*maxIter)) && (checkConv == 1)) { 
      if (std::abs(fGBest - fGBestHist(t)) < tol) { 
        fGBestHist.subvec(t+1, maxIter).fill(fGBest); t = maxIter; 
        if (COUNTER_ON) Rprintf(" The updating procedure converges. \n");
      }
    }
  }
  /* -- FINISH PSO LOOP -- */
 
 	/* -- OUTPUT -- */
  Ptr_PSO_Result->GBest = GBest;
  Ptr_PSO_Result->fGBest = fGBest;
  Ptr_PSO_Result->fGBestHist = fGBestHist;
  Ptr_PSO_Result->PBest = PBest;
  Ptr_PSO_Result->fPBest = fPBest;
  /* -- END -- */
}