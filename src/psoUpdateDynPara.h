void psoUpdateDynPara(const int LOOPID, PSO_OPTIONS PSO_OPTS[], const int iter, PSO_DYN &PSO_DYN,
											const arma::mat swarm, const arma::mat PBest, const arma::rowvec GBest,
											const arma::vec fSwarm, const arma::vec fPBest, const double fGBest)
{
  if (iter < 0) { // INITIALIZE

  	//PSO_DYN.succ_GB	= 0;
  	int w_varyfor = (int)(PSO_OPTS[LOOPID].w_var*PSO_OPTS[LOOPID].maxIter);
  	PSO_DYN.w_varyfor	= w_varyfor;
	  PSO_DYN.w_cur			= PSO_OPTS[LOOPID].w0;
	  PSO_DYN.w_dec			= (PSO_OPTS[LOOPID].w0 - PSO_OPTS[LOOPID].w1)/w_varyfor;      // Inertia weight change per iteration step

  } else { // UPDATE

		if (iter <= PSO_DYN.w_varyfor) 		PSO_DYN.w_cur = PSO_DYN.w_cur - PSO_DYN.w_dec;
  }
}

