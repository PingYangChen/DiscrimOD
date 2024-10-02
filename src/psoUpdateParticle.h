
void psoUpdateParticle(const int LOOPID, PSO_OPTIONS PSO_OPTS[], const PSO_DYN PSO_DYN,
							 				 const arma::mat PBest, const arma::rowvec GBest,
							 				 const arma::rowvec velMax, const arma::rowvec varUpper, const arma::rowvec varLower,
							 				 arma::mat &vStep, arma::mat &swarm)
{
  int dSwarm = (int)swarm.n_cols;
  int nSwarm = (int)swarm.n_rows;
	//int typePSO = PSO_OPTS[LOOPID].typePSO;
	int typePSO = 0;
	double c1 = PSO_OPTS[LOOPID].c1;
  double c2 = PSO_OPTS[LOOPID].c2;
  //double chi = PSO_OPTS[LOOPID].chi;

	arma::mat velMax_Mat = repmat(velMax, nSwarm, 1);
	arma::mat GBmat = repmat(GBest, nSwarm, 1);
	GetRNGstate();
	switch (typePSO) {
		case 0: // Linearly Decreasing Weight PSO (Shi, Y. H.	and Eberhart, R. C., 1998)
		{	// The most common one
			vStep = PSO_DYN.w_cur*vStep + c1*arma::randu(nSwarm, dSwarm) % (PBest - swarm) +
																		c2*arma::randu(nSwarm, dSwarm) % (GBmat - swarm);
			vStep = min(vStep, velMax_Mat); vStep = max(vStep, (-1)*velMax_Mat);
			swarm += vStep;
			break;
		}
	}
	PutRNGstate();
}

