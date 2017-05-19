
#include "psoHeader.h"

// RCPP FUNCTIONS
//[[Rcpp::export]]
Rcpp::List cppPSO(const int LOOPID, Rcpp::List ALG_INFO_LIST, Rcpp::List OBJ_INFO_LIST, Rcpp::List MODEL_INFO_LIST, 
                  arma::rowvec FIXEDVALUE, SEXP env, const bool IF_PARALLEL, const bool VERBOSE)
{
  //arma_rng::set_seed_random();
  /*int NCPU = omp_get_max_threads();
	if (NCPU < 3) NCPU = 2;
  omp_set_num_threads(NCPU - 1);*/

	OBJ_INFO OBJ = {}; getInfoStruct(OBJ, OBJ_INFO_LIST, env);

  int N_model = OBJ.N_model;
  MODEL_SET MODELS[N_MODEL_MAX]; getModelFunc(MODELS, MODEL_INFO_LIST, env, N_model);
	
  PSO_OPTIONS PSO_OPT[N_PSO_OPTS]; getAlgStruct(PSO_OPT, ALG_INFO_LIST);

  PSO_Result Result = {};
  
  if (VERBOSE) Rprintf("\n Calling Cpp PSO Kernel... ");
  PSO_MAIN(LOOPID, PSO_OPT, OBJ, MODELS, FIXEDVALUE, IF_PARALLEL, VERBOSE, &Result);
  if (VERBOSE) Rprintf("Done.\n");

  return List::create(Named("GBest") = wrap(Result.GBest),
                      Named("fGBest") = wrap(Result.fGBest),
                      Named("fGBestHist") = wrap(Result.fGBestHist),
                      Named("PBest") = wrap(Result.PBest),
                      Named("fPBest") = wrap(Result.fPBest));
}


//[[Rcpp::export]]
Rcpp::List cppDesignCriterion(Rcpp::List ALG_INFO_LIST, Rcpp::List OBJ_INFO_LIST, Rcpp::List MODEL_INFO_LIST,
                              arma::rowvec FIXEDVALUE, SEXP env, arma::rowvec DESIGN)
{
  OBJ_INFO OBJ = {}; getInfoStruct(OBJ, OBJ_INFO_LIST, env);

  int N_model = OBJ.N_model;
  MODEL_SET MODELS[N_model]; getModelFunc(MODELS, MODEL_INFO_LIST, env, N_model);

  PSO_OPTIONS PSO_OPT[N_PSO_OPTS]; getAlgStruct(PSO_OPT, ALG_INFO_LIST);

  arma::rowvec R_PARA;
  double val = DesignCriterion(0, PSO_OPT, OBJ, MODELS, FIXEDVALUE, DESIGN, R_PARA);

  //return val;
  return List::create(Named("val") = wrap(val),
                      Named("theta2") = wrap(R_PARA));
}
