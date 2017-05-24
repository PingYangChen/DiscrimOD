
#include "psoHeader.h"

// RCPP FUNCTIONS
//[[Rcpp::export]]
Rcpp::List cppPSO(const int LOOPID, Rcpp::List ALG_INFO_LIST, Rcpp::List OBJ_INFO_LIST, Rcpp::List MODEL_INFO_LIST, 
                  arma::rowvec FIXEDVALUE, const SEXP env, const bool IF_PARALLEL, const bool VERBOSE)
{
  //arma_rng::set_seed_random();
  /*int NCPU = omp_get_max_threads();
	if (NCPU < 3) NCPU = 2;
  omp_set_num_threads(NCPU - 1);*/

  Rcpp::EvalBase *distFunc = NULL;
  SEXP tmp = as<SEXP>(OBJ_INFO_LIST["dist_func"]);
  if (TYPEOF(tmp) == EXTPTRSXP) {   
    distFunc = new Rcpp::EvalCompiled(tmp, env); 
  } else {
    distFunc = new Rcpp::EvalStandard(tmp, env); 
  } 

	OBJ_INFO OBJ = {}; getInfoStruct(OBJ, OBJ_INFO_LIST);

  int N_model = OBJ.N_model; MODEL_SET MODELS[N_model]; 
	for (int i = 0; i < N_model; i++) {
    Rcpp::EvalBase *modelFunc = NULL;
    SEXP tmp = as<SEXP>(MODEL_INFO_LIST[i]);
    if (TYPEOF(tmp) == EXTPTRSXP) {   
      modelFunc = new Rcpp::EvalCompiled(tmp, env);
    } else {                                                
      modelFunc = new Rcpp::EvalStandard(tmp, env);
    }  
    MODELS[i].modelFunc = modelFunc;
  }

  PSO_OPTIONS PSO_OPT[N_PSO_OPTS]; getAlgStruct(PSO_OPT, ALG_INFO_LIST);

  PSO_Result Result = {};
  
  if (VERBOSE) Rprintf("\n Calling Cpp PSO Kernel... ");
  PSO_MAIN(LOOPID, PSO_OPT, OBJ, MODELS, distFunc, FIXEDVALUE, IF_PARALLEL, VERBOSE, &Result);
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
  Rcpp::EvalBase *distFunc = NULL;
  SEXP tmp = as<SEXP>(OBJ_INFO_LIST["dist_func"]);
  if (TYPEOF(tmp) == EXTPTRSXP) {   
    distFunc = new Rcpp::EvalCompiled(tmp, env);
  } else {
    distFunc = new Rcpp::EvalStandard(tmp, env); 
  } 

  OBJ_INFO OBJ = {}; getInfoStruct(OBJ, OBJ_INFO_LIST);

  int N_model = OBJ.N_model; MODEL_SET MODELS[N_model]; 
  for (int i = 0; i < N_model; i++) {
    Rcpp::EvalBase *modelFunc = NULL;
    SEXP tmp = as<SEXP>(MODEL_INFO_LIST[i]);
    if (TYPEOF(tmp) == EXTPTRSXP) {   
      modelFunc = new Rcpp::EvalCompiled(tmp, env);
    } else {                                                
      modelFunc = new Rcpp::EvalStandard(tmp, env);
    }  
    MODELS[i].modelFunc = modelFunc;
  }

  PSO_OPTIONS PSO_OPT[N_PSO_OPTS]; getAlgStruct(PSO_OPT, ALG_INFO_LIST);

  arma::mat R_PARA;
  double val = DesignCriterion(0, PSO_OPT, OBJ, MODELS, distFunc, FIXEDVALUE, DESIGN, R_PARA);

  return List::create(Named("val") = wrap(val),
                      Named("theta2") = wrap(R_PARA));
}

//[[Rcpp::export]]
List cppEquivalence(Rcpp::List ALG_INFO_LIST, Rcpp::List OBJ_INFO_LIST, Rcpp::List MODEL_INFO_LIST,
                    double GBEST_VAL, arma::mat R_PARA_SET, SEXP env, const int nGrid)
{
  Rcpp::EvalBase *distFunc = NULL;
  SEXP tmp = as<SEXP>(OBJ_INFO_LIST["dist_func"]);
  if (TYPEOF(tmp) == EXTPTRSXP) {   
    distFunc = new Rcpp::EvalCompiled(tmp, env);
  } else {
    distFunc = new Rcpp::EvalStandard(tmp, env); 
  } 

  OBJ_INFO OBJ = {}; getInfoStruct(OBJ, OBJ_INFO_LIST);

  int N_model = OBJ.N_model; MODEL_SET MODELS[N_model]; 
  for (int i = 0; i < N_model; i++) {
    Rcpp::EvalBase *modelFunc = NULL;
    SEXP tmp = as<SEXP>(MODEL_INFO_LIST[i]);
    if (TYPEOF(tmp) == EXTPTRSXP) {   
      modelFunc = new Rcpp::EvalCompiled(tmp, env);
    } else {                                                
      modelFunc = new Rcpp::EvalStandard(tmp, env);
    }  
    MODELS[i].modelFunc = modelFunc;
  }

  PSO_OPTIONS PSO_OPT[N_PSO_OPTS]; getAlgStruct(PSO_OPT, ALG_INFO_LIST);
  
  /* REVISE HERE */
  Rcpp::EvalBase *m1_func = MODELS[0].modelFunc;
  Rcpp::EvalBase *m2_func = MODELS[1].modelFunc;
  
  arma::rowvec T_PARA = OBJ.paras.submat(0, 0, 0, OBJ.dParas(0) - 1);
  arma::mat DISPVALS; 
  arma::vec xLine_1 = linspace<vec>(OBJ.dsLower(0), OBJ.dsUpper(0), nGrid); 
  arma::vec xLine_2(1, fill::zeros);
  int dSupp = OBJ.dSupp;
  if (dSupp == 1) {
    DISPVALS.set_size(1, nGrid);
    arma::rowvec eta_T, eta_R, DIV;
    eta_T = (arma::rowvec) m1_func->eval(Rcpp::wrap(xLine_1), Rcpp::wrap(T_PARA)); 
    eta_R = (arma::rowvec) m2_func->eval(Rcpp::wrap(xLine_1), Rcpp::wrap(R_PARA_SET.submat(1, 0, 1, OBJ.dParas(1) - 1))); 
    DIV = (arma::rowvec) distFunc->eval(Rcpp::wrap(eta_T), Rcpp::wrap(eta_R)); 
    DISPVALS.row(0) = DIV - GBEST_VAL;
  } else if (dSupp == 2) {
    xLine_2.reset(); xLine_2 = linspace<vec>(OBJ.dsLower(1), OBJ.dsUpper(1), nGrid);
    DISPVALS.set_size(nGrid, nGrid);
    arma::rowvec eta_T, eta_R, DIV;
    for (int i = 0; i < nGrid; i++) {
      for (int j = 0; j < nGrid; j++) {
        arma::rowvec pt(2); pt << xLine_1(i) << xLine_2(j) << endr;
        eta_T = (arma::rowvec) m1_func->eval(Rcpp::wrap(pt), Rcpp::wrap(T_PARA)); 
        eta_R = (arma::rowvec) m2_func->eval(Rcpp::wrap(pt), Rcpp::wrap(R_PARA_SET.submat(1, 0, 1, OBJ.dParas(1) - 1))); 
        DIV = (arma::rowvec) distFunc->eval(Rcpp::wrap(eta_T), Rcpp::wrap(eta_R)); 
        DISPVALS(i, j) = DIV(0) - GBEST_VAL;
      }
    }
  } 
  double MAX_DISP = DISPVALS.max();
  /* END REVISION */
  
  /* OUTPUT */  
  return List::create(Named("Grid_1") = wrap(xLine_1),
                      Named("Grid_2") = wrap(xLine_2),
                      Named("DirDeriv") = wrap(DISPVALS),
                      Named("MAX_DD") = wrap(MAX_DISP));                       
}