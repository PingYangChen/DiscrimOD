// Rcpp Header File
#include <memory>
#include <cmath>
#include <RcppArmadillo.h>
#include <R.h>
//#include <omp.h>

//using namespace Rcpp;
using namespace arma;

#ifndef Rcpp_PSO_EVALUATE_H_
#define Rcpp_PSO_EVALUATE_H_

namespace Rcpp {

  class EvalBase2 {
    public:
        EvalBase2() : neval(0) {};
        virtual Rcpp::NumericVector eval(SEXP x, SEXP p) = 0;
        virtual ~EvalBase2() {};
        //unsigned long getNbEvals() { return neval; }
    protected:
        unsigned long int neval;
  };

  class EvalStandard2 : public EvalBase2 {
    public:
        EvalStandard2(SEXP fcall_, SEXP env_) : fcall(fcall_), env(env_) {}
        Rcpp::NumericVector eval(SEXP x, SEXP p) {
          //neval++;
          return defaultfun2(x, p);
        }
        ~EvalStandard2() {};
    private:
        SEXP fcall, env;
        Rcpp::NumericVector defaultfun2(SEXP x, SEXP p) {
          //Shield<SEXP> fn(Rcpp::Rcpp_lang4(fcall, x, p, R_DotsSymbol)); // disabled by Rcpp Version 1.0.4
          Shield<SEXP> fn(::Rf_lang4(fcall, x, p, R_DotsSymbol));
          Shield<SEXP> sexp_fvec(::Rf_eval(fn, env));
          //SEXP sexp_fvec = Rcpp::Rcpp_eval(fn, env); // too slow
          Rcpp::NumericVector f_result = (Rcpp::NumericVector) Rcpp::as<Rcpp::NumericVector>(sexp_fvec);
          return f_result;
        }
    };

  typedef Rcpp::NumericVector (*funcPtr2)(SEXP, SEXP, SEXP);

  class EvalCompiled2 : public EvalBase2 {
    public:
        EvalCompiled2(Rcpp::XPtr<funcPtr2> xptr, SEXP __env) {
          funptr2 = *(xptr);
          env = __env;
        };
        EvalCompiled2(SEXP xps, SEXP __env) {
          Rcpp::XPtr<funcPtr2> xptr(xps);
          funptr2 = *(xptr);
          env = __env;
        };
        Rcpp::NumericVector eval(SEXP x, SEXP p) {
          //neval++;
          Rcpp::NumericVector f_result = funptr2(x, p, env);
          return f_result;
        }
        ~EvalCompiled2() {};
    private:
        funcPtr2 funptr2;
        SEXP env;
    };

  class EvalBase4 {
  public:
    EvalBase4() : neval(0) {};
    virtual Rcpp::NumericVector eval(SEXP xt, SEXP xr, SEXP vt, SEXP vr) = 0;
    virtual ~EvalBase4() {};
    //unsigned long getNbEvals() { return neval; }
  protected:
    unsigned long int neval;
  };

  class EvalStandard4 : public EvalBase4 {
  public:
    EvalStandard4(SEXP fcall_, SEXP env_) : fcall(fcall_), env(env_) {}
    Rcpp::NumericVector eval(SEXP xt, SEXP xr, SEXP vt, SEXP vr) {
      //neval++;
      return defaultfun4(xt, xr, vt, vr);
    }
    ~EvalStandard4() {};
  private:
    SEXP fcall, env;
    Rcpp::NumericVector defaultfun4(SEXP xt, SEXP xr, SEXP vt, SEXP vr) {
      //Shield<SEXP> fn(Rcpp::Rcpp_lang4(fcall, x, p, R_DotsSymbol)); // disabled by Rcpp Version 1.0.4
      Shield<SEXP> fn(::Rf_lang6(fcall, xt, xr, vt, vr, R_DotsSymbol));
      Shield<SEXP> sexp_fvec(::Rf_eval(fn, env));
      //SEXP sexp_fvec = Rcpp::Rcpp_eval(fn, env); // too slow
      Rcpp::NumericVector f_result = (Rcpp::NumericVector) Rcpp::as<Rcpp::NumericVector>(sexp_fvec);
      return f_result;
    }
  };

  typedef Rcpp::NumericVector (*funcPtr4)(SEXP, SEXP, SEXP, SEXP, SEXP);

  class EvalCompiled4 : public EvalBase4 {
  public:
    EvalCompiled4(Rcpp::XPtr<funcPtr4> xptr, SEXP __env) {
      funptr4 = *(xptr);
      env = __env;
    };
    EvalCompiled4(SEXP xps, SEXP __env) {
      Rcpp::XPtr<funcPtr4> xptr(xps);
      funptr4 = *(xptr);
      env = __env;
    };
    Rcpp::NumericVector eval(SEXP xt, SEXP xr, SEXP vt, SEXP vr) {
      //neval++;
      Rcpp::NumericVector f_result = funptr4(xt, xr, vt, vr, env);
      return f_result;
    }
    ~EvalCompiled4() {};
  private:
    funcPtr4 funptr4;
    SEXP env;
  };
}

#endif

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]

// DEFINE STUCTURE OF OBJECTIVE FUNCTION INFORMATION
typedef struct {
  // Optimal Design Problems
  int crit_type; // 0
  int d_type; // 0
  //int nSubj; // 1
  int dSupp; // 1
  int nSupp; // 2
  // Design Space
  arma::rowvec dsLower; // rep(0,1)
  arma::rowvec dsUpper; // rep(1,1)
  double minWt = 0.0;
  // Competing Models
  int N_PAIR;
  arma::imat MODEL_PAIR;
  arma::irowvec dParas, varParasLoc0;
  arma::mat paras, parasInit, parasUpper, parasLower, parasBdd;
  // Max-min Discrimination Design
  arma::rowvec std_vals;
} OBJ_INFO, *Ptr_OBJ_INFO;

typedef struct model_diff_func {
  Rcpp::EvalBase2* M1_FUNC;
  Rcpp::EvalBase2* M2_FUNC;
  Rcpp::EvalBase2* V1_FUNC;
  Rcpp::EvalBase2* V2_FUNC;
  Rcpp::EvalBase4* DISTFUNC;
public:
    model_diff_func() : M1_FUNC(0), M2_FUNC(0), DISTFUNC(0) {}
    model_diff_func(Rcpp::EvalBase2* m1, Rcpp::EvalBase2* m2, Rcpp::EvalBase4* d) : M1_FUNC(m1), M2_FUNC(m2), DISTFUNC(d) {}
} model_diff_func;

typedef struct lbfgs_eval {
  model_diff_func* func_input;
  arma::rowvec UPPER;
  arma::rowvec LOWER;
  arma::irowvec NBD;
  int T_VAR_PARA_L0;
  int R_VAR_PARA_L0;
  arma::rowvec T_PARA;
  arma::rowvec R_PARA;
  arma::mat DESIGN;
  arma::rowvec WT;
  double FD_DELTA;
} lbfgs_eval, *Ptr_lbfgs_eval;

typedef struct best_alpha_info {
  double CRIT_VAL;
  arma::mat DESIGN;
} best_alpha_info, *Ptr_best_alpha_info;

typedef struct inner_pso_info {
  int PAIRID;
  arma::mat DESIGN;
  arma::rowvec WT;
} inner_pso_info, *Ptr_inner_pso_info;

#define N_PSO_OPTS 100

// DEFINE STUCTURES OF PSO INFORMATION
struct PSO_OPTIONS {
	/* Copy PSO parameters to 'genCppCode_PSOParaSetting.r'
		 for generating required C++ and R codes.
		 Please DO follow the format below.
	   type1 name1; // default1
	   type2 name2; // default2
	   type3 name3; // default3
	*/
	// Basic Settings
	int nSwarm; // 64
	int dSwarm; // 2
  rowvec varUpper; // c(1,1)
  rowvec varLower; // c(0,0)
	int maxIter; // 100
  double freeRun; // 1
 	double tol; // 1e-6
 	// Basic PSO Parameters
 	double c1; // 2.05
 	double c2; // 2.05
 	double w0; // 1.2
 	double w1; // 0.2
 	double w_var; // 0.8
 	double vk; // 4
};

// DEFINE STUCTURES OF PSO PARAMETERS WHICH WILL CHANGE ITERATIVELY
typedef struct {
  // inertia weight
  int w_varyfor; // (int)(w_var*maxIter);
  double w_cur; // w0;
  double w_dec; // (w0 - w1)/w_varyfor;
} PSO_DYN, *Ptr_PSO_DYN;

typedef struct LBFGS_PARAM {
  // LBFGS OPTIONS
  int IF_INNER_LBFGS; // 1
  int LBFGS_RETRY; // 3
  int LBFGS_MAXIT; // 100
  // -- CORRECTION SETTINGS
  int LBFGS_LM; // 6
  // -- STOPPING CRITERION SETTING
  double FVAL_EPS; // 0.0   (1e-8)
  double GRAD_EPS; // 1e-5  (1e-6)
  // -- LINE SEARCH PARAMETERS
  int LINESEARCH_MAXTRIAL; // 50
  double LINESEARCH_MAX; // 1e6
  double LINESEARCH_MIN; // 1e-6
  double LINESEARCH_ARMIJO; // 1e-4
  double LINESEARCH_WOLFE; // 0.9
  double FD_DELTA; // 1e-3
} LBFGS_PARAM, *Ptr_LBFGS_PARAM;

typedef struct FED_PARAM {
  // FEDOROV OPTIONS
  int FED_MAXIT; // 100
  int FED_TRIM; // 5
  double FED_TRIM_EPS; // 1e-4
  // -- STOPPING CRITERION SETTING
  double freeRun; // 1.0
  double FED_EPS; // 0.0   (1e-8)
  // UPDATING PARAMETER
  int FED_ALPHA_GRID; // 20
} FED_PARAM, *Ptr_FED_PARAM;

typedef struct REMES_PARAM {
  // FEDOROV OPTIONS
  int REMES_MAXIT; // 100
  // -- STOPPING CRITERION SETTING
  double freeRun; // 1.0
  double REMES_EPS; // 0.0   (1e-8)
} REMES_PARAM, *Ptr_REMES_PARAM;


// DEFINE PSO RESULTS
typedef struct {
  arma::rowvec GBest;
  double fGBest;
  arma::rowvec fGBestHist;
  arma::mat PBest;
  arma::vec fPBest;
} PSO_Result, *Ptr_PSO_Result;

// DEFINE FEDOROV-WYNN RESULTS
typedef struct {
  arma::mat DESIGN;
  arma::rowvec WT;
  double F_VAL;
  arma::rowvec R_PARA;
  arma::rowvec fvalHist;
} FED_Result, *Ptr_FED_Result;

// DEFINE 1DREMES RESULTS
typedef struct {
  arma::mat DESIGN;
  arma::rowvec WT;
  arma::mat DD_DEV;
  arma::rowvec R_PARA;
  arma::rowvec CPolyVal;
} REMES_Result, *Ptr_REMES_Result;


// DECLARE FUNCTIONS
void matrixPrintf(const mat &m);
void rvecPrintf(const rowvec &v);
void getAlgStruct(PSO_OPTIONS PSO_OPT[], const Rcpp::List PSO_INFO_LIST);
void getNewtonStruct(LBFGS_PARAM &LFBGS_OPT, const Rcpp::List LBFGS_INFO_LIST);
void getInfoStruct(OBJ_INFO &OBJ, const Rcpp::List OBJ_INFO_LIST);
void getFedorovStruct(FED_PARAM &FED_OPT, const Rcpp::List FED_INFO_LIST);
void getRemesStruct(REMES_PARAM &REMES_OPT, const Rcpp::List REMES_INFO_LIST);
void PSO_MAIN(const int LOOPID, PSO_OPTIONS PSO_OPTS[], const LBFGS_PARAM LBFGS_OPTION, const OBJ_INFO OBJ, model_diff_func *MODEL_COLLECTOR[],
              void *external, const bool IF_PARALLEL, const bool COUNTER_ON, PSO_Result &PSO_Result);
void psoUpdateParticle(const int LOOPID, PSO_OPTIONS PSO_OPTS[], const PSO_DYN PSO_DYN,
                       const arma::mat PBest, const arma::rowvec GBest,
                       const arma::rowvec velMax, const arma::rowvec varUpper, const arma::rowvec varLower,
                       arma::mat &vStep, arma::mat &swarm);
void psoCheckParticle(const int LOOPID, PSO_OPTIONS PSO_OPTS[], const PSO_DYN PSO_DYN,
                      const arma::rowvec varUpper, const arma::rowvec varLower, arma::mat &swarm);
void psoUpdateDynPara(const int LOOPID, PSO_OPTIONS PSO_OPTS[], const int iter, PSO_DYN &PSO_DYN,
                      const arma::mat swarm, const arma::mat PBest, const arma::rowvec GBest,
                      const arma::vec fSwarm, const arma::vec fPBest, const double fGBest);
void psoFuncEval(const bool IF_PARALLEL, const int LOOPID, PSO_OPTIONS PSO_OPTS[], const LBFGS_PARAM LBFGS_OPTION, const OBJ_INFO OBJ,
                 const PSO_DYN PSO_DYN, model_diff_func *MODEL_COLLECTOR[], void *external, const mat swarm, vec &fSwarm);

#include "psoFuncEval.h"
#include "psoCheckParticle.h"
#include "psoUpdateParticle.h"
#include "psoUpdateDynPara.h"
#include "psoKernel.h"
#include "fedorovwynn.h"
#include "remes.h"

// BODY
void matrixPrintf(const mat &m)
{
  for (uword i = 0; i < m.n_rows; i++) {
    for (uword j = 0; j < m.n_cols; j++) Rprintf("%4.4f\t", m(i,j));
    Rprintf("\n");
  }
	Rprintf("\n\n");
}

void rvecPrintf(const rowvec &v)
{
  for (uword i = 0; i < v.n_elem; i++) Rprintf("%4.4f\t", v(i));
  Rprintf("\n\n");
}

void getInfoStruct(OBJ_INFO &OBJ, const Rcpp::List OBJ_INFO_LIST)
{
  OBJ.crit_type = as<int>(OBJ_INFO_LIST["CRIT_TYPE_NUM"]);

  OBJ.d_type = as<int>(OBJ_INFO_LIST["D_TYPE_NUM"]);

  //OBJ_INFO->nSubj     = as<int>(OBJ_INFO_LIST["nSubj"]);
  OBJ.dSupp  = as<int>(OBJ_INFO_LIST["dSupp"]);
  OBJ.nSupp  = as<int>(OBJ_INFO_LIST["nSupp"]);

  Rcpp::NumericVector dsLower_Tmp = as<NumericVector>(OBJ_INFO_LIST["dsLower"]);
  arma::rowvec dsLower(dsLower_Tmp.begin(), dsLower_Tmp.size(), false);
  OBJ.dsLower = dsLower;
  Rcpp::NumericVector dsUpper_Tmp = as<NumericVector>(OBJ_INFO_LIST["dsUpper"]);
  arma::rowvec dsUpper(dsUpper_Tmp.begin(), dsUpper_Tmp.size(), false);
  OBJ.dsUpper = dsUpper;

  OBJ.minWt = as<double>(OBJ_INFO_LIST["minWt"]);
  OBJ.N_PAIR = as<int>(OBJ_INFO_LIST["N_PAIR"]);

  Rcpp::IntegerMatrix MODEL_PAIR_Tmp   = as<IntegerMatrix>(OBJ_INFO_LIST["MODEL_PAIR"]);
  arma::imat MODEL_PAIR(MODEL_PAIR_Tmp.begin(), MODEL_PAIR_Tmp.nrow(), MODEL_PAIR_Tmp.ncol(), false);
  OBJ.MODEL_PAIR  = MODEL_PAIR;

  Rcpp::IntegerVector dParas_Tmp  = as<IntegerVector>(OBJ_INFO_LIST["dParas"]);
  arma::irowvec dParas(dParas_Tmp.begin(), dParas_Tmp.size(), false);
  OBJ.dParas  = dParas;

  Rcpp::IntegerVector vParasL_Tmp  = as<IntegerVector>(OBJ_INFO_LIST["varParasLoc0"]);
  arma::irowvec varParasLoc0(vParasL_Tmp.begin(), vParasL_Tmp.size(), false);
  OBJ.varParasLoc0  = varParasLoc0;

  Rcpp::NumericMatrix paras_Tmp   = as<NumericMatrix>(OBJ_INFO_LIST["paras"]);
  arma::mat paras(paras_Tmp.begin(), paras_Tmp.nrow(), paras_Tmp.ncol(), false);
  OBJ.paras  = paras;

  Rcpp::NumericMatrix parasInit_Tmp   = as<NumericMatrix>(OBJ_INFO_LIST["parasInit"]);
  arma::mat parasInit(parasInit_Tmp.begin(), parasInit_Tmp.nrow(), parasInit_Tmp.ncol(), false);
  OBJ.parasInit  = parasInit;

  Rcpp::NumericMatrix parasUpper_Tmp   = as<NumericMatrix>(OBJ_INFO_LIST["parasUpper"]);
  arma::mat parasUpper(parasUpper_Tmp.begin(), parasUpper_Tmp.nrow(), parasUpper_Tmp.ncol(), false);
  OBJ.parasUpper  = parasUpper;

  Rcpp::NumericMatrix parasLower_Tmp   = as<NumericMatrix>(OBJ_INFO_LIST["parasLower"]);
  arma::mat parasLower(parasLower_Tmp.begin(), parasLower_Tmp.nrow(), parasLower_Tmp.ncol(), false);
  OBJ.parasLower  = parasLower;

  Rcpp::NumericMatrix parasBdd_Tmp   = as<NumericMatrix>(OBJ_INFO_LIST["parasBdd"]);
  arma::mat parasBdd(parasBdd_Tmp.begin(), parasBdd_Tmp.nrow(), parasBdd_Tmp.ncol(), false);
  OBJ.parasBdd  = parasBdd;

  // Max-min Discrimination Design
  Rcpp::NumericVector std_vals_Tmp  = as<NumericVector>(OBJ_INFO_LIST["MaxMinStdVals"]);
  arma::rowvec std_vals(std_vals_Tmp.begin(), std_vals_Tmp.size(), false);
  OBJ.std_vals = std_vals;
}

// LBFGS OPTIONS
void getNewtonStruct(LBFGS_PARAM &LFBGS_OPT, const Rcpp::List LBFGS_INFO_LIST)
{
  LFBGS_OPT.IF_INNER_LBFGS      = as<int>(LBFGS_INFO_LIST["IF_INNER_LBFGS"]);
  LFBGS_OPT.LBFGS_RETRY         = as<int>(LBFGS_INFO_LIST["LBFGS_RETRY"]);
  LFBGS_OPT.LBFGS_MAXIT         = as<int>(LBFGS_INFO_LIST["LBFGS_MAXIT"]);
  LFBGS_OPT.LBFGS_LM            = as<int>(LBFGS_INFO_LIST["LBFGS_LM"]);
  LFBGS_OPT.FVAL_EPS            = as<double>(LBFGS_INFO_LIST["FVAL_EPS"]);
  LFBGS_OPT.GRAD_EPS            = as<double>(LBFGS_INFO_LIST["GRAD_EPS"]);
  LFBGS_OPT.LINESEARCH_MAXTRIAL = as<int>(LBFGS_INFO_LIST["LINESEARCH_MAXTRIAL"]);
  LFBGS_OPT.LINESEARCH_MAX      = as<double>(LBFGS_INFO_LIST["LINESEARCH_MAX"]);
  LFBGS_OPT.LINESEARCH_MIN      = as<double>(LBFGS_INFO_LIST["LINESEARCH_MIN"]);
  LFBGS_OPT.LINESEARCH_ARMIJO   = as<double>(LBFGS_INFO_LIST["LINESEARCH_ARMIJO"]);
  LFBGS_OPT.LINESEARCH_WOLFE    = as<double>(LBFGS_INFO_LIST["LINESEARCH_WOLFE"]);
  LFBGS_OPT.FD_DELTA            = as<double>(LBFGS_INFO_LIST["FD_DELTA"]);
}

// PSO OPTIONS
void getAlgStruct(PSO_OPTIONS PSO_OPT[], const Rcpp::List PSO_INFO_LIST)
{
  Rcpp::IntegerVector nSwarm_Tmp     = as<IntegerVector>(PSO_INFO_LIST["nSwarm"]);
  Rcpp::IntegerVector dSwarm_Tmp     = as<IntegerVector>(PSO_INFO_LIST["dSwarm"]);
  Rcpp::NumericMatrix varUpper_Tmp   = as<NumericMatrix>(PSO_INFO_LIST["varUpper"]);
  arma::mat varUpper(varUpper_Tmp.begin(), varUpper_Tmp.nrow(), varUpper_Tmp.ncol(), false);
  Rcpp::NumericMatrix varLower_Tmp   = as<NumericMatrix>(PSO_INFO_LIST["varLower"]);
  arma::mat varLower(varLower_Tmp.begin(), varLower_Tmp.nrow(), varLower_Tmp.ncol(), false);
  Rcpp::IntegerVector maxIter_Tmp    = as<IntegerVector>(PSO_INFO_LIST["maxIter"]);
  Rcpp::NumericVector freeRun_Tmp    = as<NumericVector>(PSO_INFO_LIST["freeRun"]);
  Rcpp::NumericVector tol_Tmp        = as<NumericVector>(PSO_INFO_LIST["tol"]);
  Rcpp::NumericVector c1_Tmp         = as<NumericVector>(PSO_INFO_LIST["c1"]);
  Rcpp::NumericVector c2_Tmp         = as<NumericVector>(PSO_INFO_LIST["c2"]);
  Rcpp::NumericVector w0_Tmp         = as<NumericVector>(PSO_INFO_LIST["w0"]);
  Rcpp::NumericVector w1_Tmp         = as<NumericVector>(PSO_INFO_LIST["w1"]);
  Rcpp::NumericVector w_var_Tmp      = as<NumericVector>(PSO_INFO_LIST["w_var"]);
  Rcpp::NumericVector vk_Tmp         = as<NumericVector>(PSO_INFO_LIST["vk"]);

  int N_OPTS = nSwarm_Tmp.size();

  for (int i = 0; i < N_OPTS; i++) {
    PSO_OPT[i].nSwarm     = nSwarm_Tmp[i];
    PSO_OPT[i].dSwarm     = dSwarm_Tmp[i];
    PSO_OPT[i].varUpper   = varUpper.submat(i, 0, i, dSwarm_Tmp[i] - 1);
    PSO_OPT[i].varLower   = varLower.submat(i, 0, i, dSwarm_Tmp[i] - 1);
    PSO_OPT[i].maxIter    = maxIter_Tmp[i];
    PSO_OPT[i].freeRun    = freeRun_Tmp[i];
    PSO_OPT[i].tol        = tol_Tmp[i];
    PSO_OPT[i].c1         = c1_Tmp[i];
    PSO_OPT[i].c2         = c2_Tmp[i];
    PSO_OPT[i].w0         = w0_Tmp[i];
    PSO_OPT[i].w1         = w1_Tmp[i];
    PSO_OPT[i].w_var      = w_var_Tmp[i];
    PSO_OPT[i].vk         = vk_Tmp[i];
  }
}


// FEDOROV-WYNN OPTIONS
void getFedorovStruct(FED_PARAM &FED_OPT, const Rcpp::List FED_INFO_LIST)
{
  // FEDOROV OPTIONS
  FED_OPT.FED_MAXIT    = as<int>(FED_INFO_LIST["FED_MAXIT"]);
  FED_OPT.FED_TRIM     = as<int>(FED_INFO_LIST["FED_TRIM"]);
  FED_OPT.FED_TRIM_EPS = as<double>(FED_INFO_LIST["FED_TRIM_EPS"]);
  // -- STOPPING CRITERION SETTING
  FED_OPT.freeRun      = as<double>(FED_INFO_LIST["freeRun"]);
  FED_OPT.FED_EPS      = as<double>(FED_INFO_LIST["FED_EPS"]);
  // UPDATING PARAMETER
  FED_OPT.FED_ALPHA_GRID = as<int>(FED_INFO_LIST["FED_ALPHA_GRID"]);
}

// REMES OPTIONS
void getRemesStruct(REMES_PARAM &REMES_OPT, const Rcpp::List REMES_INFO_LIST)
{
  // FEDOROV OPTIONS
  REMES_OPT.REMES_MAXIT    = as<int>(REMES_INFO_LIST["REMES_MAXIT"]);
  // -- STOPPING CRITERION SETTING
  REMES_OPT.freeRun      = as<double>(REMES_INFO_LIST["freeRun"]);
  REMES_OPT.REMES_EPS      = as<double>(REMES_INFO_LIST["REMES_EPS"]);
}


