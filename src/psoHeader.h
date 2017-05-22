// Rcpp Header File
#include <math.h>
#include <RcppArmadillo.h>

//using namespace Rcpp;
using namespace arma;

#ifndef Rcpp_PSO_EVALUATE_H_
#define Rcpp_PSO_EVALUATE_H_

namespace Rcpp {

  class EvalBase {
    public:
        EvalBase() : neval(0) {};
        virtual arma::rowvec eval(SEXP x, SEXP p) = 0;
        //unsigned long getNbEvals() { return neval; }
    protected:
        //unsigned long int neval;
        int neval;
  };

  class EvalStandard : public EvalBase {
    public:
        EvalStandard(SEXP fcall_, SEXP env_) : fcall(fcall_), env(env_) {}
        arma::rowvec eval(SEXP x, SEXP p) {
          //neval++;
          return defaultfun(x, p);
        }
    private:
        SEXP fcall, env;
        arma::rowvec defaultfun(SEXP x, SEXP p) {
          SEXP fn = ::Rf_lang4(fcall, x, p, R_DotsSymbol);
          SEXP sexp_fvec = ::Rf_eval(fn, env);
          Rcpp::NumericVector f_result_tmp = (Rcpp::NumericVector) Rcpp::as<Rcpp::NumericVector>(sexp_fvec);
          arma::rowvec f_result(f_result_tmp.begin(), f_result_tmp.size(), false);
          return f_result;
        }
    };

  typedef arma::rowvec (*funcPtr)(SEXP, SEXP, SEXP);

  class EvalCompiled : public EvalBase {
    public:
        EvalCompiled(Rcpp::XPtr<funcPtr> xptr, SEXP __env) {
          funptr = *(xptr);
          env = __env;
        };
        EvalCompiled(SEXP xps, SEXP __env) {
          Rcpp::XPtr<funcPtr> xptr(xps);
          funptr = *(xptr);
          env = __env;
        };
        arma::rowvec eval(SEXP x, SEXP p) {
          //neval++;
          arma::rowvec f_result = funptr(x, p, env);
          return f_result;
        }
    private:
        funcPtr funptr;
        SEXP env;
    };
}

#endif

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

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
  // Competing Models
  int N_model;
  arma::irowvec dParas;
  arma::mat paras, parasInit, parasUpper, parasLower, parasBdd;
  // Max-min Discrimination Design
  arma::rowvec std_vals;
} OBJ_INFO, *Ptr_OBJ_INFO;


struct MODEL_SET {
  Rcpp::EvalBase* modelFunc;
};

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
	int checkConv; // 0
	int typePSO; // 0
  double freeRun; // 0.25
 	double tol; // 1e-6
 	// Basic PSO Parameters
 	double c1; // 2.05
 	double c2; // 2.05
 	double w0; // 1.2
 	double w1; // 0.2
 	double w_var; // 0.8
 	double chi; // 0.729
 	double vk; // 4
  // Topology
  int typeTopo; // 0
  int nGroup; // 1
 	// Guarantee Convergence PSO Parameters
	int GC_S_ROOF; // 5
	int GC_F_ROOF; // 15
	double GC_RHO; // 1
	// Quantum PSO Parameters
  int Q_cen_type; // 0
	double Q_a0; // 1.7
	double Q_a1; // 0.7
	double Q_a_var; // 0.8
  // LcRiPSO
  double LcRi_L; // 0.01
};

// DEFINE STUCTURES OF PSO PARAMETERS WHICH WILL CHANGE ITERATIVELY
typedef struct {
  int succ_GB; // 0
  mat network; // 0
  // inertia weight
  int w_varyfor; // (int)(w_var*maxIter);
  double w_cur; // w0;
  double w_dec; // (w0 - w1)/w_varyfor;
  // GCPSO
  int GC_S_COUNT; // 0;
  int GC_F_COUNT; // 0;
  double GC_RHO;
  // Quantum PSO
  int Q_a_varyfor; // (int)(Q_a_var*maxIter);
  double Q_a_cur; // Q_a0;
  double Q_a_dec; // (Q_a0 - Q_a1)/Q_a_varyfor;
  // LcRiPSO
  vec LcRi_sigP;
  vec LcRi_sigG;
} PSO_DYN, *Ptr_PSO_DYN;

// DEFINE PSO RESULTS
typedef struct {
  arma::rowvec GBest;
  double fGBest;
  arma::rowvec fGBestHist;
  arma::mat PBest;
  arma::vec fPBest;
} PSO_Result, *Ptr_PSO_Result;


// DECLARE FUNCTIONS
void matrixPrintf(const mat &m);
void rvecPrintf(const rowvec &v);
void getAlgStruct(PSO_OPTIONS PSO_OPT[], const Rcpp::List &ALG_INFO_LIST);
//void getModelFunc(MODEL_SET MODELS[], const Rcpp::List MODEL_LIST, SEXP env, const int N_model);
void getInfoStruct(OBJ_INFO &OBJ, const Rcpp::List OBJ_INFO_LIST);
void PSO_MAIN(const int &LOOPID, const PSO_OPTIONS PSO_OPTS[], const OBJ_INFO &OBJ, const MODEL_SET MODELS[], Rcpp::EvalBase* distFunc, 
              const rowvec &FIXEDVALUE, const bool &IF_PARALLEL, const bool COUNTER_ON, Ptr_PSO_Result Ptr_PSO_Result);
void psoUpdateParticle(const int &LOOPID, const PSO_OPTIONS PSO_OPTS[], const PSO_DYN &PSO_DYN,
                       const arma::mat &PBest, const arma::rowvec &GBest,
                       const arma::rowvec &velMax, const arma::rowvec &varUpper, const arma::rowvec &varLower,
                       arma::mat &vStep, arma::mat &swarm);
void psoCheckParticle(const int &LOOPID, const PSO_OPTIONS PSO_OPTS[], const PSO_DYN &PSO_DYN,
                      const arma::rowvec &varUpper, const arma::rowvec &varLower, arma::mat &swarm);
void psoUpdateDynPara(const int &LOOPID, const PSO_OPTIONS PSO_OPTS[], const int iter, PSO_DYN &PSO_DYN,
                      const arma::mat &swarm, const arma::mat &PBest, const arma::rowvec &GBest,
                      const arma::vec &fSwarm, const arma::vec &fPBest, const double &fGBest);
void psoFuncEval(const bool &IF_PARALLEL, const int &LOOPID, const PSO_OPTIONS PSO_OPTS[], const OBJ_INFO &OBJ, const PSO_DYN &PSO_DYN, 
                 const MODEL_SET MODELS[], Rcpp::EvalBase* distFunc, const rowvec &FIXEDVALUE, const mat &swarm, vec &fSwarm);

#include "psoFuncEval.h"
#include "psoCheckParticle.h"
#include "psoUpdateParticle.h"
#include "psoUpdateDynPara.h"
#include "psoKernel.h"

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

  OBJ.N_model = as<int>(OBJ_INFO_LIST["N_model"]);

  Rcpp::IntegerVector dParas_Tmp  = as<IntegerVector>(OBJ_INFO_LIST["dParas"]);
  arma::irowvec dParas(dParas_Tmp.begin(), dParas_Tmp.size(), false);
  OBJ.dParas  = dParas;

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

void getAlgStruct(PSO_OPTIONS PSO_OPT[], const Rcpp::List &ALG_INFO_LIST)
{
  Rcpp::IntegerVector nSwarm_Tmp     = as<IntegerVector>(ALG_INFO_LIST["nSwarm"]);
  Rcpp::IntegerVector dSwarm_Tmp     = as<IntegerVector>(ALG_INFO_LIST["dSwarm"]);
  Rcpp::NumericMatrix varUpper_Tmp   = as<NumericMatrix>(ALG_INFO_LIST["varUpper"]);
  arma::mat varUpper(varUpper_Tmp.begin(), varUpper_Tmp.nrow(), varUpper_Tmp.ncol(), false);
  Rcpp::NumericMatrix varLower_Tmp   = as<NumericMatrix>(ALG_INFO_LIST["varLower"]);
  arma::mat varLower(varLower_Tmp.begin(), varLower_Tmp.nrow(), varLower_Tmp.ncol(), false);
  Rcpp::IntegerVector maxIter_Tmp    = as<IntegerVector>(ALG_INFO_LIST["maxIter"]);
  Rcpp::IntegerVector checkConv_Tmp  = as<IntegerVector>(ALG_INFO_LIST["checkConv"]);
  Rcpp::IntegerVector typePSO_Tmp    = as<IntegerVector>(ALG_INFO_LIST["typePSO"]);
  Rcpp::NumericVector freeRun_Tmp    = as<NumericVector>(ALG_INFO_LIST["freeRun"]);
  Rcpp::NumericVector tol_Tmp        = as<NumericVector>(ALG_INFO_LIST["tol"]);
  Rcpp::NumericVector c1_Tmp         = as<NumericVector>(ALG_INFO_LIST["c1"]);
  Rcpp::NumericVector c2_Tmp         = as<NumericVector>(ALG_INFO_LIST["c2"]);
  Rcpp::NumericVector w0_Tmp         = as<NumericVector>(ALG_INFO_LIST["w0"]);
  Rcpp::NumericVector w1_Tmp         = as<NumericVector>(ALG_INFO_LIST["w1"]);
  Rcpp::NumericVector w_var_Tmp      = as<NumericVector>(ALG_INFO_LIST["w_var"]);
  Rcpp::NumericVector chi_Tmp        = as<NumericVector>(ALG_INFO_LIST["chi"]);
  Rcpp::NumericVector vk_Tmp         = as<NumericVector>(ALG_INFO_LIST["vk"]);
  Rcpp::IntegerVector typeTopo_Tmp   = as<IntegerVector>(ALG_INFO_LIST["typeTopo"]);
  Rcpp::IntegerVector nGroup_Tmp     = as<IntegerVector>(ALG_INFO_LIST["nGroup"]);
  Rcpp::IntegerVector GC_S_ROOF_Tmp  = as<IntegerVector>(ALG_INFO_LIST["GC_S_ROOF"]);
  Rcpp::IntegerVector GC_F_ROOF_Tmp  = as<IntegerVector>(ALG_INFO_LIST["GC_F_ROOF"]);
  Rcpp::NumericVector GC_RHO_Tmp     = as<NumericVector>(ALG_INFO_LIST["GC_RHO"]);
  Rcpp::IntegerVector Q_cen_type_Tmp = as<IntegerVector>(ALG_INFO_LIST["Q_cen_type"]);
  Rcpp::NumericVector Q_a0_Tmp       = as<NumericVector>(ALG_INFO_LIST["Q_a0"]);
  Rcpp::NumericVector Q_a1_Tmp       = as<NumericVector>(ALG_INFO_LIST["Q_a1"]);
  Rcpp::NumericVector Q_a_var_Tmp    = as<NumericVector>(ALG_INFO_LIST["Q_a_var"]);
  Rcpp::NumericVector LcRi_L_Tmp     = as<NumericVector>(ALG_INFO_LIST["LcRi_L"]);

  int N_OPTS = nSwarm_Tmp.size();

  for (int i = 0; i < N_OPTS; i++) {
    PSO_OPT[i].nSwarm     = nSwarm_Tmp[i];
    PSO_OPT[i].dSwarm     = dSwarm_Tmp[i];
    PSO_OPT[i].varUpper   = varUpper.submat(i, 0, i, dSwarm_Tmp[i] - 1);
    PSO_OPT[i].varLower   = varLower.submat(i, 0, i, dSwarm_Tmp[i] - 1);
    PSO_OPT[i].maxIter    = maxIter_Tmp[i];
    PSO_OPT[i].checkConv  = checkConv_Tmp[i];
    PSO_OPT[i].typePSO    = typePSO_Tmp[i];
    PSO_OPT[i].freeRun    = freeRun_Tmp[i];
    PSO_OPT[i].tol        = tol_Tmp[i];
    PSO_OPT[i].c1         = c1_Tmp[i];
    PSO_OPT[i].c2         = c2_Tmp[i];
    PSO_OPT[i].w0         = w0_Tmp[i];
    PSO_OPT[i].w1         = w1_Tmp[i];
    PSO_OPT[i].w_var      = w_var_Tmp[i];
    PSO_OPT[i].chi        = chi_Tmp[i];
    PSO_OPT[i].vk         = vk_Tmp[i];
    PSO_OPT[i].typeTopo   = typeTopo_Tmp[i];
    PSO_OPT[i].nGroup     = nGroup_Tmp[i];
    PSO_OPT[i].GC_S_ROOF  = GC_S_ROOF_Tmp[i];
    PSO_OPT[i].GC_F_ROOF  = GC_F_ROOF_Tmp[i];
    PSO_OPT[i].GC_RHO     = GC_RHO_Tmp[i];
    PSO_OPT[i].Q_cen_type = Q_cen_type_Tmp[i];
    PSO_OPT[i].Q_a0       = Q_a0_Tmp[i];
    PSO_OPT[i].Q_a1       = Q_a1_Tmp[i];
    PSO_OPT[i].Q_a_var    = Q_a_var_Tmp[i];
    PSO_OPT[i].LcRi_L     = LcRi_L_Tmp[i];
  }
}

/*
void getTopology(const int &LOOPID, const PSO_OPTIONS PSO_OPT[]) {

  int topo = PSO_OPT[LOOPID].topology;
  switch (topo) {
    case 1: // Ring Topology
    {
      rowvec ring_set(3); uword min_id; double min_val;
      int nSwarm = PSO_OPT[LOOPID].nSwarm;
      for (i = 0; i < nSwarm; i++) {
        if (i == 0) {
          ring_set << fSwarm(nSwarm - 1) << fSwarm(i) << fSwarm(i+1) << endr;
          min_val = ring_set.min(min_id);
          if (min_val < fPBest(i)) {
            fPBest(i) = min_val;
            if (min_id == 0) PBest.row(i) = swarm.row(nSwarm - 1); else PBest.row(i) = swarm.row(i + min_id - 1);
          }
        } else if (i == (nSwarm - 1)) {
          ring_set << fSwarm(i-1) << fSwarm(i) << fSwarm(0) << endr;
          min_val = ring_set.min(min_id);
          if (min_val < fPBest(i)) {
            fPBest(i) = min_val;
            if (min_id == 2) PBest.row(i) = swarm.row(0); else PBest.row(i) = swarm.row(i + min_id - 1);
          }
        } else {
          ring_set << fSwarm(i-1) << fSwarm(i) << fSwarm(i+1) << endr;
          min_val = ring_set.min(min_id);
          if (min_val < fPBest(i)) {
            fPBest(i) = min_val; PBest.row(i) = swarm.row(i + min_id - 1);
          }
        }
      }
      break;
    }
  }
}
*/
