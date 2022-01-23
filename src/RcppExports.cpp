// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// cppPSO
Rcpp::List cppPSO(const int LOOPID, Rcpp::List PSO_INFO_LIST, Rcpp::List LBFGS_INFO_LIST, Rcpp::List OBJ_INFO_LIST, Rcpp::List MODEL_INFO_LIST, Rcpp::List EXTERNAL_LIST, const SEXP env, const bool IF_PARALLEL, const bool VERBOSE);
RcppExport SEXP _DiscrimOD_cppPSO(SEXP LOOPIDSEXP, SEXP PSO_INFO_LISTSEXP, SEXP LBFGS_INFO_LISTSEXP, SEXP OBJ_INFO_LISTSEXP, SEXP MODEL_INFO_LISTSEXP, SEXP EXTERNAL_LISTSEXP, SEXP envSEXP, SEXP IF_PARALLELSEXP, SEXP VERBOSESEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type LOOPID(LOOPIDSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type PSO_INFO_LIST(PSO_INFO_LISTSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type LBFGS_INFO_LIST(LBFGS_INFO_LISTSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type OBJ_INFO_LIST(OBJ_INFO_LISTSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type MODEL_INFO_LIST(MODEL_INFO_LISTSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type EXTERNAL_LIST(EXTERNAL_LISTSEXP);
    Rcpp::traits::input_parameter< const SEXP >::type env(envSEXP);
    Rcpp::traits::input_parameter< const bool >::type IF_PARALLEL(IF_PARALLELSEXP);
    Rcpp::traits::input_parameter< const bool >::type VERBOSE(VERBOSESEXP);
    rcpp_result_gen = Rcpp::wrap(cppPSO(LOOPID, PSO_INFO_LIST, LBFGS_INFO_LIST, OBJ_INFO_LIST, MODEL_INFO_LIST, EXTERNAL_LIST, env, IF_PARALLEL, VERBOSE));
    return rcpp_result_gen;
END_RCPP
}
// cppDesignCriterion
Rcpp::List cppDesignCriterion(Rcpp::List PSO_INFO_LIST, Rcpp::List LBFGS_INFO_LIST, Rcpp::List OBJ_INFO_LIST, Rcpp::List MODEL_INFO_LIST, Rcpp::List EXTERNAL_LIST, SEXP env, arma::rowvec DESIGN);
RcppExport SEXP _DiscrimOD_cppDesignCriterion(SEXP PSO_INFO_LISTSEXP, SEXP LBFGS_INFO_LISTSEXP, SEXP OBJ_INFO_LISTSEXP, SEXP MODEL_INFO_LISTSEXP, SEXP EXTERNAL_LISTSEXP, SEXP envSEXP, SEXP DESIGNSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type PSO_INFO_LIST(PSO_INFO_LISTSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type LBFGS_INFO_LIST(LBFGS_INFO_LISTSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type OBJ_INFO_LIST(OBJ_INFO_LISTSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type MODEL_INFO_LIST(MODEL_INFO_LISTSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type EXTERNAL_LIST(EXTERNAL_LISTSEXP);
    Rcpp::traits::input_parameter< SEXP >::type env(envSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type DESIGN(DESIGNSEXP);
    rcpp_result_gen = Rcpp::wrap(cppDesignCriterion(PSO_INFO_LIST, LBFGS_INFO_LIST, OBJ_INFO_LIST, MODEL_INFO_LIST, EXTERNAL_LIST, env, DESIGN));
    return rcpp_result_gen;
END_RCPP
}
// cppEquivalence
Rcpp::List cppEquivalence(Rcpp::List OBJ_INFO_LIST, Rcpp::List MODEL_INFO_LIST, const double GBEST_VAL, const arma::mat PARA_SET, const arma::rowvec alpha, const SEXP env, const int nGrid);
RcppExport SEXP _DiscrimOD_cppEquivalence(SEXP OBJ_INFO_LISTSEXP, SEXP MODEL_INFO_LISTSEXP, SEXP GBEST_VALSEXP, SEXP PARA_SETSEXP, SEXP alphaSEXP, SEXP envSEXP, SEXP nGridSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type OBJ_INFO_LIST(OBJ_INFO_LISTSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type MODEL_INFO_LIST(MODEL_INFO_LISTSEXP);
    Rcpp::traits::input_parameter< const double >::type GBEST_VAL(GBEST_VALSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type PARA_SET(PARA_SETSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const SEXP >::type env(envSEXP);
    Rcpp::traits::input_parameter< const int >::type nGrid(nGridSEXP);
    rcpp_result_gen = Rcpp::wrap(cppEquivalence(OBJ_INFO_LIST, MODEL_INFO_LIST, GBEST_VAL, PARA_SET, alpha, env, nGrid));
    return rcpp_result_gen;
END_RCPP
}
// cppFedorovWynn
Rcpp::List cppFedorovWynn(Rcpp::List FED_INFO_LIST, Rcpp::List LBFGS_INFO_LIST, Rcpp::List OBJ_INFO_LIST, Rcpp::List MODEL_INFO_LIST, const SEXP env, const bool VERBOSE);
RcppExport SEXP _DiscrimOD_cppFedorovWynn(SEXP FED_INFO_LISTSEXP, SEXP LBFGS_INFO_LISTSEXP, SEXP OBJ_INFO_LISTSEXP, SEXP MODEL_INFO_LISTSEXP, SEXP envSEXP, SEXP VERBOSESEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type FED_INFO_LIST(FED_INFO_LISTSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type LBFGS_INFO_LIST(LBFGS_INFO_LISTSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type OBJ_INFO_LIST(OBJ_INFO_LISTSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type MODEL_INFO_LIST(MODEL_INFO_LISTSEXP);
    Rcpp::traits::input_parameter< const SEXP >::type env(envSEXP);
    Rcpp::traits::input_parameter< const bool >::type VERBOSE(VERBOSESEXP);
    rcpp_result_gen = Rcpp::wrap(cppFedorovWynn(FED_INFO_LIST, LBFGS_INFO_LIST, OBJ_INFO_LIST, MODEL_INFO_LIST, env, VERBOSE));
    return rcpp_result_gen;
END_RCPP
}
// cppUnifApprox
Rcpp::List cppUnifApprox(Rcpp::List REMES_INFO_LIST, Rcpp::List LBFGS_INFO_LIST, Rcpp::List OBJ_INFO_LIST, Rcpp::List MODEL_INFO_LIST, const SEXP env, const bool VERBOSE);
RcppExport SEXP _DiscrimOD_cppUnifApprox(SEXP REMES_INFO_LISTSEXP, SEXP LBFGS_INFO_LISTSEXP, SEXP OBJ_INFO_LISTSEXP, SEXP MODEL_INFO_LISTSEXP, SEXP envSEXP, SEXP VERBOSESEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type REMES_INFO_LIST(REMES_INFO_LISTSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type LBFGS_INFO_LIST(LBFGS_INFO_LISTSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type OBJ_INFO_LIST(OBJ_INFO_LISTSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type MODEL_INFO_LIST(MODEL_INFO_LISTSEXP);
    Rcpp::traits::input_parameter< const SEXP >::type env(envSEXP);
    Rcpp::traits::input_parameter< const bool >::type VERBOSE(VERBOSESEXP);
    rcpp_result_gen = Rcpp::wrap(cppUnifApprox(REMES_INFO_LIST, LBFGS_INFO_LIST, OBJ_INFO_LIST, MODEL_INFO_LIST, env, VERBOSE));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_DiscrimOD_cppPSO", (DL_FUNC) &_DiscrimOD_cppPSO, 9},
    {"_DiscrimOD_cppDesignCriterion", (DL_FUNC) &_DiscrimOD_cppDesignCriterion, 7},
    {"_DiscrimOD_cppEquivalence", (DL_FUNC) &_DiscrimOD_cppEquivalence, 7},
    {"_DiscrimOD_cppFedorovWynn", (DL_FUNC) &_DiscrimOD_cppFedorovWynn, 6},
    {"_DiscrimOD_cppUnifApprox", (DL_FUNC) &_DiscrimOD_cppUnifApprox, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_DiscrimOD(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
