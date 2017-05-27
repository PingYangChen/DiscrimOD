#' Hybrid Algorithm of PSO and BFGS for Finding Optimal Discrimination Design
#'
#' Please follow the instruction.
#'
#' @param MODEL_INFO list of information of competing models. For details, run \code{emptyModelList()} and see the instruction below.
#' @param DISTANCE function. The R/C++ function of distance measure.  For T-optimal design, the function is the squared difference of two models means.
#' @param nSupp integer. The number fo support points (at least 2).
#' @param dsLower vector. The finite lower bounds of the design space. Its length should be equal to the dimension of design space.
#' @param dsUpper vector. The finite upper bounds of the design space. Its length should be equal to the dimension of design space.
#' @param crit_type string. The name of the case of the discrimination design problem. The default is 'pair_fixed_true'.
#' @param MaxMinStdVals vector. The values of demoninators in the design efficiency calculation for finding max-min discrimination design.
#' @param ALG_INFO list. PSO and BFGS options.
#' @param seed numeric. random seed.  The default is \code{NULL}.
#' @param verbose logical. If \code{TRUE}, the PSO will reports the updating progress.
#' @return An List.
#' \itemize{
#' \item{BESTDESIGN}{ the resulting design. Each row is a support point with its weight. The last column is the weights of the support points.}
#' \item{BESTVAL}{ the design criterion value of the resulting design.}
#' \item{GBESTHIST}{ a vector of design citerion values of the global best particle in the PSO search history.}
#' \item{CPUTIME}{ the computational time in seconds.}
#' }
#' @details
#' 
#' @examples
#' # Atkinson and Fedorov (1975)
#' # Two R functions of competing models are given by
#' m1 <- function(x, p) p[1] + p[2]*exp(x) + p[3]*exp(-x)
#' m2 <- function(x, p) p[1] + p[2]*x + p[3]*x^2
#' # Set the model information
#' # The nominla value in 'm1' is 4.5, -1.5, -2.0
#' # For 'm2', we set the parameter space to be [-10, 10]^3 and 
#' # the initial guess (for LBFGS) of the rival model parameter is zero vector
#' AF_para_m1 <- c(4.5, -1.5, -2.0)
#' MODEL_INFO <- list(
#'   list(model = m1, para = AF_para_m1),
#'   list(model = m2,
#'        paraLower = rep(-10, 3),
#'        paraUpper = rep(10, 3),
#'        paraInit = c(0, 0, 0))
#' )
#' # Define the R function for the distance measure
#' # Here we use T-optimal criterion
#' DISTANCE <- function(xt, xr) (xt - xr)^2
#' 
#' # Initialize PSO and BFGS options
#' ALG_INFO <- getAlgInfo(nSwarm = 32, maxIter = 100, typePSO = 0,
#'                     		LBFGS_RETRY = 2, FVAL_EPS = 0, GRAD_EPS = 1e-6,
#'                     		LINESEARCH_MAX = 1e5)
#'
#' # Run Algorithm
#' out <- DiscrimOD(MODEL_INFO, DISTANCE, nSupp = 4, dsLower = -1.0, dsUpper = 1.0, crit_type = "pair_fixed_true",
#'                  MaxMinStdVals = NULL, ALG_INFO = ALG_INFO, seed = NULL, verbose = TRUE)
#'
#' 
#' # C++ Function Input
#' library(inline)
#' m1.inc <- 'arma::rowvec m1(SEXP xx, SEXP pp){
#'   arma::rowvec x = Rcpp::as<arma::rowvec>(xx);
#'   arma::rowvec p = Rcpp::as<arma::rowvec>(pp);
#'   return p(0) + p(1)*arma::exp(x) + p(2)*arma::exp(-x);
#' }'
#' m1.body <- '
#'   typedef arma::rowvec (*funcPtr)(SEXP, SEXP);
#'   return(XPtr<funcPtr>(new funcPtr(&m1)));
#' '
#' 
#' m2.inc <- 'arma::rowvec m2(SEXP xx, SEXP pp){
#'   arma::rowvec x = Rcpp::as<arma::rowvec>(xx);
#'   arma::rowvec p = Rcpp::as<arma::rowvec>(pp);
#'   return p(0) + p(1)*x + p(2)*(x%x);
#' }'
#' m2.body <- '
#'   typedef arma::rowvec (*funcPtr)(SEXP, SEXP);
#'   return(XPtr<funcPtr>(new funcPtr(&m2)));
#' '
#' 
#' m1_Cpp <- cxxfunction(signature(), body = m1.body, inc = m1.inc,
#'                       plugin = "RcppArmadillo")
#' m2_Cpp <- cxxfunction(signature(), body = m2.body, inc = m2.inc,
#'                       plugin = "RcppArmadillo")
#' 
#' MODEL_INFO_Cpp <- list(
#'   list(model = m1_Cpp(), para = AF_para_m1),
#'   list(model = m2_Cpp(), paraLower = rep(-10, 3),
#'                          paraUpper = rep(10, 3),
#'                          paraInit = c(0,0,0))
#' )
#' 
#' dist.inc <- 'arma::rowvec t_optimal(SEXP xt, SEXP xr){
#'   arma::rowvec val_t = Rcpp::as<arma::rowvec>(xt);
#'   arma::rowvec val_r = Rcpp::as<arma::rowvec>(xr);
#'   return (val_t - val_r)%(val_t - val_r);
#' }'
#' 
#' dist.body <- '
#'   typedef arma::rowvec (*funcPtr)(SEXP, SEXP);
#'   return(XPtr<funcPtr>(new funcPtr(&t_optimal)));
#' '
#' 
#' DISTANCE_Cpp <- cxxfunction(signature(), body = dist.body, inc = dist.inc,
#'                             plugin = "RcppArmadillo")
#' 
#' # Run Algorithm with C++ functions
#' out <- DiscrimOD(MODEL_INFO_Cpp, DISTANCE_Cpp(), nSupp = 4, dsLower = -1.0, dsUpper = 1.0, crit_type = "pair_fixed_true",
#'                  MaxMinStdVals = NULL, ALG_INFO = ALG_INFO, seed = NULL, verbose = TRUE)
#'
#' @name DiscrimOD
#' @rdname DiscrimOD
#' @export
DiscrimOD <- function(MODEL_INFO, DISTANCE, nSupp, dsLower, dsUpper, crit_type = "pair_fixed_true", MaxMinStdVals = NULL, 
											ALG_INFO = NULL, seed = NULL, verbose = TRUE, environment, ...) {

  assert_that(nSupp >= 2L, all(is.finite(dsLower)), all(is.finite(dsUpper)),
              length(dsLower) == length(dsUpper), all(dsUpper > dsLower))

	MODEL_LIST <- lapply(1:length(MODEL_INFO), function(k) MODEL_INFO[[k]]$model)

	if (is.null(MaxMinStdVals)) MaxMinStdVals <- 0
	D_INFO <- getDesignInfo(D_TYPE = "approx", MODEL_INFO = MODEL_INFO, dist_func = DISTANCE, 
                          crit_type = crit_type, MaxMinStdVals = MaxMinStdVals,
                          dSupp = length(dsLower), nSupp = nSupp, dsLower = dsLower, dsUpper = dsUpper)

	if (is.null(ALG_INFO)) {
		ALG_INFO <- getAlgInfo()
		if (verbose) cat(paste0("Use the default settings for PSO. See 'getAlgInfo()'.\n"))
	}

	if (!hasArg(environment)) environment <- new.env()

	# Adjust ALG_INFO according to D_INFO
	swarmSetting <- algInfoUpdate(D_INFO)
	ALG_INFO$varUpper <- swarmSetting$UB
	ALG_INFO$varLower <- swarmSetting$LB
	ALG_INFO$dSwarm	<- ncol(swarmSetting$UB)

	# Start
	set.seed(seed)
	cputime <- system.time(
		psoOut <- cppPSO(0, ALG_INFO, D_INFO, MODEL_LIST, 0, environment, FALSE, verbose)
	)[3]

	if (verbose) cat(paste0("CPU time: ", round(cputime, 2), " seconds.\n"))

	BESTDESIGN <- designM2V(psoOut$GBest, D_INFO)

	list(BESTDESIGN = BESTDESIGN, BESTVAL = -psoOut$fGBest, GBESTHIST = -psoOut$fGBestHist, 
	     CPUTIME = cputime)
}

#' Calculate Design Criterion Values
#'
#' Please follow the instruction.
#'
#' @param DESIGN1 matrix. The approximate design.  
#' @param MODEL_INFO list of information of competing models. For details, run \code{emptyModelList()} and see the instruction below.
#' @param DISTANCE function. The R/C++ function of distance measure.  For T-optimal design, the function is the squared difference of two models means.
#' @param dsLower vector. The finite lower bounds of the design space. Its length should be equal to the dimension of design space.
#' @param dsUpper vector. The finite upper bounds of the design space. Its length should be equal to the dimension of design space.
#' @param crit_type string. The name of the case of the discrimination design problem. The default is 'pair_fixed_true'.
#' @param MaxMinStdVals vector. The values of demoninators in the design efficiency calculation for finding max-min discrimination design.
#' @param ALG_INFO list. PSO and BFGS options.
#' @return An List.
#' \itemize{
#' \item{cri_val}{ the design criterion value of the resulting design in the PSO result or the given design.}
#' \item{theta2}{ a matrix of the parameters of each input model that result the design criterion value.}
#' }
#' @name designCriterion
#' @rdname designCriterion
#' @export
designCriterion <- function(DESIGN1, MODEL_INFO, DISTANCE, dsLower, dsUpper, crit_type = "pair_fixed_true", MaxMinStdVals = NULL, 
														ALG_INFO = NULL, environment, ...) {

	assert_that(all(is.finite(dsLower)), all(is.finite(dsUpper)),
              length(dsLower) == length(dsUpper), all(dsUpper > dsLower))

	nSupp <- nrow(DESIGN1)
	MODEL_LIST <- lapply(1:length(MODEL_INFO), function(k) MODEL_INFO[[k]]$model)

	if (is.null(MaxMinStdVals)) MaxMinStdVals <- 0
	D_INFO <- getDesignInfo(D_TYPE = "approx", MODEL_INFO = MODEL_INFO, dist_func = DISTANCE, 
                          crit_type = crit_type, MaxMinStdVals = MaxMinStdVals,
                          dSupp = length(dsLower), nSupp = nSupp, dsLower = dsLower, dsUpper = dsUpper)

	if (is.null(ALG_INFO)) {
		ALG_INFO <- getAlgInfo()
	}

	if (!hasArg(environment)) environment <- new.env()

	# Adjust ALG_INFO according to D_INFO
	swarmSetting <- algInfoUpdate(D_INFO)
	ALG_INFO$varUpper <- swarmSetting$UB
	ALG_INFO$varLower <- swarmSetting$LB
	ALG_INFO$dSwarm	<- ncol(swarmSetting$UB)

	DESIGN1_M <- designV2M(DESIGN1, D_INFO)

	cri_1 <- cppDesignCriterion(ALG_INFO, D_INFO, MODEL_LIST, 0, environment, DESIGN1_M)

	rownames(cri_1$theta2) <- paste0("model_", 1:length(MODEL_INFO))
	
  return(list(cri_val = -cri_1$val, theta2 = cri_1$theta2))
}

#' Equivalence Theorem
#'
#' Please follow the instruction.
#'
#' @param DESIGN1 matrix.
#' @param PSO_RESULT list.
#' @param ngrid integer. The default is 100.
#' @param IFPLOT logical.
#' @param MODEL_INFO list of information of competing models. For details, run \code{emptyModelList()} and see the instruction below.
#' @param DISTANCE function. 
#' @param dsLower vector. The finite lower bounds of the design space. Its length should be equal to the dimension of design space.
#' @param dsUpper vector. The finite upper bounds of the design space. Its length should be equal to the dimension of design space.
#' @param crit_type string. The name of the case of the discrimination design problem. The default is 'pair_fixed_true'.
#' @param MaxMinStdVals vector. The values of demoninators in the design efficiency calculation for finding max-min discrimination design.
#' @param ALG_INFO list. PSO and BFGS options.
#' @return An List.
#' \itemize{
#' \item{Grid_1}{ a vector of grid points on the design space.  For 2-D case, this vector is the grid of the first dimension.}
#' \item{Grid_2}{ a vector of grid points on the second dimension of the 2-D design space.  For 1_d case, this output can be ignored.}
#' \item{DirDeriv}{ a vector (1-D case) or a matrix (2-D case) of the directional derivative function values on the grid points.}
#' \item{MAX_DD}{ the maximal value of the directional derivative function in \code{DirDeriv}.}
#' \item{alpha}{ a vector of weights of efficiency values that satisfies the equivalence theorem for max-min optimal discrimination design.}
#' }
#' @name equivalence
#' @rdname equivalence
#' @export
equivalence <- function(DESIGN = NULL, PSO_RESULT = NULL, ngrid = 100, IFPLOT = FALSE,
												MODEL_INFO, DISTANCE, dsLower, dsUpper, crit_type = "pair_fixed_true", 
												MaxMinStdVals = NULL, ALG_INFO = NULL, environment, ...) {

	assert_that(all(is.finite(dsLower)), all(is.finite(dsUpper)),
              length(dsLower) == length(dsUpper), all(dsUpper > dsLower))

	if (is.null(DESIGN)) DESIGN <- PSO_RESULT$BESTDESIGN

	nSupp <- nrow(DESIGN)
	dSupp <- ncol(DESIGN) - 1

	MODEL_LIST <- lapply(1:length(MODEL_INFO), function(k) MODEL_INFO[[k]]$model)
	
	if (!hasArg(environment)) environment <- new.env()

	if (is.null(MaxMinStdVals)) MaxMinStdVals <- 0
	D_INFO <- getDesignInfo(D_TYPE = "approx", MODEL_INFO = MODEL_INFO, dist_func = DISTANCE, 
                          crit_type = crit_type, MaxMinStdVals = MaxMinStdVals,
                          dSupp = length(dsLower), nSupp = nSupp, dsLower = dsLower, dsUpper = dsUpper)

	# Adjust ALG_INFO according to D_INFO
	swarmSetting <- algInfoUpdate(D_INFO)
	ALG_INFO$varUpper <- swarmSetting$UB
	ALG_INFO$varLower <- swarmSetting$LB
	ALG_INFO$dSwarm	<- ncol(swarmSetting$UB)

	DESIGN_M <- designV2M(DESIGN, D_INFO)
	
	CRIT_VAL <- cppDesignCriterion(ALG_INFO, D_INFO, MODEL_LIST, 0, environment, DESIGN_M)
	PARA_SET <- CRIT_VAL$theta2

	ALPHA <- 0
	if (crit_type == "maxmin_fixed_true") {
		#message("Looking for best weight...")
		# Find the weight vector first
		ALPHA_INFO <- getDesignInfo(D_TYPE = "maxmin_eqv_wt", MODEL_INFO = MODEL_INFO, dist_func = DISTANCE, 
                          	 		crit_type = crit_type, MaxMinStdVals = MaxMinStdVals,
                             		dSupp = length(dsLower), nSupp = nSupp, dsLower = dsLower, dsUpper = dsUpper)
		ALPHA_INFO$paras <- PARA_SET
		ALPHA_ALG_INFO <- getAlgInfo(nSwarm = 32, maxIter = 100)
		swarmSetting <- algInfoUpdate(ALPHA_INFO)
		ALPHA_ALG_INFO$varUpper <- swarmSetting$UB
		ALPHA_ALG_INFO$varLower <- swarmSetting$LB
		ALPHA_ALG_INFO$dSwarm	<- ncol(swarmSetting$UB)
		FIXEDVALUE <- c(DESIGN_M[1:(nSupp*dSupp)], PSO_RESULT$BESTVAL)
		psoOut <- cppPSO(0, ALPHA_ALG_INFO, ALPHA_INFO, MODEL_LIST, FIXEDVALUE, environment, FALSE, FALSE)
		ALPHA <- designM2V(psoOut$GBest, ALPHA_INFO)
	}

	equiv <- cppEquivalence(ALG_INFO, D_INFO, MODEL_LIST, -CRIT_VAL$val, PARA_SET, ALPHA, environment, ngrid)
	
	if (crit_type == "maxmin_fixed_true") { equiv$alpha <- ALPHA }

	return(equiv)
}

#' Create An Empty Model List
#'
#' Run \code{emptyModelList()} for instruction of model settings.
#'
#' @param N_model Integer. The numebr of competing models including the true model. The default value is 2.
#' @return An List.
#' \itemize{
#' \item{True}{a list of informaiton of true model.}
#' \itemize{
#' \item{model}{the model formula. Input the R function or xPtr function pointer generated by \code{cxxfunction} in the \code{inline} package.}	
#' \item{para}{the nominal values of parameeters. Input a vector of parameter values.}	
#' }
#' \item{Rival\code{k}}{a vector of grid points on the second dimension of the 2-D design space.  For 1_d case, this output can be ignored.}
#' \itemize{
#' \item{model}{ the model formula. Input the R function or xPtr function pointer generated by \code{cxxfunction} in the \code{inline} package.}	
#' \item{paraLower}{ the lower bound of teh parameters in the k-th rival model. Input a vector of lower bound values.}	
#' \item{paraUpper}{ the upper bound of teh parameters in the k-th rival model. Input a vector of upper bound values.}	
#' \item{paraInit}{ the initial guess of parameters for the LBFGS algorithm. Input a vector of parameter values.}	
#' }
#' }
#' @examples
#' emptyModelList()  # True model and one rival model
#' emptyModelList(3) # True model and two rival models
#' emptyModelList(5) # True model and four rival models
#' @name emptyModelList
#' @rdname emptyModelList
#' @export
emptyModelList <- function(N_model = 2) {
	out <- lapply(1:N_model, function(k) {
		if (k == 1) list(model = 'R or C++ Function', para = 'Nominal Values of Parameters in True Model')
		else list(model = 'R or C++ Function', 
							paraLower = paste0('Lower Bound of Parameters in Rival', k-1), 
							paraUpper = paste0('Upper Bound of Parameters in Rival', k-1),
							paraInit = paste0('Initial Guess of Parameters for LBFGS Algorithm in Rival', k-1))
	})
	names(out) <- c("True", paste0("Rival", 1:(N_model-1)))
	return(out)
}
