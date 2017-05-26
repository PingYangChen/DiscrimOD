#' Hybrid Algorithm of PSO and BFGS for Finding Optimal Discrimination Design
#'
#' Please follow the instruction.
#'
#' @param MODEL_INFO list of information of competing models. For details, run \code{emptyModelList()} and see the instruction below.
#' @param DISTANCE function. 
#' @param nSupp integer. The number fo support points (at least 2).
#' @param dsLower vector. The finite lower bounds of the design space. Its length should be equal to the dimension of design space.
#' @param dsUpper vector. The finite upper bounds of the design space. Its length should be equal to the dimension of design space.
#' @param crit_type string. The name of the case of the discrimination design problem. The default is 'pair_fixed_true'.
#' @param MaxMinStdVals vector. The values of demoninators in the design efficiency calculation for finding max-min discrimination design.
#' @param ALG_INFO list. PSO and BFGS options.
#' @param seed numeric. random seed.  The default is \code{NULL}.
#' @param verbose logical. If \code{TRUE}, the PSO will reports the updating progress.
#' @return An List.
#' @examples
#' # Atkinson and Fedorov (1975)
#' # Two models are given by
#'
#' # Initialize PSO and BFGS options
#'
#' # Run Algorithm
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
#' @param DESIGN1 matrix.
#' @param MODEL_INFO list of information of competing models. For details, run \code{emptyModelList()} and see the instruction below.
#' @param DISTANCE function. 
#' @param dsLower vector. The finite lower bounds of the design space. Its length should be equal to the dimension of design space.
#' @param dsUpper vector. The finite upper bounds of the design space. Its length should be equal to the dimension of design space.
#' @param crit_type string. The name of the case of the discrimination design problem. The default is 'pair_fixed_true'.
#' @param MaxMinStdVals vector. The values of demoninators in the design efficiency calculation for finding max-min discrimination design.
#' @param ALG_INFO list. PSO and BFGS options.
#' @return An List.
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
#' @param N_model Integer. 
#' @return An List.
#' @name emptyModelList
#' @rdname emptyModelList
#' @export
emptyModelList <- function(N_model = 2) {
	out <- lapply(1:N_model, function(k) {
		if (k == 1) list(model = 'R/C++ Function', para = 'Nominal Values of Parameters in True Model')
		else list(model = 'R/C++ Function', 
							paraLower = paste0('Lower Bound of Parameters in Rival', k-1), 
							paraUpper = paste0('Upper Bound of Parameters in Rival', k-1))
	})
	names(out) <- c("True", paste0("Rival", 1:(N_model-1)))
	return(out)
}
