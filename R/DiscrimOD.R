#' PSO for Optimal Discrimination Design
#'
#' Please follow the instruction.
#'
#' @param D_INFO list. Optimal design problem.
#' @param ALG_INFO list. PSO options.
#' @param seed numeric. random seed.  The default is \code{NULL}.
#' @param verbose logical. If \code{TRUE}, the PSO will reports the updating progress.
#' @return An List.
#' @examples
#' # Michaelis-Menten (MM) model with i.i.d. observations
# the parameter vector only consists of parameters in the mean function,
# that is, tmPara = c(v, km).
#' mm_mean <- function(N, x, p) {
#'   f <- matrix(0, N, 1)
#'   for (i in 1:N) f[i,1] <- p[1]*x[i,1]/(p[2]+x[i,1])
#'   return(f)
#' }
#' iid_cor <- function(N, x, p) {
#'   # in fact, for this case, diag(N) would be OK
#'   s <- matrix(0, N, N)
#'   for (i in 1:N) s[i,i] <- 1.0
#'   return(s)
#' }
#' # This MM model has 1 covariate (x > 0) and 2 parameters (v > 0, km > 0).
#' # Thus, 'dSupp = 1' and the lengthes of 'dsLower' and 'dsUpper' should be 1.
#' # Consider a 2-point D-optimal design, the setting is
#' D_INFO <- getObjInfo(CRIT_TYPE = "D", nSupp = 2,
#'                      dSupp = 1, dsLower = 0, dsUpper = 1,
#'                      tmPara = c(1., 2.), tmPOI = c(1, 2),
#'                      tmParaUpper = c(Inf, Inf), tmParaLower = c(0, 0),
#'                      mean_func = mm_mean, corr_func = iid_cor)
#'
#' # Initialize PSO options
#' ALG_INFO <- getAlgInfo(nSwarm = 64, maxIter = 50, typePSO = 0)
#'
#' # Run PSO
#' RESULT <- rOptimalDesignPSO(D_INFO, ALG_INFO, verbose = TRUE)
#' # Check the 2-point D-optimal design
#' RESULT$GBEST
#'
#' # Michaelis-Menten model with autocorrelated observations
#' # Now, the parameter vector consists of parameters in the mean function
#' # and the correlation function.
#' # That is, tmPara = c(v, km, lambda, sigmaSq),
#' # where the last two parameters are correlation parameters.
#' atuo_cor <- function(N, x, p) {
#'   s <- matrix(0, N, N)
#'   for (i in 1:N) {
#'     for (j in 1:i) {
#'       s[i,j] <- p[4]*(p[3]^abs(x[i,1]-x[j,1]))
#'       if (j < i) s[j,i] <- s[i,j]
#'     }
#'   }
#'   return(s)
#' }
#' # This MM model has 1 covariate (x > 0) and 4 parameters,
#' # (0 < v, 0 < km, , 0 < lambda < 1, 0 < sigmaSq).
#' # For covariate, 'dSupp = 1' and the lengthes of 'dsLower' and 'dsUpper' should be 1.
#' # Consider a 3-point D-optimal design.
#' # Let nominal values tmPara = c(1., 1.2, 0.5, 1.0).
#' # Suppose the parameters of interest are v and km, the setting is
#' D_INFO <- getObjInfo(CRIT_TYPE = "D", nSupp = 3,
#'                      dSupp = 1, dsLower = 0, dsUpper = 1,
#'                      tmPara = c(1., 1.2, 0.5, 1.), tmPOI = c(1, 2),
#'                      tmParaUpper = c(Inf, Inf, 1., Inf), tmParaLower = c(0, 0, 0, 0),
#'                      mean_func = mm_mean, corr_func = atuo_cor)
#' # Run PSO
#' RESULT <- rOptimalDesignPSO(D_INFO, ALG_INFO, verbose = TRUE)
#' RESULT$GBEST
#' @name DiscrimOD
#' @rdname DiscrimOD
#' @export
DiscrimOD <- function(MODEL_INFO, DISTANCE, nSupp, dsLower, dsUpper, MaxMinStdVals = NULL, ALG_INFO = NULL, seed = NULL, verbose = TRUE, ...) {

  assert_that(nSupp >= 2L, all(is.finite(dsLower)), all(is.finite(dsUpper)),
              length(dsLower) == length(dsUpper), all(dsUpper > dsLower))

	MODEL_LIST <- lapply(1:length(MODEL_INFO), function(k) MODEL_INFO[[k]]$model)

	if (is.null(MaxMinStdVals)) MaxMinStdVals <- 0
	D_INFO <- getDesignInfo(D_TYPE = "approx", MODEL_INFO = MODEL_INFO, dist_func = DISTANCE, 
                          crit_type = 0, MaxMinStdVals = MaxMinStdVals,
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
#' @param DESIGN1 string. "exact"
#' @param DESIGN2 integer. Number of support points.
#' @param CRIT_TYPE string.
#' @param D_INFO list. Optimal design problem.
#' @param ALG_INFO list. PSO options.
#' @return An List.
#' @name designCriterion
#' @rdname designCriterion
#' @export
designCriterion <- function(DESIGN1, MODEL_INFO, DISTANCE, dsLower, dsUpper, MaxMinStdVals = NULL, ALG_INFO = NULL, ...) {

	assert_that(all(is.finite(dsLower)), all(is.finite(dsUpper)),
              length(dsLower) == length(dsUpper), all(dsUpper > dsLower))

	nSupp <- nrow(DESIGN1)
	MODEL_LIST <- lapply(1:length(MODEL_INFO), function(k) MODEL_INFO[[k]]$model)

	if (is.null(MaxMinStdVals)) MaxMinStdVals <- 0
	D_INFO <- getDesignInfo(D_TYPE = "approx", MODEL_INFO = MODEL_INFO, dist_func = DISTANCE, 
                          crit_type = 0, MaxMinStdVals = MaxMinStdVals,
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

	DESIGN1_M <- designV2M(DESIGN1, D_INFO)

	cri_1 <- cppDesignCriterion(ALG_INFO, D_INFO, MODEL_LIST, 0, environment, DESIGN1_M)

  return(list(cri_val = -cri_1$val, theta2 = cri_1$theta2[-1,]))
}

#' Equivalence Theorem
#' @name equivalence
#' @rdname equivalence
#' @export
equivalence <- function(DESIGN = NULL, PSO_RESULT = NULL, ngrid = 100, IFPLOT = FALSE,
												MODEL_INFO, DISTANCE, dsLower, dsUpper, MaxMinStdVals = NULL, ALG_INFO = NULL, ...) {

	if (is.null(DESIGN)) DESIGN <- PSO_RESULT$BESTDESIGN

	assert_that(all(is.finite(dsLower)), all(is.finite(dsUpper)),
              length(dsLower) == length(dsUpper), all(dsUpper > dsLower))

	nSupp <- nrow(DESIGN)
	MODEL_LIST <- lapply(1:length(MODEL_INFO), function(k) MODEL_INFO[[k]]$model)

	if (is.null(MaxMinStdVals)) MaxMinStdVals <- 0
	D_INFO <- getDesignInfo(D_TYPE = "approx", MODEL_INFO = MODEL_INFO, dist_func = DISTANCE, 
                          crit_type = 0, MaxMinStdVals = MaxMinStdVals,
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

	DESIGN_M <- designV2M(DESIGN, D_INFO)
	
	CRIT_VAL <- cppDesignCriterion(ALG_INFO, D_INFO, MODEL_LIST, 0, environment, DESIGN_M)
	
	PARA_SET <- D_INFO$parasInit
	PARA_SET[2,] <- CRIT_VAL$theta2[2,]
	
	equiv <- cppEquivalence(ALG_INFO, D_INFO, MODEL_LIST, -CRIT_VAL$val, PARA_SET, environment, ngrid)

	return(equiv)
}

#' Create An Empty Model List
#'
#' Please follow the instruction.
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
