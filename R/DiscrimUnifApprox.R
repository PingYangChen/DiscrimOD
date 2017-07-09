#' Hybrid Algorithm of PSO and L-BFGS for Finding Optimal Discrimination Design
#'
#' Implement the hybridized PSO algorithm to find the optimal discrimination designs
#' for 2 or more than 2 competing models
#'
#'
#' @references Atkinson, A. C. and Fedorov, V. V. (1975a). The design of experiments for discriminating between two rival models. Biometrika, 62(1):57-70.
#' @references Lopez-Fidalgo, J., Tommasi, C., and Trandafir, P. C. (2007). An optimal experimental design criterion for discriminating between non-normal models. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 69(2):231-242.
#' @name DiscrimUnifApproxT
#' @rdname DiscrimUnifApproxT
#' @export
DiscrimUnifApproxT <- function(MODEL_INFO, nSupp, dsLower, dsUpper, REMES_MAXIT = 50, REMES_FreeRun = 1.0, REMES_EPS = 1e-3, 
													 		LBFGS_INFO = NULL, seed = NULL, verbose = TRUE, environment, ...) {

  stopifnot(all(is.finite(dsLower)), all(is.finite(dsUpper)), length(dsLower) == 1, 
            length(dsLower) == length(dsUpper), all(dsUpper > dsLower),
            all(names(FED_INFO) == names(getFEDInfo())),
            all(names(LBFGS_INFO) == names(getLBFGSInfo())))

	MODEL_LIST <- lapply(1:length(MODEL_INFO), function(k) MODEL_INFO[[k]]$model)

  dSupp <- length(dsLower)
  DISTANCE <- function(xt, xr) (xt - xr)^2

	D_INFO <- getDesignInfo(D_TYPE = "approx", MODEL_INFO = MODEL_INFO, dist_func = DISTANCE,
                          crit_type = "pair_fixed_true", MaxMinStdVals = 0,
                          dSupp = length(dsLower), nSupp = nSupp, dsLower = dsLower, dsUpper = dsUpper)

	stopifnot(D_INFO$nSupp > D_INFO$dParas[2])
	
	REMES_INFO <- list(REMES_MAXIT = REMES_MAXIT, freeRun = REMES_FreeRun, REMES_EPS = REMES_EPS)
	
	if (is.null(LBFGS_INFO)) {
		LBFGS_INFO <- getLBFGSInfo()
		if (verbose) message(paste0("Use the default settings for LBFGS. See '?getLBFGSInfo'."))
	} else {
		if (!LBFGS_INFO$IF_INNER_LBFGS) stop("L-BFGS is needed in the Remes algorithm.")
	}

	if (!hasArg(environment)) environment <- new.env()

	PSO_INFO <- getPSOInfo(nSwarm = c(32, 32), maxIter = c(100, 100))
	# Adjust PSO_INFO according to D_INFO
	swarmSetting <- algInfoUpdate(D_INFO)
	PSO_INFO$varUpper <- matrix(swarmSetting$UB, length(PSO_INFO$nSwarm), ncol(swarmSetting$UB), byrow = TRUE)
	PSO_INFO$varLower <- matrix(swarmSetting$LB, length(PSO_INFO$nSwarm), ncol(swarmSetting$UB), byrow = TRUE)
	PSO_INFO$dSwarm <- rep(ncol(swarmSetting$UB), length(PSO_INFO$nSwarm))
	
	# Find Initial Guess for parameters of rival model
	UNIFDESIGN_M <- designV2M(cbind(seq(dsLower, dsUpper, length = nSupp), 1/nSupp), D_INFO)
	iniGuess <- cppDesignCriterion(PSO_INFO, LBFGS_INFO, D_INFO, MODEL_LIST, 0, environment, UNIFDESIGN_M)
	D_INFO$parasInit <- iniGuess$theta2

	# Start
	set.seed(seed)
	cputime <- system.time(
		remesOut <- cppUnifApprox(REMES_INFO, LBFGS_INFO, D_INFO, MODEL_LIST, environment, verbose)
	)[3]

	if (verbose) message(paste0("CPU time: ", round(cputime, 2), " seconds."))

	BESTDESIGN <- cbind(remesOut$DESIGN, t(remesOut$WT))
	dimnames(BESTDESIGN) <- list(paste0("obs_", 1:nrow(BESTDESIGN)), 
															 c(paste0("dim_", 1:(ncol(BESTDESIGN) - 1)), "weight"))

	BESTVAL <- "NOT A DESIGN"
	if (all(remesOut$WT >= 0) & all(remesOut$WT <= 1)) {
		BESTDESIGN_M <- designV2M(BESTDESIGN, D_INFO)
		tmp <- cppDesignCriterion(PSO_INFO, LBFGS_INFO, D_INFO, MODEL_LIST, 0, environment, BESTDESIGN_M)
		BESTVAL <- -tmp$val
	} 

	list(BESTDESIGN = BESTDESIGN, BESTVAL = BESTVAL, CPoly = remesOut$CPolyVal, CPUTIME = cputime)
}
