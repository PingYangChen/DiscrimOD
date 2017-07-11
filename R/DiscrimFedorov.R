#' Hybrid Algorithm of PSO and L-BFGS for Finding Optimal Discrimination Design
#'
#' Implement the hybridized PSO algorithm to find the optimal discrimination designs
#' for 2 or more than 2 competing models
#'
#'
#' @references Atkinson, A. C. and Fedorov, V. V. (1975a). The design of experiments for discriminating between two rival models. Biometrika, 62(1):57-70.
#' @references Lopez-Fidalgo, J., Tommasi, C., and Trandafir, P. C. (2007). An optimal experimental design criterion for discriminating between non-normal models. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 69(2):231-242.
#' @name DiscrimFedWynn
#' @rdname DiscrimFedWynn
#' @export
DiscrimFedWynn <- function(MODEL_INFO, DISTANCE, dsLower, dsUpper, 
													 FED_INFO = NULL, LBFGS_INFO = NULL, seed = NULL, verbose = TRUE, environment, ...) {

  stopifnot(all(is.finite(dsLower)), all(is.finite(dsUpper)),
            length(dsLower) == length(dsUpper), all(dsUpper > dsLower),
            all(names(FED_INFO) == names(getFEDInfo())),
            all(names(LBFGS_INFO) == names(getLBFGSInfo())))

	MODEL_LIST <- lapply(1:length(MODEL_INFO), function(k) MODEL_INFO[[k]]$model)

  dSupp <- length(dsLower)

	D_INFO <- getDesignInfo(D_TYPE = "approx", MODEL_INFO = MODEL_INFO, dist_func = DISTANCE,
                          crit_type = "pair_fixed_true", MaxMinStdVals = 0,
                          dSupp = length(dsLower), nSupp = 2, dsLower = dsLower, dsUpper = dsUpper)

	D_INFO$nSupp <- D_INFO$dParas[2] + 1
	nSupp <- D_INFO$nSupp

	if (is.null(FED_INFO)) {
		FED_INFO <- getFEDInfo()
		if (verbose) message(paste0("Use the default settings for Fedorov-Wynn Algorithm. See '?getFEDInfo'."))
	}

	if (is.null(LBFGS_INFO)) {
		LBFGS_INFO <- getLBFGSInfo()
		if (verbose) message(paste0("Use the default settings for LBFGS. See '?getLBFGSInfo'."))
	} else {
		if (!LBFGS_INFO$IF_INNER_LBFGS) stop("L-BFGS is needed in Fedorov-Wynn algorithm.")
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
		fedOut <- cppFedorovWynn(FED_INFO, LBFGS_INFO, D_INFO, MODEL_LIST, environment, verbose)
	)[3]

	if (verbose) message(paste0("CPU time: ", round(cputime, 2), " seconds."))

	nonzero_wt <- which(fedOut$WT[1,] > 1e-12)
	BESTDESIGN <- cbind(fedOut$DESIGN[nonzero_wt,], fedOut$WT[nonzero_wt])
	BESTDESIGN <- BESTDESIGN[order(BESTDESIGN[,1]),]
	dimnames(BESTDESIGN) <- list(paste0("obs_", 1:nrow(BESTDESIGN)), 
															 c(paste0("dim_", 1:(ncol(BESTDESIGN) - 1)), "weight"))
	
	list(BESTDESIGN = BESTDESIGN, BESTVAL = fedOut$F_VAL, FVALHIST = fedOut$fvalHist,
			 R_PARA = fedOut$R_PARA, CPUTIME = cputime)
}
