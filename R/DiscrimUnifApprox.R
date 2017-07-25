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
#' @importFrom limSolve lsei
DiscrimUnifApproxT <- function(MODEL_INFO, nSupp, dsLower, dsUpper, REMES_MAXIT = 50, REMES_FreeRun = 1.0, REMES_EPS = 1e-3, 
													 		LBFGS_INFO = NULL, seed = NULL, verbose = TRUE, environment, ...) {

  stopifnot(all(is.finite(dsLower)), all(is.finite(dsUpper)), length(dsLower) == 1, 
            length(dsLower) == length(dsUpper), all(dsUpper > dsLower),
            all(names(LBFGS_INFO) == names(getLBFGSInfo())),
            nSupp > 1)

	MODEL_LIST <- lapply(1:length(MODEL_INFO), function(k) MODEL_INFO[[k]]$model)

  dSupp <- length(dsLower)
  DISTANCE <- function(xt, xr) (xt - xr)^2

	D_INFO <- getDesignInfo(D_TYPE = "approx", MODEL_INFO = MODEL_INFO, dist_func = DISTANCE,
                          crit_type = "pair_fixed_true", MaxMinStdVals = 0,
                          dSupp = length(dsLower), nSupp = nSupp, dsLower = dsLower, dsUpper = dsUpper)
	
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

	wout <- tryCatch(
		{
			limSolve::lsei(A = remesOut$DD_DEV, B = rep(0, nrow(remesOut$DD_DEV)), 
										 E = t(rep(1, ncol(remesOut$DD_DEV))), F = 1, 
										 G = rbind(diag(ncol(remesOut$DD_DEV)), -diag(ncol(remesOut$DD_DEV))), 
										 H = c(rep(0, ncol(remesOut$DD_DEV)), rep(-1, ncol(remesOut$DD_DEV))))
		},
		error = function(cond) {
			message("weight cannot be found by limSolve::lsei")
			message("orignal error message:")
			message(cond)
			return(-1)
		},
		warning = function(cond) {
			message("these is a warning message from limSolve::lsei:")
			message(cond)
			return(-1)
		}
	)

	if (!is.list(wout)) {
		OUTPUT <- "CANNOT FIND A DESIGN DUE TO NO SOLUTION IN limSolve::lsei"	
	} else {
		BESTDESIGN <- cbind(remesOut$DESIGN, wout$X)
		dimnames(BESTDESIGN) <- list(paste0("obs_", 1:nrow(BESTDESIGN)), 
																 c(paste0("dim_", 1:(ncol(BESTDESIGN) - 1)), "weight"))
		BESTDESIGN_M <- designV2M(BESTDESIGN, D_INFO)
		tmp <- cppDesignCriterion(PSO_INFO, LBFGS_INFO, D_INFO, MODEL_LIST, 0, environment, BESTDESIGN_M)
		BESTVAL <- -tmp$val
		OUTPUT <- list(BESTDESIGN = BESTDESIGN, BESTVAL = BESTVAL, CPoly = remesOut$CPolyVal, CPUTIME = cputime)
	}
	OUTPUT
}
