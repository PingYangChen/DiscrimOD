#' Generation function of PSO parameter settings
#'
#' Create a list with PSO parameters for optimal discrimination design search.
#'
#' @param nSwarm A integer number of swarm size in PSO algorithm.
#' @param maxIter A integer number of maximal PSO iterations.
#' @param checkConv A logical value which controls whether PSO checks the stopping criterion during updating procedure.
#' Specify \code{TRUE} for PSO to compute the stopping criterion \eqn{|f'-f|<\varepsilon}
#' where \eqn{f'} and \eqn{f} are the objective function values in the previous and current iterations, respectively.
#' The default is \code{FALSE}.
#' @param freeRun A number between \eqn{[0,1]} that controls the percentage of PSO iterations which are free from examining the 
#' stopping criterion, \eqn{|f'-f|<\varepsilon}
#' where \eqn{f'} and \eqn{f} are the objective function values in the previous and current iterations, respectively.
#' The default is 1.00 implying the PSO will completely ignore the stopping criterion. 
#' Otherwise, the PSO checks the stopping criterion after free iterations.
#' @param tol A small value for the tolerance, \eqn{\varepsilon}, in the stopping criterion.
#' If \code{checkConv = TRUE}, the default is \code{1e-6}. Otherwise, this value would not affect the algorithm.
# @param typePSO integer. The type of PSO. In this package, we have the following types:
# \itemize{
# \item{0}{ Linearly Decreasing Weight PSO (Shi, Y. H.	and Eberhart, R. C., 1998)}
# \item{1}{ GCPSO (van den Bergh, F. and	Engelbrecht, A. P., 2002)}
# \item{2}{ Quantum PSO (Sun, J., Feng, B. and Xu, W., 2004)}
# \item{3}{ LcRiPSO (Bonyadi, M. R., Michalewicz, Z., 2014)}
# }
#' @param c1 The value of cognitive parameter in PSO updating procedure. The default is 2.05.
#' @param c2 The value of social parameter in PSO updating procedure. The default is 2.05.
#' @param w0 The value of starting inertia weight in PSO updating procedure. The default is 1.2.
#' @param w1 The value of ending inertia weight in PSO updating procedure. The default is 0.2.
#' @param w_var A number between \eqn{[0,1]} that controls the percentage of iterations that PSO linearly decreases the inertia weight
#' from \code{w0} to \code{w1}. The default is 0.8.
#' @param vk The value of velocity clamping parameter. The default is 4.
#' @return A list of PSO parameter settings.
#' @examples
#' # Get default settings with specified swarm size and maximal number of iterations.
#' PSO_INFO <- getPSOInfo(nSwarm = 32, maxIter = 100)
#'
#' # If wanted to disable L-BFGS for the inner optimization loop and
#' # use NestedPSO algorithm (Chen et al., 2015), we need the options
#' # for the two-layer PSO: c(outer loop option, inner loop option)
#' NESTEDPSO_INFO <- getPSOInfo(nSwarm = c(16, 32), maxIter = c(100, 200))
#' # Also, disable the L-BFGS algorithm
#' LBFGS_NOTRUN <- getLBFGSInfo(IF_INNER_LBFGS = FALSE)
#'
#' @references Chen, R.-B., Chang, S.-P., Wang, W., Tung, H.-C., and Wong, W. K. (2015). Minimax optimal designs via particle swarm optimization methods. Statistics and Computing, 25(5):975-988.
#' @name getPSOInfo
#' @rdname getPSOInfo
#' @export
getPSOInfo <- function(nSwarm = 32, maxIter = 100,
  #typePSO = 0, #dSwarm = NULL, varUpper = NULL, varLower = NULL, checkConv = 0, 
  freeRun = 1.0, tol = 1e-6, c1 = 2.05, c2 = 2.05,
  w0 = 1.2, w1 = 0.2, w_var = 0.8, vk = 4 #, chi = NULL,
  #typeTopo = NULL, nGroup = NULL, GC_S_ROOF = 5, GC_F_ROOF = 15, GC_RHO = 1,
  #Q_cen_type = 1, Q_a0 = 1.7, Q_a1 = 0.7, Q_a_var = 0.8, LcRi_L = 0.01,
  ) {

  stopifnot(length(nSwarm) == length(maxIter))

  nLoop <- length(nSwarm)
  #if (length(dSwarm))     dSwarm     <- numeric(nLoop)
  #if (length(varUpper))   varUpper   <- matrix(0, nLoop)
  #if (length(varLower))   varLower   <- matrix(0, nLoop)
  #if (length(maxIter))    maxIter    <- rep(100   , nLoop)
  #if (length(checkConv) < nLoop)  checkConv  <- rep(checkConv, nLoop)
  #if (length(typePSO)) < nLoop)    typePSO    <- rep(0     , nLoop)
  if (length(freeRun) < nLoop)    freeRun    <- rep(freeRun, nLoop)
  if (length(tol) < nLoop)        tol        <- rep(tol, nLoop)
  if (length(c1) < nLoop)         c1         <- rep(c1, nLoop)
  if (length(c2) < nLoop)         c2         <- rep(c2, nLoop)
  if (length(w0) < nLoop)         w0         <- rep(w0, nLoop)
  if (length(w1) < nLoop)         w1         <- rep(w1, nLoop)
  if (length(w_var) < nLoop)      w_var      <- rep(w_var, nLoop)
  #if (length(chi) < nLoop)        chi        <- rep(0.729 , nLoop)
  if (length(vk) < nLoop)         vk         <- rep(vk, nLoop)
  #if (length(typeTopo) < nLoop)   typeTopo   <- rep(0     , nLoop)
  #if (length(nGroup) < nLoop)     nGroup     <- rep(1     , nLoop)
  #if (length(GC_S_ROOF) < nLoop)  GC_S_ROOF  <- rep(5     , nLoop)
  #if (length(GC_F_ROOF) < nLoop)  GC_F_ROOF  <- rep(15    , nLoop)
  #if (length(GC_RHO) < nLoop)     GC_RHO     <- rep(1     , nLoop)
  #if (length(Q_cen_type) < nLoop) Q_cen_type <- rep(1     , nLoop)
  #if (length(Q_a0) < nLoop)       Q_a0       <- rep(1.7   , nLoop)
  #if (length(Q_a1) < nLoop)       Q_a1       <- rep(0.7   , nLoop)
  #if (length(Q_a_var) < nLoop)    Q_a_var    <- rep(0.8   , nLoop)
  #if (length(LcRi_L) < nLoop)     LcRi_L     <- rep(0.01  , nLoop)

  list(nSwarm = nSwarm, dSwarm = "autogen", varUpper = "autogen", varLower = "autogen", maxIter = maxIter, #typePSO = typePSO, checkConv = checkConv, 
    freeRun = freeRun, tol = tol, c1 = c1, c2 = c2, w0 = w0, w1 = w1, w_var = w_var, #chi = chi,
    vk = vk #, #typeTopo = typeTopo, nGroup = nGroup, GC_S_ROOF = GC_S_ROOF, GC_F_ROOF = GC_F_ROOF,
    #GC_RHO = GC_RHO, Q_cen_type = Q_cen_type, Q_a0 = Q_a0, Q_a1 = Q_a1, Q_a_var = Q_a_var,
    #LcRi_L = LcRi_L,
    )
}

#' Generation function of L-BFGS parameter settings
#'
#' Create a list with L-BFGS parameters for optimal discrimination design search.
#'
#' @param IF_INNER_LBFGS The logical input \code{TRUE}/\code{FALSE} to turn on/off L-BFGS algorithm for the
#' inner optimization problem (minimizing the distance among parameter space).
#' If \code{IF_INNER_LBFGS = FALSE}, then the NestedPSO algorithm in Chen et al. (2015) is used for optimal discrimination design search.
#' @param LBFGS_RETRY The integer number of total trials of L-BFGS in computing the minimal distance
#' with randomly generated initial values. The default is 1.
#' @param LBFGS_MAXIT The integer number of maximal iteration of L-BFGS algorithm.
#' The default is 0 and L-BFGS stops when it converges.
#' @param LBFGS_LM The integer number of corrections to approximate the inverse hessian matrix.
#' The default is 6.
#' @param FVAL_EPS The tolerance value, \eqn{\varepsilon_f}, of the stopping criterion on objective function value,
#' \eqn{|f'-f|/f<\varepsilon_f} where \eqn{f'} and \eqn{f} are the objective function values
#' in the previous and current iterations, respectively.
#' The default is 0 and L-BFGS ignores this stopping criterion.
#' @param GRAD_EPS The tolerance value, \eqn{\varepsilon_g}, of the stopping criterion on gradients,
#' \eqn{\|g\|/\|x\|<\varepsilon_g} where \eqn{g} is the gradient, \eqn{x} is the current position
#' and \eqn{\|a\|=\sqrt{a^\top a}}. The default is \code{1e-5}.
#' @param LINESEARCH_MAXTRIAL The integer number of maximal trial of More-Thuente line search routine. The default is 20.
#' @param LINESEARCH_MAX The maximal step size in the line search routine. The default is 1e20.
#' @param LINESEARCH_MIN The minimal step size in the line search routine. The default is 1e-20.
#' @param LINESEARCH_ARMIJO A parameter to control the accuracy of the line search routine. The default value
#' is \code{1e-4}. This parameter should be greater than zero and smaller than 0.5.
#' @param LINESEARCH_WOLFE A coefficient for the Wolfe condition in the line search routine. The default
#' value is 0.9. This parameter should be greater than \code{LINESEARCH_ARMIJO} parameter and smaller than 1.0.
#' @param FD_DELTA A small value for the gap of finite difference method for computing the gradient. The default is \code{1e-3}.
#' @return The list of L-BFGS parameter settings.
#' @examples
#' # Get default settings with 2 repeatedly trails for L-BFGS algorithm.
#' LBFGS_INFO <- getLBFGSInfo(LBFGS_RETRY = 2)
#'
#' # If wanted to disable L-BFGS for the inner optimization loop and
#' # use NestedPSO algorithm (Chen et al., 2015), we need the options
#' # for the two-layer PSO: c(outer loop option, inner loop option)
#' NESTEDPSO_INFO <- getPSOInfo(nSwarm = c(16, 32), maxIter = c(100, 200))
#' # Also, disable the L-BFGS algorithm
#' LBFGS_NOTRUN <- getLBFGSInfo(IF_INNER_LBFGS = FALSE)
#'
#' @name getLBFGSInfo
#' @rdname getLBFGSInfo
#' @export
getLBFGSInfo <- function(IF_INNER_LBFGS = TRUE, LBFGS_RETRY = 1, LBFGS_MAXIT = 0, LBFGS_LM = 6, FVAL_EPS = 0, GRAD_EPS = 1e-5,
                         LINESEARCH_MAXTRIAL = 20, LINESEARCH_MAX = 1e20, LINESEARCH_MIN = 1e-20,
                         LINESEARCH_ARMIJO = 1e-4, LINESEARCH_WOLFE = 0.9, FD_DELTA = 1e-3) {

  list(IF_INNER_LBFGS = as.integer(IF_INNER_LBFGS), LBFGS_RETRY = LBFGS_RETRY, LBFGS_MAXIT = LBFGS_MAXIT, LBFGS_LM = LBFGS_LM,
       FVAL_EPS = FVAL_EPS, GRAD_EPS = GRAD_EPS,
       LINESEARCH_MAXTRIAL = LINESEARCH_MAXTRIAL, LINESEARCH_MAX = LINESEARCH_MAX, LINESEARCH_MIN = LINESEARCH_MIN,
       LINESEARCH_ARMIJO = LINESEARCH_ARMIJO, LINESEARCH_WOLFE = LINESEARCH_WOLFE, FD_DELTA = FD_DELTA)
}

#' Generation function of Fedorov-Wynn algorithm parameter settings
#'
#' Create a list with Fedorov-Wynn algorithm parameters for optimal discrimination design search.
#'
#' @param freeRun A number between \eqn{[0,1]} that controls the percentage of updating iterations which are free from examining the 
#' stopping criterion, \eqn{|f'-f|<\varepsilon}
#' where \eqn{f'} and \eqn{f} are the objective function values in the previous and current iterations, respectively.
#' The default is 1.00 implying the algorithm will completely ignore the stopping criterion. 
#' Otherwise, the algorithm checks the stopping criterion after free iterations.
#' @return The list of L-BFGS parameter settings.
#' @examples
#' # Get default settings for Fedorov-Wynn algorithm.
#' FED_INFO <- getFEDInfo(FED_MAXIT = 200)
#'
#' @name getFEDInfo 
#' @rdname getFEDInfo
#' @export
getFEDInfo <- function(FED_MAXIT = 200, FED_TRIM = 5, FED_TRIM_EPS = 1e-4, freeRun = 0.5, FED_EPS = 1e-6, FED_ALPHA_GRID = 20) {
   
  list(FED_MAXIT = FED_MAXIT, FED_TRIM = FED_TRIM, FED_TRIM_EPS = FED_TRIM_EPS, 
       freeRun = freeRun, FED_EPS = FED_EPS, FED_ALPHA_GRID = FED_ALPHA_GRID)
}

#' @export
getDesignInfo <- function(D_TYPE = "approx", MODEL_INFO = NULL, dist_func = NULL,
                          crit_type = "pair_fixed_true", MaxMinStdVals = NULL,
                          dSupp = 1L, nSupp = 2L, dsLower = NULL, dsUpper = NULL) {

  dParas <- sapply(1:length(MODEL_INFO), function(k) if (k == 1) length(MODEL_INFO[[k]]$para) else length(MODEL_INFO[[k]]$paraUpper))

  paras <- parasInit <- parasUpper <- parasLower <- parasBdd <- matrix(0, length(MODEL_INFO), max(dParas))
  for (k in 1:length(MODEL_INFO)) {
    if (k == 1) { paras[k,] <- c(MODEL_INFO[[k]]$para, rep(0, max(dParas) - dParas[k])) } else {

      parasUpper[k,] <- c(ifelse(is.finite(MODEL_INFO[[k]]$paraUpper), MODEL_INFO[[k]]$paraUpper, 0), rep(0, max(dParas) - dParas[k]))
      parasLower[k,] <- c(ifelse(is.finite(MODEL_INFO[[k]]$paraLower), MODEL_INFO[[k]]$paraLower, 0), rep(0, max(dParas) - dParas[k]))

      tmp <- is.finite(MODEL_INFO[[k]]$paraLower) + 10*is.finite(MODEL_INFO[[k]]$paraUpper)
      tmp2 <- ifelse(tmp == 0, 0, ifelse(tmp == 1, 1, ifelse(tmp == 10, 3, 2)))
      parasBdd[k,] <- c(tmp2, rep(0, max(dParas) - dParas[k]))
      parasInit[k,] <- runif(max(dParas), as.vector(parasLower[k,]), as.vector(parasUpper[k,]))
    }
  }

  CRIT_TYPE_NUM <- ifelse(crit_type == "pair_fixed_true", 0,
                      ifelse(crit_type == "maxmin_fixed_true", 1, 2))

  if (D_TYPE == "maxmin_eqv_wt") { D_TYPE_NUM <- 1001 } else { D_TYPE_NUM <- 1 }

  N_PAIR <- length(MODEL_INFO) - 1
  MODEL_PAIR <- cbind(0, 1:N_PAIR)

  return(list(D_TYPE = D_TYPE, D_TYPE_NUM = D_TYPE_NUM, dist_func = dist_func,
              CRIT_TYPE_NUM = CRIT_TYPE_NUM,
              dSupp = dSupp, nSupp = nSupp, dsLower = dsLower, dsUpper = dsUpper,
              N_PAIR = N_PAIR, MODEL_PAIR = MODEL_PAIR, dParas = dParas, paras = paras, parasInit = parasInit,
              parasUpper = parasUpper, parasLower = parasLower, parasBdd = parasBdd,
              MaxMinStdVals = MaxMinStdVals))
}
