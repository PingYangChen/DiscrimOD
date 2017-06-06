#' Initialize PSO options
#'
#' Please follow the instruction.
#'
#' @param nSwarm integer. Swarm size.
#' @param maxIter integer. Maximal number of iterations.
#' @param checkConv logical. Specify \code{TRUE} if one wants PSO to stop when
#' satisfyning the stopping criterion \eqn{|f'-f|<\varepsilon}
#' where \eqn{f'} and \eqn{f} are the objective function values in the previous and current iterations, respectively.
#' The default is \code{FALSE}.
#' @param freeRun numeric. The persentage of iterations which are free from examing the stopping criterion.
#' The default is 0.25.
#' @param tol numeric. The value of \eqn{\varepsilon} in the stopping criterion, \eqn{|f'-f|<\varepsilon}.
#' The default is \code{1e-6}.
# @param typePSO integer. The type of PSO. In this package, we have the following types:
# \itemize{
# \item{0}{ Linearly Decreasing Weight PSO (Shi, Y. H.	and Eberhart, R. C., 1998)}
# \item{1}{ GCPSO (van den Bergh, F. and	Engelbrecht, A. P., 2002)}
# \item{2}{ Quantum PSO (Sun, J., Feng, B. and Xu, W., 2004)}
# \item{3}{ LcRiPSO (Bonyadi, M. R., Michalewicz, Z., 2014)}
# }
#' @param c1 numeric. The cognitive parameter in the PSO update procedure. The default is 2.05.
#' @param c2 numeric. The social parameter in the PSO update procedure. The default is 2.05.
#' @param w0 numeric. The starting inertia weight in the PSO update procedure. The default is 0.95.
#' @param w1 numeric. The ending inertia weight in the PSO update procedure. The default is 0.2.
#' @param w_var numeric. The persentage of iterations which linearly decrease the inertia weight
#' from \code{w0} to \code{w1}. The default is 0.8.
#' @param vk numeric. The velocity clamping parameter. The default is 4.
#' @param LBFGS_RETRY integer. The total times of trials in calculating the minimal distance
#' with randomly generated initial starters. The default is 3.
#' @param LBFGS_MAXIT integer. The maximal number of iterations of the L-BFGS algorithm.
#' The default is 100.
#' @param LBFGS_LM integer. The number of corrections to approximate the inverse hessian matrix.
#' The default is 6.
#' @param FVAL_EPS numeric. The value of \eqn{\varepsilon_f} in the stopping criterion,
#' \eqn{|f'-f|/f<\varepsilon_f} where \eqn{f'} and \eqn{f} are the objective function values
#' in the previous and current iterations, respectively.
#' The default is 0 which means not using this stopping criteiron.
#' @param GRAD_EPS numeric. The value of \eqn{\varepsilon_g} in the stopping criterion,
#' \eqn{\|g\|/\|x\|<\varepsilon_g} where \eqn{g} is the gradient vector, \eqn{x} is the current position
#' and \eqn{\|a\|=\sqrt{a^\top a}}. The default is 1e-5.
#' @param LINESEARCH_MAXTRIAL integer The maximal trial of the backtrackign line search. The default is 50.
#' @param LINESEARCH_MAX numeric. The initial step size in the line search. The default is 1.
#' @param LINESEARCH_MIN numeric. The minimal step size in the line search. The default is 1e-9.
#' @param LINESEARCH_RHO numeric. The persentage of the decreasement on the step size in each
#' iteration of the line search. The default is 0.618.
#' @param LINESEARCH_ARMIJO numeric. The parameter to
#' control the accuracy of the backtracking line search algorithm.
#' The default is \code{1e-4}. The value should be positive and smaller than 0.5.
#' @return The list of algorithm setting parameters as the input of \code{DiscrimOD}.
#' @examples
#' # Get default settings with specified swarm size and maximal number of iterations.
#' algInfo <- getPSOInfo(nSwarm = 32, maxIter = 100)
#'
#' @name getPSOInfo
#' @rdname getPSOInfo
#' @export
getPSOInfo <- function(nSwarm = 32, maxIter = 100,
  #typePSO = 0, #dSwarm = NULL, varUpper = NULL, varLower = NULL,
  checkConv = 0, freeRun = 0.25, tol = 1e-6, c1 = 2.05, c2 = 2.05,
  w0 = 0.95, w1 = 0.2, w_var = 0.8, vk = 4 #, chi = NULL,
  #typeTopo = NULL, nGroup = NULL, GC_S_ROOF = 5, GC_F_ROOF = 15, GC_RHO = 1,
  #Q_cen_type = 1, Q_a0 = 1.7, Q_a1 = 0.7, Q_a_var = 0.8, LcRi_L = 0.01,
  ) {

  stopifnot(length(nSwarm) == length(maxIter))

  nLoop <- length(nSwarm)
  #if (length(dSwarm))     dSwarm     <- numeric(nLoop)
  #if (length(varUpper))   varUpper   <- matrix(0, nLoop)
  #if (length(varLower))   varLower   <- matrix(0, nLoop)
  #if (length(maxIter))    maxIter    <- rep(100   , nLoop)
  if (length(checkConv) < nLoop)  checkConv  <- rep(0     , nLoop)
  #if (length(typePSO)) < nLoop)    typePSO    <- rep(0     , nLoop)
  if (length(freeRun) < nLoop)    freeRun    <- rep(0.25  , nLoop)
  if (length(tol) < nLoop)        tol        <- rep(1e-6  , nLoop)
  if (length(c1) < nLoop)         c1         <- rep(2.05  , nLoop)
  if (length(c2) < nLoop)         c2         <- rep(2.05  , nLoop)
  if (length(w0) < nLoop)         w0         <- rep(0.95  , nLoop)
  if (length(w1) < nLoop)         w1         <- rep(0.2   , nLoop)
  if (length(w_var) < nLoop)      w_var      <- rep(0.8   , nLoop)
  #if (length(chi) < nLoop)        chi        <- rep(0.729 , nLoop)
  if (length(vk) < nLoop)         vk         <- rep(4     , nLoop)
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

  list(nSwarm = nSwarm, dSwarm = "autogen", varUpper = "autogen", varLower = "autogen", maxIter = maxIter, #typePSO = typePSO,
    checkConv = checkConv, freeRun = freeRun, tol = tol, c1 = c1, c2 = c2, w0 = w0, w1 = w1, w_var = w_var, #chi = chi,
    vk = vk #, #typeTopo = typeTopo, nGroup = nGroup, GC_S_ROOF = GC_S_ROOF, GC_F_ROOF = GC_F_ROOF,
    #GC_RHO = GC_RHO, Q_cen_type = Q_cen_type, Q_a0 = Q_a0, Q_a1 = Q_a1, Q_a_var = Q_a_var,
    #LcRi_L = LcRi_L,
    )
}

#' Initialize LBFGS options
#'
#' Please follow the instruction.
#'
#' @param IF_INNER_LBFGS logical. 
#' @param LBFGS_RETRY integer. The total times of trials in calculating the minimal distance
#' with randomly generated initial starters. The default is 3.
#' @param LBFGS_MAXIT integer. The maximal number of iterations of the L-BFGS algorithm.
#' The default is 100.
#' @param LBFGS_LM integer. The number of corrections to approximate the inverse hessian matrix.
#' The default is 6.
#' @param FVAL_EPS numeric. The value of \eqn{\varepsilon_f} in the stopping criterion,
#' \eqn{|f'-f|/f<\varepsilon_f} where \eqn{f'} and \eqn{f} are the objective function values
#' in the previous and current iterations, respectively.
#' The default is 0 which means not using this stopping criteiron.
#' @param GRAD_EPS numeric. The value of \eqn{\varepsilon_g} in the stopping criterion,
#' \eqn{\|g\|/\|x\|<\varepsilon_g} where \eqn{g} is the gradient vector, \eqn{x} is the current position
#' and \eqn{\|a\|=\sqrt{a^\top a}}. The default is 1e-5.
#' @param LINESEARCH_MAXTRIAL integer The maximal trial of the backtrackign line search. The default is 50.
#' @param LINESEARCH_MAX numeric. The initial step size in the line search. The default is 1.
#' @param LINESEARCH_MIN numeric. The minimal step size in the line search. The default is 1e-9.
#' @param LINESEARCH_ARMIJO numeric. The parameter to
#' control the accuracy of the backtracking line search algorithm.
#' The default is \code{1e-4}. The value should be positive and smaller than 0.5.
#' @param LINESEARCH_WOLFE numeric. The persentage of the decreasement on the step size in each
#' iteration of the line search. The default is 0.618.
#' @return The list of algorithm setting parameters as the input of \code{DiscrimOD}.
#' @examples
#' lbfgsInfo <- getLBFGSInfo(LBFGS_RETRY = 2)
#'
#' @name getLBFGSInfo
#' @rdname getLBFGSInfo
#' @export
getLBFGSInfo <- function(IF_INNER_LBFGS = TRUE, LBFGS_RETRY = 1, LBFGS_MAXIT = 0, LBFGS_LM = 6, FVAL_EPS = 0, GRAD_EPS = 1e-5,
                         LINESEARCH_MAXTRIAL = 20, LINESEARCH_MAX = 1e20, LINESEARCH_MIN = 1e-20, 
                         LINESEARCH_ARMIJO = 1e-4, LINESEARCH_WOLFE = 0.9) {

  list(IF_INNER_LBFGS = as.integer(IF_INNER_LBFGS), LBFGS_RETRY = LBFGS_RETRY, LBFGS_MAXIT = LBFGS_MAXIT, LBFGS_LM = LBFGS_LM, 
       FVAL_EPS = FVAL_EPS, GRAD_EPS = GRAD_EPS,
       LINESEARCH_MAXTRIAL = LINESEARCH_MAXTRIAL, LINESEARCH_MAX = LINESEARCH_MAX, LINESEARCH_MIN = LINESEARCH_MIN,
       LINESEARCH_ARMIJO = LINESEARCH_ARMIJO, LINESEARCH_WOLFE = LINESEARCH_WOLFE)
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
