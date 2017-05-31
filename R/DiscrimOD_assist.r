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
#' algInfo <- getAlgInfo(nSwarm = 32, maxIter = 100)
#'
#' # Set the L-BFGS algorithm to repeat 4 times.
#' algInfo <- getAlgInfo(nSwarm = 32, maxIter = 100, LBFGS_RETRY = 4)
#' @name getAlgInfo
#' @rdname getAlgInfo
#' @export
getAlgInfo <- function(nSwarm = 32, maxIter = 100,
  #typePSO = NULL, #dSwarm = NULL, varUpper = NULL, varLower = NULL,
  checkConv = NULL, freeRun = NULL, tol = NULL, c1 = NULL, c2 = NULL,
  w0 = NULL, w1 = NULL, w_var = NULL, vk = NULL, #chi = NULL,
  #typeTopo = NULL, nGroup = NULL, GC_S_ROOF = NULL, GC_F_ROOF = NULL, GC_RHO = NULL,
  #Q_cen_type = NULL, Q_a0 = NULL, Q_a1 = NULL, Q_a_var = NULL, LcRi_L = NULL,
  LBFGS_RETRY = NULL, LBFGS_MAXIT = NULL, LBFGS_LM = NULL, FVAL_EPS = NULL, GRAD_EPS = NULL,
  LINESEARCH_MAXTRIAL = NULL, LINESEARCH_MAX = NULL, LINESEARCH_MIN = NULL, LINESEARCH_RHO = NULL,
  LINESEARCH_ARMIJO = NULL) {#, LINESEARCH_WOLFE = NULL) {

  nLoop <- length(nSwarm)
  #if(is.null(dSwarm))     dSwarm     <- NULL
  #if(is.null(varUpper))   varUpper   <- NULL
  #if(is.null(varLower))   varLower   <- NULL
  #if(is.null(maxIter))    maxIter    <- rep(100   , nLoop)
  if(is.null(checkConv))  checkConv  <- rep(0     , nLoop)
  #if(is.null(typePSO))    typePSO    <- rep(0     , nLoop)
  if(is.null(freeRun))    freeRun    <- rep(0.25  , nLoop)
  if(is.null(tol))        tol        <- rep(1e-6  , nLoop)
  if(is.null(c1))         c1         <- rep(2.05  , nLoop)
  if(is.null(c2))         c2         <- rep(2.05  , nLoop)
  if(is.null(w0))         w0         <- rep(0.95  , nLoop)
  if(is.null(w1))         w1         <- rep(0.2   , nLoop)
  if(is.null(w_var))      w_var      <- rep(0.8   , nLoop)
  #if(is.null(chi))        chi        <- rep(0.729 , nLoop)
  if(is.null(vk))         vk         <- rep(4     , nLoop)
  #if(is.null(typeTopo))   typeTopo   <- rep(0     , nLoop)
  #if(is.null(nGroup))     nGroup     <- rep(1     , nLoop)
  #if(is.null(GC_S_ROOF))  GC_S_ROOF  <- rep(5     , nLoop)
  #if(is.null(GC_F_ROOF))  GC_F_ROOF  <- rep(15    , nLoop)
  #if(is.null(GC_RHO))     GC_RHO     <- rep(1     , nLoop)
  #if(is.null(Q_cen_type)) Q_cen_type <- rep(1     , nLoop)
  #if(is.null(Q_a0))       Q_a0       <- rep(1.7   , nLoop)
  #if(is.null(Q_a1))       Q_a1       <- rep(0.7   , nLoop)
  #if(is.null(Q_a_var))    Q_a_var    <- rep(0.8   , nLoop)
  #if(is.null(LcRi_L))     LcRi_L     <- rep(0.01  , nLoop)

  if(is.null(LBFGS_RETRY))          LBFGS_RETRY         <- rep(3    , nLoop)
  if(is.null(LBFGS_MAXIT))          LBFGS_MAXIT         <- rep(100  , nLoop)
  if(is.null(LBFGS_LM))             LBFGS_LM            <- rep(6    , nLoop)
  if(is.null(FVAL_EPS))             FVAL_EPS            <- rep(0    , nLoop)
  if(is.null(GRAD_EPS))             GRAD_EPS            <- rep(1e-10, nLoop)
  if(is.null(LINESEARCH_MAXTRIAL))  LINESEARCH_MAXTRIAL <- rep(50   , nLoop)
  if(is.null(LINESEARCH_MAX))       LINESEARCH_MAX      <- rep(1    , nLoop)
  if(is.null(LINESEARCH_MIN))       LINESEARCH_MIN      <- rep(1e-6 , nLoop)
  if(is.null(LINESEARCH_RHO))       LINESEARCH_RHO      <- rep(0.618, nLoop)
  if(is.null(LINESEARCH_ARMIJO))    LINESEARCH_ARMIJO   <- rep(1e-4 , nLoop)
  #if(is.null(LINESEARCH_WOLFE))     LINESEARCH_WOLFE    <- rep(0.9 , nLoop)

  list(nSwarm = nSwarm, dSwarm = "autogen", varUpper = "autogen", varLower = "autogen", maxIter = maxIter, #typePSO = typePSO,
    checkConv = checkConv, freeRun = freeRun, tol = tol, c1 = c1, c2 = c2, w0 = w0, w1 = w1, w_var = w_var, #chi = chi,
    vk = vk, #typeTopo = typeTopo, nGroup = nGroup, GC_S_ROOF = GC_S_ROOF, GC_F_ROOF = GC_F_ROOF,
    #GC_RHO = GC_RHO, Q_cen_type = Q_cen_type, Q_a0 = Q_a0, Q_a1 = Q_a1, Q_a_var = Q_a_var,
    #LcRi_L = LcRi_L,
    LBFGS_RETRY = LBFGS_RETRY, LBFGS_MAXIT = LBFGS_MAXIT, LBFGS_LM = LBFGS_LM, FVAL_EPS = FVAL_EPS, GRAD_EPS = GRAD_EPS,
    LINESEARCH_MAXTRIAL = LINESEARCH_MAXTRIAL, LINESEARCH_MAX = LINESEARCH_MAX, LINESEARCH_MIN = LINESEARCH_MIN,
    LINESEARCH_RHO = LINESEARCH_RHO, LINESEARCH_ARMIJO = LINESEARCH_ARMIJO)#, LINESEARCH_WOLFE = LINESEARCH_WOLFE)
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
      parasInit[k,] <- c(MODEL_INFO[[k]]$paraInit, rep(0, max(dParas) - dParas[k]))
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
