#' Initialize PSO options
#'
#' Please follow the instruction.
#'
#' @param nSwarm integer. Swarm size.
#' @param maxIter integer. Maximal number of iterations.
#' @param checkConv logical. The stopping criterion.
#' @param freeRun numeric. The persentage of iterations which is free from examing the stopping criterion.
#' @param tol numeric. Convergence test. The default is \code{1e-6}.
#' @param typePSO integer, PSO type.
#' \itemize{
#' \item{0}{Linearly Decreasing Weight PSO (Shi, Y. H.	and Eberhart, R. C., 1998)}
#' \item{1}{GCPSO (van den Bergh, F. and	Engelbrecht, A. P., 2002)}
#' \item{2}{Quantum PSO (Sun, J., Feng, B. and Xu, W., 2004)}
#' \item{3}{LcRiPSO (Bonyadi, M. R., Michalewicz, Z., 2014)}
#' }
#' @param c1 numeric. Default is 2.
#' @param c2 numeric. Default is 2.
#' @return An List.
#' @name getAlgInfo
#' @rdname getAlgInfo
#' @export
getAlgInfo <- function(nSwarm = 64, typePSO = NULL, #dSwarm = NULL, varUpper = NULL, varLower = NULL,
  maxIter = NULL, checkConv = NULL, freeRun = NULL, tol = NULL, c1 = NULL, c2 = NULL,
  w0 = NULL, w1 = NULL, w_var = NULL, chi = NULL, vk = NULL,
  typeTopo = NULL, nGroup = NULL, GC_S_ROOF = NULL, GC_F_ROOF = NULL, GC_RHO = NULL,
  Q_cen_type = NULL, Q_a0 = NULL, Q_a1 = NULL, Q_a_var = NULL, LcRi_L = NULL,
  LBFGS_RETRY = NULL, LBFGS_MAXIT = NULL, LBFGS_LM = NULL, FVAL_EPS = NULL, GRAD_EPS = NULL, 
  LINESEARCH_MAXTRIAL = NULL, LINESEARCH_MAX = NULL, LINESEARCH_MIN = NULL,
  LINESEARCH_C = NULL, LINESEARCH_TAU = NULL) {


  nLoop <- length(nSwarm)
  #if(is.null(dSwarm))     dSwarm     <- NULL
  #if(is.null(varUpper))   varUpper   <- NULL
  #if(is.null(varLower))   varLower   <- NULL
  if(is.null(maxIter))    maxIter    <- rep(100   , nLoop)
  if(is.null(checkConv))  checkConv  <- rep(0     , nLoop)
  if(is.null(typePSO))    typePSO    <- rep(0     , nLoop)
  if(is.null(freeRun))    freeRun    <- rep(0.25  , nLoop)
  if(is.null(tol))        tol        <- rep(1e-6  , nLoop)
  if(is.null(c1))         c1         <- rep(2.05  , nLoop)
  if(is.null(c2))         c2         <- rep(2.05  , nLoop)
  if(is.null(w0))         w0         <- rep(1.2   , nLoop)
  if(is.null(w1))         w1         <- rep(0.2   , nLoop)
  if(is.null(w_var))      w_var      <- rep(0.8   , nLoop)
  if(is.null(chi))        chi        <- rep(0.729 , nLoop)
  if(is.null(vk))         vk         <- rep(4     , nLoop)
  if(is.null(typeTopo))   typeTopo   <- rep(0     , nLoop)
  if(is.null(nGroup))     nGroup     <- rep(1     , nLoop)
  if(is.null(GC_S_ROOF))  GC_S_ROOF  <- rep(5     , nLoop)
  if(is.null(GC_F_ROOF))  GC_F_ROOF  <- rep(15    , nLoop)
  if(is.null(GC_RHO))     GC_RHO     <- rep(1     , nLoop)
  if(is.null(Q_cen_type)) Q_cen_type <- rep(1     , nLoop)
  if(is.null(Q_a0))       Q_a0       <- rep(1.7   , nLoop)
  if(is.null(Q_a1))       Q_a1       <- rep(0.7   , nLoop)
  if(is.null(Q_a_var))    Q_a_var    <- rep(0.8   , nLoop)
  if(is.null(LcRi_L))     LcRi_L     <- rep(0.01  , nLoop)
  
  if(is.null(LBFGS_RETRY))          LBFGS_RETRY         <- rep(3   , nLoop)
  if(is.null(LBFGS_MAXIT))          LBFGS_MAXIT         <- rep(100 , nLoop)
  if(is.null(LBFGS_LM))             LBFGS_LM            <- rep(6   , nLoop)
  if(is.null(FVAL_EPS))             FVAL_EPS            <- rep(1e-8, nLoop)
  if(is.null(GRAD_EPS))             GRAD_EPS            <- rep(1e-5, nLoop)
  if(is.null(LINESEARCH_MAXTRIAL))  LINESEARCH_MAXTRIAL <- rep(50  , nLoop)
  if(is.null(LINESEARCH_MAX))       LINESEARCH_MAX      <- rep(1e6 , nLoop)
  if(is.null(LINESEARCH_MIN))       LINESEARCH_MIN      <- rep(1e-6, nLoop)
  if(is.null(LINESEARCH_C))         LINESEARCH_C        <- rep(1e-4, nLoop)
  if(is.null(LINESEARCH_TAU))       LINESEARCH_TAU      <- rep(0.25, nLoop)

  list(nSwarm = nSwarm, dSwarm = "autogen", varUpper = "autogen", varLower = "autogen", maxIter = maxIter,
    checkConv = checkConv, typePSO = typePSO, freeRun = freeRun, tol = tol, c1 = c1,
    c2 = c2, w0 = w0, w1 = w1, w_var = w_var, chi = chi,
    vk = vk, typeTopo = typeTopo, nGroup = nGroup, GC_S_ROOF = GC_S_ROOF, GC_F_ROOF = GC_F_ROOF,
    GC_RHO = GC_RHO, Q_cen_type = Q_cen_type, Q_a0 = Q_a0, Q_a1 = Q_a1, Q_a_var = Q_a_var,
    LcRi_L = LcRi_L,  LBFGS_RETRY = LBFGS_RETRY, LBFGS_MAXIT = LBFGS_MAXIT, LBFGS_LM = LBFGS_LM, FVAL_EPS = FVAL_EPS, GRAD_EPS = GRAD_EPS, 
    LINESEARCH_MAXTRIAL = LINESEARCH_MAXTRIAL, LINESEARCH_MAX = LINESEARCH_MAX, LINESEARCH_MIN = LINESEARCH_MIN,
    LINESEARCH_C = LINESEARCH_C, LINESEARCH_TAU = LINESEARCH_TAU)
}

#' Initialize for information of the optimal design problem
#'
#' Please follow the instruction.
#'
#' @param D_TYPE string. "exact"
#' @param nSupp integer. Number of support points.
#' @return An List.
#' @name getDesignInfo
#' @rdname getDesignInfo
#' @export
getDesignInfo <- function(D_TYPE = "approx", MODEL_INFO = NULL, dist_func = NULL, 
                          crit_type = 0, MaxMinStdVals = NULL,
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

  N_model <- length(MODEL_INFO)
  D_TYPE_NUM <- 1
  # Fit C++ indexing

  return(list(D_TYPE = D_TYPE, D_TYPE_NUM = D_TYPE_NUM, dist_func = dist_func,
              CRIT_TYPE_NUM = 0,
              dSupp = dSupp, nSupp = nSupp, dsLower = dsLower, dsUpper = dsUpper,
              N_model = N_model, dParas = dParas, paras = paras, parasInit = parasInit, 
              parasUpper = parasUpper, parasLower = parasLower, parasBdd = parasBdd, 
              MaxMinStdVals = MaxMinStdVals))
}
