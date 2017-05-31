## ----LOAD_PKG, cache=TRUE,eval=TRUE,echo=TRUE----------------------------
library(DiscrimOD)

## ----AF1975b_MODEL, cache=T----------------------------------------------
m1 <- function(x, p) p[1] + p[2]*exp(x) + p[3]*exp(-x)
m2 <- function(x, p) p[1] + p[2]*x + p[3]*x^2
m3 <- function(x, p) p[1] + p[2]*sin(0.5*pi*x) + p[3]*cos(0.5*pi*x) + p[4]*sin(pi*x)

## ----call_empty, cache=T-------------------------------------------------
emptyModelList(N_model = 2)

## ----AF1975b_CASE_PAIR, cache=T------------------------------------------
AF_para_m1 <- c(4.5, -1.5, -2)
# The first pair: \eta_1 vs \eta_2
MODEL_INFO_12 <- list(
  # The first list should be the true model
  list(model = m1, para = AF_para_m1),
  # Then the rival models are listed accordingly
  list(model = m2, 
       paraLower = rep(-10, 3), paraUpper = rep(10, 3),
       paraInit = c(0, 0, 0))
)
# The second pair: \eta_1 vs \eta_3
MODEL_INFO_13 <- list(
  list(model = m1, para = AF_para_m1),
  list(model = m3,
       paraLower = rep(-10, 4), paraUpper = rep(10, 4),
       paraInit = c(0, 0, 0, 0))
)

## ----T_OPTIMAL, cache=T--------------------------------------------------
# xt is the mean values of the true model
# xr is the mean values of the rival model
DISTANCE <- function(xt, xr) (xt - xr)^2

## ----ALG_SETTING, cache=T------------------------------------------------
# getAlgInfo() # Call for seeing the default settings
ALG_INFO <- getAlgInfo(nSwarm = 16, maxIter = 200,        # PSO basic settings
                       LBFGS_RETRY = 3, GRAD_EPS = 1e-10) # L-BFGS basic settings

## ----AF1975b_RESULT_PAIR_12, cache=T-------------------------------------
# Run PSO-QN Algorithm for discriminating \eta_1 and \eta_2
DISC_12 <- DiscrimOD(MODEL_INFO = MODEL_INFO_12, DISTANCE = DISTANCE, 
                     crit_type = "pair_fixed_true",
                     nSupp = 4, dsLower = -1.0, dsUpper = 1.0, 
                     ALG_INFO = ALG_INFO, verbose = FALSE)
round(DISC_12$BESTDESIGN, 3)
DISC_12$BESTVAL
DISC_12$CPUTIME

## ----AF1975b_EQUV_PAIR_12, cache=T, fig.align='center', fig.height = 6, fig.width = 6, fig.cap='Directional derivative function of T-optimal design in the case DISC_12.'----
EQUV_12 <- equivalence(PSO_RESULT = DISC_12, MODEL_INFO = MODEL_INFO_12, 
                       DISTANCE = DISTANCE, crit_type = "pair_fixed_true",
                       dsLower = -1, dsUpper = 1, ngrid = 100)
# Draw the directional derivative curve
plot(EQUV_12$Grid_1, EQUV_12$DirDeriv, type = "l", col = "blue", 
     xlab = "x", ylab = "Directional Derivative"); abline(h = 0)
points(DISC_12$BESTDESIGN[,1], rep(0, nrow(DISC_12$BESTDESIGN)), pch = 19)

