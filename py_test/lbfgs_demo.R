

AF_para_m1 <- c(4.5, -1.5, -2)
optD <- cbind(c(-1,-.669,.144,.957), c(.253,.428,.247,.072))

opjFun <- function(para, optD) {
  val <- 0
  for (i in 1:nrow(optD)) {
    y1 <- MODEL_INFO_Cpp[[1]]$model(optD[i,1], c(4.5, -1.5, -2))
    y2 <- MODEL_INFO_Cpp[[2]]$model(optD[i,1], para)
    div_val <- DISTANCE_Cpp(y1, y2)
    val <- val + optD[i,2]*div_val
  }
  val
}

opjFunGrad <- function(para, optD) {
  gr <- numeric(length(para))
  for (i in 1:length(para)) {
    para_f <- para_b <- para
    para_f[i] <- para_f[i] + 1e-3
    para_b[i] <- para_b[i] - 1e-3
    fv <- opjFun(para_f, optD)
    bv <- opjFun(para_b, optD)
    gr[i] <- (fv - bv)/2e-3
  }
  gr
}



maxit <- 100

par0 <- rep(0, 3)
eyeM <- H0 <- diag(length(par0))
f0 <- opjFun(par0, optD)
g0 <- opjFunGrad(par0, optD)

H1 <- g1 <- NULL
par1 <- s1 <- y1 <- matrix(0, length(par0), maxit)
f1 <- rho <- numeric(maxit)

m <- 6
for (it in 1:maxit) {

  if (it > m) {
    q0 <- g0
    ai <- numeric(m)
    loc <- 1
    for (i in (it-1):(it-m)) {
      ai[loc] <- rho[i]*(t(s1[,i]) %*% q0)
      q0 <- q0 - ai[loc]*y1[,i]
      loc <- loc + 1
    }
    ai <- rev(ai)
    loc <- 1
    for (i in (it-m):(it-1)) {
      bi <- rho[i]*(t(y1[,i]) %*% q0)
      q0 <- q0 + s1[,i]*(ai[loc] - bi)
      loc <- loc + 1
    }
    d <- -q0
  } else {
    d <- -H0 %*% g0
  }

  std <- -0.5*(t(g0) %*% d)
  alpha <- 1e6
  counter <- 0; FLAG <- 1
  while ((counter < 20) & FLAG) {
    counter <- counter + 1;
    fv <- opjFun(par0, optD) - opjFun(par0 + alpha*d, optD)
    if (fv < alpha*std) {
      alpha <- 0.5*alpha
    } else {
      FLAG <- 0
    }
  }
  alpha <- max(alpha, 1e-20)

  par1[,it] <- par0 + alpha*d
  g1 <- opjFunGrad(par1[,it], optD)
  s1[,it] <- par1[,it] - par0
  y1[,it] <- g1 - g0

  rho[it] <- 1/sum(s1[,it]*y1[,it])
  H1 <- (eyeM - rho[it]*(y1[,it] %*% t(s1[,it])))  %*%
   H0 %*% (eyeM - rho[it]*(s1[,it] %*% t(y1[,it]))) + rho[it]*(s1[,it] %*% t(s1[,it]))

  par0 <- par1[,it]
  f1[it] <- opjFun(par0, optD)
  g0 <- g1
  H0 <- H1
}






