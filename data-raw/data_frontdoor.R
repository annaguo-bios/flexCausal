## Simulation Front-door model ==
set.seed(7)

generate_data <- function(n,parA = c(0.3,0.2), parU=c(1,1,1,0), parM = c(-1,1,1,0), parY = c(1, 1, 1, 0), sd.U=1, sd.Y=1){ # change the parM to c(-1,1,1,0)

  X <- runif(n, 0, 1) # p(X)

  A <- rbinom(n, 1, (parA[1] + parA[2]*X)) # p(A|X)

  U <- parU[1] + parU[2]*A + parU[3]*X + parU[4]*A*X + rnorm(n,0,sd.U) # p(U|A,X)

  M <- rbinom(n,1,plogis(parM[1] + parM[2]*A + parM[3]*X + parM[4]*A*X)) # p(M|A,X)

  Y <- parY[1]*U + parY[2]*M + parY[3]*X + parY[4]*M*X + rnorm(n, 0, sd.Y) # p(Y|U,M,X)

  data <- data.frame(X=X, U=U, A=A, M=M, Y=Y)

  # propensity score
  ps <- A*(parA[1] + parA[2]*X)+(1-A)*(1-(parA[1] + parA[2]*X))

  # mediator density ratio: p(M|a,X)/p(M|A,X)
  m.ratio.a1 <- dbinom(M,1,plogis(parM[1] + parM[2]*1 + parM[3]*X + parM[4]*1*X))/dbinom(M,1,plogis(parM[1] + parM[2]*A + parM[3]*X + parM[4]*A*X))
  m.ratio.a0 <- dbinom(M,1,plogis(parM[1] + parM[2]*0 + parM[3]*X + parM[4]*0*X))/dbinom(M,1,plogis(parM[1] + parM[2]*A + parM[3]*X + parM[4]*A*X))


  return(list(data = data,
              parA=parA,
              parU=parU,
              parM=parM,
              parY=parY,
              sd.U=sd.U,
              sd.Y=sd.Y,
              ps=ps,
              m.ratio.a1=m.ratio.a1,
              m.ratio.a0=m.ratio.a0))
}


data_frontdoor <- generate_data(2000)$data

usethis::use_data(data_frontdoor, overwrite = T)
