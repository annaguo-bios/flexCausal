library(here)
library(MASS)
## Example (a) ==

set.seed(7)
## generate continuous outcome Y, continuous mediator M, single measured covariate X====
generate_data <- function(n,parU1= 1, parU2= -4, parA = c(0.3,0.2,1), parM = c(1,1,1,0), parL = c(0.5,2,1,1, -1), parY = c(1, 1, 1, -1), sd.M=1, sd.U1=1, sd.U2=1, sd.L =1, sd.Y=1){

  U1 <- rnorm(n,parU1,sd.U1)

  U2 <- rnorm(n,parU2,sd.U2)

  X <- runif(n, 0, 1) # p(X)

  A <- rbinom(n, 1, plogis(parA[1] + parA[2]*X+ parA[3]*U1)) # p(A|X)

  #A <- a

  M <- parM[1] + parM[2]*A + parM[3]*X + parM[4]*A*X + rnorm(n,0,sd.M) # p(M|A,X)

  L <- parL[1] + parL[2]*M + parL[3]*X + parL[4]*U1 + parL[5]*U2 + rnorm(n,0,sd.L) # p(L|A,X)

  Y <- parY[1] + parY[2]*L + parY[3]*X + parY[4]*U2 + rnorm(n, 0, sd.Y) # p(Y|U,M,X)

  data <- data.frame(X=X, U1=U1,U2=U2, A=A, M=M, L=L, Y=Y)


  return(list(data = data,
              parA=parA,
              parU1=parU1,
              parU2=parU2,
              parM=parM,
              parL=parL,
              parY=parY,
              sd.U1=sd.U1,
              sd.U2=sd.U2,
              sd.L=sd.L,
              sd.Y=sd.Y,
              sd.M=sd.M))
}

a.true.A1 <- mean(generate_data(50000)$data[,"Y"])


data.a <- generate_data(2000)$data
save(list= c("data.a", "a.true.A1"), file = "/Users/apple/Library/CloudStorage/Dropbox/primal-fixability/code/ADMGtmle/example-a.RData")



## Example (b) ====

set.seed(7)
## generate continuous outcome Y, continuous mediator M, single measured covariate X====
generate_data <- function(n,parU1= 1, parU2= -4, parA = c(0.3,0.2,1), parM = c(1,1,1,-1), parL = c(0.5,2,1,1), parY = c(1, 1, 1, -1,2), sd.M=1, sd.U1=1, sd.U2=1, sd.L =1, sd.Y=1){

  U1 <- rnorm(n,parU1,sd.U1)

  U2 <- rnorm(n,parU2,sd.U2)

  X <- runif(n, 0, 1) # p(X)

  A <- rbinom(n, 1, plogis(parA[1] + parA[2]*X+ parA[3]*U1)) # p(A|X)

  #A <- a

  M <- parM[1] + parM[2]*A + parM[3]*X + parM[4]*U2 + rnorm(n,0,sd.M) # p(M|A,X)

  L <- parL[1] + parL[2]*M + parL[3]*X + parL[4]*U1 + rnorm(n,0,sd.L) # p(L|A,X)

  Y <- parY[1] + parY[2]*L + parY[3]*X + parY[4]*U2 + parY[5]*A + rnorm(n, 0, sd.Y) # p(Y|U,M,X)

  data <- data.frame(X=X, U1=U1,U2=U2, A=A, M=M, L=L, Y=Y)


  return(list(data = data,
              parA=parA,
              parU1=parU1,
              parU2=parU2,
              parM=parM,
              parL=parL,
              parY=parY,
              sd.U1=sd.U1,
              sd.U2=sd.U2,
              sd.L=sd.L,
              sd.Y=sd.Y,
              sd.M=sd.M))
}

b.true.A1 <- mean(generate_data(50000)$data[,"Y"])


data.b <- generate_data(2000)$data
save(list= c("data.b", "b.true.A1"), file = "/Users/apple/Library/CloudStorage/Dropbox/primal-fixability/code/ADMGtmle/example-b.RData")



## Front-door ==
set.seed(7)

## generate continuous outcome Y, continuous mediator M, single measured covariate X====
generate_data <- function(n,parA = c(1,2), parU=c(1,1,1,0), parM = c(1,1,1,0), parY = c(1, 1, 1, 0), sd.M=1, sd.U=1, sd.Y=1){

  X <- runif(n, 0, 1) # p(X)

  A <- rbinom(n, 1, plogis(parA[1] + parA[2]*X)) # p(A|X)

 # A <- a

  U <- parU[1] + parU[2]*A + parU[3]*X + parU[4]*A*X + rnorm(n,0,sd.U) # p(U|A,X)

  M <- parM[1] + parM[2]*A + parM[3]*X + parM[4]*A*X + rnorm(n,0,sd.M) # p(M|A,X)

  Y <- parY[1]*U + parY[2]*M + parY[3]*X + parY[4]*M*X + rnorm(n, 0, sd.Y) # p(Y|U,M,X)

  data <- data.frame(X=X, U=U, A=A, M=M, Y=Y)

  # propensity score
  ps <- A*(parA[1] + parA[2]*X)+(1-A)*(1-(parA[1] + parA[2]*X))

  # mediator density ratio: p(M|a,X)/p(M|A,X)
  m.ratio.a1 <- dnorm(M,parM[1] + parM[2]*1 + parM[3]*X + parM[4]*1*X,sd.M)/dnorm(M,parM[1] + parM[2]*A + parM[3]*X + parM[4]*A*X,sd.M)
  m.ratio.a0 <- dnorm(M,parM[1] + parM[2]*0 + parM[3]*X + parM[4]*0*X,sd.M)/dnorm(M,parM[1] + parM[2]*A + parM[3]*X + parM[4]*A*X,sd.M)

  return(list(data = data,
              parA=parA,
              parU=parU,
              parM=parM,
              parY=parY,
              sd.U=sd.U,
              sd.Y=sd.Y,
              sd.M=sd.M,
              ps=ps,
              m.ratio.a1=m.ratio.a1,
              m.ratio.a0=m.ratio.a0))
}


frontdoor.true.A1 <- mean(generate_data(50000)$data[,"Y"])


data.frontdoor <- generate_data(2000)$data
save(list= c("data.frontdoor", "frontdoor.true.A1"), file = "/Users/apple/Library/CloudStorage/Dropbox/primal-fixability/code/ADMGtmle/example-frontdoor.RData")
