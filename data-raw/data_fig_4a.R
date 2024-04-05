## Simulation Figure4(a) ==
set.seed(7)

generate_data <- function(n,parA = c(1,1), parU=c(1,1,1,0), parM = matrix(c(1, 1, 1, 0,-1,-0.5,2,0), nrow = 2,byrow = T), parL= c(1,1,1, 1, 1) ,parY = c(1, 1, 1, 1, 1, 1), sd.U=1, sd.M= matrix(c(2, 1, 1, 3), nrow = 2), sd.L=1, sd.Y=1){

  X <- runif(n, 0, 1) # p(X)

  A <- rbinom(n, 1, plogis(parA[1] + parA[2]*X)) # p(A|X)

  U <- parU[1] + parU[2]*A + parU[3]*X + parU[4]*A*X + rnorm(n,0,sd.U) # p(U|A,X)

  M <- cbind(parM[1,1] + parM[1,2]*A + parM[1,3]*X + parM[1,4]*A*X,
             parM[2,1] + parM[2,2]*A + parM[2,3]*X + parM[2,4]*A*X)+ mvrnorm(n , mu =c(0,0) , Sigma = sd.M) # p(M|A,X)

  L <- parL[1] + parL[2]*A + parL[3]*M[,1] + parL[4]*M[,2] + parL[5]*X + rnorm(n,0,sd.L) # p(L|A,X)

  Y <- parY[1] + parY[2]*L  + parY[3]*M[,1] + parY[4]*M[,2] + parY[5]*X + parY[6]*U  + rnorm(n, 0, sd.Y) # p(Y|U,M,X)

  data <- data.frame(X=X, U=U, A=A, M=M, L=L, Y=Y)

  # propensity score
  ps <- A*(plogis(parA[1] + parA[2]*X))+(1-A)*(1-(plogis(parA[1] + parA[2]*X)))

  return(list(data = data,
              parA=parA,
              parU=parU,
              parM=parM,
              parL=parL,
              parY=parY,
              sd.U=sd.U,
              sd.M=sd.M,
              sd.L=sd.L,
              sd.Y=sd.Y,
              ps=ps))
}


data_fig_4a <- generate_data(2000)$data

usethis::use_data(data_fig_4a, overwrite = T)

