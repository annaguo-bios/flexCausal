
set.seed(7)
## Simulation Figure4(b) ==

generate_data <- function(n,parA = c(1,1), parU1=c(1,1,1,0),parU2=c(1,1,1,1,1), parM = matrix(c(1, 1, 1, 0,-1,-0.5,2,0), nrow = 2,byrow = T), parL= c(1,1,1, 1, 1) ,parY = c(1, 1, 1, 1, 1, 1), sd.U1=1, sd.U2=1, sd.M= matrix(c(2, 1, 1, 3), nrow = 2), sd.L=1, sd.Y=1){

  X <- runif(n, 0, 1) # p(X)

  A <- rbinom(n, 1, plogis(parA[1] + parA[2]*X)) # p(A|X)

  U1 <- parU1[1] + parU1[2]*A + parU1[3]*X + parU1[4]*A*X + rnorm(n,0,sd.U1) # p(U1|A,X)

  M <- cbind(parM[1,1] + parM[1,2]*A + parM[1,3]*X + parM[1,4]*A*X,
             parM[2,1] + parM[2,2]*A + parM[2,3]*X + parM[2,4]*A*X)+ mvrnorm(n , mu =c(0,0) , Sigma = sd.M) # p(M|A,X)

  U2 <- parU2[1] + parU2[2]*M[,1] + parU2[3]*M[,2] + parU2[4]*A + parU2[5]*X + rnorm(n,0,sd.U2) # p(U2|A,X,M)

  L <- parL[1] + parL[2]*M[,1] + parL[3]*M[,2] + parL[4]*X + parL[5]*U1 + rnorm(n,0,sd.L) # p(L|M,X,U1)

  Y <- parY[1] + parY[2]*L  + parY[3]*A + parY[4]*X + parY[5]*U2  + rnorm(n, 0, sd.Y) # p(Y|L,A,X,U2)

  data <- data.frame(X=X, U1=U1, U2=U2, A=A, M=M, L=L, Y=Y)

  # propensity score
  ps <- A*(plogis(parA[1] + parA[2]*X))+(1-A)*(1-(plogis(parA[1] + parA[2]*X)))

  return(list(data = data,
              parA=parA,
              parU1=parU1,
              parU2=parU2,
              parM=parM,
              parL=parL,
              parY=parY,
              sd.U1=sd.U1,
              sd.U2=sd.U2,
              sd.M=sd.M,
              sd.L=sd.L,
              sd.Y=sd.Y,
              ps=ps))
}


data_fig_4b <- generate_data(2000)$data

usethis::use_data(data_fig_4b, overwrite = T)
