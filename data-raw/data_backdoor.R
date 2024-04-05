## Simulation Figure4(a) ==
set.seed(7)

generate_data <- function(n,parA = c(1,1), parY = c(1, 1, 1, 0), sd.Y=1){

  X <- runif(n, 0, 1) # p(X)

  A <- rbinom(n, 1, plogis(parA[1] + parA[2]*X)) # p(A|X)

  Y <- parY[1] + parY[2]*X + parY[3]*A + parY[4]*A*X  + rnorm(n, 0, sd.Y) # p(Y|U,M,X)

  data <- data.frame(X=X, A=A, Y=Y)

  # propensity score
  ps <- A*(plogis(parA[1] + parA[2]*X))+(1-A)*(1-(plogis(parA[1] + parA[2]*X)))

  return(list(data = data,
              parA=parA,
              parY=parY,
              sd.Y=sd.Y,
              ps=ps))
}


data_backdoor <- generate_data(2000)$data

usethis::use_data(data_backdoor, overwrite = T)
