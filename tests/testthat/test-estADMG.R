library(testthat)
library(flexCausal)

test_that("estADMG runs on simple backdoor example", {
  set.seed(1)

  n <- 200
  X <- rnorm(n)
  A <- rbinom(n, 1, plogis(X))
  Y <- 2 * A + X + rnorm(n)

  dat <- data.frame(A = A, Y = Y, X = X)

  graph <- make.graph(vertices=c('A','Y','X'),
                        bi_edges=NULL,
                        di_edges=list(c('X','A'), c('X','Y'), c('A','Y')))

  res <- estADMG(
    data = dat,
    graph = graph,
    a = c(1,0),
    treatment = "A",
    outcome = "Y"
  )

  expect_true(is.list(res))
  expect_true("TMLE" %in% names(res))
  expect_true("Onestep" %in% names(res))
  expect_true(is.numeric(res$TMLE$ACE))
})


test_that("estADMG works with different variable names", {
  set.seed(1)

  n <- 200
  X1 <- rnorm(n)
  Z <- rbinom(n, 1, plogis(X1))
  W <- 3 * Z + X1 + rnorm(n)

  dat <- data.frame(Z = Z, W = W, X1 = X1)

  graph <- make.graph(vertices=c('Z','W','X1'),
                      bi_edges=NULL,
                      di_edges=list(c('X1','Z'), c('X1','W'), c('Z','W')))

  res <- estADMG(
    data = dat,
    graph = graph,
    a = c(1,0),
    treatment = "Z",
    outcome = "W"
  )

  expect_true(is.numeric(res$TMLE$ACE))
  expect_true(is.numeric(res$Onestep$ACE))
})



test_that("estADMG works with frontdoor model", {
  set.seed(1)

  n <- 200
  X <- runif(n, 0, 1) # p(X)
  A <- rbinom(n, 1, (0.3 + 0.2*X)) # p(A|X)
  U <- 1 + 1*A + 1*X + 0*A*X + rnorm(n,0,1) # p(U|A,X)
  M <- rbinom(n,1,plogis(-1 + 1*A + 1*X + 0*A*X)) # p(M|A,X)
  Y <- 1*U + 1*M + 1*X + 0*M*X + rnorm(n, 0, 1) # p(Y|U,M,X)

  dat <- data.frame(X = X, A = A, M = M, Y=Y)

  graph <- make.graph(vertices=c('A','Y','X','M'),
                      bi_edges=list(c('A','Y')),
                      di_edges=list(c('X','A'), c('X','Y'), c('X','M'), c('A','M'), c('M','Y'))
                      )

  res <- estADMG(
    data = dat,
    graph = graph,
    a = c(1,0),
    treatment = "A",
    outcome = "Y"
  )

  expect_true(is.numeric(res$TMLE$ACE))
  expect_true(is.numeric(res$Onestep$ACE))
})


test_that("estADMG works with single intervention level a and output E(Y(a))", {
  set.seed(1)

  n <- 200
  X <- runif(n, 0, 1) # p(X)
  A <- rbinom(n, 1, (0.3 + 0.2*X)) # p(A|X)
  U <- 1 + 1*A + 1*X + 0*A*X + rnorm(n,0,1) # p(U|A,X)
  M <- rbinom(n,1,plogis(-1 + 1*A + 1*X + 0*A*X)) # p(M|A,X)
  Y <- 1*U + 1*M + 1*X + 0*M*X + rnorm(n, 0, 1) # p(Y|U,M,X)

  dat <- data.frame(X = X, A = A, M = M, Y=Y)

  graph <- make.graph(vertices=c('A','Y','X','M'),
                      bi_edges=list(c('A','Y')),
                      di_edges=list(c('X','A'), c('X','Y'), c('X','M'), c('A','M'), c('M','Y'))
  )

  res <- estADMG(
    data = dat,
    graph = graph,
    a = 1,
    treatment = "A",
    outcome = "Y"
  )

  expect_true(is.numeric(res$TMLE$EYa))
  expect_true(is.numeric(res$Onestep$EYa))
})

