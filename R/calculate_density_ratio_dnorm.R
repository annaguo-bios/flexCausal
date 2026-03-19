#' Calculate the density ratio of M given its Markov pillow at two treatment levels.
#'
#' Computes the ratio of two conditional densities of M evaluated at treatment
#' levels \eqn{a_0} and \eqn{1 - a_0}. Let \eqn{mp(M) \setminus A} denote the
#' Markov pillow of M excluding the treatment A. The density ratio is defined as
#' \deqn{\frac{p(M \mid mp(M) \setminus A,\, A = a_0)}{p(M \mid mp(M) \setminus A,\, A = 1 - a_0)}.}
#' The conditional density of M is modelled as follows: Gaussian for univariate
#' or multivariate continuous M (via linear regression), Bernoulli for binary M
#' (via logistic regression), and multivariate Gaussian for multivariate
#' continuous M (via separate linear regressions with a shared residual
#' covariance). Multivariate variables with binary elements are not supported.
#'
#' @param a0 Numeric. The reference treatment level; must be 0 or 1. The
#'   density ratio is computed as \eqn{p(M \mid A = a_0) / p(M \mid A = 1 - a_0)}.
#' @param M A character string naming the variable for which the density ratio
#'   is computed. May refer to a univariate or multivariate vertex as defined in
#'   \code{graph}.
#' @param graph A graph object created by \code{\link{make.graph}}.
#' @param treatment A character string naming the binary treatment variable A in
#'   \code{data}.
#' @param data A data frame containing all variables in the graph.
#' @param formula An optional named list of regression formulas, where each name
#'   is a variable name and the corresponding value is the formula to use for
#'   that variable's regression on its Markov pillow. Variables not included in
#'   this list are regressed using all Markov pillow variables as predictors.
#'   See \code{\link{estADMG}} for details.
#'
#' @return A numeric vector of length \code{nrow(data)} containing the
#'   density ratio for each observation. Returns a vector of ones if the
#'   treatment A is not in the Markov pillow of M (i.e., the ratio is
#'   identically 1).
#'
#' @note Currently only supports binary treatment coded as 0/1. Multivariate
#'   variables with binary elements are not supported.
#'
#' @keywords internal
#' @importFrom mvtnorm dmvnorm
calculate_density_ratio_dnorm <- function(a0, M, graph, treatment, data, formula=NULL){ # A is a vector, M and X are data frame

  # This function only allow M to be either univariate binary / continuous
  # or multivariate continuous.
  # It does not allow M to be a multivariate variable and has binary elements.


  if (!(a0 %in% c(0, 1))) stop("`a0` must be 0 or 1 for binary treatment.")

  # extract elements from graph
  multivariate.variables <- graph$multivariate.variables

  num_rows <- nrow(data)

  mp <- f.markov_pillow(graph, M, treatment=treatment) # Markov pillow for M
  mp <- replace.vector(mp, multivariate.variables) # replace multivariate vertices with their elements

  # if treatment is not in the Markov pillow of M, return 1
  if (!(treatment %in% mp)){
    return(rep(1, nrow(data)))
  }


  # prediction data
  data.a0 <- data[, mp, drop=F] %>% mutate(!!treatment := a0)

  data.a1 <- data[, mp, drop=F] %>% mutate(!!treatment := (1-a0))

  # for which variables, user specified formula for regression
  if (!is.null(formula)){
    formula.variables <- names(formula)
  }else{
    formula.variables <- c()
  }


  if (M %in% names(multivariate.variables)){ ## if M is multivariate ##

    num_columns <- length(multivariate.variables[[M]]) # number of variables under vertex M
    variables <- multivariate.variables[[M]]

    for (var in variables){

      if (all(data[, var] %in% c(0,1))) {

        stop("This function doesn't support multivariate variables with binary elements.")

      }

    } # end of for loop over variables

    ## Fit regression model to each component of M ##

    # model prediction
    predict.M.a0 <- list() # store the prediction results for each variable in M under A=a0
    predict.M.a1 <- list() # store the prediction results for each variable in M under A=a1

    # currently only support linear models

    # Initialize vectors to store errors
    model_errors <- matrix(NA, nrow = num_rows, ncol = num_columns)

    for (i in 1:num_columns){ # loop over each variable in M

      if (variables[i] %in% formula.variables){ # if the user specified formula for this variable
        model <- lm(as.formula(formula[[variables[i]]]), data=data[, c(variables[i],mp)])

      }else{ # if the user didn't specify formula for this variable
        model <- lm(data[, variables[i]] ~ . , data=data[, mp])
      }

      model_errors[, i] <- data[, variables[i]] - predict(model) # errors of the model

      predict.M.a0[[variables[i]]] <- predict(model, newdata=data.a0) # store the prediction result for A=a0
      predict.M.a1[[variables[i]]] <- predict(model, newdata=data.a1) # store the prediction result for A=a1

    } # end of for loop over variables

    # compute the variance-covariance matrix of the errors
    varcov <- cov(data.frame(model_errors))

    # Define a function for the ratio calculation
    f.m.ratio.a <- function(j) {

      mean.a0 <- sapply(predict.M.a0, `[`, j) # mean of the normal distribution for A=a0
      mean.a1 <- sapply(predict.M.a1, `[`, j) # mean of the normal distribution for A=a1

      dmvnorm(
        x = data[j, variables],
        mean = mean.a0, # coeff*mp
        sigma = varcov
      ) / dmvnorm(
        x = data[j, variables],
        mean = mean.a1, # coeff*mp
        sigma = varcov
      )
    }

    # Apply the function to each row of M
    ratio <- sapply(1:num_rows, f.m.ratio.a)




  }else{ ## if M is univariate ##

    if (all(data[, M] %in% c(0,1))) { # binary variable

      # For binary columns, use glm

      if (M %in% formula.variables){ # if the user specified formula for this variable
        model <- glm(formula[[M]], data=data[, c(mp,M)], family = binomial())
      }else{
        model <- glm(data[, M] ~ . , data=data[, mp], family = binomial())
      }

      # calculate the density ratio
      ratio <- predict(model, newdata=data.a0, type="response") / predict(model, newdata=data.a1, type="response")

    } else { # continuous variable

      # For continuous columns, use lm
      var <- data[, M]

      if (M %in% formula.variables){ # if the user specified formula for this variable
        model <- lm(as.formula(formula[[M]]), data=data[, c(mp,M)])

      }else{
        model <- lm( var ~ . , data=data[, mp])}

      # model prediction
      predict.M.a0 <- predict(model, newdata=data.a0) # store the prediction result for A=a0
      predict.M.a1 <- predict(model, newdata=data.a1) # store the prediction result for A=a1

      # Store model errors
      model_errors <- data[, M] - predict(model)

      ratio <- dnorm(as.vector(data[,M]), mean = predict.M.a0, sd = sd(model_errors))/
        dnorm(as.vector(data[,M]), mean = predict.M.a1, sd = sd(model_errors))

    } # end of if-else

  }

  return(ratio)

}

