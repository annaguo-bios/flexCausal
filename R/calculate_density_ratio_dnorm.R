#' @title Calculate the density ratio of M given its Markov pillow at treatment level a0 verse a1.
#' @description This function calculates the density ratio of two conditional densities. Let \eqn{mp(M)-A} denote the Markov pillow of M excluding A. The density ratio is calculated as \eqn{p(M|mp(M)-A, A=a0)/p(M|mp(M)-A, A=a1)}.
#' If M is univariate or multivariate continuous variable, it assume \eqn{p(M|mp(M)} is normally distributed. If M is a binary variable, it assume \eqn{p(M|mp(M)} is Bernoulli distributed.
#' @param a0 A numeric value indicating the treatment level. The density ratio is calculated at this level verse level 1-a0.
#' @param M A character string indicating the name of the variable, for which the density ratio is calculated.
#' @param graph A graph object created by \code{\link{make.graph}}.
#' @param treatment A character string indicating the name of the treatment variable A.
#' @param data A data frame containing the data.
#' @param formula A regression formula specifying the relationship between M or L and its Markov pillow. See the description for `ADMGtmle()` for more details.
#' @keywords graph density ratio
#' @return A vector of density ratios for each row of the data.
#' @export
#' @importFrom mvtnorm dmvnorm
#' @examples
#' graph <- make.graph(vertices=c('A','M','L','Y','X'),
#' bi_edges=list(c('A','Y')),
#' di_edges=list(c('X','A'), c('X','M'), c('X','L'),
#' c('X','Y'), c('M','Y'), c('A','M'), c('A','L'), c('M','L'), c('L','Y')),
#' multivariate.variables = list(M=c('M.1','M.2')))
#' ratio <- calculate_density_ratio_dnorm(a0=0, "M", graph, treatment="A", data=data_example_a)
#' head(ratio)
#'
calculate_density_ratio_dnorm <- function(a0, M, graph, treatment, data, formula=NULL){ # A is a vector, M and X are data frame

  # This function only allow M to be either univariate binary / continuous
  # or multivariate continuous.
  # It does not allow M to be a multivariate variable and has binary elements.


  # extract elements from graph
  multivariate.variables <- graph$multivariate.variables

  num_rows <- nrow(data)

  mp <- f.markov_pillow(graph, M, treatment=treatment) # Markov pillow for M
  mp <- replace.vector(mp, multivariate.variables) # replace multivariate vertices with their elements

  # if treatment is not in the Markov pillow of M, return 1
  if (!(treatment %in% mp)){
    return(1)
  }


  # prediction data
  data.a0 <- data[, mp] %>% mutate(!!treatment := a0)

  data.a1 <- data[, mp] %>% mutate(!!treatment := (1-a0))

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
        model <- glm(formula[[M]], data=data[, c(mp,M)])
      }else{
        model <- glm(data[, M] ~ . , data=data[, mp])
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

