#' Estimate the average counterfactual outcome E(Y(a)) via backdoor adjustment.
#'
#' Estimates the average counterfactual outcome \eqn{E\{Y(a)\}} using the Markov
#' pillow of the treatment variable as the adjustment set. Four estimators are
#' returned: G-computation, inverse probability weighting (IPW), one-step
#' corrected plug-in (AIPW), and targeted maximum likelihood estimation (TMLE).
#'
#' @param a Numeric scalar or length-two numeric vector specifying the treatment
#'   level(s) of interest. The treatment must be coded as 0/1. If a scalar,
#'   the function returns \eqn{E\{Y(a)\}}. If a length-two vector
#'   \code{c(a1, a0)}, the function returns the contrast
#'   \eqn{E\{Y(a1)\} - E\{Y(a0)\}}.
#' @param data A data frame containing all variables listed in \code{vertices}.
#' @param vertices A character vector of variable names in the causal graph.
#' @param di_edges A list of length-two character vectors specifying directed
#'   edges. For example, \code{di_edges = list(c('A', 'B'))} encodes a directed
#'   edge from A to B. Ignored if \code{graph} is provided.
#' @param bi_edges A list of length-two character vectors specifying bidirected
#'   edges. For example, \code{bi_edges = list(c('A', 'B'))} encodes a
#'   bidirected edge between A and B. Ignored if \code{graph} is provided.
#' @param multivariate.variables A named list mapping compound variable names to
#'   their component column names in \code{data}. For example,
#'   \code{list(M = c('M1', 'M2'))} indicates that M is bivariate with columns
#'   M1 and M2 in \code{data}. Ignored if \code{graph} is provided.
#' @param graph A graph object created by \code{\link{make.graph}}. If supplied,
#'   the arguments \code{vertices}, \code{di_edges}, \code{bi_edges}, and
#'   \code{multivariate.variables} are ignored. The user need only specify
#'   either \code{graph} or those four arguments.
#' @param treatment A character string naming the treatment variable in
#'   \code{data}. Must be coded as 0/1.
#' @param outcome A character string naming the outcome variable in \code{data}.
#' @param superlearner.Y Logical. If \code{TRUE}, the outcome regression is
#'   estimated via \code{\link[SuperLearner]{SuperLearner}} using the library
#'   specified in \code{lib.Y}. If \code{FALSE} (default), a GLM with formula
#'   \code{formulaY} is used.
#' @param superlearner.A Logical. If \code{TRUE}, the propensity score is
#'   estimated via \code{\link[SuperLearner]{SuperLearner}} using the library
#'   specified in \code{lib.A}. If \code{FALSE} (default), a logistic GLM with
#'   formula \code{formulaA} is used.
#' @param crossfit Logical. If \code{TRUE}, cross-fitting with \code{K} folds
#'   is applied when using SuperLearner. Ignored if both \code{superlearner.Y}
#'   and \code{superlearner.A} are \code{FALSE}.
#' @param K A positive integer specifying the number of folds for cross-fitting.
#'   Used only when \code{crossfit = TRUE}. Default is 5.
#' @param lib.Y A character vector specifying the SuperLearner library for
#'   outcome regression. Used only when \code{superlearner.Y = TRUE} or
#'   \code{crossfit = TRUE}. Default is
#'   \code{c("SL.glm", "SL.earth", "SL.ranger", "SL.mean")}.
#' @param lib.A A character vector specifying the SuperLearner library for
#'   propensity score estimation. Used only when \code{superlearner.A = TRUE}
#'   or \code{crossfit = TRUE}. Default is
#'   \code{c("SL.glm", "SL.earth", "SL.ranger", "SL.mean")}.
#' @param formulaY A formula or character string for the outcome regression of Y
#'   on its Markov pillow. Used only when \code{superlearner.Y = FALSE}.
#'   Default is \code{"Y ~ ."}.
#' @param formulaA A formula or character string for the propensity score
#'   regression of A on its Markov pillow. Used only when
#'   \code{superlearner.A = FALSE}. Default is \code{"A ~ ."}.
#' @param linkY_binary A character string specifying the link function for
#'   outcome regression when Y is binary and \code{superlearner.Y = FALSE}.
#'   Default is \code{"logit"}.
#' @param linkA A character string specifying the link function for propensity
#'   score regression when \code{superlearner.A = FALSE}. Default is
#'   \code{"logit"}.
#' @param truncate_lower Numeric. Propensity score values below this threshold
#'   are clipped to \code{truncate_lower}. Default is \code{0} (no clipping).
#' @param truncate_upper Numeric. Propensity score values above this threshold
#'   are clipped to \code{truncate_upper}. Default is \code{1} (no clipping).
#'
#' @return A named list with four components:
#' \describe{
#'   \item{\code{TMLE}}{Results from the TMLE estimator, a list containing:
#'     \describe{
#'       \item{\code{estimated_psi}}{Point estimate of \eqn{E\{Y(a)\}}.}
#'       \item{\code{lower.ci}}{Lower bound of the 95\% confidence interval.}
#'       \item{\code{upper.ci}}{Upper bound of the 95\% confidence interval.}
#'       \item{\code{EIF}}{Estimated efficient influence function evaluated at the observed data.}
#'       \item{\code{p.a.mpA}}{Estimated propensity score \eqn{P(A=a \mid mp(A))}.}
#'       \item{\code{mu.Y_a}}{Targeted outcome regression estimate \eqn{E[Y \mid A=a, mp(A)]}.}
#'     }
#'   }
#'   \item{\code{Onestep}}{Results from the one-step (AIPW) estimator, with the
#'     same fields as \code{TMLE}.}
#'   \item{\code{gcomp}}{Results from G-computation, a list containing:
#'     \describe{
#'       \item{\code{estimated_psi}}{Point estimate of \eqn{E\{Y(a)\}}.}
#'       \item{\code{mu.Y_a}}{Estimated outcome regression \eqn{E[Y \mid A=a, mp(A)]}.}
#'     }
#'   }
#'   \item{\code{ipw}}{Results from IPW, a list containing:
#'     \describe{
#'       \item{\code{estimated_psi}}{Point estimate of \eqn{E\{Y(a)\}}.}
#'       \item{\code{p.a.mpA}}{Estimated propensity score \eqn{P(A=a \mid mp(A))}.}
#'     }
#'   }
#' }
#' @keywords internal
#' @noRd
backdoor.TMLE.a <- function(a=NULL,data=NULL,vertices=NULL, di_edges=NULL, bi_edges=NULL, treatment=NULL, outcome=NULL, multivariate.variables=NULL, graph=NULL,
                            superlearner.Y=F, # whether run superlearner for outcome regression
                            superlearner.A=F, # whether run superlearner for propensity score
                            crossfit=F, K=5,
                            lib.Y = c("SL.glm","SL.earth","SL.ranger","SL.mean"), # superlearner library for outcome regression
                            lib.A = c("SL.glm","SL.earth","SL.ranger","SL.mean"), # superlearner library for propensity score
                            formulaY="Y ~ .", formulaA="A ~ .", # regression formula for outcome regression and propensity score if superlearner is not used
                            linkY_binary="logit", linkA="logit", # link function for outcome regression and propensity score if superlearner is not used
                            truncate_lower=0, truncate_upper=1){

  n <- nrow(data)

  a0 <- a
  a1 <- 1-a

  # extract graph components
  if (!is.null(graph)){

    vertices <- graph$vertices
    di_edges <- graph$di_edges
    bi_edges <- graph$bi_edges
    multivariate.variables <- graph$multivariate.variables

  }else{

    # make a graph object
    graph <- make.graph(vertices=vertices, bi_edges=bi_edges, di_edges=di_edges, multivariate.variables=multivariate.variables)

  }

  # Variables
  A <- data[,treatment] # treatment
  Y <- data[,outcome] # outcome



  ################################################
  ############### OUTCOME REGRESSION #############
  ################################################

  #### Fit nuisance models ####

  # Find Markov pillow of the treatment
  mpA <- f.markov_pillow(graph, treatment, treatment) # Markov pillow for A
  mpA <- replace.vector(mpA, multivariate.variables) # replace vertices with it's components if vertices are multivariate

  # prepare dataset for regression and prediction
  dat_Y <- data[,c(mpA,treatment), drop = F] # extract data for Markov pillow for outcome
  dat_Y.a0 <- dat_Y %>% mutate(!!treatment := a0) # set treatment to a0
  dat_Y.a1 <- dat_Y %>% mutate(!!treatment := a1) # set treatment to a1

  if (crossfit==T){ #### cross fitting + super learner #####

    fit.family <- if(all(Y %in% c(0,1))){binomial()}else{gaussian()} # family for super learner depending on whether Y is binary or continuous

    or_fit <- CV.SuperLearner(Y=Y, X=dat_Y, family = fit.family, V = K, SL.library = lib.Y, control = list(saveFitLibrary=T),saveAll = T)

    mu.Y_a1 <- unlist(lapply(1:K, function(x) predict(or_fit$AllSL[[x]], newdata=dat_Y.a1[or_fit$folds[[x]],])[[1]] %>% as.vector()))[order(unlist(lapply(1:K, function(x) or_fit$folds[[x]])))]
    mu.Y_a0 <- unlist(lapply(1:K, function(x) predict(or_fit$AllSL[[x]], newdata=dat_Y.a0[or_fit$folds[[x]],])[[1]] %>% as.vector()))[order(unlist(lapply(1:K, function(x) or_fit$folds[[x]])))]


  } else if (superlearner.Y==T){ #### super learner #####

    fit.family <- if(all(Y %in% c(0,1))){binomial()}else{gaussian()} # family for super learner depending on whether Y is binary or continuous

    or_fit <- SuperLearner(Y=Y, X=dat_Y, family = fit.family, SL.library = lib.Y)

    mu.Y_a1 <- predict(or_fit, newdata=dat_Y.a1)[[1]] %>% as.vector()
    mu.Y_a0 <- predict(or_fit, newdata=dat_Y.a0)[[1]] %>% as.vector()

  } else { #### simple linear regression with user input regression formula: default="Y ~ ." ####

    fit.family <- if(all(Y %in% c(0,1))){binomial()}else{gaussian()} # family for super learner depending on whether Y is binary or continuous

    or_fit <- glm(as.formula(formulaY), data=dat_Y, family = fit.family)

    mu.Y_a1 <- predict(or_fit, newdata=dat_Y.a1)
    mu.Y_a0 <- predict(or_fit, newdata=dat_Y.a0)

    if(all(Y %in% c(0,1))){
      mu.Y_a1 <- plogis(mu.Y_a1)
      mu.Y_a0 <- plogis(mu.Y_a0)
    }


  }


  ################################################
  ############### PROPENSITY SCORE ###############
  ################################################

  #### Fit nuisance models ####

  # if mpA is an empty set, then we just get the propensity score as the mean of the treatment
  if (length(mpA)==0){

    p.A1.mpA <- rep(mean(A==1), nrow(data))

  }else{

    # prepare dataset for regression and prediction
    dat_mpA <- data[,mpA, drop = F] # extract data for Markov pillow for outcome

    if (crossfit==T){ #### cross fitting + super learner #####

      # fit model
      ps_fit <- CV.SuperLearner(Y=A, X=dat_mpA, family = binomial(), V = K, SL.library = lib.A, control = list(saveFitLibrary=T),saveAll = T)

      # make prediction: p(A=1|mp(A))
      p.A1.mpA <- ps_fit$SL.predict

    } else if (superlearner.A==T){ #### super learner #####

      # fit model
      ps_fit <- SuperLearner(Y=A, X=dat_mpA, family = binomial(), SL.library = lib.A)

      # make prediction: p(A=1|mp(A))
      p.A1.mpA <- predict(ps_fit, type = "response")[[1]] %>% as.vector()

    } else { # without weak overlapping issue. Run A~X via logistic regression

      # fit model
      ps_fit <- glm(as.formula(formulaA), data=dat_mpA,  family = binomial())

      # make prediction: p(A=1|mp(A))
      p.A1.mpA <- predict(ps_fit, type = "response")  # p(A=1|X)

    }
  }

  # apply truncation to propensity score to deal with weak overlap.
  # truncated propensity score within the user specified range of [truncate_lower, truncate_upper]: default=[0,1]
  p.a0.mpA <- a0*p.A1.mpA + (1-a0)*(1-p.A1.mpA) # p(A=a1|mp(A))

  p.a0.mpA[p.a0.mpA < truncate_lower] <- truncate_lower
  p.a0.mpA[p.a0.mpA > truncate_upper] <- truncate_upper


  ## AIPW
  aipw <- mean((A==a0)*(Y-mu.Y_a0)/p.a0.mpA+mu.Y_a0)

  ## EIF
  EIF <- (A==a0)*(Y-mu.Y_a0)/p.a0.mpA+mu.Y_a0 - aipw


  # confidence interval
  lower.ci <- aipw-1.96*sqrt(mean(EIF^2)/nrow(data))
  upper.ci <- aipw+1.96*sqrt(mean(EIF^2)/nrow(data))

  aipw <- list(estimated_psi=aipw, # estimated parameter
               lower.ci=lower.ci, # lower bound of 95% CI
               upper.ci=upper.ci, # upper bound of 95% CI
               EIF=EIF, # E(Dstar) for Y|M,A,X and M|A,X, and A|X
               p.a.mpA = p.a0.mpA, # estimated E[A=a1|mp(A)]
               mu.Y_a = mu.Y_a0 # estimated E[Y|A=a0,mp(A)]
  )

  ## G-comp
  gcomp <- mean(mu.Y_a0)

  # # confidence interval
  # lower.ci <- gcomp-1.96*sqrt(mean(EIF^2)/nrow(data))
  # upper.ci <- gcomp+1.96*sqrt(mean(EIF^2)/nrow(data))

  gcomp <- list(estimated_psi=gcomp, # estimated parameter
                # lower.ci=lower.ci, # lower bound of 95% CI
                # upper.ci=upper.ci, # upper bound of 95% CI
                # EIF=EIF, # E(Dstar) for Y|M,A,X and M|A,X, and A|X
                mu.Y_a = mu.Y_a0 # estimated E[Y|A=a0,mp(A)]
  )

  ## IPW
  ipw <- mean((A==a0)*Y/p.a0.mpA)

  # confidence interval
  # lower.ci <- ipw-1.96*sqrt(mean(EIF^2)/nrow(data))
  # upper.ci <- ipw+1.96*sqrt(mean(EIF^2)/nrow(data))

  ipw <- list(estimated_psi=ipw, # estimated parameter
              # lower.ci=lower.ci, # lower bound of 95% CI
              # upper.ci=upper.ci, # upper bound of 95% CI
              # EIF=EIF, # E(Dstar) for Y|M,A,X and M|A,X, and A|X
              p.a.mpA = p.a0.mpA # estimated E[A=a1|mp(A)]
  )

  ## TMLE
  if(all(Y %in% c(0,1))){

    or_model <- glm(Y ~ offset(qlogis(mu.Y_a0))+((A==a0)*1/p.a0.mpA)-1, family=binomial(), start=0)

    eps.Y = coef(or_model)

    mu.Y_a0 <- plogis(qlogis(mu.Y_a0)+eps.Y*(A==a0)*1/p.a0.mpA)

  }else{

    or_model <- glm(Y ~ offset(mu.Y_a0)+1,weights = (A==a0)*1/p.a0.mpA)
    coef_Y <- coef(or_model)

    mu.Y_a0 <- mu.Y_a0+coef_Y

  }

  tmle.est <- mean((A==a0)*(Y-mu.Y_a0)/p.a0.mpA+mu.Y_a0)

  ## EIF
  EIF <- (A==a0)*(Y-mu.Y_a0)/p.a0.mpA+mu.Y_a0 - tmle.est


  # confidence interval
  lower.ci <- tmle.est-1.96*sqrt(mean(EIF^2)/nrow(data))
  upper.ci <- tmle.est+1.96*sqrt(mean(EIF^2)/nrow(data))

  tmle.out <- list(estimated_psi=tmle.est, # estimated parameter
                   lower.ci=lower.ci, # lower bound of 95% CI
                   upper.ci=upper.ci, # upper bound of 95% CI
                   EIF=EIF, # E(Dstar) for Y|M,A,X and M|A,X, and A|X
                   p.a.mpA = p.a0.mpA, # estimated E[A=a1|mp(A)]
                   mu.Y_a = mu.Y_a0 # estimated E[Y|A=a0,mp(A)]
  )


  return(list(TMLE=tmle.out, Onestep=aipw, gcomp=gcomp, ipw=ipw))








}
