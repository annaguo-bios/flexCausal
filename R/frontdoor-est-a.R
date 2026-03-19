#' Estimate the average counterfactual outcome E(Y(a)) under a
#' nonparametrically saturated (NPS) model.
#'
#' Estimates \eqn{E\{Y(a)\}} for a treatment \eqn{A} and outcome \eqn{Y}
#' connected by a DAG that may contain hidden variables (encoded via bidirected
#' edges). The graph vertices are partitioned into sets C (pre-treatment), L
#' (in the district of A, post-treatment), and M (post-treatment, outside the
#' district of A). The efficient influence function (EIF) is constructed from
#' the outcome regression, propensity score, and density ratios
#' \eqn{p(v \mid mp(v))|_{a_0} / p(v \mid mp(v))|_{a_1}} for \eqn{v \in L \cup M}.
#' Two estimators are returned: a one-step corrected plug-in (AIPW) estimator
#' and a targeted maximum likelihood estimator (TMLE).
#'
#' @param a Numeric scalar or length-two numeric vector specifying the treatment
#'   level(s) of interest. The treatment must be coded as 0/1. If a scalar,
#'   the function returns \eqn{E\{Y(a)\}}. If a length-two vector
#'   \code{c(a1, a0)}, the function returns \eqn{E\{Y(a1)\} - E\{Y(a0)\}}.
#' @param data A data frame containing all variables listed in \code{vertices}.
#' @param vertices A character vector of variable names in the causal graph.
#' @param di_edges A list of length-two character vectors specifying directed
#'   edges. For example, \code{list(c('A','B'))} encodes A -> B.
#'   Ignored if \code{graph} is provided.
#' @param bi_edges A list of length-two character vectors specifying bidirected
#'   edges. For example, \code{list(c('A','B'))} encodes A <-> B.
#'   Ignored if \code{graph} is provided.
#' @param multivariate.variables A named list mapping compound vertex names to
#'   their column names in \code{data}. For example,
#'   \code{list(M = c('M.1','M.2'))} indicates M is bivariate with columns M.1
#'   and M.2. Ignored if \code{graph} is provided.
#' @param graph A graph object created by \code{\link{make.graph}}. If
#'   supplied, \code{vertices}, \code{di_edges}, \code{bi_edges}, and
#'   \code{multivariate.variables} are ignored.
#' @param treatment A character string naming the binary (0/1) treatment
#'   variable in \code{data}.
#' @param outcome A character string naming the outcome variable in \code{data}.
#' @param superlearner.seq Logical. If \code{TRUE}, SuperLearner is used for
#'   sequential regression of intermediate variables. If \code{FALSE} (default),
#'   linear regression is used.
#' @param superlearner.Y Logical. If \code{TRUE}, SuperLearner is used for
#'   outcome regression. If \code{FALSE} (default), a GLM with
#'   \code{formulaY} is used.
#' @param superlearner.A Logical. If \code{TRUE}, SuperLearner is used for
#'   propensity score estimation. If \code{FALSE} (default), a logistic GLM
#'   with \code{formulaA} is used.
#' @param superlearner.M Logical. If \code{TRUE}, SuperLearner is used when
#'   estimating the density ratios for variables in M via the Bayes method
#'   (\code{ratio.method.M = "bayes"}). Ignored otherwise.
#' @param superlearner.L Logical. If \code{TRUE}, SuperLearner is used when
#'   estimating the density ratios for variables in L via the Bayes method
#'   (\code{ratio.method.L = "bayes"}). Ignored otherwise.
#' @param crossfit Logical. If \code{TRUE}, cross-fitting with \code{K} folds
#'   is applied to all SuperLearner fits. Ignored if all SuperLearner flags are
#'   \code{FALSE}.
#' @param K A positive integer specifying the number of cross-fitting folds.
#'   Used only when \code{crossfit = TRUE}. Default is 5.
#' @param ratio.method.L A character string specifying the method for
#'   estimating density ratios \eqn{p(L|mp(L))|_{a_0}/p(L|mp(L))|_{a_1}}.
#'   Options are:
#'   \describe{
#'     \item{\code{"bayes"}}{(Default) Rewrites the ratio via Bayes' rule as
#'       \eqn{[p(A=a_0|L,mp(L))/p(A=a_1|L,mp(L))] / [p(A=a_0|mp(L))/p(A=a_1|mp(L))]}
#'       and estimates each factor via logistic regression or SuperLearner.}
#'     \item{\code{"dnorm"}}{Assumes \eqn{L | mp(L), A} is Gaussian (continuous L)
#'       or Bernoulli (binary L), estimated via linear or logistic regression
#'       with linear terms only.}
#'     \item{\code{"densratio"}}{Uses the \code{\link[densratio]{densratio}}
#'       package. Supports only numeric/integer variables and is computationally
#'       expensive; not recommended for graphs with many variables.}
#'   }
#' @param ratio.method.M A character string specifying the method for
#'   estimating density ratios \eqn{p(M|mp(M))|_{a_0}/p(M|mp(M))|_{a_1}}.
#'   Same three options as \code{ratio.method.L}. Default is \code{"bayes"}.
#' @param dnorm.formula.L An optional named list of regression formulas for
#'   variables in L, used when \code{ratio.method.L = "dnorm"}. Names are
#'   variable names; values are formula strings. Variables omitted from the
#'   list are regressed on all Markov pillow variables. For multivariate L,
#'   specify one formula per component, e.g.
#'   \code{list(L.1 = "L.1 ~ A + X", L.2 = "L.2 ~ A + X + I(M^2)")}.
#' @param dnorm.formula.M An optional named list of regression formulas for
#'   variables in M, used when \code{ratio.method.M = "dnorm"}. Same structure
#'   as \code{dnorm.formula.L}.
#' @param lib.seq A character vector specifying the SuperLearner library for
#'   sequential regression. Used only when \code{superlearner.seq = TRUE} or
#'   \code{crossfit = TRUE}.
#' @param lib.L A character vector specifying the SuperLearner library for
#'   density ratio estimation for L via the Bayes method.
#' @param lib.M A character vector specifying the SuperLearner library for
#'   density ratio estimation for M via the Bayes method.
#' @param lib.Y A character vector specifying the SuperLearner library for
#'   outcome regression.
#' @param lib.A A character vector specifying the SuperLearner library for
#'   propensity score estimation.
#' @param formulaY A formula or character string for the outcome regression of
#'   Y on its Markov pillow. Used only when \code{superlearner.Y = FALSE}.
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
#' @param cvg.criteria Numeric. TMLE convergence threshold: the iterative
#'   update stops when \eqn{|\text{mean}(D^*)| <} \code{cvg.criteria}.
#'   Default is \code{0.01}.
#' @param n.iter A positive integer specifying the maximum number of TMLE
#'   iterations. Default is 500.
#' @param truncate_lower Numeric. Propensity score values below this threshold
#'   are clipped to \code{truncate_lower}. Default is \code{0} (no clipping).
#' @param truncate_upper Numeric. Propensity score values above this threshold
#'   are clipped to \code{truncate_upper}. Default is \code{1} (no clipping).
#' @param zerodiv.avoid Numeric. Values of density ratios or propensity scores
#'   below this threshold are clipped to \code{zerodiv.avoid} to prevent
#'   division by zero. Default is \code{0} (no clipping).
#'
#' @return A named list with two components, \code{Onestep} and \code{TMLE},
#'   each containing:
#' \describe{
#'   \item{\code{estimated_psi}}{Point estimate of \eqn{E\{Y(a)\}}.}
#'   \item{\code{lower.ci}}{Lower bound of the 95\% confidence interval.}
#'   \item{\code{upper.ci}}{Upper bound of the 95\% confidence interval.}
#'   \item{\code{EIF}}{Estimated efficient influence function, a numeric vector
#'     of length \code{nrow(data)}.}
#'   \item{\code{EIF.Y}}{EIF component from the outcome regression \eqn{E[Y|mp(Y)]}.}
#'   \item{\code{EIF.A}}{EIF component from the propensity score \eqn{p(A|mp(A))}.}
#'   \item{\code{EIF.v}}{Summed EIF components from sequential regressions for
#'     variables between A and Y in topological order.}
#'   \item{\code{p.a1.mpA}}{Estimated \eqn{P(A = 1-a \mid mp(A))}.}
#'   \item{\code{mu.next.A}}{Estimated \eqn{E[v \mid mp(v)]} for the vertex
#'     immediately after A in topological order, evaluated at \eqn{A = a}.}
#' }
#' The \code{TMLE} component additionally contains:
#' \describe{
#'   \item{\code{EDstar}}{Final value of \eqn{\text{mean}(D^*)} at convergence.
#'     Should be close to 0 if converged.}
#'   \item{\code{EDstar.record}}{Numeric vector recording \eqn{\text{mean}(D^*)}
#'     at each iteration.}
#'   \item{\code{iter}}{Number of iterations performed.}
#' }
#' @noRd
#' @keywords internal
#' @importFrom dplyr %>% mutate select
#' @importFrom MASS mvrnorm
#' @importFrom SuperLearner CV.SuperLearner SuperLearner
#' @importFrom mvtnorm dmvnorm
#' @importFrom densratio densratio
#' @importFrom utils combn
#' @importFrom stats rnorm runif rbinom dnorm dbinom binomial gaussian predict glm as.formula qlogis plogis lm coef cov sd
#'
#'
NPS.TMLE.a <- function(a=NULL,data=NULL,vertices=NULL, di_edges=NULL, bi_edges=NULL, treatment=NULL, outcome=NULL, multivariate.variables=NULL, graph=NULL,
                     superlearner.seq = F, # whether run superlearner for sequential regression
                     superlearner.Y=F, # whether run superlearner for outcome regression
                     superlearner.A=F, # whether run superlearner for propensity score
                     superlearner.M=F, # whether run superlearner for estimating densratio for M using bayes method
                     superlearner.L=F, # whether run superlearner for estimating densratio for L using bayes method
                     crossfit=F, K=5,
                     ratio.method.L="bayes", # method for estimating the density ratio associated with L
                     ratio.method.M="bayes", # method for estimating the density ratio associated with M
                     dnorm.formula.M=NULL, # formula for regression of M when the ratio.method.M is 'dnorm'
                     dnorm.formula.L=NULL, # formula for regression of L when the ratio.method.L is 'dnorm'
                     lib.seq = c("SL.glm","SL.earth","SL.ranger","SL.mean"), # superlearner library for sequential regression
                     lib.L = c("SL.glm","SL.earth","SL.ranger","SL.mean"), # superlearner library for density ratio estimation via bayes rule for variables in L
                     lib.M = c("SL.glm","SL.earth","SL.ranger","SL.mean"), # superlearner library for density ratio estimation via bayes rule for variables in M
                     lib.Y = c("SL.glm","SL.earth","SL.ranger","SL.mean"), # superlearner library for outcome regression
                     lib.A = c("SL.glm","SL.earth","SL.ranger","SL.mean"), # superlearner library for propensity score
                     formulaY="Y ~ .", formulaA="A ~ .", # regression formula for outcome regression and propensity score if superlearner is not used
                     linkY_binary="logit", linkA="logit", # link function for outcome regression and propensity score if superlearner is not used
                     n.iter=500, cvg.criteria=0.01,
                     truncate_lower=0, truncate_upper=1, zerodiv.avoid=0){

  # attach(data, warn.conflicts=FALSE)

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

  # return topological ordering
  tau <- f.top_order(graph, treatment)
  tau.df <- data.frame(tau=tau, order = 1:length(tau))

  # Get set C, M, L
  setCML <- CML(graph, treatment) # get set C, M, L

  C <- setCML$C # everything comes before the treatment following topological order tau

  L <- setCML$L # variables within the district of treatment and comes after the treatment (including the treatment itself) following topological order tau
  # L <- setdiff(L, outcome) # remove outcome from L if it's there

  M <- setCML$M # everything else
  # M <- setdiff(M,outcome) # remove outcome from M if it's there

  # re-order vertices according to their topological order in tau

  C <- rerank(C, tau) # re-order vertices in C according to their topological order in tau

  L <- rerank(L, tau) # re-order vertices in L according to their topological order in tau
  L.removedA <- L[L!=treatment] # remove treatment from L

  M <- rerank(M, tau) # re-order vertices in M according to their topological order in tau

  # Variables
  A <- data[,treatment] # treatment
  Y <- data[,outcome] # outcome


  are_same <- function(vec1, vec2) {
    identical(sort(vec1), sort(vec2))
  }

  #///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////#
  ##///////////////////////////// STEP1: TMLE initialization for sequential regression based estimator///////////////////////##
  #///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////#

  ################################################
  ############### OUTCOME REGRESSION #############
  ################################################

  #### Fit nuisance models ####

  # Find Markov pillow of outcome
  mpY <- f.markov_pillow(graph, node=outcome, treatment=treatment) # Markov pillow for outcome
  mpY <- replace.vector(mpY, multivariate.variables) # replace vertices with it's components if vertices are multivariate

  # prepare dataset for regression and prediction
  dat_mpY <- data[,mpY, drop=F] # extract data for Markov pillow for outcome
  dat_mpY.a0 <- dat_mpY %>% mutate(!!treatment := a0) # set treatment to a0
  dat_mpY.a1 <- dat_mpY %>% mutate(!!treatment := a1) # set treatment to a1

  if (crossfit==T){ #### cross fitting + super learner #####

    fit.family <- if(all(Y %in% c(0,1))){binomial(linkY_binary)}else{gaussian()} # family for super learner depending on whether Y is binary or continuous

    or_fit <- CV.SuperLearner(Y=Y, X=dat_mpY, family = fit.family, V = K, SL.library = lib.Y, control = list(saveFitLibrary=T),saveAll = T)

    mu.Y_a1 <- unlist(lapply(1:K, function(x) predict(or_fit$AllSL[[x]], newdata=dat_mpY.a1[or_fit$folds[[x]],])[[1]] %>% as.vector()))[order(unlist(lapply(1:K, function(x) or_fit$folds[[x]])))]
    mu.Y_a0 <- unlist(lapply(1:K, function(x) predict(or_fit$AllSL[[x]], newdata=dat_mpY.a0[or_fit$folds[[x]],])[[1]] %>% as.vector()))[order(unlist(lapply(1:K, function(x) or_fit$folds[[x]])))]

    assign(paste0("mu.",outcome,"_a0"), mu.Y_a0)
    assign(paste0("mu.",outcome,"_a1"), mu.Y_a1)

  } else if (superlearner.Y==T){ #### super learner #####

    fit.family <- if(all(Y %in% c(0,1))){binomial(linkY_binary)}else{gaussian()} # family for super learner depending on whether Y is binary or continuous

    or_fit <- SuperLearner(Y=Y, X=dat_mpY, family = fit.family, SL.library = lib.Y)

    mu.Y_a1 <- predict(or_fit, newdata=dat_mpY.a1)[[1]] %>% as.vector()
    mu.Y_a0 <- predict(or_fit, newdata=dat_mpY.a0)[[1]] %>% as.vector()

    assign(paste0("mu.",outcome,"_a0"), mu.Y_a0)
    assign(paste0("mu.",outcome,"_a1"), mu.Y_a1)

  } else { #### simple linear regression with user input regression formula: default="Y ~ ." ####

    fit.family <- if(all(Y %in% c(0,1))){binomial(linkY_binary)}else{gaussian()} # family for super learner depending on whether Y is binary or continuous

    or_fit <- glm(as.formula(formulaY), data=dat_mpY, family = fit.family)

    mu.Y_a1 <- predict(or_fit, newdata=dat_mpY.a1, type = "response")
    mu.Y_a0 <- predict(or_fit, newdata=dat_mpY.a0, type = "response")

    assign(paste0("mu.",outcome,"_a0"), mu.Y_a0)
    assign(paste0("mu.",outcome,"_a1"), mu.Y_a1)

  }


  ################################################
  ############### PROPENSITY SCORE ###############
  ################################################

  #### Fit nuisance models ####

  # Find Markov pillow of treatment
  mpA <- f.markov_pillow(graph, node=treatment, treatment=treatment) # Markov pillow for treatment
  mpA <- replace.vector(mpA, multivariate.variables) # replace vertices with it's components if vertices are multivariate

  # if mpA is an empty set, then we just get the propensity score as the mean of the treatment
  if (length(mpA)==0){
    p.A1.mpA <- rep(mean(A==1), nrow(data))

  }else{

    # prepare dataset for regression and prediction
    dat_mpA <- data[,mpA, drop = F] # extract data for Markov pillow for outcome

    if (crossfit==T){ #### cross fitting + super learner #####

      # fit model
      ps_fit <- CV.SuperLearner(Y=A, X=dat_mpA, family = binomial(linkA), V = K, SL.library = lib.A, control = list(saveFitLibrary=T),saveAll = T)

      # make prediction: p(A=1|mp(A))
      p.A1.mpA <- ps_fit$SL.predict

    } else if (superlearner.A==T){ #### super learner #####

      # fit model
      ps_fit <- SuperLearner(Y=A, X=dat_mpA, family = binomial(linkA), SL.library = lib.A)

      # make prediction: p(A=1|mp(A))
      p.A1.mpA <- predict(ps_fit, type = "response")[[1]] %>% as.vector()

    } else { # without weak overlapping issue. Run A~X via logistic regression

        # fit model
        ps_fit <- glm(as.formula(formulaA), data=dat_mpA,  family = binomial(linkA))

        # make prediction: p(A=1|mp(A))
        p.A1.mpA <- predict(ps_fit, type = "response")  # p(A=1|X)

      }
    }

 # end of if (length(mpA)==0)

  # apply truncation to propensity score to deal with weak overlap.
  # truncated propensity score within the user specified range of [truncate_lower, truncate_upper]: default=[0,1]
  p.A1.mpA[p.A1.mpA < truncate_lower] <- truncate_lower
  p.A1.mpA[p.A1.mpA > truncate_upper] <- truncate_upper

  p.a1.mpA <- a1*p.A1.mpA + (1-a1)*(1-p.A1.mpA) # p(A=a1|mp(A))
  p.a0.mpA <-1-p.a1.mpA # p(A=a0|mp(A))

  # avoid zero indivision error
  p.a0.mpA[p.a0.mpA<zerodiv.avoid] <- zerodiv.avoid
  p.a1.mpA[p.a1.mpA<zerodiv.avoid] <- zerodiv.avoid


  assign(paste0("densratio_",treatment), p.a0.mpA/p.a1.mpA) # density ratio regarding the treatment p(A|mp(A))|_{a_0}/p(A|mp(A))|_{a_1}




  ################################################
  ############### DENSITY RATIO.   ###############
  ################################################

  ###### DENSITY RATIO ASSOCIATED WITH L ######

  if (ratio.method.L=="densratio"){ ################### METHOD 2A: densratio method  ###################

    # Error3: densratio method doesn't support factor variables
    if (!all(sapply(replace.vector(unique(c(L.removedA, unlist(lapply(1:length(L.removedA),function(i) f.markov_pillow(graph, L.removedA, treatment))))), multivariate.variables), function(var) is.numeric(data[,var]) | is.integer(data[,var])))){
      stop("Error in estimating density ratios associated with variables in L: densratio method only support numeric/integer variables, try bayes method instead.")
      }


    L.use.densratioA <- c()

    # if M,A,X only consists numeric/integer variables: apply density ratio estimation
    for (v in L.removedA){ ## Iterate over each variable in L\A

      # used for estimating numerator of the density ratio
      dat_v.num.a0 <- data[data[[treatment]] == a0, replace.vector(c(v, f.markov_pillow(graph, v, treatment)), multivariate.variables)] # select rows where A=a0
      dat_v.num.a1 <- data[data[[treatment]] == a1, replace.vector(c(v, f.markov_pillow(graph, v, treatment)), multivariate.variables)] # select rows where A=a1


      # since for v in L, the ratio used is p(L|mp(L))|_{a_1}/p(L|mp(L))|_{a_0}. Therefore, if we calculate the ratio p(L|mp(L))|_{a_0}/p(L|mp(L))|_{a_1} and make it be divided by 1,
      # it may result in large values. Therefore, we calculate the ratio p(L|mp(L))|_{a_1}/p(L|mp(L))|_{a_0} and assign it to 1/ratio to make the ratio estimation more stable.
      densratio.v.num <- densratio(dat_v.num.a1, dat_v.num.a0) # calculate the ratio of a1/a0

      ratio.num <- densratio.v.num$compute_density_ratio(data[, replace.vector(c(v, f.markov_pillow(graph, v, treatment)), multivariate.variables)])

      ratio.num[ratio.num<zerodiv.avoid] <- zerodiv.avoid # control for very small numerator values

      ## Estimate the denomator via bays rule to avoid zero division error
      # used for estimating the denominator of the density ratio
      dat_bayes.v <- data[, setdiff(replace.vector(f.markov_pillow(graph, v, treatment), multivariate.variables), treatment), drop=F]
      names(dat_bayes.v)[names(dat_bayes.v)=="Y"] <- "outcome" # Super Learner get confused of the regresor contains a variable named Y. Thus rename it to outcome

      if (are_same(setdiff(replace.vector(f.markov_pillow(graph, v, treatment), multivariate.variables) ,treatment), replace.vector(f.markov_pillow(graph,treatment, treatment), multivariate.variables))){

        ratio.den <- 1/get(paste0("densratio_",treatment))

        assign(paste0("densratio.densratio_",v), ratio.num)

        L.use.densratioA <- c(L.use.densratioA, v)

      }else{



        if (crossfit==T){

          bayes_fit <- CV.SuperLearner(Y=A, X=dat_bayes.v, family = binomial(), V = K, SL.library = lib.L, control = list(saveFitLibrary=T),saveAll = T)

          # p(A=1|mp(v)\A,v)
          p.A1.mpv <- bayes_fit$SL.predict

          #p(v=a0|mp(v)\A,v)
          p.a0.mpv <- a0*p.A1.mpv+(1-a0)*(1-p.A1.mpv)

          #p(v=a0|mp(v)\A,v)
          p.a1.mpv <- 1-p.a0.mpv

          p.a1.mpv[p.a1.mpv<zerodiv.avoid] <- zerodiv.avoid # added to avoid INF
          p.a0.mpv[p.a0.mpv<zerodiv.avoid] <- zerodiv.avoid # added to avoid INF

          # save the density ratio: p(a1|mp(V))/p(a0|mp(V))
          ratio.den <- p.a1.mpv/p.a0.mpv


        }else if (superlearner.L==T){

          bayes_fit <- SuperLearner(Y=A, X=dat_bayes.v, family = binomial(), SL.library = lib.L)

          # p(A=1|mp(v)\A,v)
          p.A1.mpv <- predict(bayes_fit, type = "response")[[1]] %>% as.vector()  # p(A=1|X)

          #p(v=a0|mp(v)\A,v)
          p.a0.mpv <- a0*p.A1.mpv+(1-a0)*(1-p.A1.mpv)

          #p(v=a0|mp(v)\A,v)
          p.a1.mpv <- 1-p.a0.mpv

          p.a1.mpv[p.a1.mpv<zerodiv.avoid] <- zerodiv.avoid # added to avoid INF
          p.a0.mpv[p.a0.mpv<zerodiv.avoid] <- zerodiv.avoid # added to avoid INF

          # save the density ratio: p(a1|mp(V))/p(a0|mp(V))
          ratio.den <- p.a1.mpv/p.a0.mpv

        } else {

          # estimate density ratio using bayes rule
          bayes_fit <- glm(A ~ ., data=dat_bayes.v, family = binomial())

          # p(A=1|mp(v)\A,v)
          p.A1.mpv <- predict(bayes_fit, type = "response")

          #p(v=a0|mp(v)\A,v)
          p.a0.mpv <- a0*p.A1.mpv+(1-a0)*(1-p.A1.mpv)

          #p(v=a0|mp(v)\A,v)
          p.a1.mpv <- 1-p.a0.mpv

          p.a1.mpv[p.a1.mpv<zerodiv.avoid] <- zerodiv.avoid # added to avoid INF
          p.a0.mpv[p.a0.mpv<zerodiv.avoid] <- zerodiv.avoid # added to avoid INF

          # save the density ratio: p(a1|mp(V))/p(a0|mp(V))
          ratio.den <- p.a1.mpv/p.a0.mpv

        }





      }


      ratio <- ratio.num/ratio.den # p(L|mp(L))|_{a_1}/p(L|mp(L))|_{a_0}

      assign(paste0("densratio_",v), 1/ratio) # but assign ratio a0/a1
    }


  } else if (ratio.method.L=="dnorm"){ ################### METHOD 2B: dnorm method  ###################

    if (!all(sapply(replace.vector(L.removedA, multivariate.variables), function(var) is.numeric(data[,var]) | is.integer(data[,var]) | length(unique(data[,var]))==2))){
      stop("Error in estimating density ratios associated with variables in L: dnorm method only support continuous or binary variables, try bayes method instead.")
      }


    # if M,A,X only consists numeric/integer variables: apply density ratio estimation
    for (v in L.removedA){ ## Iterate over each variable in L\A

      ratio <- calculate_density_ratio_dnorm(a0=a0, v , graph, treatment=treatment, data=data, formula=dnorm.formula.L) # p(L|mp(L))|_{a_0}/p(L|mp(L))|_{a_1}

      assign(paste0("densratio_",v), ratio)
    }


  }else if (ratio.method.L=="bayes"){ ################### METHOD 2C: Bayes method ###################

    L.use.densratioA <- c()

    for (v in L.removedA){ ## Iterate over each variable in L\A

      #### Prepare data for regression and prediction ####
      dat_bayes.v <- data[,setdiff( replace.vector(c(v, f.markov_pillow(graph, v, treatment)), multivariate.variables) ,treatment), drop=F] # contains variable v + Markov pillow of v - treatment
      names(dat_bayes.v)[names(dat_bayes.v)=="Y"] <- "outcome" # Super Learner get confused of the regresor contains a variable named Y. Thus rename it to outcome

      dat_bayes.v_v <- data[,setdiff( replace.vector(f.markov_pillow(graph, v, treatment), multivariate.variables) ,treatment), drop=F] # contains variable Markov pillow of v - treatment

      #### Fit nuisance models ####

      if (crossfit==T){

        bayes_fit <- CV.SuperLearner(Y=A, X=dat_bayes.v, family = binomial(), V = K, SL.library = lib.L, control = list(saveFitLibrary=T),saveAll = T)

        # p(A=1|mp(v)\A,v)
        p.A1.mpv <- bayes_fit$SL.predict

        #p(v=a0|mp(v)\A,v)
        p.a0.mpv <- a0*p.A1.mpv+(1-a0)*(1-p.A1.mpv)

        #p(v=a0|mp(v)\A,v)
        p.a1.mpv <- 1-p.a0.mpv

        p.a1.mpv[p.a1.mpv<zerodiv.avoid] <- zerodiv.avoid # added to avoid INF
        p.a0.mpv[p.a0.mpv<zerodiv.avoid] <- zerodiv.avoid # added to avoid INF


        if (are_same(setdiff(replace.vector(f.markov_pillow(graph, v, treatment), multivariate.variables) ,treatment), replace.vector(f.markov_pillow(graph,treatment, treatment), multivariate.variables))){

          # save the density ratio: p(v|mp(V))|_{a0}/p(v|mp(V))|_{a1}
          assign(paste0("densratio_",v), {p.a0.mpv/p.a1.mpv}/get(paste0("densratio_",treatment)))

          # save the ratio of p(A|mp(v)\A,v)|_{a0}/p(A|mp(v)\A,v)|_{a1}
          # such that we can come back to update the density ratio of v once we update the densratioA
          assign(paste0("bayes.densratio_",v), {p.a0.mpv/p.a1.mpv})

          L.use.densratioA <- c(L.use.densratioA, v)

        }else{

          bayes_fit_v <- CV.SuperLearner(Y=A, X=dat_bayes.v_v, family = binomial(), V = K, SL.library = lib.L, control = list(saveFitLibrary=T),saveAll = T)

          # p(A=1|mp(v)\A)
          p.A1.mpv_v <- bayes_fit_v$SL.predict

          #p(v=a0|mp(v)\A)
          p.a0.mpv_v <- a0*p.A1.mpv_v+(1-a0)*(1-p.A1.mpv_v)

          #p(v=a0|mp(v)\A)
          p.a1.mpv_v <- 1-p.a0.mpv_v

          p.a1.mpv_v[p.a1.mpv_v<zerodiv.avoid] <- zerodiv.avoid # added to avoid INF
          p.a0.mpv_v[p.a0.mpv_v<zerodiv.avoid] <- zerodiv.avoid # added to avoid INF

          # save the density ratio: p(v|mp(V))|_{a0}/p(v|mp(V))|_{a1}
          assign(paste0("densratio_",v), {p.a0.mpv/p.a1.mpv}/{p.a0.mpv_v/p.a1.mpv_v})

        }





      }else if (superlearner.L==T){

        bayes_fit <- SuperLearner(Y=A, X=dat_bayes.v, family = binomial(), SL.library = lib.L)

        # p(A=1|mp(v)\A,v)
        p.A1.mpv <- predict(bayes_fit, type = "response")[[1]] %>% as.vector()  # p(A=1|X)

        #p(v=a0|mp(v)\A,v)
        p.a0.mpv <- a0*p.A1.mpv+(1-a0)*(1-p.A1.mpv)

        #p(v=a0|mp(v)\A,v)
        p.a1.mpv <- 1-p.a0.mpv

        p.a1.mpv[p.a1.mpv<zerodiv.avoid] <- zerodiv.avoid # added to avoid INF
        p.a0.mpv[p.a0.mpv<zerodiv.avoid] <- zerodiv.avoid # added to avoid INF

        if (are_same(setdiff(replace.vector(f.markov_pillow(graph, v, treatment), multivariate.variables) ,treatment), replace.vector(f.markov_pillow(graph,treatment, treatment), multivariate.variables))){

          # save the density ratio: p(v|mp(V))|_{a0}/p(v|mp(V))|_{a1}
          assign(paste0("densratio_",v), {p.a0.mpv/p.a1.mpv}/get(paste0("densratio_",treatment)))

          # save the ratio of p(A|mp(v)\A,v)|_{a0}/p(A|mp(v)\A,v)|_{a1}
          # such that we can come back to update the density ratio of v once we update the densratioA
          assign(paste0("bayes.densratio_",v), {p.a0.mpv/p.a1.mpv})

          L.use.densratioA <- c(L.use.densratioA, v)

        }else{

          bayes_fit_v <- SuperLearner(Y=A, X=dat_bayes.v_v, family = binomial(), SL.library = lib.L)

          # p(A=1|mp(v)\A)
          p.A1.mpv_v <- predict(bayes_fit_v, type = "response")[[1]] %>% as.vector()  # p(A=1|X)

          #p(v=a0|mp(v)\A)
          p.a0.mpv_v <- a0*p.A1.mpv_v+(1-a0)*(1-p.A1.mpv_v)

          #p(v=a0|mp(v)\A)
          p.a1.mpv_v <- 1-p.a0.mpv_v

          p.a1.mpv_v[p.a1.mpv_v<zerodiv.avoid] <- zerodiv.avoid # added to avoid INF
          p.a0.mpv_v[p.a0.mpv_v<zerodiv.avoid] <- zerodiv.avoid # added to avoid INF

          # save the density ratio: p(v|mp(V))|_{a0}/p(v|mp(V))|_{a1}
          assign(paste0("densratio_",v), {p.a0.mpv/p.a1.mpv}/{p.a0.mpv_v/p.a1.mpv_v})

        }

      } else { # Linear model

        # estimate density ratio using bayes rule
        bayes_fit <- glm(A ~ ., data=dat_bayes.v, family = binomial())

        # p(A=1|mp(v)\A,v)
        p.A1.mpv <- predict(bayes_fit, type = "response")

        #p(v=a0|mp(v)\A,v)
        p.a0.mpv <- a0*p.A1.mpv+(1-a0)*(1-p.A1.mpv)

        #p(v=a0|mp(v)\A,v)
        p.a1.mpv <- 1-p.a0.mpv

        p.a1.mpv[p.a1.mpv<zerodiv.avoid] <- zerodiv.avoid # added to avoid INF
        p.a0.mpv[p.a0.mpv<zerodiv.avoid] <- zerodiv.avoid # added to avoid INF


        if (are_same(setdiff(replace.vector(f.markov_pillow(graph, v, treatment), multivariate.variables) ,treatment), replace.vector(f.markov_pillow(graph,treatment, treatment), multivariate.variables))){

          # save the density ratio: p(v|mp(V))|_{a0}/p(v|mp(V))|_{a1}
          assign(paste0("densratio_",v), {p.a0.mpv/p.a1.mpv}/get(paste0("densratio_",treatment)))

          # save the ratio of p(A|mp(v)\A,v)|_{a0}/p(A|mp(v)\A,v)|_{a1}
          # such that we can come back to update the density ratio of v once we update the densratioA
          assign(paste0("bayes.densratio_",v), {p.a0.mpv/p.a1.mpv})

          L.use.densratioA <- c(L.use.densratioA, v)

        }else{

          bayes_fit_v <- glm(A ~ ., data=dat_bayes.v_v, family = binomial())

          # p(A=1|mp(v)\A)
          p.A1.mpv_v <- predict(bayes_fit_v, type = "response")

          #p(v=a0|mp(v)\A)
          p.a0.mpv_v <- a0*p.A1.mpv_v+(1-a0)*(1-p.A1.mpv_v)

          #p(v=a0|mp(v)\A)
          p.a1.mpv_v <- 1-p.a0.mpv_v

          p.a1.mpv_v[p.a1.mpv_v<zerodiv.avoid] <- zerodiv.avoid # added to avoid INF
          p.a0.mpv_v[p.a0.mpv_v<zerodiv.avoid] <- zerodiv.avoid # added to avoid INF

          # save the density ratio: p(v|mp(V))|_{a0}/p(v|mp(V))|_{a1}
          assign(paste0("densratio_",v), {p.a0.mpv/p.a1.mpv}/{p.a0.mpv_v/p.a1.mpv_v})

        }



      }


    } ## Iterate over each variable in L\A


  } else {

    stop("Invalid ratio.method.L input.")

  }

  # Get all the densratio vectors based on their names
  densratio.vectors.L <- mget(c(paste0("densratio_",L.removedA),paste0("densratio_",treatment)))

  # Create a data frame using all the vectors
  densratio.L <- data.frame(densratio.vectors.L)

  ###### DENSITY RATIO ASSOCIATED WITH M ######


  # Find vertices in set M, that the Markov pillow of which include A and don't include A
  # For the vertices whose Markov pillow include A, we will use either the "densratio" or "bayes" method to estimate the density ratio p(M|mp(M))|_{a0}/p(M|mp(M))|_{a1}
  # For the vertices whose Markov pillow don't include A, density ratio p(M|mp(M))|_{a0}/p(M|mp(M))|_{a1} is 1

  M.mpM.includeA <- c()
  M.mpM.excludeA <- c()

  for (v in M){ ## Iterate over each variable in M

    if (treatment %in% f.markov_pillow(graph, v, treatment)){ # vertices whose Markov pillow include A

      M.mpM.includeA <- c(M.mpM.includeA, v)


    } else { # vertices whose Markov pillow don't include A

      M.mpM.excludeA <- c(M.mpM.excludeA, v)

    }

  } ## Iterate over each variable in M

  # assign ratio=1 for vertices in M.mpM.excludeA
  for (m in M.mpM.excludeA) { assign(paste0("densratio_", m), 1) }

  if (ratio.method.M=="densratio"){ ################### METHOD 2A: densratio method  ###################

    # Error: densratio method doesn't support factor variables

    if (!all(sapply(replace.vector(unique(c(M.mpM.includeA, unlist(lapply(1:length(M.mpM.includeA), function(i) f.markov_pillow(graph, M.mpM.includeA[i], treatment))))), multivariate.variables), function(var) is.numeric(data[,var]) | is.integer(data[,var])))){
      stop("Error in estimating density ratios associated with variables in M: densratio method only support numeric/integer variables, try bayes method instead.")
      }


    M.use.densratioA <- c()


    # if M and mpi(M) only consists numeric/integer variables: apply density ratio estimation
    for (v in M.mpM.includeA){

      # used for estimating numerator of the density ratio
      dat_v.num.a0 <- data[data[[treatment]] == a0, replace.vector(c(v, f.markov_pillow(graph, v, treatment)), multivariate.variables)] # select rows where A=a0
      dat_v.num.a1 <- data[data[[treatment]] == a1, replace.vector(c(v, f.markov_pillow(graph, v, treatment)), multivariate.variables)] # select rows where A=a1

      densratio.v.num <- densratio(dat_v.num.a0, dat_v.num.a1)

      ratio.num <- densratio.v.num$compute_density_ratio(data[, replace.vector(c(v, f.markov_pillow(graph, v, treatment)), multivariate.variables)])


      ## Estimator the denomator via bayes rule to avoid zero division
      # used for estimating the denominator of the density ratio
      dat_bayes.v <- data[, setdiff(replace.vector(f.markov_pillow(graph, v, treatment), multivariate.variables), treatment), drop=F]
      names(dat_bayes.v)[names(dat_bayes.v)=="Y"] <- "outcome" # Super Learner get confused of the regresor contains a variable named Y. Thus rename it to outcome

      if (are_same(setdiff(replace.vector(f.markov_pillow(graph, v, treatment), multivariate.variables) ,treatment), replace.vector(f.markov_pillow(graph,treatment, treatment), multivariate.variables))){

        ratio.den <- get(paste0("densratio_",treatment))

        assign(paste0("densratio.densratio_",v), ratio.num)

        M.use.densratioA <- c(M.use.densratioA, v)

      }else{



        if (crossfit==T){

          bayes_fit <- CV.SuperLearner(Y=A, X=dat_bayes.v, family = binomial(), V = K, SL.library = lib.M, control = list(saveFitLibrary=T),saveAll = T)

          # p(A=1|mp(v)\A,v)
          p.A1.mpv <- bayes_fit$SL.predict

          #p(v=a0|mp(v)\A,v)
          p.a0.mpv <- a0*p.A1.mpv+(1-a0)*(1-p.A1.mpv)

          #p(v=a0|mp(v)\A,v)
          p.a1.mpv <- 1-p.a0.mpv

          p.a1.mpv[p.a1.mpv<zerodiv.avoid] <- zerodiv.avoid # added to avoid INF
          p.a0.mpv[p.a0.mpv<zerodiv.avoid] <- zerodiv.avoid # added to avoid INF

          # save the density ratio: p(a0|mp(V))/p(a1|mp(V))
          ratio.den <- p.a0.mpv/p.a1.mpv


        }else if (superlearner.M==T){

          bayes_fit <- SuperLearner(Y=A, X=dat_bayes.v, family = binomial(), SL.library = lib.M)

          # p(A=1|mp(v)\A,v)
          p.A1.mpv <- predict(bayes_fit, type = "response")[[1]] %>% as.vector()  # p(A=1|X)

          #p(v=a0|mp(v)\A,v)
          p.a0.mpv <- a0*p.A1.mpv+(1-a0)*(1-p.A1.mpv)

          #p(v=a0|mp(v)\A,v)
          p.a1.mpv <- 1-p.a0.mpv

          p.a1.mpv[p.a1.mpv<zerodiv.avoid] <- zerodiv.avoid # added to avoid INF
          p.a0.mpv[p.a0.mpv<zerodiv.avoid] <- zerodiv.avoid # added to avoid INF

          # save the density ratio: p(a0|mp(V))/p(a1|mp(V))
          ratio.den <- p.a0.mpv/p.a1.mpv

        } else {

          # estimate density ratio using bayes rule
          bayes_fit <- glm(A ~ ., data=dat_bayes.v, family = binomial())

          # p(A=1|mp(v)\A,v)
          p.A1.mpv <- predict(bayes_fit, type = "response")

          #p(v=a0|mp(v)\A,v)
          p.a0.mpv <- a0*p.A1.mpv+(1-a0)*(1-p.A1.mpv)

          #p(v=a0|mp(v)\A,v)
          p.a1.mpv <- 1-p.a0.mpv

          p.a1.mpv[p.a1.mpv<zerodiv.avoid] <- zerodiv.avoid # added to avoid INF
          p.a0.mpv[p.a0.mpv<zerodiv.avoid] <- zerodiv.avoid # added to avoid INF

          # save the density ratio: p(a0|mp(V))/p(a1|mp(V))
          ratio.den <- p.a0.mpv/p.a1.mpv

        }





      }


      ratio <- ratio.num/ratio.den # p(M|mp(M))|_{a_0}/p(M|mp(M))|_{a_1}

      assign(paste0("densratio_",v), ratio)
    }


  } else if (ratio.method.M=="dnorm"){ ################### METHOD 2B: dnorm method  ###################


    if (!all(sapply(replace.vector(M.mpM.includeA, multivariate.variables), function(var) is.numeric(data[,var]) | is.integer(data[,var]) | length(unique(data[,var]))==2 ))){

      stop("Error in estimating density ratios associated with variables in M: dnorm method only support continuous or binary variables, try bayes method instead.")
      }


    # if M consists of continuous or binary variables: apply density ratio estimation via dnorm
    for (v in M.mpM.includeA){ ## Iterate over each variable in L\A

      ratio <- calculate_density_ratio_dnorm(a0=a0, v , graph, treatment=treatment, data=data, formula=dnorm.formula.M) # p(M|mp(M))|_{a_0}/p(M|mp(M))|_{a_1}

      assign(paste0("densratio_",v), ratio)
    }


  }else if (ratio.method.M=="bayes"){ ################### METHOD 2C: Bayes method ###################

    M.use.densratioA <- c()

    for (v in M.mpM.includeA){ ## Iterate over each variable in M

      #### Prepare data for regression and prediction ####
      dat_bayes.v <- data[,setdiff(replace.vector(c(v, f.markov_pillow(graph, v, treatment)), multivariate.variables) ,treatment), drop=F] # contains variable v + Markov pillow of v - treatment
      names(dat_bayes.v)[names(dat_bayes.v)=="Y"] <- "outcome" # Super Learner get confused of the regresor contains a variable named Y. Thus rename it to outcome

      dat_bayes.v_v <- data[,setdiff(replace.vector(f.markov_pillow(graph, v, treatment), multivariate.variables) ,treatment), drop=F] # contains variable Markov pillow of v - treatment

      #### Fit nuisance models ####

      if (crossfit==T){

        # fit p(A|mp(v)\A,v)
        bayes_fit <- CV.SuperLearner(Y=A, X=dat_bayes.v, family = binomial(), V = K, SL.library = lib.M, control = list(saveFitLibrary=T),saveAll = T)

        # p(A=1|mp(v)\A,v)
        p.A1.mpv <- bayes_fit$SL.predict

        #p(A=a0|mp(v)\A,v)
        p.a0.mpv <- a0*p.A1.mpv+(1-a0)*(1-p.A1.mpv)

        #p(A=a1|mp(v)\A,v)
        p.a1.mpv <- 1-p.a0.mpv

        p.a1.mpv[p.a1.mpv<zerodiv.avoid] <- zerodiv.avoid # added to avoid INF
        p.a0.mpv[p.a0.mpv<zerodiv.avoid] <- zerodiv.avoid # added to avoid INF

        if (are_same(setdiff(replace.vector(f.markov_pillow(graph, v, treatment), multivariate.variables) ,treatment), replace.vector(f.markov_pillow(graph,treatment, treatment), multivariate.variables))){

          # save the density ratio: p(v|mp(V))|_{a0}/p(v|mp(V))|_{a1}
          assign(paste0("densratio_",v), {p.a0.mpv/p.a1.mpv}/get(paste0("densratio_",treatment)))

          # save the ratio of p(A|mp(v)\A,v)|_{a0}/p(A|mp(v)\A,v)|_{a1}
          # such that we can come back to update the density ratio of v once we update the densratioA
          assign(paste0("bayes.densratio_",v), {p.a0.mpv/p.a1.mpv})

          M.use.densratioA <- c(M.use.densratioA, v)

        }else{

          bayes_fit_v <- CV.SuperLearner(Y=A, X=dat_bayes.v_v, family = binomial(), V = K, SL.library = lib.M, control = list(saveFitLibrary=T),saveAll = T)

          # p(A=1|mp(v)\A)
          p.A1.mpv_v <- bayes_fit_v$SL.predict  # p(A=1|X)

          #p(v=a0|mp(v)\A)
          p.a0.mpv_v <- a0*p.A1.mpv_v+(1-a0)*(1-p.A1.mpv_v)

          #p(v=a0|mp(v)\A)
          p.a1.mpv_v <- 1-p.a0.mpv_v

          p.a1.mpv_v[p.a1.mpv_v<zerodiv.avoid] <- zerodiv.avoid # added to avoid INF
          p.a0.mpv_v[p.a0.mpv_v<zerodiv.avoid] <- zerodiv.avoid # added to avoid INF

          # save the density ratio: p(v|mp(V))|_{a0}/p(v|mp(V))|_{a1}
          assign(paste0("densratio_",v), {p.a0.mpv/p.a1.mpv}/{p.a0.mpv_v/p.a1.mpv_v})

        }

      }else if (superlearner.M==T){

        # fit p(A|mp(v)\A,v)
        bayes_fit <- SuperLearner(Y=A, X=dat_bayes.v, family = binomial(), SL.library = lib.M)

        # p(A=1|mp(v)\A,v)
        p.A1.mpv <- predict(bayes_fit, type = "response")[[1]] %>% as.vector()  # p(A=1|X)

        #p(A=a0|mp(v)\A,v)
        p.a0.mpv <- a0*p.A1.mpv+(1-a0)*(1-p.A1.mpv)

        #p(A=a1|mp(v)\A,v)
        p.a1.mpv <- 1-p.a0.mpv

        p.a1.mpv[p.a1.mpv<zerodiv.avoid] <- zerodiv.avoid # added to avoid INF
        p.a0.mpv[p.a0.mpv<zerodiv.avoid] <- zerodiv.avoid # added to avoid INF



        if (are_same(setdiff(replace.vector(f.markov_pillow(graph, v, treatment), multivariate.variables) ,treatment), replace.vector(f.markov_pillow(graph,treatment, treatment), multivariate.variables))){

          # save the density ratio: p(v|mp(V))|_{a0}/p(v|mp(V))|_{a1}
          assign(paste0("densratio_",v), {p.a0.mpv/p.a1.mpv}/get(paste0("densratio_",treatment)))

          # save the ratio of p(A|mp(v)\A,v)|_{a0}/p(A|mp(v)\A,v)|_{a1}
          # such that we can come back to update the density ratio of v once we update the densratioA
          assign(paste0("bayes.densratio_",v), {p.a0.mpv/p.a1.mpv})

          M.use.densratioA <- c(M.use.densratioA, v)

        }else{

          bayes_fit_v <- SuperLearner(Y=A, X=dat_bayes.v_v, family = binomial(), SL.library = lib.M)

          # p(A=1|mp(v)\A)
          p.A1.mpv_v <- predict(bayes_fit_v, type = "response")[[1]] %>% as.vector()  # p(A=1|X)

          #p(v=a0|mp(v)\A)
          p.a0.mpv_v <- a0*p.A1.mpv_v+(1-a0)*(1-p.A1.mpv_v)

          #p(v=a0|mp(v)\A)
          p.a1.mpv_v <- 1-p.a0.mpv_v

          p.a1.mpv_v[p.a1.mpv_v<zerodiv.avoid] <- zerodiv.avoid # added to avoid INF
          p.a0.mpv_v[p.a0.mpv_v<zerodiv.avoid] <- zerodiv.avoid # added to avoid INF

          # save the density ratio: p(v|mp(V))|_{a0}/p(v|mp(V))|_{a1}
          assign(paste0("densratio_",v), {p.a0.mpv/p.a1.mpv}/{p.a0.mpv_v/p.a1.mpv_v})

        }

      } else {

        # estimate density ratio using bayes rule
        bayes_fit <- glm(A ~ ., data=dat_bayes.v, family = binomial())


        # p(A=1|mp(v)\A,v)
        p.A1.mpv <- predict(bayes_fit, type = "response")


        #p(v=a0|mp(v)\A,v)
        p.a0.mpv <- a0*p.A1.mpv+(1-a0)*(1-p.A1.mpv)


        #p(v=a0|mp(v)\A,v)
        p.a1.mpv <- 1-p.a0.mpv

        p.a1.mpv[p.a1.mpv<zerodiv.avoid] <- zerodiv.avoid # added to avoid INF
        p.a0.mpv[p.a0.mpv<zerodiv.avoid] <- zerodiv.avoid # added to avoid INF



        if (are_same(setdiff(replace.vector(f.markov_pillow(graph, v, treatment), multivariate.variables) ,treatment), replace.vector(f.markov_pillow(graph,treatment, treatment), multivariate.variables))){

          # save the density ratio: p(v|mp(V))|_{a0}/p(v|mp(V))|_{a1}
          assign(paste0("densratio_",v), {p.a0.mpv/p.a1.mpv}/get(paste0("densratio_",treatment)))

          # save the ratio of p(A|mp(v)\A,v)|_{a0}/p(A|mp(v)\A,v)|_{a1}
          # such that we can come back to update the density ratio of v once we update the densratioA
          assign(paste0("bayes.densratio_",v), {p.a0.mpv/p.a1.mpv})

          M.use.densratioA <- c(M.use.densratioA, v)

        }else{

          bayes_fit_v <- glm(A ~ ., data=dat_bayes.v_v, family = binomial())

          # p(A=1|mp(v)\A)
          p.A1.mpv_v <- predict(bayes_fit_v, type = "response")

          #p(v=a0|mp(v)\A)
          p.a0.mpv_v <- a0*p.A1.mpv_v+(1-a0)*(1-p.A1.mpv_v)

          #p(v=a0|mp(v)\A)
          p.a1.mpv_v <- 1-p.a0.mpv_v

          p.a1.mpv_v[p.a1.mpv_v<zerodiv.avoid] <- zerodiv.avoid # added to avoid INF
          p.a0.mpv_v[p.a0.mpv_v<zerodiv.avoid] <- zerodiv.avoid # added to avoid INF

          # save the density ratio: p(v|mp(V))|_{a0}/p(v|mp(V))|_{a1}
          assign(paste0("densratio_",v), {p.a0.mpv/p.a1.mpv}/{p.a0.mpv_v/p.a1.mpv_v})

        }





      }


    } ## Iterate over each variable in M


  } else {

    stop("Invalid ratio.method.M input.")

  }

  # Get the densratio vectors based on their names
  densratio.vectors.M <- mget(paste0("densratio_",M))

  # Create a data frame using all the vectors
  densratio.M <- data.frame(densratio.vectors.M)

  # fit sequential regression for v and save the predictions as mu.v_a1 and mu.v_a0 in the parent environment of the function
  outer_env <- environment()

  fit_seq_reg <- function(v) {

    # Markov pillow and data prep
    mpv <- replace.vector(f.markov_pillow(graph, v, treatment), multivariate.variables)
    dat_mpv <- data[, mpv, drop=F]

    next.v <- tau.df$tau[tau.df$order == tau.df$order[tau.df$tau==v] + 1]
    next.mu <- if (next.v %in% L) get(paste0("mu.", next.v, "_a1"), envir=outer_env) else
      get(paste0("mu.", next.v, "_a0"), envir=outer_env)
    next.mu.transform <- if (all(Y %in% c(0,1))) qlogis(next.mu) else next.mu

    if (treatment %in% mpv) {
      dat_mpv.a0 <- dat_mpv %>% mutate(!!treatment := a0)
      dat_mpv.a1 <- dat_mpv %>% mutate(!!treatment := a1)
    }

    # Fit model
    if (crossfit) {
      v_fit <- CV.SuperLearner(Y=next.mu.transform, X=dat_mpv, family=gaussian(),
                               V=K, SL.library=lib.seq,
                               control=list(saveFitLibrary=T), saveAll=T)
      if (treatment %in% mpv) {
        reg_a1 <- unlist(lapply(1:K, function(x) predict(v_fit$AllSL[[x]], newdata=dat_mpv.a1[v_fit$folds[[x]],])[[1]] %>% as.vector()))[order(unlist(lapply(1:K, function(x) v_fit$folds[[x]])))]
        reg_a0 <- unlist(lapply(1:K, function(x) predict(v_fit$AllSL[[x]], newdata=dat_mpv.a0[v_fit$folds[[x]],])[[1]] %>% as.vector()))[order(unlist(lapply(1:K, function(x) v_fit$folds[[x]])))]
      } else {
        reg_a1 <- reg_a0 <- v_fit$SL.predict
      }

    } else if (superlearner.seq) {
      v_fit <- SuperLearner(Y=next.mu.transform, X=dat_mpv, family=gaussian(), SL.library=lib.seq)
      if (treatment %in% mpv) {
        reg_a1 <- predict(v_fit, newdata=dat_mpv.a1)[[1]] %>% as.vector()
        reg_a0 <- predict(v_fit, newdata=dat_mpv.a0)[[1]] %>% as.vector()
      } else {
        reg_a1 <- reg_a0 <- predict(v_fit)[[1]] %>% as.vector()
      }

    } else {
      v_fit <- lm(next.mu.transform ~ ., data=dat_mpv)
      if (treatment %in% mpv) {
        reg_a1 <- predict(v_fit, newdata=dat_mpv.a1)
        reg_a0 <- predict(v_fit, newdata=dat_mpv.a0)
      } else {
        reg_a1 <- reg_a0 <- predict(v_fit)
      }
    }

    # Assign results — back-transform if Y is binary
    assign(paste0("mu.", v, "_a1"),
           if (all(Y %in% c(0,1))) plogis(reg_a1) else reg_a1,
           envir=outer_env)
    assign(paste0("mu.", v, "_a0"),
           if (all(Y %in% c(0,1))) plogis(reg_a0) else reg_a0,
           envir=outer_env)
  }





  #///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////#
  #/////////////////////////////////////////// STEP2: One-step estimator /////////////////////////////////////////////////////#
  #///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////#

  ############################################################
  ## SEQUENTIAL REGRESSION FOR ESETIMATING mu(Z_i, a_{Z_i})
  ############################################################

  # the sequential regression will be perform for v that locates between A and Y accroding to topological order tau
  vertices.between.AY <- tau[{which(tau==treatment)+1}:{which(tau==outcome)-1}]

  # sequential regression
  for (v in rev(vertices.between.AY)) fit_seq_reg(v)

  ############################################################
  ## END OF SEQUENTIAL REGRESSION FOR ESETIMATING mu(Z_i, a_{Z_i})
  ############################################################

  ######################
  # EIF calculations
  ######################

  ## EIF for Y|mp(Y):
  #if Y in L: I(A=a1)*I(M < Y){p(M|mp(M))|_{a0}/p(M|mp(M))|_{a1}}*(Y-E(Y|mp(Y))|_{a1})
  #if Y in M: I(A=a0)*I(L < Y){p(L|mp(L))|_{a0}/p(:|mp(L))|_{a1}}*(Y-E(Y|mp(Y))|_{a0})

  # density ratio before Y
  selected.M <- M[sapply(M, function(m) tau.df$order[tau.df$tau==m] < tau.df$order[tau.df$tau==outcome])] # M precedes Y
  selected.L <- L[sapply(L, function(l) tau.df$order[tau.df$tau==l] < tau.df$order[tau.df$tau==outcome])] # L precedes Y

  f.M_preY <- Reduce(`*`, densratio.M[,paste0("densratio_",selected.M), drop=F]) # Mi precede Y = the set M
  f.L_preY <- Reduce(`*`, densratio.L[,paste0("densratio_",selected.L), drop=F]) # Li precede Y = the set L


  EIF.Y <- if(outcome %in% L){(A==a1)*f.M_preY*(Y-mu.Y_a1)}else{(A==a0)*1/f.L_preY*(Y-mu.Y_a0)}

  ## EIF for Z|mp(Z)

  for (v in rev(vertices.between.AY)){ ## iterate over all vertices between A and Y

    # select M and L that precede v
    selected.M <- M[sapply(M, function(m) tau.df$order[tau.df$tau==m] < tau.df$order[tau.df$tau==v])] # M precedes Z
    selected.L <- L[sapply(L, function(l) tau.df$order[tau.df$tau==l] < tau.df$order[tau.df$tau==v])] # L precedes Z

    # vertex that right after v according to tau
    next.v <- tau.df$tau[tau.df$order=={tau.df$order[tau.df$tau==v]+1}]

    # mu(next.v, a_{next.v})
    next.mu <- if(next.v %in% L){ get(paste0("mu.",next.v,"_a1"))}else{get(paste0("mu.",next.v,"_a0"))}

    EIF.v <- if(v %in% L){ # v in L

      # product of the selected variables density ratio
      f.M_prev <- Reduce(`*`, densratio.M[,paste0("densratio_",selected.M), drop=F]) # Mi precede v

      # EIF for v|mp(v)
      (A==a1)*f.M_prev*( next.mu - get(paste0("mu.",v,"_a1")) )

    }else{ # v in M

      # product of the selected variables density ratio
      f.L_prev <- Reduce(`*`, densratio.L[,paste0("densratio_",selected.L), drop=F]) # Li precede v

      # EIF for v|mp(v)
      (A==a0)*1/f.L_prev*( next.mu - get(paste0("mu.",v,"_a0")) )
    }

    assign(paste0("EIF.",v), EIF.v)

  } ## End of iteration over all vertices between A and Y

  ## EIF for A=a1|mp(A)

  # vertex that right after A according to tau
  next.A <- tau.df$tau[tau.df$order=={tau.df$order[tau.df$tau==treatment]+1}]

  EIF.A <- {(A==a1) - p.a1.mpA}*get(paste0("mu.",next.A,"_a0"))


  ######################
  # estimate E[Y(a)]
  ######################

  # estimated psi
  estimated_psi = mean( EIF.Y + rowSums(as.data.frame(mget(paste0("EIF.",vertices.between.AY)))) +  EIF.A + p.a1.mpA*get(paste0("mu.",next.A,"_a0")) + (A==a0)*Y )

  # EIF
  EIF <- EIF.Y + # EIF of Y|mp(Y)
    rowSums(as.data.frame(mget(paste0("EIF.",vertices.between.AY)))) + # EIF of v|mp(v) for v between A and Y
    EIF.A + # EIF of A|mp(A)
    p.a1.mpA*get(paste0("mu.",next.A,"_a0")) + # EIF of mp(A)
    (A==a0)*Y -
    estimated_psi


  # confidence interval
  lower.ci <- estimated_psi-1.96*sqrt(mean(EIF^2)/nrow(data))
  upper.ci <- estimated_psi+1.96*sqrt(mean(EIF^2)/nrow(data))


  onestep.out <- list(estimated_psi=estimated_psi, # estimated parameter
                      lower.ci=lower.ci, # lower bound of 95% CI
                      upper.ci=upper.ci, # upper bound of 95% CI
                      EIF=EIF, # E(Dstar) for Y|M,A,X and M|A,X, and A|X
                      EIF.Y=EIF.Y, # EIF of Y|mp(Y)
                      EIF.A=EIF.A, # EIF of A|mp(A)
                      EIF.v = rowSums(as.data.frame(mget(paste0("EIF.",vertices.between.AY)))), # EIF of v|mp(v) for v between A and Y
                      p.a1.mpA = p.a1.mpA, # estimated E[A=a1|mp(A)]
                      mu.next.A = get(paste0("mu.",next.A,"_a0")) # estimated E[v|mp(v)] for v that comes right after A
  )









  #///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////#
  #/////////////////////////////////////////// STEP2: Sequential regression based TMLE ///////////////////////////////////////#
  #///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////#

  # initialize EDstar: the mean of the EIF for A, Y, and v
  EDstar <- 10 # random large number

  # record EDstar over iterations
  EDstar.record <- c()

  # initialize iteration counter
  iter <- 0



  while(abs(EDstar) > cvg.criteria & iter < n.iter){


    ######################
    # update p(A=a1|mp(A))
    ######################

    # clever coefficient for propensity score: mu(Z_1,a0)
    clevercoef.A <- get(paste0("mu.",next.A,"_a0"))

    # derive epsA
    ind <- A==a1

    ps_model <- glm(
      ind ~ offset(qlogis(p.a1.mpA))+ clevercoef.A -1, family=binomial(), start=0
    )

    eps.A <- coef(ps_model)

    # update p(A=a1|mp(A))
    p.a1.mpA <- plogis(qlogis(p.a1.mpA)+eps.A*(clevercoef.A))

    p.a0.mpA <-1-p.a1.mpA # p(A=a0|mp(A))

    # avoid zero indivision error
    p.a0.mpA[p.a0.mpA<zerodiv.avoid] <- zerodiv.avoid
    p.a1.mpA[p.a1.mpA<zerodiv.avoid] <- zerodiv.avoid

    # update density ratio of A
    assign(paste0("densratio_",treatment), p.a0.mpA/p.a1.mpA) # density ratio regarding the treatment p(A|mp(A))|_{a_0}/p(A|mp(A))|_{a_1}


    # update density ratio for v in L if the ratio.method.L = "bayes" OR "densratio"
    if (ratio.method.L == "bayes"){ for (v in L.use.densratioA){assign(paste0("densratio_",v), get(paste0("bayes.densratio_",v))/get(paste0("densratio_",treatment)))} }
    if (ratio.method.L == "densratio"){ for (v in L.use.densratioA){assign(paste0("densratio_",v), 1/(get(paste0("densratio.densratio_",v))*get(paste0("densratio_",treatment))))} }


    # Update the density ratio data frame for L
    densratio.vectors.L <- mget(paste0("densratio_",L))
    densratio.L <- data.frame(densratio.vectors.L)

    # update density ratio for v in M if the ratio.method.M = "bayes" OR "densratio"
    if (ratio.method.M == "bayes"){ for (v in M.use.densratioA){assign(paste0("densratio_",v), get(paste0("bayes.densratio_",v))/get(paste0("densratio_",treatment)))} }
    if (ratio.method.M == "densratio"){ for (v in M.use.densratioA){assign(paste0("densratio_",v), get(paste0("densratio.densratio_",v))/get(paste0("densratio_",treatment))) } }

      # Update the density ratio data frame for M
      densratio.vectors.M <- mget(paste0("densratio_",M))
      densratio.M <- data.frame(densratio.vectors.M)



      ######################
      # update E[Y|mp(Y)]
      ######################

      # density ratio before Y
      selected.M <- M[sapply(M, function(m) tau.df$order[tau.df$tau==m] < tau.df$order[tau.df$tau==outcome])] # M precedes Y
      selected.L <- L[sapply(L, function(l) tau.df$order[tau.df$tau==l] < tau.df$order[tau.df$tau==outcome])] # L precedes Y

      f.M_preY <- Reduce(`*`, densratio.M[,paste0("densratio_",selected.M), drop=F]) # Mi precede Y = the set M
      f.L_preY <- Reduce(`*`, densratio.L[,paste0("densratio_",selected.L), drop=F]) # Li precede Y = the set L

      # clever coefficient for outcome regression: E(Y|mp(Y))
      weight.Y <- if(outcome %in% L){(A==a1)*f.M_preY}else{(A==a0)*1/f.L_preY}

      offset.Y <- if(outcome %in% L){mu.Y_a1}else{mu.Y_a0}

      if (all(Y %in% c(0,1))){ # binary Y

        # one iteration
        or_model <- glm(
          Y ~ offset(offset.Y)+weight.Y-1, family=binomial(), start=0
        )

        eps.Y = coef(or_model)

        # updated outcome regression
        # E[Y|mp(Y)]|_{a1} is updated if Y in L
        # E[Y|mp(Y)]|_{a0} is updated if Y in M
        if(outcome %in% L){mu.Y_a1 <- plogis(qlogis(offset.Y)+eps.Y*weight.Y)}else{mu.Y_a0 <- plogis(qlogis(offset.Y)+eps.Y*weight.Y)}


      } else { # continuous Y

        # one iteration
        or_model <- lm(
          Y ~ offset(offset.Y)+1, weights = weight.Y
        )

        eps.Y <- coef(or_model)

        # updated outcome regression
        # E[Y|mp(Y)]|_{a1} is updated if Y in L
        # E[Y|mp(Y)]|_{a0} is updated if Y in M
        if(outcome %in% L){mu.Y_a1 <- offset.Y+eps.Y}else{mu.Y_a0 <- offset.Y+eps.Y}
      }

      # update EIF of Y
      EIF.Y <- if(outcome %in% L){(A==a1)*f.M_preY*(Y-mu.Y_a1)}else{(A==a0)*1/f.L_preY*(Y-mu.Y_a0)} # if Y in M


      ######################
      # update mu(v,a_v)
      ######################

      for (v in rev(vertices.between.AY)){ ## iterative over the vertices between A and Y to update mu(v,a_v)

        # mu(v,a_v) estimates
        fit_seq_reg(v)

        ## Perform TMLE update for mu(v,a_v)
        # offset
        if (all(Y %in% c(0,1))){
          # L(mu(v)) = mu(next.v)*log(mu(v)(eps.v))+(1-mu(next.v))*log(1-mu(v)(eps.v)), where mu(v)(eps.v)=expit{logit(mu(v))+eps.v*weight}
          offset.v <- if(v %in% L){qlogis(get(paste0("mu.",v,"_a1")))}else{qlogis(get(paste0("mu.",v,"_a0")))}

        }else{
          # L(mu(v)) = weight*{mu(next.v)-mu(v)(eps.v)}^2, where mu(v)(eps.v)=mu(v)+eps.v
          offset.v <- if(v %in% L){get(paste0("mu.",v,"_a1"))}else{get(paste0("mu.",v,"_a0"))}
          }

        # weight for regression
        # select M and L that precede v
        selected.M <- M[sapply(M, function(m) tau.df$order[tau.df$tau==m] < tau.df$order[tau.df$tau==v])] # M precedes Z
        selected.L <- L[sapply(L, function(l) tau.df$order[tau.df$tau==l] < tau.df$order[tau.df$tau==v])] # L precedes Z


        weight.v <- if(v %in% L){
          # product of the selected variables density ratio
          f.M_prev <- Reduce(`*`, densratio.M[,paste0("densratio_",selected.M), drop=F]) # Mi precede Z

          # weight
          (A==a1)*f.M_prev
          }else{
          # product of the selected variables density ratio
          f.L_prev <- Reduce(`*`, densratio.L[,paste0("densratio_",selected.L),drop=F]) # Mi precede Z

          # weight
          (A==a0)*1/f.L_prev
          }


        # compute next.mu for the fluctuation step
        next.v  <- tau.df$tau[tau.df$order == tau.df$order[tau.df$tau==v] + 1]
        next.mu <- if (next.v %in% L) get(paste0("mu.", next.v, "_a1")) else
          get(paste0("mu.", next.v, "_a0"))

        # one iteration update
        if (all(Y %in% c(0,1))){
          v_model <- glm(next.mu ~ offset(offset.v)+weight.v-1, family=binomial(link="logit")) # fit logistic regression if Y is binary, this is like a generalized logistic regression

        }else{
          v_model <- lm(next.mu ~ offset(offset.v)+1, weights = weight.v) # fit linear regression if Y is continuous
        }

        # the optimized submodel index
        eps.v <- coef(v_model)

        # updated outcome regression
        # E[Y|mp(Y)]|_{a1} is updated if v in L
        # E[Y|mp(Y)]|_{a0} is updated if v in M
        if(v %in% L){

          assign(paste0("mu.",v,"_a1"), if(all(Y %in% c(0,1))){plogis(offset.v+eps.v*weight.v)}else{offset.v+eps.v}) # transform back to probability scale if Y is binary

        }else{

          assign(paste0("mu.",v,"_a0"), if(all(Y %in% c(0,1))){plogis(offset.v+eps.v*weight.v)}else{offset.v+eps.v}) # transform back to probability scale if Y is binary

        }

        EIF.v <- if(v %in% L){ weight.v*( next.mu - get(paste0("mu.",v,"_a1"))) }else{ weight.v*( next.mu - get(paste0("mu.",v,"_a0")))}

        assign(paste0("EIF.",v), EIF.v) # if Y in M

      } ## End of update of mu(v,a_v) and EIF.v


      # update EIF for A
      EIF.A <- ((A==a1) - p.a1.mpA)*get(paste0("mu.",next.A,"_a0"))

      # update stoping criteria
      EDstar <- mean(EIF.A) + mean(EIF.Y) + mean(rowSums(as.data.frame(mget(paste0("EIF.",vertices.between.AY)))))


      # update iteration counter
      iter <- iter + 1

      # record EDstar
      EDstar.record <- c(EDstar.record, EDstar)



    } ## End of while loop



  ######################
  # estimate E[Y(a)]
  ######################

  # estimated psi
  estimated_psi = mean( EIF.Y + rowSums(as.data.frame(mget(paste0("EIF.",vertices.between.AY)))) +  EIF.A + p.a1.mpA*get(paste0("mu.",next.A,"_a0")) + (A==a0)*Y )

  # EIF
  EIF <- EIF.Y + # EIF of Y|mp(Y)
    rowSums(as.data.frame(mget(paste0("EIF.",vertices.between.AY)))) + # EIF of v|mp(v) for v between A and Y
    EIF.A + # EIF of A|mp(A)
    p.a1.mpA*get(paste0("mu.",next.A,"_a0")) + # EIF of mp(A)
    (A==a0)*Y -
    estimated_psi


  # confidence interval
  lower.ci <- estimated_psi-1.96*sqrt(mean(EIF^2)/nrow(data))
  upper.ci <- estimated_psi+1.96*sqrt(mean(EIF^2)/nrow(data))


  tmle.out <- list(estimated_psi=estimated_psi, # estimated parameter
                   lower.ci=lower.ci, # lower bound of 95% CI
                   upper.ci=upper.ci, # upper bound of 95% CI
                   EIF=EIF, # EIF
                   EIF.Y=EIF.Y, # EIF of Y|mp(Y)
                   EIF.A=EIF.A, # EIF of A|mp(A)
                   EIF.v = rowSums(as.data.frame(mget(paste0("EIF.",vertices.between.AY)))), # EIF of v|mp(v) for v between A and Y
                   p.a1.mpA = p.a1.mpA, # estimated E[A=a1|mp(A)]
                   mu.next.A = get(paste0("mu.",next.A,"_a0")), # estimated E[v|mp(v)] for v that comes right after A
                   EDstar = EDstar, # stopping criteria, the mean of EIF of A, Y, and v between A and Y
                   iter = iter, # number of iterations to achieve convergence
                   EDstar.record = EDstar.record,
                   densratio.M=densratio.M,
                   densratio.L=densratio.L,
                   mu.Y_a1=mu.Y_a1,
                   mu.Y_a0=mu.Y_a0) # record of EDstar



  return(list(Onestep=onestep.out, TMLE=tmle.out))

}

