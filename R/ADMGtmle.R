#' Estimate the average causal effect (ACE) or the average counterfactual outcome E(Y(a)) from observational data under a
#' DAG with hidden variables.
#'
#' The main user-facing function of the package. Given a causal graph specified
#' as an acyclic directed mixed graph (ADMG), this function automatically
#' determines the identifiability status of the treatment effect and dispatches
#' to the appropriate estimator:
#' \itemize{
#'   \item If the treatment is \emph{fixable} (i.e., backdoor-adjustable),
#'     estimation proceeds via \code{.call_backdoor}, returning G-computation,
#'     IPW, one-step (AIPW), and TMLE estimators.
#'   \item If the treatment is \emph{primal fixable} (extended front-door
#'     functional), estimation proceeds via \code{.call_nps}, returning
#'     one-step and TMLE estimators.
#'   \item If the treatment is neither fixable nor primal fixable, the function
#'     stops with an error.
#' }
#' A message is also printed indicating whether the graph is nonparametrically
#' saturated, in which case the returned estimators are semiparametrically
#' efficient.
#'
#' @param a Numeric scalar or length-two numeric vector specifying the treatment
#'   level(s) of interest. The treatment must be coded as 0/1. If a scalar,
#'   the function returns \eqn{E\{Y(a)\}}. If a length-two vector
#'   \code{c(a1, a0)}, the function returns the contrast
#'   \eqn{E\{Y(a1)\} - E\{Y(a0)\}}.
#' @param data A data frame containing all variables listed in \code{vertices}.
#' @param vertices A character vector of variable names in the causal graph.
#'   Ignored if \code{graph} is provided.
#' @param di_edges A list of length-two character vectors specifying directed
#'   edges. For example, \code{list(c('A', 'B'))} encodes A -> B.
#'   Ignored if \code{graph} is provided.
#' @param bi_edges A list of length-two character vectors specifying bidirected
#'   edges. For example, \code{list(c('A', 'B'))} encodes A <-> B.
#'   Ignored if \code{graph} is provided.
#' @param multivariate.variables A named list mapping compound vertex names to
#'   their column names in \code{data}. For example,
#'   \code{list(M = c('M.1', 'M.2'))} indicates M is bivariate with columns
#'   M.1 and M.2. Ignored if \code{graph} is provided.
#' @param graph A graph object created by \code{\link{make.graph}}. If
#'   supplied, \code{vertices}, \code{di_edges}, \code{bi_edges}, and
#'   \code{multivariate.variables} are ignored.
#' @param treatment A character string naming the binary (0/1) treatment
#'   variable in \code{data}.
#' @param outcome A character string naming the outcome variable in \code{data}.
#' @param superlearner.seq Logical. If \code{TRUE}, SuperLearner is used for
#'   sequential regression of intermediate variables (primal fixable case only).
#'   Default is \code{FALSE}.
#' @param superlearner.Y Logical. If \code{TRUE}, SuperLearner is used for
#'   outcome regression. Default is \code{FALSE}.
#' @param superlearner.A Logical. If \code{TRUE}, SuperLearner is used for
#'   propensity score estimation. Default is \code{FALSE}.
#' @param superlearner.M Logical. If \code{TRUE}, SuperLearner is used for
#'   density ratio estimation for variables in M via the Bayes method
#'   (primal fixable case only). Default is \code{FALSE}.
#' @param superlearner.L Logical. If \code{TRUE}, SuperLearner is used for
#'   density ratio estimation for variables in L via the Bayes method
#'   (primal fixable case only). Default is \code{FALSE}.
#' @param crossfit Logical. If \code{TRUE}, cross-fitting with \code{K} folds
#'   is applied to all SuperLearner fits. Default is \code{FALSE}.
#' @param K A positive integer specifying the number of cross-fitting folds.
#'   Used only when \code{crossfit = TRUE}. Default is 5.
#' @param ratio.method.L A character string specifying the method for
#'   estimating density ratios for variables in L (primal fixable case only).
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
#'   estimating density ratios for variables in M (primal fixable case only).
#'   Same options as \code{ratio.method.L}. Default is \code{"bayes"}.
#' @param dnorm.formula.L An optional named list of regression formulas for
#'   variables in L, used when \code{ratio.method.L = "dnorm"}. Names are
#'   variable names; values are formula strings. Variables omitted from the
#'   list are regressed on all Markov pillow variables. For multivariate L,
#'   specify one formula per component, e.g.
#'   \code{list(L.1 = "L.1 ~ A + X", L.2 = "L.2 ~ A + X + I(M^2)")}.
#' @param dnorm.formula.M An optional named list of regression formulas for
#'   variables in M, used when \code{ratio.method.M = "dnorm"}. Same structure
#'   as \code{dnorm.formula.L}.
#' @param lib.seq SuperLearner library for sequential regression.
#'   Default is \code{c("SL.glm", "SL.earth", "SL.ranger", "SL.mean")}.
#' @param lib.L SuperLearner library for density ratio estimation for L.
#'   Default is \code{c("SL.glm", "SL.earth", "SL.ranger", "SL.mean")}.
#' @param lib.M SuperLearner library for density ratio estimation for M.
#'   Default is \code{c("SL.glm", "SL.earth", "SL.ranger", "SL.mean")}.
#' @param lib.Y SuperLearner library for outcome regression.
#'   Default is \code{c("SL.glm", "SL.earth", "SL.ranger", "SL.mean")}.
#' @param lib.A SuperLearner library for propensity score estimation.
#'   Default is \code{c("SL.glm", "SL.earth", "SL.ranger", "SL.mean")}.
#' @param formulaY A formula or character string for outcome regression of Y on
#'   its Markov pillow. Used only when \code{superlearner.Y = FALSE}.
#'   Default is \code{"Y ~ ."}.
#' @param formulaA A formula or character string for propensity score regression
#'   of A on its Markov pillow. Used only when \code{superlearner.A = FALSE}.
#'   Default is \code{"A ~ ."}.
#' @param linkY_binary A character string specifying the link function for
#'   outcome regression when Y is binary and \code{superlearner.Y = FALSE}.
#'   Default is \code{"logit"}.
#' @param linkA A character string specifying the link function for propensity
#'   score regression when \code{superlearner.A = FALSE}.
#'   Default is \code{"logit"}.
#' @param cvg.criteria Numeric. TMLE convergence threshold. The iterative
#'   update stops when \eqn{|\text{mean}(D^*)| <} \code{cvg.criteria}.
#'   Default is \code{0.01}.
#' @param n.iter Maximum number of TMLE iterations. Default is 500.
#' @param truncate_lower Numeric. Propensity score values below this threshold
#'   are clipped. Default is \code{0} (no clipping).
#' @param truncate_upper Numeric. Propensity score values above this threshold
#'   are clipped. Default is \code{1} (no clipping).
#' @param zerodiv.avoid Numeric. Density ratio or propensity score values below
#'   this threshold are clipped to prevent division by zero.
#'   Default is \code{0} (no clipping).
#'
#' @return The return structure depends on the identifiability path:
#' \describe{
#'   \item{Fixable (backdoor)}{A named list with components \code{TMLE},
#'     \code{Onestep}, \code{IPW}, and \code{Gcomp}, plus per-treatment-level
#'     sub-lists.}
#'   \item{Primal fixable (front-door)}{A named list with components
#'     \code{TMLE} and \code{Onestep}.}
#' }
#'
#' @examples
#' # Fixable graph: simple backdoor adjustment
#' test <- estADMG(
#'   a = 1,
#'   data = data_backdoor,
#'   vertices = c('A', 'Y', 'X'),
#'   di_edges = list(c('X', 'A'), c('X', 'Y'), c('A', 'Y')),
#'   treatment = 'A',
#'   outcome = 'Y'
#' )
#'
#' # Primal fixable graph: extended front-door functional
#' test <- estADMG(
#'   a = 1,
#'   data = data_example_a,
#'   vertices = c('A', 'M', 'L', 'Y', 'X'),
#'   bi_edges = list(c('A', 'Y')),
#'   di_edges = list(c('X', 'A'), c('X', 'M'), c('X', 'L'),
#'                   c('X', 'Y'), c('M', 'Y'), c('A', 'M'),
#'                   c('A', 'L'), c('M', 'L'), c('L', 'Y')),
#'   treatment = 'A',
#'   outcome = 'Y',
#'   multivariate.variables = list(M = c('M.1', 'M.2'))
#' )
#'
#' # ACE estimation E(Y(1)) - E(Y(0))
#' test <- estADMG(
#'   a = c(1, 0),
#'   data = data_example_a,
#'   vertices = c('A', 'M', 'L', 'Y', 'X'),
#'   bi_edges = list(c('A', 'Y')),
#'   di_edges = list(c('X', 'A'), c('X', 'M'), c('X', 'L'),
#'                   c('X', 'Y'), c('M', 'Y'), c('A', 'M'),
#'                   c('A', 'L'), c('M', 'L'), c('L', 'Y')),
#'   treatment = 'A',
#'   outcome = 'Y',
#'   multivariate.variables = list(M = c('M.1', 'M.2'))
#' )
#' @importFrom rlang :=
#' @importFrom dplyr %>% mutate select
#' @importFrom MASS mvrnorm
#' @importFrom SuperLearner CV.SuperLearner SuperLearner
#' @importFrom mvtnorm dmvnorm
#' @importFrom densratio densratio
#' @importFrom utils combn
#' @importFrom stats rnorm runif rbinom dnorm dbinom binomial gaussian predict glm as.formula qlogis plogis lm coef cov sd
#' @export
#'
#'
estADMG <- function(a=NULL,data=NULL,vertices=NULL, di_edges=NULL, bi_edges=NULL, treatment=NULL, outcome=NULL, multivariate.variables=NULL, graph=NULL,
                 superlearner.seq = F, # whether run superlearner for sequential regression
                 superlearner.Y=F, # whether run superlearner for outcome regression
                 superlearner.A=F, # whether run superlearner for propensity score
                 superlearner.M=F, # whether run superlearner for estimating densratio for M using bayes method
                 superlearner.L=F, # whether run superlearner for estimating densratio for L using bayes method
                 crossfit=F, K=5,
                 ratio.method.L="bayes", # method for estimating the density ratio associated with L
                 ratio.method.M="bayes", # method for estimating the density ratio associated with M
                 dnorm.formula.L=NULL, # formula for the density ratio associated with L
                 dnorm.formula.M=NULL, # formula for the density ratio associated with M
                 lib.seq = c("SL.glm","SL.earth","SL.ranger","SL.mean"), # superlearner library for sequential regression
                 lib.L = c("SL.glm","SL.earth","SL.ranger","SL.mean"), # superlearner library for density ratio estimation via bayes rule for variables in L
                 lib.M = c("SL.glm","SL.earth","SL.ranger","SL.mean"), # superlearner library for density ratio estimation via bayes rule for variables in M
                 lib.Y = c("SL.glm","SL.earth","SL.ranger","SL.mean"), # superlearner library for outcome regression
                 lib.A = c("SL.glm","SL.earth","SL.ranger","SL.mean"), # superlearner library for propensity score
                 formulaY="Y ~ .", formulaA="A ~ .", # regression formula for outcome regression and propensity score if superlearner is not used
                 linkY_binary="logit", linkA="logit", # link function for outcome regression and propensity score if superlearner is not used
                 n.iter=500, cvg.criteria=0.01,
                 truncate_lower=0, truncate_upper=1, zerodiv.avoid=0){


  # make a graph object if it's not provided
  if (is.null(graph)){ graph <- make.graph(vertices=vertices, bi_edges=bi_edges, di_edges=di_edges, multivariate.variables=multivariate.variables)}


  #####################################################
  # Fixable vs Primal fixable
  #####################################################

  ## if the graph is fixable==> backdoor adjustment functional
  if (suppressMessages(is.fix(graph, treatment))){
    message("The treatment is fixable. Estimation provided via backdoor adjustment.\n")

    np.out <- .call_backdoor(a = a, data = data, vertices = vertices,
                             di_edges = di_edges, bi_edges = bi_edges, treatment = treatment, outcome = outcome,
                             multivariate.variables = multivariate.variables, graph = graph,

                             superlearner.Y = superlearner.Y, # whether run superlearner for outcome regression
                             superlearner.A = superlearner.A, # whether run superlearner for propensity score

                             crossfit = crossfit, K = K,

                             lib.Y = lib.Y, # superlearner library for outcome regression
                             lib.A = lib.A, # superlearner library for propensity score

                             formulaY = formulaY, formulaA = formulaA, # regression formula for outcome regression and propensity score if superlearner is not used
                             linkY_binary = linkY_binary, linkA = linkA, # link function for outcome regression and propensity score if superlearner is not used

                             truncate_lower = truncate_lower, truncate_upper = truncate_upper)

    return(np.out)

  }

  ## if the graph is not fixable nor primal fixable==> stop with an error
  if (!suppressMessages(is.p.fix(graph, treatment))){
    stop("The treatment is not fixable nor primal-fixable. The treatment effect may or may not be identified. Further investigation is needed.\n")
  }

  ## if the graph is not fixable but primal fixable==> extended front-door functional
  message("The treatment is not fixable but is primal fixable. Estimation provided via extended front-door functional.")
  np.out <- .call_nps(a = a, data = data, vertices = vertices,
                      di_edges = di_edges, bi_edges = bi_edges, treatment = treatment, outcome = outcome,
                      multivariate.variables = multivariate.variables, graph = graph,

                      superlearner.seq = superlearner.seq, # whether run superlearner for sequential regression
                      superlearner.Y = superlearner.Y, # whether run superlearner for outcome regression
                      superlearner.A = superlearner.A, # whether run superlearner for propensity score
                      superlearner.M = superlearner.M, # whether run superlearner for estimating densratio for M using bayes method
                      superlearner.L = superlearner.L, # whether run superlearner for estimating densratio for L using bayes method

                      crossfit = crossfit, K = K,

                      ratio.method.L = ratio.method.L, # method for estimating the density ratio associated with L
                      ratio.method.M = ratio.method.M, # method for estimating the density ratio associated with M

                      dnorm.formula.L = dnorm.formula.L, # formula for the density ratio associated with L
                      dnorm.formula.M = dnorm.formula.M, # formula for the density ratio associated with M

                      lib.seq = lib.seq, # superlearner library for sequential regression
                      lib.L = lib.L, # superlearner library for density ratio estimation via bayes rule for variables in L
                      lib.M = lib.M, # superlearner library for density ratio estimation via bayes rule for variables in M
                      lib.Y = lib.Y, # superlearner library for outcome regression
                      lib.A = lib.A, # superlearner library for propensity score

                      formulaY = formulaY, formulaA = formulaA, # regression formula for outcome regression and propensity score if superlearner is not used
                      linkY_binary = linkY_binary, linkA = linkA, # link function for outcome regression and propensity score if superlearner is not used

                      n.iter = n.iter, cvg.criteria = cvg.criteria,
                      truncate_lower = truncate_lower, truncate_upper = truncate_upper, zerodiv.avoid=zerodiv.avoid)


  #####################################################
  # NPS and semi-parametric model
  #####################################################

  if (suppressMessages(is.np.saturated(graph))){
    message("The graph is nonparametrically saturated. Results from the one-step estimator and TMLE are provided, which are in theory the most efficient estimators.\n")

  }else{
    message("The graph is NOT nonparametrically saturated. Note that there may be more efficient estimators.\n")

  }

  return(np.out)

} # end of estADMG function


