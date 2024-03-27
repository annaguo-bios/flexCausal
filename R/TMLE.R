## usethis namespace: start
#' @importFrom dplyr %>%
#' @importFrom SuperLearner SuperLearner
#' @importFrom np npcdensbw npcdens
#' @importFrom densratio densratio
#' @importFrom stats coef cov lm glm as.formula binomial gaussian integrate optimize plogis predict qlogis
#' @importFrom mvtnorm dmvnorm
## usethis namespace: end
NULL

## The main function that execute all TMLE estimators ====
#' Function for average counterfactual outcome E(Y(a)) and average causal effect (ACE) estimation. Most of the arguments are passed directly to TMLE.all.
#' @param a Treatment level at which the average counterfactual outcome is computed
#' @param data A Dataframe contains treatment, mediators, outcome, and measured confounders
#' @param treatment Variable name for the unvariate binary treatment
#' @param mediators Variable name for the continuous univariate mediator
#' @param outcome Variable name for the continuous univariate outcome
#' @param covariates Variable name for the measured confounders
#' @param onestep A logical indicator determines whether one-step estimation is executed. When 'onestep=T', the one-step estimation result is provided. Conversely, if 'onestep=F', the result is withheld.
#' @param mediator.method When M is univariate binary, regression is adopted for estimating \eqn{p(M|A,X)}. Otherwise, four methods for mediator density estimation is provided, namely "bayes", "densratio", "dnorm", and "np".
#' The "bayes" method estimates the density ratio \eqn{p(M|A,X)/p(M|a,X)} by rewriting it as \eqn{(p(a|M,X)/p(A|M,X))/(p(a|X)/p(A|X))}, where p(A|M,X) is then estimated via regression. TMLE estimators using 'bayes' method avoids updating the mediator density p(M|A,X).
#' The "densratio" method estimates \eqn{p(M|A,X)/p(M|a,X)} by rewriting it as \eqn{(p(M,a,X)/p(M,A,X))/(p(a|X)/p(A|X))}, where \eqn{p(M,A,X)/p(M,a,X)} is estimated with the \link[densratio]{densratio} function,
#' and \eqn{p(a|X)/p(A|X)} is estimated via regression. "densratio" method is only applicable if all mediators are numeric or integer valued. TMLE estimators using 'densratio' method avoids updating the mediator density p(M|A,X).
#' The "dnorm" method estimates the mediator density ratio \eqn{p(M|A,X)/p(M|a,X)} assuming p(M|A,X) follows a normal distribution. TMLE estimators using 'dnorm' method avoids updating the mediator density p(M|A,X).
#' TMLE estimators using 'np' method directly update the mediator density p(M|A,X). When np.dnorm=F, the "np" method estimates \eqn{p(M|a,X)/p(M|A,X)} by estimating \eqn{p(M|A,X)} with the \link[np]{npcdens} function.
#' When np.dnorm=T, the "np" method estimates \eqn{p(M|a,X)/p(M|A,X)} assuming normal distribution. Due to the computational burden, "np" method is only available for univariate continuous mediator.
#' @param superlearner A logical indicator determines whether SuperLearner via the \link[SuperLearner]{SuperLearner} function is adopted for estimating the outcome regression, mediator density, and the propensity score.
#' SuperLearner is the preferred option in cases where complex relationships among variables exist, potentially leading to model misspecification issues when using simple linear models.
#' @param crossfit A logical indicator determines whether SuperLearner+Cross-fitting is adopted for estimating the outcome regression, mediator density, and the propensity score.
#' @param K A integer indicating the number of folds for cross-fitting, the default is 5.
#' @param lib Library of algorithms for SuperLearner.
#' @param n.iter The maximum number of iterations performed when iteratively updating the mediator density and propensity score.
#' @param eps A logical indicator determines the stopping criteria used when iteratively updating the mediator density and propensity score. The default is 'eps=T'.
#' When 'eps=T', \eqn{\sqrt{\epsilon_2^2+\epsilon_3^2}} is used, where \eqn{\epsilon_2} and \eqn{\epsilon_3} are the index of the sub-models for mediator density and propensity score. When 'eps=F', \eqn{max(|\Phi_M|,|\Phi_A|)} is used,
#' where \eqn{\Phi_M} and \eqn{\Phi_A} are the mapping of efficient influcence function (EIF) into the tangent space of \eqn{M|A,X} and \eqn{A|X}. In general, adoption of 'eps=F' results in better convergence while it takes longer time.
#' Conversely, adoption of 'eps=T' usually requires less time but can result in algorithm divergence.
#' @param cvg.criteria A numerical value representing the convergence criteria when iteratively updating the mediator density and propensity score.
#'  The default value is 0.01, meaning update stops when stopping criteria < 0.01. The stopping criteria is chosen by the eps.
#' @param formulaY Regression formula for the outcome regression of Y on M, A, X. The default is 'Y ~ 1+ M + A + X'.
#' @param linkY_binary The link function used for the logistic regression of Y on M, A, X when Y is binary. The default is the 'logit' link.
#' @param formulaA Regression formula for the propensity score regression of A on X. The default is 'A ~ 1 + X'.
#' @param linkA The link function used for the logistic regression of A on X. The default is the 'logit' link.
#' @param formulaM Regression formula for the mediator density regression of M on A and X. The default is 'M ~ 1 + A + X'. This parameter is only needed when M is a univariate binary mediator.
#' @param linkM_binary The link function used for the logistic regression of M on A and X. The default is the 'logit' link. This parameter is only needed when M is a univariate binary mediator.
#' @param formula_bayes Regression formula for the regression of A on M and X. The default is 'A ~ 1 + M + X'. This parameter is only needed when mediator.method="bayes".
#' @param link_bayes The link function used for the logistic regression of A on M and X. The default is 'logit' link. This parameter is only needed when mediator.method="bayes".
#' @param truncate_lower A numeric variable, setting lower bound for the truncated propensity score. The default is 0.
#' @param truncate_upper A numeric variable, setting upper bound for the truncated propensity score. The default is 1.
#' @param np.dnorm A logic variable. If np.dnorm=T, p(M|A,X) is directly estimated assuming normal distribution. If np.dnorm=F, p(M|A,X) is directly estimated using the \link[np]{npcdens} function.
#' @return Function outputs a list containing TMLE results (and Onestep results if 'onestep=T' is specified). When 'a=c(1,0)', function also outputs corresponding results on \eqn{E(Y^1)} and \eqn{E(Y^1)}:
#' \describe{
#'       \item{\code{ATE}}{The estimated Average Causal Effect: \eqn{E(Y^1)-E(Y^0)}}
#'       \item{\code{estimated_psi}}{The estimated parameter of interest: \eqn{E(Y^a)}}
#'       \item{\code{lower.ci}}{Lower bound of the 95% confidence interval for \code{ATE} or \code{estimated_psi}}
#'       \item{\code{upper.ci}}{Upper bound of the 95% confidence interval for \code{ATE} or \code{estimated_psi}}
#'       \item{\code{theta_x}}{\eqn{\int E(Y|M,A,X)p(M|A=a,X)p(A|X) dM dA}}
#'       \item{\code{p.m1.aX}}{\eqn{\int p(M=1|A=a,X)}}
#'       \item{\code{p.a1.X}}{\eqn{p(A=1|X)}}
#'       \item{\code{or_pred}}{\eqn{E(Y|M,A,X)}}
#'       \item{\code{EIF}}{The estimated efficient influence function evaluated at the observed data}
#'       \item{\code{EDstar}}{A vector of the mapping of \code{EIF} in the tangent space of \eqn{Y|M,A,X}; \eqn{M|A,X}; and \eqn{A|X}.}
#'       \item{\code{EDstar_M.vec}}{A vector containing the average value of the mapping of EIF in tangent space \eqn{M|A,X} over iterations. This is useful for checking the convergence behavior of the mediator density. It's expected to be close to 0 when convergence is achieved.}
#'       \item{\code{EDstar_ps.vec}}{A vector containing the average value of the mapping of EIF in tangent space \eqn{A|X} over iterations. This is useful for checking the convergence behavior of the propensity score. It's expected to be close to 0 when convergence is achieved.}
#'       \item{\code{eps2_vec}}{A vector containing the index for submodels of the mediator density over iterations. This is useful for checking the convergence behavior of the mediator density. It's expected to be close to 0 when convergence is achieved.}
#'       \item{\code{eps3_vec}}{A vector containing the index for submodels of the propensity score over iterations. This is useful for checking the convergence behavior of the propensity score. It's expected to be close to 0 when convergence is achieved.}
#'       \item{\code{iter}}{Number of iterations where convergence is achieved for the iterative update of the mediator density and propensity score.}}
#' @examples
#' # ATE estimation. For binary outcome Y and binary mediator M.
#' # Data is generated under p(A=1|X) = 0.3 + 0.2X. Therefore, setting link for propensity score to be "identity".
#' TMLE(a=c(1,0),data=binaryY_binaryM,treatment="A", mediators="M", outcome="Y", covariates="X", onestep=T, linkA="identity")
#'
#' # ATE estimation. For continuous outcome Y and binary mediator M
#' # Data is generated under p(A=1|X) = 0.3 + 0.2X. Therefore, setting link for propensity score to be "identity".
#' TMLE(a=c(1,0),data=continuousY_binaryM,treatment="A", mediators="M", outcome="Y", covariates="X", onestep=T, linkA="identity")
#'
#' # ATE estimation. For continuous outcome Y and binary mediator M.
#' # Data is generated under p(A=1|X) = 0.001 + 0.998X. And X~Uniform(0,1). Therefore, this dataset suffers from weak overlapping. Below we apply truncation to the propensity score to truncate it between (0.001, 0.999).
#' TMLE(a=c(1,0),data=continuousY_binaryM_weakoverlap, treatment="A", mediators="M", outcome="Y", covariates="X", onestep=T, linkA="identity", truncate_lower=0.001, truncate_upper=0.999)
#'
#' # ATE estimation. For continuous outcome Y and univariate continuous mediator M. Using 'np' method for mediator density estimation. Setting np.dnorm=F, so that mediator density is estimated via the np function.
#' # Data is generated under p(A=1|X) = 0.3 + 0.2X. Therefore, setting link for propensity score to be "identity".
#' TMLE(a=c(1,0),data=continuousY_continuousM, treatment="A", mediators="M", outcome="Y", covariates="X", onestep=T, linkA="identity", mediator.method="np", np.dnorm=F)
#'
#' # ATE estimation. For continuous outcome Y and univariate continuous mediator M. Using 'np' method for mediator density estimation. Setting np.dnorm=T, so that mediator density is estimated assuming normal distribution.
#' # Data is generated under p(A=1|X) = 0.3 + 0.2X. Therefore, setting link for propensity score to be "identity".
#' TMLE(a=c(1,0),data=continuousY_continuousM, treatment="A", mediators="M", outcome="Y", covariates="X", onestep=T, linkA="identity", mediator.method="np", np.dnorm=T)
#'
#' # ATE estimation. For continuous outcome Y and univariate continuous mediator M. Using 'densratio' method for mediator density ratio estimation.
#' # Data is generated under p(A=1|X) = 0.001 + 0.998X. And X~Uniform(0,1). Therefore, this dataset suffers from weak overlapping. Below we apply truncation to the propensity score to truncate it between (0.001, 0.999).
#' TMLE(a=c(1,0),data=continuousY_continuousM_weakoverlap, treatment="A", mediators="M", outcome="Y", covariates="X", onestep=T, linkA="identity", mediator.method="densratio", truncate_lower=0.001, truncate_upper=0.999)
#'
#' # ATE estimation. For continuous outcome Y and bivariate mediator M. Using 'densratio' method for mediator density ratio estimation.
#' # Data is generated under p(A=1|X) = 0.3 + 0.2X. Therefore, setting link for propensity score to be "identity".
#' TMLE(a=c(1,0),data=continuousY_bivariateM,treatment="A", mediators=c("M.2","M.2"), outcome="Y", covariates="X", onestep=T, linkA="identity", mediator.method="densratio")
#'
#' # ATE estimation. For continuous outcome Y and bivariate mediator M. Using 'bayes' method for mediator density ratio estimation.
#' # Data is generated under p(A=1|X) = 0.3 + 0.2X. Therefore, setting link for propensity score to be "identity".
#' TMLE(a=c(1,0),data=continuousY_bivariateM,treatment="A", mediators=c("M.2","M.2"), outcome="Y", covariates="X", onestep=T, linkA="identity", mediator.method="bayes")
#'
#' # ATE estimation. For continuous outcome Y and bivariate mediator M. Using 'dnorm' method for mediator density ratio estimation.
#' # Data is generated under p(A=1|X) = 0.3 + 0.2X. Therefore, setting link for propensity score to be "identity".
#' TMLE(a=c(1,0),data=continuousY_bivariateM,treatment="A", mediators=c("M.2","M.2"), outcome="Y", covariates="X", onestep=T, linkA="identity", mediator.method="dnorm")
#'
#' # ATE estimation. For continuous outcome Y and quadrivariate mediator M. Using 'densratio' method for mediator density ratio estimation.
#' # Data is generated under p(A=1|X) = 0.3 + 0.2X. Therefore, setting link for propensity score to be "identity".
#' TMLE(a=c(1,0),data=continuousY_bivariateM,treatment="A", mediators=c("M.2","M.2","M.3","M.4"), outcome="Y", covariates="X", onestep=T, linkA="identity", mediator.method="densratio")
#'
#' # ATE estimation. For continuous outcome Y and quadrivariate mediator M. Using 'bayes' method for mediator density ratio estimation.
#' # Data is generated under p(A=1|X) = 0.3 + 0.2X. Therefore, setting link for propensity score to be "identity".
#' TMLE(a=c(1,0),data=continuousY_bivariateM,treatment="A", mediators=c("M.2","M.2","M.3","M.4"), outcome="Y", covariates="X", onestep=T, linkA="identity", mediator.method="bayes")
#'
#' # ATE estimation. For continuous outcome Y and quadrivariate mediator M. Using 'dnorm' method for mediator density ratio estimation.
#' # Data is generated under p(A=1|X) = 0.3 + 0.2X. Therefore, setting link for propensity score to be "identity".
#' TMLE(a=c(1,0),data=continuousY_bivariateM,treatment="A", mediators=c("M.2","M.2","M.3","M.4"), outcome="Y", covariates="X", onestep=T, linkA="identity", mediator.method="dnorm")
#'
#' @import dplyr MASS densratio SuperLearner stats
#' @export
#'
#'
TMLE <- function(a,data,treatment, mediators, outcome, covariates,
                 onestep=T, mediator.method="bayes", np.dnorm=T, superlearner=F,crossfit=F,K=5,
                 lib = c("SL.glm","SL.earth","SL.ranger","SL.mean"), n.iter=500, eps=T, cvg.criteria=0.01,
                 formulaY="Y ~ .", formulaA="A ~ .", formulaM="M~.", linkY_binary="logit", linkA="logit", linkM_binary="logit",
                 formula_bayes="A ~ .",link_bayes="logit",
                 truncate_lower=0, truncate_upper=1){
#
#   # sample size
#
#   n <- nrow(data)
#
#   if (is.vector(a) & length(a)>2){ ## Invalid input ==
#
#     print("Invalid input. Enter a=c(1,0) for Average Causal Effect estimation. Enter a=1 or a=0 for average counterfactual outcome estimation at the specified treatment level.")
#
#   }else if (is.vector(a) & length(a)==2){ ## ATE estimate ==
#
#     ## TMLE estimator
#
#     out.a1 <- TMLE.all(a=a[1],data=data,treatment=treatment, mediators=mediators, outcome=outcome, covariates=covariates,
#                        onestep=onestep, mediator.method=mediator.method, superlearner=superlearner,crossfit=crossfit,K=K,
#                        lib = lib, n.iter=n.iter, eps=eps, cvg.criteria=cvg.criteria,
#                        formulaY=formulaY, formulaA=formulaA, formulaM=formulaM, linkY_binary=linkY_binary, linkA=linkA, linkM_binary=linkM_binary,
#                        formula_bayes=formula_bayes,link_bayes=link_bayes,
#                        truncate_lower=truncate_lower, truncate_upper=truncate_upper)
#
#     out.a0 <- TMLE.all(a=a[2],data=data,treatment=treatment, mediators=mediators, outcome=outcome, covariates=covariates,
#                        onestep=onestep, mediator.method=mediator.method, superlearner=superlearner,crossfit=crossfit,K=K,
#                        lib = lib, n.iter=n.iter, eps=eps, cvg.criteria=cvg.criteria,
#                        formulaY=formulaY, formulaA=formulaA, formulaM=formulaM, linkY_binary=linkY_binary, linkA=linkA, linkM_binary=linkM_binary,
#                        formula_bayes=formula_bayes,link_bayes=link_bayes,
#                        truncate_lower=truncate_lower, truncate_upper=truncate_upper)
#
#     # run TMLE
#     tmle_output_Y1 <- out.a1$TMLE
#     tmle_output_Y0 <- out.a0$TMLE
#
#     # estimate E[Y(1)], E[Y(0)], and ATE
#     hat_E.Y1 = tmle_output_Y1$estimated_psi
#     hat_E.Y0 = tmle_output_Y0$estimated_psi
#     hat_ATE = hat_E.Y1 - hat_E.Y0
#
#     # lower CI
#     lower.ci_ATE = hat_ATE - 1.96*sqrt(mean((tmle_output_Y1$EIF-tmle_output_Y0$EIF)^2)/n)
#
#     # upper CI
#     upper.ci_ATE = hat_ATE + 1.96*sqrt(mean((tmle_output_Y1$EIF-tmle_output_Y0$EIF)^2)/n)
#
#     tmle.out <- list(ATE=hat_ATE, # estimated parameter
#                      lower.ci=lower.ci_ATE, # lower bound of 95% CI
#                      upper.ci=upper.ci_ATE, # upper bound of 95% CI
#                      EIF=tmle_output_Y1$EIF-tmle_output_Y0$EIF # EIF
#                      )
#
#     if (onestep==T){
#
#       # run TMLE
#       onestep_output_Y1 <- out.a1$Onestep
#       onestep_output_Y0 <- out.a0$Onestep
#
#       # estimate E[Y(1)], E[Y(0)], and ATE
#       hat_E.Y1 = onestep_output_Y1$estimated_psi
#       hat_E.Y0 = onestep_output_Y0$estimated_psi
#       hat_ATE = hat_E.Y1 - hat_E.Y0
#
#       # lower CI
#       lower.ci_ATE = hat_ATE - 1.96*sqrt(mean((onestep_output_Y1$EIF-onestep_output_Y0$EIF)^2)/n)
#
#       # upper CI
#       upper.ci_ATE = hat_ATE + 1.96*sqrt(mean((onestep_output_Y1$EIF-onestep_output_Y0$EIF)^2)/n)
#
#       onestep.out <- list(ATE=hat_ATE, # estimated parameter
#                        lower.ci=lower.ci_ATE, # lower bound of 95% CI
#                        upper.ci=upper.ci_ATE, # upper bound of 95% CI
#                        EIF=onestep_output_Y1$EIF-onestep_output_Y0$EIF # EIF
#       )
#
#       cat(paste0("TMLE estimated ACE: ",round(tmle.out$ATE,2),"; 95% CI: (",round(tmle.out$lower.ci,2),", ",round(tmle.out$upper.ci,2),") \n","Onestep estimated ACE: ",round(onestep.out$ATE,2),"; 95% CI: (",round(onestep.out$lower.ci,2),", ",round(onestep.out$upper.ci,2),")"))
#
#       return(list(TMLE=tmle.out,Onestep=onestep.out, TMLE.Y1=tmle_output_Y1, TMLE.Y0 = tmle_output_Y0, Onestep.Y1=onestep_output_Y1, Onestep.Y0=onestep_output_Y0))
#
#     }else {
#
#       cat(paste0("TMLE estimated ATE: ",round(tmle.out$ATE,2),"; 95% CI: (",round(tmle.out$lower.ci,2),", ",round(tmle.out$upper.ci,2),")"))
#
#       return(list(TMLE=tmle.out,TMLE.Y1=tmle_output_Y1, TMLE.Y0 = tmle_output_Y0))
#
#       }
#
#   }else if (a==1) { ## E(Y^1) estimate ==
#
#     out.a1 <- TMLE.all(a=1,data=data,treatment=treatment, mediators=mediators, outcome=outcome, covariates=covariates,
#                        onestep=onestep, mediator.method=mediator.method, superlearner=superlearner,crossfit=crossfit,K=K,
#                        lib = lib, n.iter=n.iter, eps=eps, cvg.criteria=cvg.criteria,
#                        formulaY=formulaY, formulaA=formulaA, formulaM=formulaM, linkY_binary=linkY_binary, linkA=linkA, linkM_binary=linkM_binary,
#                        formula_bayes=formula_bayes,link_bayes=link_bayes,
#                        truncate_lower=truncate_lower, truncate_upper=truncate_upper)
#     return(out.a1)
#
#   }else if (a==0){ ## E(Y^0) estimate ==
#
#     out.a0 <- TMLE.all(a=0,data=data,treatment=treatment, mediators=mediators, outcome=outcome, covariates=covariates,
#                        onestep=onestep, mediator.method=mediator.method, superlearner=superlearner,crossfit=crossfit,K=K,
#                        lib = lib, n.iter=n.iter, eps=eps, cvg.criteria=cvg.criteria,
#                        formulaY=formulaY, formulaA=formulaA, formulaM=formulaM, linkY_binary=linkY_binary, linkA=linkA, linkM_binary=linkM_binary,
#                        formula_bayes=formula_bayes,link_bayes=link_bayes,
#                        truncate_lower=truncate_lower, truncate_upper=truncate_upper)
#
#     return(out.a0)
#
#   }

}


