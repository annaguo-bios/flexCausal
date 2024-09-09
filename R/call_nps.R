.call_nps <- function(a, data, vertices,
                      di_edges, bi_edges, treatment, outcome,
                      multivariate.variables, graph,
                      superlearner.seq, superlearner.Y, superlearner.A, superlearner.M, superlearner.L,
                      crossfit, K,
                      ratio.method.L, ratio.method.M,
                      dnorm.formula.L, dnorm.formula.M,
                      lib.seq, lib.L, lib.M, lib.Y, lib.A,
                      formulaY, formulaA,
                      linkY_binary, linkA,
                      n.iter, cvg.criteria,
                      truncate_lower, truncate_upper, zerodiv.avoid
){

  n <- nrow(data)

  ############################################################
  # Perform estimation with the NPS method
  ############################################################

  if (is.vector(a) & length(a)>2){ ## Invalid input ==

    stop("Invalid input. Enter a=c(vaule1,value2) for average causal effect estimation: (Y(a=value1)) - E(Y(a=value2)). Enter a=value1 for average counterfactual outcome estimation at the specified treatment level value1.")

  }else if (is.vector(a) & length(a)==2){ ## ATE estimate ==

    ## TMLE estimator

    out.a1 <- NPS.TMLE.a(a = a[1], data = data, vertices = vertices, di_edges = di_edges, bi_edges = bi_edges, treatment = treatment, outcome = outcome, multivariate.variables = multivariate.variables, graph = graph,

                         superlearner.seq = superlearner.seq, # whether run superlearner for sequential regression
                         superlearner.Y = superlearner.Y, # whether run superlearner for outcome regression
                         superlearner.A = superlearner.A, # whether run superlearner for propensity score
                         superlearner.M = superlearner.M, # whether run superlearner for estimating densratio for M using bayes method
                         superlearner.L = superlearner.L, # whether run superlearner for estimating densratio for L using bayes method

                         crossfit = crossfit, K = K,

                         ratio.method.L = ratio.method.L, # method for estimating the density ratio associated with M
                         ratio.method.M = ratio.method.M, # method for estimating the density ratio associated with L

                         dnorm.formula.L = dnorm.formula.L, # formula for estimating the density ratio associated with M
                         dnorm.formula.M = dnorm.formula.M, # formula for estimating the density ratio associated with L

                         lib.seq = lib.seq, # superlearner library for sequential regression
                         lib.L = lib.L, # superlearner library for density ratio estimation via bayes rule for variables in L
                         lib.M = lib.M, # superlearner library for density ratio estimation via bayes rule for variables in M
                         lib.Y = lib.Y, # superlearner library for outcome regression
                         lib.A = lib.A, # superlearner library for propensity score

                         formulaY = formulaY, formulaA = formulaA, # regression formula for outcome regression and propensity score if superlearner is not used
                         linkY_binary = linkY_binary, linkA = linkA, # link function for outcome regression and propensity score if superlearner is not used

                         n.iter = n.iter, cvg.criteria = cvg.criteria,
                         truncate_lower = truncate_lower, truncate_upper = truncate_upper,zerodiv.avoid=zerodiv.avoid)

    out.a0 <- NPS.TMLE.a(a = a[2], data = data, vertices = vertices, di_edges = di_edges, bi_edges = bi_edges, treatment = treatment, outcome = outcome, multivariate.variables = multivariate.variables, graph = graph,

                         superlearner.seq = superlearner.seq, # whether run superlearner for sequential regression
                         superlearner.Y = superlearner.Y, # whether run superlearner for outcome regression
                         superlearner.A = superlearner.A, # whether run superlearner for propensity score
                         superlearner.M = superlearner.M, # whether run superlearner for estimating densratio for M using bayes method
                         superlearner.L = superlearner.L, # whether run superlearner for estimating densratio for L using bayes method

                         crossfit = crossfit, K = K,

                         ratio.method.L = ratio.method.L, # method for estimating the density ratio associated with M
                         ratio.method.M = ratio.method.M, # method for estimating the density ratio associated with L

                         dnorm.formula.L = dnorm.formula.L, # formula for estimating the density ratio associated with M
                         dnorm.formula.M = dnorm.formula.M, # formula for estimating the density ratio associated with L

                         lib.seq = lib.seq, # superlearner library for sequential regression
                         lib.L = lib.L, # superlearner library for density ratio estimation via bayes rule for variables in L
                         lib.M = lib.M, # superlearner library for density ratio estimation via bayes rule for variables in M
                         lib.Y = lib.Y, # superlearner library for outcome regression
                         lib.A = lib.A, # superlearner library for propensity score

                         formulaY = formulaY, formulaA = formulaA, # regression formula for outcome regression and propensity score if superlearner is not used
                         linkY_binary = linkY_binary, linkA = linkA, # link function for outcome regression and propensity score if superlearner is not used

                         n.iter = n.iter, cvg.criteria = cvg.criteria,
                         truncate_lower = truncate_lower, truncate_upper = truncate_upper,zerodiv.avoid=zerodiv.avoid)

    ############################ TMLE ############################
    # run TMLE
    tmle_output_Y1 <- out.a1$TMLE
    tmle_output_Y0 <- out.a0$TMLE

    # estimate E[Y(1)], E[Y(0)], and ATE
    hat_E.Y1 = tmle_output_Y1$estimated_psi
    hat_E.Y0 = tmle_output_Y0$estimated_psi
    hat_ATE = hat_E.Y1 - hat_E.Y0

    # lower CI
    lower.ci_ATE = hat_ATE - 1.96*sqrt(mean((tmle_output_Y1$EIF-tmle_output_Y0$EIF)^2)/n)

    # upper CI
    upper.ci_ATE = hat_ATE + 1.96*sqrt(mean((tmle_output_Y1$EIF-tmle_output_Y0$EIF)^2)/n)

    tmle.out <- list(ATE=hat_ATE, # estimated parameter
                     lower.ci=lower.ci_ATE, # lower bound of 95% CI
                     upper.ci=upper.ci_ATE, # upper bound of 95% CI
                     EIF=tmle_output_Y1$EIF-tmle_output_Y0$EIF # EIF
    )


    ############################ onestep ############################
    # run onestep
    onestep_output_Y1 <- out.a1$Onestep
    onestep_output_Y0 <- out.a0$Onestep

    # estimate E[Y(1)], E[Y(0)], and ATE
    hat_E.Y1 = onestep_output_Y1$estimated_psi
    hat_E.Y0 = onestep_output_Y0$estimated_psi
    hat_ATE = hat_E.Y1 - hat_E.Y0

    # lower CI
    lower.ci_ATE = hat_ATE - 1.96*sqrt(mean((onestep_output_Y1$EIF-onestep_output_Y0$EIF)^2)/n)

    # upper CI
    upper.ci_ATE = hat_ATE + 1.96*sqrt(mean((onestep_output_Y1$EIF-onestep_output_Y0$EIF)^2)/n)

    onestep.out <- list(ATE=hat_ATE, # estimated parameter
                        lower.ci=lower.ci_ATE, # lower bound of 95% CI
                        upper.ci=upper.ci_ATE, # upper bound of 95% CI
                        EIF=onestep_output_Y1$EIF-onestep_output_Y0$EIF # EIF
    )

    cat(paste0("TMLE estimated ACE: ",round(tmle.out$ATE,2),"; 95% CI: (",round(tmle.out$lower.ci,2),", ",round(tmle.out$upper.ci,2),") \n","Onestep estimated ACE: ",round(onestep.out$ATE,2),"; 95% CI: (",round(onestep.out$lower.ci,2),", ",round(onestep.out$upper.ci,2),")"))

    np.out <- list(TMLE=tmle.out,Onestep=onestep.out, TMLE.Y1=tmle_output_Y1, TMLE.Y0 = tmle_output_Y0, Onestep.Y1=onestep_output_Y1, Onestep.Y0=onestep_output_Y0)


  }else if (length(a)==1) { ## E(Y^1) estimate ==

    out.a <- NPS.TMLE.a(a = a, data = data, vertices = vertices,
                        di_edges = di_edges, bi_edges = bi_edges, treatment = treatment, outcome = outcome,
                        multivariate.variables = multivariate.variables, graph = graph,

                        superlearner.seq = superlearner.seq, # whether run superlearner for sequential regression
                        superlearner.Y = superlearner.Y, # whether run superlearner for outcome regression
                        superlearner.A = superlearner.A, # whether run superlearner for propensity score
                        superlearner.M = superlearner.M, # whether run superlearner for estimating densratio for M using bayes method
                        superlearner.L = superlearner.L, # whether run superlearner for estimating densratio for L using bayes method

                        crossfit = crossfit, K = K,

                        ratio.method.L = ratio.method.L, # method for estimating the density ratio associated with M
                        ratio.method.M = ratio.method.M, # method for estimating the density ratio associated with L

                        lib.seq = lib.seq, # superlearner library for sequential regression
                        lib.L = lib.L, # superlearner library for density ratio estimation via bayes rule for variables in L
                        lib.M = lib.M, # superlearner library for density ratio estimation via bayes rule for variables in M
                        lib.Y = lib.Y, # superlearner library for outcome regression
                        lib.A = lib.A, # superlearner library for propensity score

                        formulaY = formulaY, formulaA = formulaA, # regression formula for outcome regression and propensity score if superlearner is not used
                        linkY_binary = linkY_binary, linkA = linkA, # link function for outcome regression and propensity score if superlearner is not used

                        n.iter = n.iter, cvg.criteria = cvg.criteria,
                        truncate_lower = truncate_lower, truncate_upper = truncate_upper, zerodiv.avoid=zerodiv.avoid)
    np.out <- out.a

  } # end of if else condition for testing the length of a

  return(np.out)
}
