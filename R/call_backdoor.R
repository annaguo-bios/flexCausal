.call_backdoor <- function(a, data, vertices,
                           di_edges, bi_edges, treatment, outcome,
                           multivariate.variables, graph,
                           superlearner.Y, superlearner.A,
                           crossfit, K,
                           lib.Y, lib.A,
                           formulaY, formulaA,
                           linkY_binary, linkA,
                           truncate_lower, truncate_upper){

  n <- nrow(data)

  if (is.vector(a) & length(a)>2){ ## Invalid input ==

    print("Invalid input. Enter a=c(vaule1,value2) for average causal effect estimation: (Y(a=value1)) - E(Y(a=value2)). Enter a=value1 for average counterfactual outcome estimation at the specified treatment level value1.")

  }else if (is.vector(a) & length(a)==2){ ## ATE estimate ==

    ## TMLE estimator

    out.a1 <- backdoor.TMLE.a(a = a[1], data = data, vertices = vertices, di_edges = di_edges, bi_edges = bi_edges, treatment = treatment, outcome = outcome, multivariate.variables = multivariate.variables, graph = graph,

                              superlearner.Y = superlearner.Y, # whether run superlearner for outcome regression
                              superlearner.A = superlearner.A, # whether run superlearner for propensity score

                              crossfit = crossfit, K = K,

                              lib.Y = lib.Y, # superlearner library for outcome regression
                              lib.A = lib.A, # superlearner library for propensity score

                              formulaY = formulaY, formulaA = formulaA, # regression formula for outcome regression and propensity score if superlearner is not used
                              linkY_binary = linkY_binary, linkA = linkA, # link function for outcome regression and propensity score if superlearner is not used

                              truncate_lower = truncate_lower, truncate_upper = truncate_upper)

    out.a0 <- backdoor.TMLE.a(a = a[2], data = data, vertices = vertices, di_edges = di_edges, bi_edges = bi_edges, treatment = treatment, outcome = outcome, multivariate.variables = multivariate.variables, graph = graph,

                              superlearner.Y = superlearner.Y, # whether run superlearner for outcome regression
                              superlearner.A = superlearner.A, # whether run superlearner for propensity score

                              crossfit = crossfit, K = K,

                              lib.Y = lib.Y, # superlearner library for outcome regression
                              lib.A = lib.A, # superlearner library for propensity score

                              formulaY = formulaY, formulaA = formulaA, # regression formula for outcome regression and propensity score if superlearner is not used
                              linkY_binary = linkY_binary, linkA = linkA, # link function for outcome regression and propensity score if superlearner is not used

                              truncate_lower = truncate_lower, truncate_upper = truncate_upper)

    ############################ aipw ############################
    # run aipw
    aipw_output_Y1 <- out.a1$Onestep
    aipw_output_Y0 <- out.a0$Onestep

    # estimate E[Y(1)], E[Y(0)], and ATE
    hat_E.Y1 = aipw_output_Y1$estimated_psi
    hat_E.Y0 = aipw_output_Y0$estimated_psi
    hat_ATE = hat_E.Y1 - hat_E.Y0

    # lower CI
    lower.ci_ATE = hat_ATE - 1.96*sqrt(mean((aipw_output_Y1$EIF-aipw_output_Y0$EIF)^2)/n)

    # upper CI
    upper.ci_ATE = hat_ATE + 1.96*sqrt(mean((aipw_output_Y1$EIF-aipw_output_Y0$EIF)^2)/n)

    aipw.out <- list(ATE=hat_ATE, # estimated parameter
                     lower.ci=lower.ci_ATE, # lower bound of 95% CI
                     upper.ci=upper.ci_ATE, # upper bound of 95% CI
                     EIF=aipw_output_Y1$EIF-aipw_output_Y0$EIF # EIF
    )



    ############################ tmle ############################
    # run aipw
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


    ############################ gcomp ############################
    # run gcomp
    gcomp_output_Y1 <- out.a1$gcomp
    gcomp_output_Y0 <- out.a0$gcomp

    # estimate E[Y(1)], E[Y(0)], and ATE
    hat_E.Y1 = gcomp_output_Y1$estimated_psi
    hat_E.Y0 = gcomp_output_Y0$estimated_psi
    hat_ATE = hat_E.Y1 - hat_E.Y0

    # # lower CI
    # lower.ci_ATE = hat_ATE - 1.96*sqrt(mean((gcomp_output_Y1$EIF-gcomp_output_Y0$EIF)^2)/n)
    #
    # # upper CI
    # upper.ci_ATE = hat_ATE + 1.96*sqrt(mean((gcomp_output_Y1$EIF-gcomp_output_Y0$EIF)^2)/n)

    gcomp.out <- list(ATE=hat_ATE, # estimated parameter
                      # lower.ci=lower.ci_ATE, # lower bound of 95% CI
                      # upper.ci=upper.ci_ATE, # upper bound of 95% CI
                      EIF=gcomp_output_Y1$EIF-gcomp_output_Y0$EIF # EIF
    )

    ############################ ipw ############################
    # run ipw
    ipw_output_Y1 <- out.a1$ipw
    ipw_output_Y0 <- out.a0$ipw

    # estimate E[Y(1)], E[Y(0)], and ATE
    hat_E.Y1 = ipw_output_Y1$estimated_psi
    hat_E.Y0 = ipw_output_Y0$estimated_psi
    hat_ATE = hat_E.Y1 - hat_E.Y0

    # # lower CI
    # lower.ci_ATE = hat_ATE - 1.96*sqrt(mean((ipw_output_Y1$EIF-ipw_output_Y0$EIF)^2)/n)
    #
    # # upper CI
    # upper.ci_ATE = hat_ATE + 1.96*sqrt(mean((ipw_output_Y1$EIF-ipw_output_Y0$EIF)^2)/n)

    ipw.out <- list(ATE=hat_ATE, # estimated parameter
                    # lower.ci=lower.ci_ATE, # lower bound of 95% CI
                    # upper.ci=upper.ci_ATE, # upper bound of 95% CI
                    EIF=ipw_output_Y1$EIF-ipw_output_Y0$EIF # EIF
    )

    cat(paste0("Onestep estimated ACE: ",round(aipw.out$ATE,2),"; 95% CI: (",round(aipw.out$lower.ci,2),", ",round(aipw.out$upper.ci,2),") \n",
               "TMLE estimated ACE: ",round(tmle.out$ATE,2),"; 95% CI: (",round(tmle.out$lower.ci,2),", ",round(tmle.out$upper.ci,2),") \n",
               "IPW estimated ACE: ",round(ipw.out$ATE,2),"; 95% CI needs to be calculated via bootstrap \n",
               "G-comp estimated ACE: ",round(gcomp.out$ATE,2),"; 95% CI needs to be calculated via bootstrap"))

    np.out <- list(TMLE=tmle.out, Onestep=aipw.out, IPW=ipw.out, Gcomp=gcomp.out,
                   TMLE.Y1 = tmle_output_Y1, TMLE.Y0 = tmle_output_Y0,
                   Onestep.Y1 = aipw_output_Y1, Onestep.Y0 = aipw_output_Y0,
                   gcomp.Y1 = gcomp_output_Y1,
                   gcomp.Y0 = gcomp_output_Y0, ipw.Y1 = ipw_output_Y1, ipw.Y0 = ipw_output_Y0)


  }else if (length(a)==1) { ## E(Y^1) estimate ==

    out.a <- NPS.TMLE.a(a = a, data = data, vertices = vertices,
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

    np.out <- out.a

  } # end of if else condition for testing the length of a


  return(np.out)

}
