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

  if (length(a) == 0) stop("`a` must be a numeric scalar or length-two vector.")
  if (length(a) > 2)  stop("Invalid input. Enter a = c(value1, value2) for ACE estimation or a = value1 for E(Y(a)) estimation.")

  # call estimation function
  call_nps <- function(a_val) {
    NPS.TMLE.a(a = a_val, data = data, vertices = vertices, di_edges = di_edges, bi_edges = bi_edges, treatment = treatment, outcome = outcome, multivariate.variables = multivariate.variables, graph = graph,
                    superlearner.seq = superlearner.seq, # whether run superlearner for sequential regression
                    superlearner.Y = superlearner.Y, # whether run superlearner for outcome regression
                    superlearner.A = superlearner.A, # whether run superlearner for propensity score
                    superlearner.M = superlearner.M, # whether run superlearner for estimating densratio for M using bayes method
                    superlearner.L = superlearner.L, # whether run superlearner for estimating densratio for L using bayes method
                    crossfit = crossfit, K = K,
                    ratio.method.L = ratio.method.L, # method for estimating the density ratio associated with L
                    ratio.method.M = ratio.method.M, # method for estimating the density ratio associated with M
                    dnorm.formula.L = dnorm.formula.L, # formula for estimating the density ratio associated with L
                    dnorm.formula.M = dnorm.formula.M, # formula for estimating the density ratio associated with M
                    lib.seq = lib.seq, # superlearner library for sequential regression
                    lib.L = lib.L, # superlearner library for density ratio estimation via bayes rule for variables in L
                    lib.M = lib.M, # superlearner library for density ratio estimation via bayes rule for variables in M
                    lib.Y = lib.Y, # superlearner library for outcome regression
                    lib.A = lib.A, # superlearner library for propensity score
                    formulaY = formulaY, formulaA = formulaA, # regression formula for outcome regression and propensity score if superlearner is not used
                    linkY_binary = linkY_binary, linkA = linkA, # link function for outcome regression and propensity score if superlearner is not used
                    n.iter = n.iter, cvg.criteria = cvg.criteria,
                    truncate_lower = truncate_lower, truncate_upper = truncate_upper,zerodiv.avoid=zerodiv.avoid)}

  # organize output
  make_ate_out <- function(out1, out0 = NULL) {
    if (is.null(out0)) {
      list(EYa     = out1$estimated_psi,
           lower.ci = out1$lower.ci,
           upper.ci = out1$upper.ci,
           EIF      = out1$EIF)

    } else {
      ate <- out1$estimated_psi - out0$estimated_psi
      eif <- out1$EIF - out0$EIF
      list(ACE      = ate,
           lower.ci = ate - 1.96 * sqrt(mean(eif^2) / n),
           upper.ci = ate + 1.96 * sqrt(mean(eif^2) / n),
           EIF      = eif)

    }
  }

  ## ACE estimate ==
  if (length(a)==2){

    out.a1 <- call_nps(a[1])
    out.a0 <- call_nps(a[2])

    tmle.out  <- make_ate_out(out.a1$TMLE,    out.a0$TMLE)
    aipw.out  <- make_ate_out(out.a1$Onestep, out.a0$Onestep)

    message(paste0("Onestep estimated ACE: ",round(aipw.out$ACE,2),"; 95% CI: (",round(aipw.out$lower.ci,2),", ",round(aipw.out$upper.ci,2),") \n",
               "TMLE estimated ACE: ",round(tmle.out$ACE,2),"; 95% CI: (",round(tmle.out$lower.ci,2),", ",round(tmle.out$upper.ci,2),")"))

    np.out <- list(TMLE=tmle.out,Onestep=aipw.out, TMLE.Y1=out.a1$TMLE, TMLE.Y0 = out.a0$TMLE, Onestep.Y1=out.a1$Onestep, Onestep.Y0=out.a0$Onestep)

    ## E(Y^1) estimate ==
  }else if (length(a)==1) {

    out.a <- call_nps(a)

    tmle.out  <- make_ate_out(out.a$TMLE,    NULL)
    aipw.out  <- make_ate_out(out.a$Onestep, NULL)

    message(paste0("Onestep estimated E(Y(a)): ",  round(aipw.out$EYa, 2),  "; 95% CI: (", round(aipw.out$lower.ci, 2),  ", ", round(aipw.out$upper.ci, 2),  ") \n",
               "TMLE estimated E(Y(a)): ",     round(tmle.out$EYa, 2),  "; 95% CI: (", round(tmle.out$lower.ci, 2),  ", ", round(tmle.out$upper.ci, 2),  ")"))

    np.out <- list(TMLE=tmle.out,Onestep=aipw.out, TMLE.Ya=out.a$TMLE, Onestep.Ya=out.a$Onestep)

  } # end of if else condition for testing the length of a

  return(np.out)
}
