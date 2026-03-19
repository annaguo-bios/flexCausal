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

  if (length(a) == 0) stop("`a` must be a numeric scalar or length-two vector.")
  if (length(a) > 2)  stop("Invalid input. Enter a = c(value1, value2) for ACE estimation or a = value1 for E(Y(a)) estimation.")

  # call estimation function
  call_one <- function(a_val) {
    backdoor.TMLE.a(a = a_val, data = data, vertices = vertices,
                    di_edges = di_edges, bi_edges = bi_edges,
                    treatment = treatment, outcome = outcome,
                    multivariate.variables = multivariate.variables,
                    graph = graph, superlearner.Y = superlearner.Y,
                    superlearner.A = superlearner.A, crossfit = crossfit,
                    K = K, lib.Y = lib.Y, lib.A = lib.A,
                    formulaY = formulaY, formulaA = formulaA,
                    linkY_binary = linkY_binary, linkA = linkA,
                    truncate_lower = truncate_lower,
                    truncate_upper = truncate_upper)}


  # organize output
  make_ate_out <- function(out1, out0 = NULL, has_eif = TRUE) {
    if (is.null(out0)) {
      if (has_eif) {
        list(EYa     = out1$estimated_psi,
             lower.ci = out1$lower.ci,
             upper.ci = out1$upper.ci,
             EIF      = out1$EIF)
      } else {
        list(EYa = out1$estimated_psi)
      }
    } else {
      ate <- out1$estimated_psi - out0$estimated_psi
      if (has_eif) {
        eif <- out1$EIF - out0$EIF
        list(ACE      = ate,
             lower.ci = ate - 1.96 * sqrt(mean(eif^2) / n),
             upper.ci = ate + 1.96 * sqrt(mean(eif^2) / n),
             EIF      = eif)
      } else {
        list(ACE = ate)
      }
    }
  }


    if (length(a)==2){ ## ACE estimate ==

      ## TMLE estimator


      out.a1 <- call_one(a[1])
      out.a0 <- call_one(a[2])

      tmle.out  <- make_ate_out(out.a1$TMLE,    out.a0$TMLE)
      aipw.out  <- make_ate_out(out.a1$Onestep, out.a0$Onestep)
      gcomp.out <- make_ate_out(out.a1$gcomp,   out.a0$gcomp,  has_eif = FALSE)
      ipw.out   <- make_ate_out(out.a1$ipw,     out.a0$ipw,    has_eif = FALSE)

      message(paste0("Onestep estimated ACE: ",round(aipw.out$ACE,2),"; 95% CI: (",round(aipw.out$lower.ci,2),", ",round(aipw.out$upper.ci,2),") \n",
                 "TMLE estimated ACE: ",round(tmle.out$ACE,2),"; 95% CI: (",round(tmle.out$lower.ci,2),", ",round(tmle.out$upper.ci,2),") \n",
                 "IPW estimated ACE: ",round(ipw.out$ACE,2),"; 95% CI needs to be calculated via bootstrap \n",
                 "G-comp estimated ACE: ",round(gcomp.out$ACE,2),"; 95% CI needs to be calculated via bootstrap"))

      np.out <- list(TMLE=tmle.out, Onestep=aipw.out, IPW=ipw.out, Gcomp=gcomp.out,
                     TMLE.Y1 = out.a1$TMLE, TMLE.Y0 = out.a0$TMLE,
                     Onestep.Y1 = out.a1$Onestep, Onestep.Y0 = out.a0$Onestep,
                     gcomp.Y1 = out.a1$gcomp, gcomp.Y0 = out.a0$gcomp,
                     ipw.Y1 = out.a1$ipw, ipw.Y0 = out.a0$ipw)


    }else if (length(a)==1) { ## E(Y^1) estimate ==

      out.a <- call_one(a)

      tmle.out  <- make_ate_out(out.a$TMLE,    NULL, has_eif = TRUE)
      aipw.out  <- make_ate_out(out.a$Onestep, NULL, has_eif = TRUE)
      gcomp.out <- make_ate_out(out.a$gcomp,   NULL,  has_eif = FALSE)
      ipw.out   <- make_ate_out(out.a$ipw,     NULL,    has_eif = FALSE)

      message(paste0("Onestep estimated E(Y(a)): ",  round(aipw.out$EYa, 2),  "; 95% CI: (", round(aipw.out$lower.ci, 2),  ", ", round(aipw.out$upper.ci, 2),  ") \n",
                 "TMLE estimated E(Y(a)): ",     round(tmle.out$EYa, 2),  "; 95% CI: (", round(tmle.out$lower.ci, 2),  ", ", round(tmle.out$upper.ci, 2),  ") \n",
                 "IPW estimated E(Y(a)): ",      round(ipw.out$EYa, 2),   "; 95% CI needs to be calculated via bootstrap \n",
                 "G-comp estimated E(Y(a)): ",   round(gcomp.out$EYa, 2), "; 95% CI needs to be calculated via bootstrap"))

      np.out <- list(TMLE=tmle.out, Onestep=aipw.out, IPW=ipw.out, Gcomp=gcomp.out,
                     TMLE.Ya = out.a$TMLE, Onestep.Ya = out.a$Onestep, gcomp.Ya = out.a$gcomp, ipw.Ya = out.a$ipw)

    } # end of if else condition for testing the length of a




  return(np.out)

}
