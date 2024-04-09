.force_weight <- function(fit, model, lib){

  if (length(lib)==1){ # only force weight to be one if there is only one algorithm in lib




    if (model=="cv"){ # code for force weight when fit is get from CV.SupearLearner

      K <- length(fit$coef) # number of folds

      fit$coef <- rep(1, K)

      for (i in 1:K){

        v_fit$AllSL[[i]]$coef <- 1

      }

    }else if (model=="sl"){ # code for force weight when fit is get from SupearLearner

      fit$coef <- 1

    }




  } # only force weight to be one if there is only one algorithm in lib



  return(fit)

}

