.force_weight <- function(fit, model, lib){

  if (length(lib)==1){ # only force weight to be one if there is only one algorithm in lib




    if (model=="cv"){ # code for force weight when fit is get from CV.SupearLearner

      fit$coef <- rep(1, length(fit$coef))

    }else if (model=="sl"){ # code for force weight when fit is get from SupearLearner

      fit$coef <- 1

    }




  } # only force weight to be one if there is only one algorithm in lib



  return(fit)

}

