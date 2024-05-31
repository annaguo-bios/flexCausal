.force_weight <- function(fit, model, lib){

  if (length(lib)==1){ # only force weight to be one if there is only one algorithm in lib




    if (model=="cv"){ # code for force weight when fit is get from CV.SupearLearner

      K <- length(fit$coef) # number of folds

      fit$coef <- rep(1, K)

      for (i in 1:K){

        fit$AllSL[[i]]$coef <- 1

      }

    }else if (model=="sl"){ # code for force weight when fit is get from SupearLearner

      fit$coef <- 1

    }




  } # only force weight to be one if there is only one algorithm in lib



  if (length(lib)>1){ # if there are more than one learner, and the coefficient for all learners are 0




    if (model=="cv"){ # code for force weight when fit is get from CV.SupearLearner

      K <- nrow(fit$coef) # number of folds
      n.learners <- ncol(fit$coef) # number of learners


      for (i in 1:K){

        if (sum(fit$AllSL[[i]]$coef)==0){ # if the coefficient for all learners are 0, then force equal weights to all learners

          fit$AllSL[[i]]$coef <- rep(1/n.learners,n.learners)
          fit$coef[i,] <- rep(1/n.learners,n.learners)

        }



      }

    }else if (model=="sl" & sum(fit$coef)==0){ # code for force weight when fit is get from SupearLearner

      n.learners <- length(lib) # number of learners

      fit$coef <- rep(1/n.learners,n.learners)

    }




  } # only force weight to be one if there is only one algorithm in lib



  return(fit)

}

