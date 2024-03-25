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
  ps_model <- glm(
    (A==a1) ~ offset(qlogis(p.a1.mpA))+ clevercoef.A -1, family=binomial(), start=0
  )
  
  eps.A <- coef(ps_model)
  
  # update p(A=a1|mp(A))
  p.a1.mpA <- plogis(qlogis(p.a1.mpA)+eps.A*(clevercoef.A))
  
  # update density ratio of A
  assign("densratio_A", (1-p.a1.mpA)/p.a1.mpA) # density ratio regarding the treatment p(A|mp(A))|_{a_0}/p(A|mp(A))|_{a_1}
  
  # update density ratio for v in L if the ratio.method.L = "bayes"
  if (ratio.method.L == "bayes"){ for (v in L.removedA){assign(paste0("densratio_",v), get(paste0("bayes.densratio_",v))/densratio_A)} }
  
  # Update the density ratio data frame for L
  densratio.vectors.L <- mget(paste0("densratio_",L))
  densratio.L <- data.frame(densratio.vectors.L)
  
  # update density ratio for v in M if the ratio.method.M = "bayes"
  if (ratio.method.M == "bayes"){ for (v in M.mpM.excludeA){assign(paste0("densratio_",v), get(paste0("bayes.densratio_",v))/densratio.A)} }
  
  # Update the density ratio data frame for M
  densratio.vectors.M <- mget(paste0("densratio_",M))
  densratio.M <- data.frame(densratio.vectors.M)
  
  
  
  ######################
  # update E[Y|mp(Y)]
  ######################
  
  # density ratio before Y
  f.M_preY <- Reduce(`*`, densratio.M) # Mi precede Y = the set M
  f.L_preY <- Reduce(`*`, densratio.L) # Li precede Y = the set L
  
  # clever coefficient for outcome regression: E(Y|mp(Y))
  weight.Y <- if(outcome %in% L){(A==a1)*f.M_preY}else{(A==a0)*f.L_preY}
  
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
  EIF.Y <- if(outcome %in% L){(A==a1)*f.M_preY*(Y-mu.Y_a1)}else{(A==a0)*f.L_preY*(Y-mu.Y_a0)} # if Y in M
  
  ######################
  # update mu(v,a_v)
  ######################
  
  for (v in rev(vertices.between.AY)){ ## iterative over the vertices between A and Y to update mu(v,a_v)
    
    #### Update initial estimate of the sequential regression  ###
    
    # vertex that right after Z according to tau
    next.v <- tau.df$tau[tau.df$order=={tau.df$order[tau.df$tau==v]+1}]
    
    next.mu <- if(next.v %in% L){ get(paste0("mu.",next.v,"_a1"))}else{get(paste0("mu.",next.v,"_a0"))}
    
    next.mu.transform <- if(all(Y %in% c(0,1))){qlogis(next.mu)}else{next.mu} # if Y is binary, logit transform the sequential regression first
    
    
    if (crossfit==T){ #### cross fitting + super learner #####
      
      v_fit <- CV.SuperLearner(Y=next.mu.transform, X=dat_mpv, family = gaussian(), V = K, SL.library = lib.seq, control = list(saveFitLibrary=T),saveAll = T)
      
      
      #################### prediction: A in mp(v) vs A NOT in mp(v) ######################
      if (treatment %in% mpv){
        
        reg_a1 <- unlist(lapply(1:K, function(x) predict(v_fit$AllSL[[x]], newdata=dat_mpv.a1[v_fit$folds[[x]],])[[1]] %>% as.vector()))[order(unlist(lapply(1:K, function(x) v_fit$folds[[x]])))]
        reg_a0 <- unlist(lapply(1:K, function(x) predict(v_fit$AllSL[[x]], newdata=dat_mpv.a0[v_fit$folds[[x]],])[[1]] %>% as.vector()))[order(unlist(lapply(1:K, function(x) v_fit$folds[[x]])))]
        
        # assign prediction as mu.v_a1 and mu.v_a0
        assign(paste0("mu.",v,"_a1"), if(all(Y %in% c(0,1))){plogis(reg_a1)}else{reg_a1}) # transform back to probability scale if Y is binary
        assign(paste0("mu.",v,"_a0"), if(all(Y %in% c(0,1))){plogis(reg_a0)}else{reg_a0}) # transform back to probability scale if Y is binary
        
      }else{
        
        # assign prediction as mu.v_a1 and mu.v_a0, where mu.v_a1 = mu.v_a0 = mu.v
        assign(paste0("mu.",v,"_a1"), if(all(Y %in% c(0,1))){plogis(v_fit$SL.predict)}else{v_fit$SL.predict}) # transform back to probability scale if Y is binary
        assign(paste0("mu.",v,"_a0"), if(all(Y %in% c(0,1))){plogis(v_fit$SL.predict)}else{v_fit$SL.predict}) # transform back to probability scale if Y is binary
      }
      
      
      
    } else if (superlearner.seq==T){ #### super learner #####
      
      v_fit <- SuperLearner(Y=next.mu.transform, X=dat_mpv, family = gaussian(), SL.library = lib.seq)
      
      
      #################### prediction: A in mp(v) vs A NOT in mp(v) ######################
      if (treatment %in% mpv){
        
        reg_a1 <- predict(v_fit, newdata=dat_mpv.a1)[[1]] %>% as.vector()
        reg_a0 <- predict(v_fit, newdata=dat_mpv.a0)[[1]] %>% as.vector()
        
        # assign prediction as mu.v_a1 and mu.v_a0
        assign(paste0("mu.",v,"_a1"), if(all(Y %in% c(0,1))){plogis(reg_a1)}else{reg_a1}) # transform back to probability scale if Y is binary
        assign(paste0("mu.",v,"_a0"), if(all(Y %in% c(0,1))){plogis(reg_a0)}else{reg_a0}) # transform back to probability scale if Y is binary)
        
      }else{
        
        # assign prediction as mu.v_a1 and mu.v_a0, where mu.v_a1 = mu.v_a0 = mu.v
        assign(paste0("mu.",v,"_a1"), if(all(Y %in% c(0,1))){plogis(predict(v_fit)[[1]] %>% as.vector())}else{predict(v_fit)[[1]] %>% as.vector()}) # transform back to probability scale if Y is binary
        assign(paste0("mu.",v,"_a0"), if(all(Y %in% c(0,1))){plogis(predict(v_fit)[[1]] %>% as.vector())}else{predict(v_fit)[[1]] %>% as.vector()}) # transform back to probability scale if Y is binary
      }
      
      
      
    }else { #### linear regression #####
      
      v_fit <- lm(next.mu.transform ~ ., data=dat_mpv) # fit linear regression/ logistic regression depending on type of v
      
      #################### prediction: A in mp(v) vs A NOT in mp(v) ######################
      if (treatment %in% mpv){
        
        reg_a1 <- predict(v_fit, newdata=dat_mpv.a1)
        reg_a0 <- predict(v_fit, newdata=dat_mpv.a0)
        
        # assign prediction as mu.v_a1 and mu.v_a0
        assign(paste0("mu.",v,"_a1"), if(all(Y %in% c(0,1))){plogis(reg_a1)}else{reg_a1}) # transform back to probability scale if Y is binary
        assign(paste0("mu.",v,"_a0"), if(all(Y %in% c(0,1))){plogis(reg_a0)}else{reg_a0}) # transform back to probability scale if Y is binary
        
      }else{
        
        # assign prediction as mu.v_a1 and mu.v_a0, where mu.v_a1 = mu.v_a0 = mu.v
        assign(paste0("mu.",v,"_a1"), if(all(Y %in% c(0,1))){plogis(predict(v_fit))}else{predict(v_fit)}) # transform back to probability scale if Y is binary
        assign(paste0("mu.",v,"_a0"), if(all(Y %in% c(0,1))){plogis(predict(v_fit))}else{predict(v_fit)}) # transform back to probability scale if Y is binary
        
      }
      
    } ## End of update of initial estimates of mu(v,a_v)
    
    ##########################################################################################################
    
    
    
    ## Perform TMLE update for mu(v,a_v)
    
    # offset
    if (all(Y %in% c(0,1))){ # L(mu(v)) = mu(next.v)*log(mu(v)(eps.v))+(1-mu(next.v))*log(1-mu(v)(eps.v)), where mu(v)(eps.v)=expit{logit(mu(v))+eps.v*weight}
      
      offset.v <- if(v %in% L){qlogis(get(paste0("mu.",v,"_a1")))}else{qlogis(get(paste0("mu.",v,"_a0")))}
      
    }else{ # L(mu(v)) = weight*{mu(next.v)-mu(v)(eps.v)}^2, where mu(v)(eps.v)=mu(v)+eps.v
      offset.v <- if(v %in% L){get(paste0("mu.",v,"_a1"))}else{get(paste0("mu.",v,"_a0"))}}
    
    
    # weight for regression
    # select M and L that precede v
    selected.M <- M[sapply(M, function(m) tau.df$order[tau.df$tau==m] < tau.df$order[tau.df$tau==v])] # M precedes Z
    selected.L <- L[sapply(L, function(l) tau.df$order[tau.df$tau==l] < tau.df$order[tau.df$tau==v])] # L precedes Z
    
    
    weight.v <- if(v %in% L){ 
      # product of the selected variables density ratio
      f.M_prev <- Reduce(`*`, densratio.M[,paste0("densratio_",selected.M), drop=F]) # Mi precede Z 
      
      # weight
      (A==a1)*f.M_prev }else{
        
        # product of the selected variables density ratio
        f.L_prev <- Reduce(`*`, densratio.L[,paste0("densratio_",selected.L),drop=F]) # Mi precede Z 
        
        # weight
        (A==a0)*f.L_prev }
    
    
    # one iteration update
    if (all(Y %in% c(0,1))){
      v_model <- lm(next.mu ~ offset(offset.v)+1, weights = weight.v) # fit linear regression if Y is continuous
    }else{
      v_model <- glm(next.mu ~ offset(offset.v)+weight.v-1, family=binomial(link="logit")) # fit logistic regression if Y is binary, this is like a generalized logistic regression
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
  mean((A==a0)*Y) - 
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
                 EDstar.record = EDstar.record # record of EDstar
)
