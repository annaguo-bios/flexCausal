######################################################################################################################################
####################################### SEQUENTIAL REGRESSION FOR ESETIMATING mu(Z_i, a_{Z_i}) #######################################
######################################################################################################################################

# the sequential regression will be perform for v that locates between A and Y accroding to topological order tau
vertices.between.AY <- tau[{which(tau==treatment)+1}:{which(tau==outcome)-1}]

# sequential regression
for (v in rev(vertices.between.AY)){ ## iterate through vertices between A and Y according to topological order tau
  
  # perform sequential regressions for v with larger order first: rev()
  
  ### Prepare Markov pillow and dataset for regression and prediction ####
  
  # Find Markov pillow of v
  mpv <- markov_pillow(vertices, di_edges, bi_edges, v) # Markov pillow for outcome
  
  # prepare dataset for regression and prediction
  dat_mpv <- data[,mpv] # extract data for Markov pillow for v
  
  if (treatment %in% mpv){ # we only need to consider evaluate regression at specific level of A if treatment is in Markov pillow for v
    
    # set treatment to a0
    dat_mpv.a0 <- dat_mpv %>% mutate(!!treatment := a0) 
    
    # set treatment to a1 
    dat_mpv.a1 <- dat_mpv %>% mutate(!!treatment := a1) }
  
  # vertex that right after Z according to tau
  next.v <- tau.df$tau[tau.df$order=={tau.df$order[tau.df$tau==v]+1}]
  
  next.mu <- if(next.v %in% L){ get(paste0("mu.",next.v,"_a1")) }else{ get(paste0("mu.",next.v,"_a0")) } # mu is evaluated at a0 for next.v in M and at a1 for next.v in L
  
  next.mu.transform <- if(all(Y %in% c(0,1))){qlogis(next.mu)}else{next.mu} # if Y is binary, logit transform the sequential regression first
  
  ### End of preparation ####
  
  if (crossfit==T){ #### cross fitting + super learner #####
    
    v_fit <- CV.SuperLearner(Y=next.mu.transform, X=dat_mpv, family = gaussian(), V = K, SL.library = lib.seq, control = list(saveFitLibrary=T), saveAll = T)
    

    ######## prediction: A in mp(v) vs A NOT in mp(v) ########
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
    
    #########################################################
    
    
  } else if (superlearner.seq==T){ #### super learner #####
    
    v_fit <- SuperLearner(Y=next.mu.transform, X=dat_mpv, family = gaussian(), SL.library = lib.seq)
    
    ######## prediction: A in mp(v) vs A NOT in mp(v) ########
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
    #########################################################
    
    
    
  }else { #### linear regression #####
    
    v_fit <- lm(next.mu.transform ~ ., data=dat_mpv) # fit linear regression/ logistic regression depending on type of v
    
    ######## prediction: A in mp(v) vs A NOT in mp(v)########
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
    ##########################################################
  }
  
}  ## End of iteration over all vertics between A and Y

######################################################################################################################################
####################################### END OF SEQUENTIAL REGRESSION FOR ESETIMATING mu(Z_i, a_{Z_i}) ################################
######################################################################################################################################

######################
# EIF calculations
######################

## EIF for Y|mp(Y): 
#if Y in L: I(A=a1)*I(M < Y){p(M|mp(M))|_{a0}/p(M|mp(M))|_{a1}}*(Y-E(Y|mp(Y))|_{a1})
#if Y in M: I(A=a0)*I(L < Y){p(L|mp(L))|_{a0}/p(:|mp(L))|_{a1}}*(Y-E(Y|mp(Y))|_{a0})

f.M_preY <- Reduce(`*`, densratio.M) # Mi precede Y = the set M
f.L_preY <- Reduce(`*`, densratio.L) # Li precede Y = the set L


EIF.Y <- if(outcome %in% L){
  
  (A==a1)*f.M_preY*(Y-mu.Y_a1)}else{ # if Y in L
    
    (A==a0)*f.L_preY*(Y-mu.Y_a0)} # if Y in M

## EIF for Z|mp(Z)

for (v in rev(vertices.between.AY)){ ## iterate over all vertices between A and Y
  
  # select M and L that precede v
  selected.M <- M[sapply(M, function(m) tau.df$order[tau.df$tau==m] < tau.df$order[tau.df$tau==v])] # M precedes Z
  selected.L <- L[sapply(L, function(l) tau.df$order[tau.df$tau==l] < tau.df$order[tau.df$tau==v])] # L precedes Z
  
  # vertex that right after v according to tau
  next.v <- tau.df$tau[tau.df$order=={tau.df$order[tau.df$tau==v]+1}]
  
  # mu(next.v, a_{next.v})
  next.mu <- if(next.v %in% L){ get(paste0("mu.",next.v,"_a1"))}else{get(paste0("mu.",next.v,"_a0"))}
  
  EIF.v <- if(v %in% L){ # v in L
    
    # product of the selected variables density ratio
    f.M_prev <- Reduce(`*`, densratio.M[,paste0("densratio_",selected.M), drop=F]) # Mi precede v
    
    # EIF for v|mp(v)
    (A==a1)*f.M_prev*( next.mu - get(paste0("mu.",v,"_a1")) )
    
  }else{ # v in M
    
    # product of the selected variables density ratio
    f.L_prev <- Reduce(`*`, densratio.L[,paste0("densratio_",selected.L), drop=F]) # Li precede v 
    
    # EIF for v|mp(v)
    (A==a0)*f.L_prev*( next.mu - get(paste0("mu.",v,"_a0")) )
  }
  
  assign(paste0("EIF.",v), EIF.v) 
  
} ## End of iteration over all vertices between A and Y

## EIF for A=a1|mp(A)

# vertex that right after A according to tau
next.A <- tau.df$tau[tau.df$order=={tau.df$order[tau.df$tau==treatment]+1}]

EIF.A <- {(A==a1) - p.a1.mpA}*get(paste0("mu.",next.A,"_a0"))


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


onestep.out <- list(estimated_psi=estimated_psi, # estimated parameter
                    lower.ci=lower.ci, # lower bound of 95% CI
                    upper.ci=upper.ci, # upper bound of 95% CI
                    EIF=EIF # E(Dstar) for Y|M,A,X and M|A,X, and A|X
)

