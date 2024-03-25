################################################
############### OUTCOME REGRESSION #############
################################################

#### Fit nuisance models ####

# Find Markov pillow of outcome
mpY <- markov_pillow(vertices, di_edges, bi_edges, outcome) # Markov pillow for outcome

# prepare dataset for regression and prediction
dat_mpY <- data[,mpY] # extract data for Markov pillow for outcome
dat_mpY.a0 <- dat_mpY %>% mutate(!!treatment := a0) # set treatment to a0
dat_mpY.a1 <- dat_mpY %>% mutate(!!treatment := a1) # set treatment to a1

if (crossfit==T){ #### cross fitting + super learner #####
  
  fit.family <- if(all(Y %in% c(0,1))){binomial()}else{gaussian()} # family for super learner depending on whether Y is binary or continuous
  
  or_fit <- CV.SuperLearner(Y=Y, X=dat_mpY, family = fit.family, V = K, SL.library = lib.Y, control = list(saveFitLibrary=T),saveAll = T)
  
  mu.Y_a1 <- unlist(lapply(1:K, function(x) predict(or_fit$AllSL[[x]], newdata=dat_mpY.a1[or_fit$folds[[x]],])[[1]] %>% as.vector()))[order(unlist(lapply(1:K, function(x) or_fit$folds[[x]])))]
  mu.Y_a0 <- unlist(lapply(1:K, function(x) predict(or_fit$AllSL[[x]], newdata=dat_mpY.a0[or_fit$folds[[x]],])[[1]] %>% as.vector()))[order(unlist(lapply(1:K, function(x) or_fit$folds[[x]])))]

    
} else if (superlearner.Y==T){ #### super learner #####
  
  fit.family <- if(all(Y %in% c(0,1))){binomial()}else{gaussian()} # family for super learner depending on whether Y is binary or continuous
    
  or_fit <- SuperLearner(Y=Y, X=dat_mpY, family = fit.family, SL.library = lib.Y)
  
  mu.Y_a1 <- predict(or_fit, newdata=dat_mpY.a1)[[1]] %>% as.vector()
  mu.Y_a0 <- predict(or_fit, newdata=dat_mpY.a0)[[1]] %>% as.vector()

  
} else { #### simple linear regression with user input regression formula: default="Y ~ ." ####
    
  fit.family <- if(all(Y %in% c(0,1))){binomial()}else{gaussian()} # family for super learner depending on whether Y is binary or continuous
  
  or_fit <- glm(as.formula(formulaY), data=dat_mpY, family = fit.family)
  
  mu.Y_a1 <- predict(or_fit, newdata=dat_mpY.a1)
  mu.Y_a0 <- predict(or_fit, newdata=dat_mpY.a0)

  
}


################################################
############### PROPENSITY SCORE ###############
################################################

#### Fit nuisance models ####

# Find Markov pillow of treatment
mpA <- markov_pillow(vertices, di_edges, bi_edges, treatment) # Markov pillow for treatment

# prepare dataset for regression and prediction
dat_mpA <- data[,mpA, drop = F] # extract data for Markov pillow for outcome

if (crossfit==T){ #### cross fitting + super learner #####
  
  # fit model
  ps_fit <- CV.SuperLearner(Y=A, X=dat_mpA, family = binomial(), V = K, SL.library = lib.A, control = list(saveFitLibrary=T),saveAll = T)
  
  # make prediction: p(A=1|mp(A))
  p.A1.mpA <- ps_fit$SL.predict
  
} else if (superlearner.A==T){ #### super learner #####
  
  # fit model
  ps_fit <- SuperLearner(Y=A, X=dat_mpA, family = binomial(), SL.library = lib.A)
  
  # make prediction: p(A=1|mp(A))
  p.A1.mpA <- predict(ps_fit, type = "response")[[1]] %>% as.vector()
  
} else { #### simple linear regression with user input regression formula: default="A ~ ." ####
  
  if (truncate_lower!=0 | truncate_upper!=1){ # under weak overlapping issue, it's more stable to run A~X via linear regression
    
    # fit model
    ps_fit <- lm(as.formula(formulaA), data=dat_mpA)
    
    # make prediction: p(A=1|mp(A))
    p.A1.mpA <- predict(ps_fit)
    
  } else { # without weak overlapping issue. Run A~X via logistic regression
    
    # fit model
    ps_fit <- glm(as.formula(formulaA), data=dat_mpA,  family = binomial())
    
    # make prediction: p(A=1|mp(A))
    p.A1.mpA <- predict(ps_fit, type = "response")  # p(A=1|X)
    
  }
}

# apply truncation to propensity score to deal with weak overlap. 
# truncated propensity score within the user specified range of [truncate_lower, truncate_upper]: default=[0,1]
p.A1.mpA[p.A1.mpA < truncate_lower] <- truncate_lower
p.A1.mpA[p.A1.mpA > truncate_upper] <- truncate_upper

p.a1.mpA <- a1*p.A1.mpA + (1-a1)*(1-p.A1.mpA) # p(A=a1|mp(A))

assign("densratio_A", (1-p.a1.mpA)/p.a1.mpA) # density ratio regarding the treatment p(A|mp(A))|_{a_0}/p(A|mp(A))|_{a_1}




################################################
############### DENSITY RATIO.   ###############
################################################

###### DENSITY RATIO ASSOCIATED WITH L ######

if (ratio.method.L=="densratio"){ ################### METHOD 2A: densratio method  ###################
  
  # Error3: densratio method doesn't support factor variables
  
  if (!all(sapply(c(L.removedA, markov_pillow(vertices, di_edges, bi_edges, L.removedA)), function(var) is.numeric(data[,var]) | is.integer(data[,var])))){
    
    print("Error in estimating density ratios associated with variables in L: densratio method only support numeric/integer variables, try bayes method instead.")
    
    stop() }
  
  
  # if M,A,X only consists numeric/integer variables: apply density ratio estimation
  for (v in L.removedA){ ## Iterate over each variable in L\A
    
    dat_v.a0 <- data[data[[treatment]] == a0,c(v, markov_pillow(vertices, di_edges, bi_edges, v))] # select rows where A=a0
    dat_v.a1 <- data[data[[treatment]] == a1,c(v, markov_pillow(vertices, di_edges, bi_edges, v))] # select rows where A=a1
    
    densratio.v <- densratio(dat_v.a0, dat_v.a1)
    
    ratio <- densratio.v$compute_density_ratio(data[,c(v,markov_pillow(vertices, di_edges, bi_edges, v))]) # p(L|mp(L))|_{a_0}/p(L|mp(L))|_{a_1}
    
    assign(paste0("densratio_",v), ratio)
  }
  
  
  
} else if (ratio.method.L=="bayes"){ ################### METHOD 2B: Bayes method ###################
  
  for (v in L.removedA){ ## Iterate over each variable in L\A
    
    #### Prepare data for regression and prediction ####
    dat_bayes.v <- data[,setdiff(c(v,markov_pillow(vertices, di_edges, bi_edges, v)),treatment), drop=F] # contains variable v + Markov pillow of v - treatment
    
    #### Fit nuisance models ####
    
    if (crossfit==T){
      
      bayes_fit <- CV.SuperLearner(Y=A, X=dat_bayes.v, family = binomial(), V = K, SL.library = lib.L, control = list(saveFitLibrary=T),saveAll = T)
      
      # p(A=1|mp(v)\A,v)
      p.A1.mpv <- bayes_fit$SL.predict
      
      #p(v=a0|mp(v)\A,v)
      p.a0.mpv <- a0*p.A1.mpv+(1-a0)*(1-p.A1.mpv)
      
      #p(v=a0|mp(v)\A,v)
      p.a1.mpv <- 1-p.a0.mpv
      
      # save the density ratio: p(v|mp(V))|_{a0}/p(v|mp(V))|_{a1}
      assign(paste0("densratio_",v), {p.a0.mpv/p.a1.mpv}/densratio_A)
      
      # save the ratio of p(A|mp(v)\A,v)|_{a0}/p(A|mp(v)\A,v)|_{a1}
      # such that we can come back to update the density ratio of v once we update the densratioA
      assign(paste0("bayes.densratio_",v), {p.a0.mpv/p.a1.mpv})
      
    }else if (superlearner.L==T){
      
      bayes_fit <- SuperLearner(Y=A, X=dat_bayes.v, family = binomial(), SL.library = lib.L)
      
      # p(A=1|mp(v)\A,v)
      p.A1.mpv <- predict(bayes_fit, type = "response")[[1]] %>% as.vector()  # p(A=1|X)
      
      #p(v=a0|mp(v)\A,v)
      p.a0.mpv <- a0*p.A1.mpv+(1-a0)*(1-p.A1.mpv)
      
      #p(v=a0|mp(v)\A,v)
      p.a1.mpv <- 1-p.a0.mpv
      
      # save the density ratio: p(v|mp(V))|_{a0}/p(v|mp(V))|_{a1}
      assign(paste0("densratio_",v), {p.a0.mpv/p.a1.mpv}/densratio_A)
      
      # save the ratio of p(A|mp(v)\A,v)|_{a0}/p(A|mp(v)\A,v)|_{a1}
      # such that we can come back to update the density ratio of v once we update the densratioA
      assign(paste0("bayes.densratio_",v), {p.a0.mpv/p.a1.mpv})
      
    } else {
      
      # estimate density ratio using bayes rule
      bays_fit <- glm(A ~ ., data=dat_bayes.v, family = binomial())
      
      # p(A=1|mp(v)\A,v)
      p.A1.mpv <- predict(bays_fit, type = "response")
      
      #p(v=a0|mp(v)\A,v)
      p.a0.mpv <- a0*p.A1.mpv+(1-a0)*(1-p.A1.mpv)
      
      #p(v=a0|mp(v)\A,v)
      p.a1.mpv <- 1-p.a0.mpv
      
      # save the density ratio: p(v|mp(V))|_{a0}/p(v|mp(V))|_{a1}
      assign(paste0("densratio_",v), {p.a0.mpv/p.a1.mpv}/densratio_A)
      
      # save the ratio of p(A|mp(v)\A,v)|_{a0}/p(A|mp(v)\A,v)|_{a1}
      # such that we can come back to update the density ratio of v once we update the densratioA
      assign(paste0("bayes.densratio_",v), {p.a0.mpv/p.a1.mpv})
      
    }
    
    
  } ## Iterate over each variable in L\A
  
  
} else {
  
  print("Invalid ratio.method.L input.")
  
}

# Get all the densratio vectors based on their names
densratio.vectors.L <- mget(c(paste0("densratio_",L.removedA),"densratio_A"))

# Create a data frame using all the vectors
densratio.L <- data.frame(densratio.vectors.L)

###### DENSITY RATIO ASSOCIATED WITH M ######


# Find vertices in set M, that the Markov pillow of which include A and don't include A
# For the vertices whose Markov pillow include A, we will use either the "densratio" or "bayes" method to estimate the density ratio p(M|mp(M))|_{a0}/p(M|mp(M))|_{a1}
# For the vertices whose Markov pillow don't include A, density ratio p(M|mp(M))|_{a0}/p(M|mp(M))|_{a1} is 1

M.mpM.includeA <- c()
M.mpM.excludeA <- c()

for (v in M){ ## Iterate over each variable in M
  
  if (treatment %in% markov_pillow(vertices, di_edges, bi_edges, v)){ # vertices whose Markov pillow include A
    M.mpM.includeA <- c(M.mpM.includeA, v)
    
  } else { # vertices whose Markov pillow don't include A
    
    M.mpM.excludeA <- c(M.mpM.excludeA, v)
    
  }
  
} ## Iterate over each variable in M

# assign ratio=1 for vertices in M.mpM.excludeA
for (m in M.mpM.excludeA) { assign(paste("densratio_", m), 1) }

if (ratio.method.M=="densratio"){ ################### METHOD 2A: densratio method  ###################
  
  # Error: densratio method doesn't support factor variables
  
  if (!all(sapply(c(M.mpM.includeA, markov_pillow(vertices, di_edges, bi_edges, M.mpM.includeA)), function(var) is.numeric(data[,var]) | is.integer(data[,var])))){
    
    print("Error in estimating density ratios associated with variables in M: densratio method only support numeric/integer variables, try bayes method instead.")
    
    stop() }
  
  
  # if M and mpi(M) only consists numeric/integer variables: apply density ratio estimation
  for (v in M.mpM.includeA){
    
    dat_v.a0 <- data[data[[treatment]] == a0,c(v, markov_pillow(vertices, di_edges, bi_edges, v))] # select rows where A=a0
    dat_v.a1 <- data[data[[treatment]] == a1,c(v, markov_pillow(vertices, di_edges, bi_edges, v))] # select rows where A=a1
    
    densratio.v <- densratio(dat_v.a0, dat_v.a1)
    
    ratio <- densratio.v$compute_density_ratio(data[,c(v,markov_pillow(vertices, di_edges, bi_edges, v))]) # p(M|mp(M))|_{a_0}/p(M|mp(M))|_{a_1}
    
    assign(paste0("densratio_",v), ratio)
  }
  
  
} else if (ratio.method.M=="bayes"){ ################### METHOD 2B: Bayes method ###################
  
  for (v in M.mpM.includeA){ ## Iterate over each variable in M
    
    #### Prepare data for regression and prediction ####
    dat_bayes.v <- data[,setdiff(c(v,markov_pillow(vertices, di_edges, bi_edges, v)),treatment), drop=F] # contains variable v + Markov pillow of v - treatment
    
    #### Fit nuisance models ####
    
    if (crossfit==T){
      
      # fit p(A|mp(v)\A,v)
      bayes_fit <- CV.SuperLearner(Y=A, X=dat_bayes.v, family = binomial(), V = K, SL.library = lib.M, control = list(saveFitLibrary=T),saveAll = T)
      
      # p(A=1|mp(v)\A,v)
      p.A1.mpv <- bayes_fit$SL.predict
      
      #p(A=a0|mp(v)\A,v)
      p.a0.mpv <- a0*p.A1.mpv+(1-a0)*(1-p.A1.mpv)
      
      #p(A=a1|mp(v)\A,v)
      p.a1.mpv <- 1-p.a0.mpv
      
      # save the density ratio: p(v|mp(V))|_{a0}/p(v|mp(V))|_{a1}
      assign(paste0("densratio_",v), {p.a0.mpv/p.a1.mpv}/densratio_A)
      
      # save the ratio of p(A|mp(v)\A,v)|_{a0}/p(A|mp(v)\A,v)|_{a1}
      # such that we can come back to update the density ratio of v once we update the densratioA
      assign(paste0("bayes.densratio_",v), {p.a0.mpv/p.a1.mpv})
      
    }else if (superlearner.M==T){
      
      # fit p(A|mp(v)\A,v)
      bayes_fit <- SuperLearner(Y=A, X=dat_bayes.v, family = binomial(), SL.library = lib.L)
      
      # p(A=1|mp(v)\A,v)
      p.A1.mpv <- predict(bayes_fit, type = "response")[[1]] %>% as.vector()  # p(A=1|X)
      
      #p(A=a0|mp(v)\A,v)
      p.a0.mpv <- a0*p.A1.mpv+(1-a0)*(1-p.A1.mpv)
      
      #p(A=a1|mp(v)\A,v)
      p.a1.mpv <- 1-p.a0.mpv
      
      # save the density ratio: p(v|mp(V))|_{a0}/p(v|mp(V))|_{a1}
      assign(paste0("densratio_",v), {p.a0.mpv/p.a1.mpv}/densratio_A)
      
      # save the ratio of p(A|mp(v)\A,v)|_{a0}/p(A|mp(v)\A,v)|_{a1}
      # such that we can come back to update the density ratio of v once we update the densratioA
      assign(paste0("bayes.densratio_",v), {p.a0.mpv/p.a1.mpv})
      
    } else {
      
      # estimate density ratio using bayes rule
      bays_fit <- glm(A ~ ., data=dat_bayes.v, family = binomial())
      
      # p(A=1|mp(v)\A,v)
      p.a1.mpv <- predict(bays_fit, type = "response")
      
      #p(v=a0|mp(v)\A,v)
      p.a0.mpv <- a0*p.a1.mpv+(1-a0)*(1-p.a1.mpv)
      
      #p(v=a0|mp(v)\A,v)
      p.a1.mpv <- 1-p.a0.mpv
      
      # save the density ratio: p(v|mp(V))|_{a0}/p(v|mp(V))|_{a1}
      assign(paste0("densratio_",v), {p.a0.mpv/p.a1.mpv}/densratio_A)
      
      # save the ratio of p(A|mp(v)\A,v)|_{a0}/p(A|mp(v)\A,v)|_{a1}
      # such that we can come back to update the density ratio of v once we update the densratioA
      assign(paste0("bayes.densratio_",v), {p.a0.mpv/p.a1.mpv})
      
    }
    
    
  } ## Iterate over each variable in M
  
  
} else {
  
  print("Invalid ratio.method.M input.")
  
  stop()
  
}

# Get the densratio vectors based on their names
densratio.vectors.M <- mget(paste0("densratio_",M))

# Create a data frame using all the vectors
densratio.M <- data.frame(densratio.vectors.M)