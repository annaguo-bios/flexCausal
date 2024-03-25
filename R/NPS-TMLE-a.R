NPS.TMLE.a <- function(a,data,vertices, di_edges, bi_edges, treatment, outcome,
                     onestep=T, 
                     superlearner.seq = F, # whether run superlearner for sequential regression
                     superlearner.Y=F, # whether run superlearner for outcome regression
                     superlearner.A=F, # whether run superlearner for propensity score
                     superlearner.M=F, # whether run superlearner for estimating densratio for M using bayes method
                     superlearner.L=F, # whether run superlearner for estimating densratio for L using bayes method
                     crossfit=F, K=5,
                     ratio.method.L="bayes", # method for estimating the density ratio associated with M
                     ratio.method.M="bayes", # method for estimating the density ratio associated with L
                     lib.seq = c("SL.glm","SL.earth","SL.ranger","SL.mean"), # superlearner library for sequential regression
                     lib.L = c("SL.glm","SL.earth","SL.ranger","SL.mean"), # superlearner library for density ratio estimation via bayes rule for variables in L
                     lib.M = c("SL.glm","SL.earth","SL.ranger","SL.mean"), # superlearner library for density ratio estimation via bayes rule for variables in M
                     lib.Y = c("SL.glm","SL.earth","SL.ranger","SL.mean"), # superlearner library for outcome regression
                     lib.A = c("SL.glm","SL.earth","SL.ranger","SL.mean"), # superlearner library for propensity score
                     formulaY="Y ~ .", formulaA="A ~ .", # regression formula for outcome regression and propensity score if superlearner is not used
                     linkY_binary="logit", linkA="logit", # link function for outcome regression and propensity score if superlearner is not used
                     n.iter=500, cvg.criteria=0.01,
                     truncate_lower=0, truncate_upper=1){
  
  # attach(data, warn.conflicts=FALSE)
  
  n <- nrow(data)
  
  a0 <- a
  a1 <- 1-a
  
  # return topological ordering
  tau <- top_order(vertices, di_edges, bi_edges)
  tau.df <- data.frame(tau=tau, order = 1:length(tau))
  
  # Get set C, M, L
  setCML <- CML(vertices, di_edges, bi_edges, treatment) # get set C, M, L
  rm(order) # this function from python is required for running the top_order but it contradicts with order in R, so remove it
  
  C <- setCML$C # everything comes before the treatment following topolofical order tau
  
  L <- setCML$L # variables within the district of treatment and comes after the treatment (including the treatment itself) following topolofical order tau
  
  M <- setCML$M # everything else
  
  # re-order vertices according to their topological order in tau
  
  C <- rerank(C, tau) # re-order vertices in C according to their topological order in tau
  
  L <- rerank(L, tau) # re-order vertices in L according to their topological order in tau
  L.removedA <- L[L!=treatment] # remove treatment from L
  
  M <- rerank(M, tau) # re-order vertices in M according to their topological order in tau
  
  # Variables
  A <- data[,treatment] # treatment
  Y <- data[,outcome] # outcome
  
  
  ##################################################################
  ## TMLE initialization for sequential regression based estimator
  ##################################################################
  
  source("NPS-01-initial-nuisance-estimate.R")
  
  ##################################################################
  #################### One-step estimator ##########################
  ##################################################################

  source("NPS-02-onestep-estimator.R")
  
  ##################################################################
  #################### Sequential regression based TMLE ############
  ##################################################################
  
  source("NPS-03-tmle-estimator.R")
  
  return(list(TMLE=tmle.out,Onestep=onestep.out))
  
}

