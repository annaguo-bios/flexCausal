## Test the package with front-door


# call python code
setwd("/Users/apple/Library/CloudStorage/Dropbox/primal-fixability/code/ADMGtmle")


library(reticulate)
library(dplyr)
library(SuperLearner)
library(densratio)
library(itertools)


use_condaenv("phd")

source_python("pre-functions-ananke.py")

# rm(order) # this function from python contradicts with order in R, so remove it

# load functions
source("/Users/apple/Library/CloudStorage/Dropbox/primal-fixability/code/ADMGtmle/sets-CML.R")
source("/Users/apple/Library/CloudStorage/Dropbox/primal-fixability/code/ADMGtmle/rerank.R")



# define graph for example a
load("/Users/apple/Library/CloudStorage/Dropbox/primal-fixability/code/ADMGtmle/example-a.RData")
vertices = c('A','M','Y','X','L')
di_edges = list(c('A', 'M'), c('M', 'L'), c('L', 'Y'), c('X','A'),c('X', 'M'), c('X','L'), c('X', 'Y'))
bi_edges = list(c('A','L'),c('L','Y'))

# define graph for example b
load("/Users/apple/Library/CloudStorage/Dropbox/primal-fixability/code/ADMGtmle/example-b.RData")
vertices = c('A','M','Y','X','L')
di_edges = list(c('A', 'M'),c('A','Y')  , c('M', 'L'), c('L', 'Y'), c('X','A'),c('X', 'M'), c('X','L'), c('X', 'Y'))
bi_edges = list(c('A','L'),c('M','Y'))

# define graph for example frontdoor
load("/Users/apple/Library/CloudStorage/Dropbox/primal-fixability/code/ADMGtmle/example-frontdoor.RData")
vertices = c('A','M','Y','X')
di_edges = list(c('A', 'M'), c('M', 'Y'), c('X','A'),c('X', 'M'), c('X', 'Y'))
bi_edges = list(c('A','Y'))

treatment='A'
outcome='Y'


a=1
data=data.b
vertices=vertices
di_edges=di_edges
bi_edges=bi_edges
treatment=treatment
outcome=outcome
onestep=T
superlearner.seq = F # whether run superlearner for sequential regression
superlearner.Y=F # whether run superlearner for outcome regression
superlearner.A=F # whether run superlearner for propensity score
superlearner.M=F # whether run superlearner for estimating densratio for M using bayes method
superlearner.L=F # whether run superlearner for estimating densratio for L using bayes method
crossfit=F
K=5
ratio.method.L="bayes" # method for estimating the density ratio associated with M
ratio.method.M="bayes" # method for estimating the density ratio associated with L
lib.seq = c("SL.glm","SL.earth","SL.ranger","SL.mean") # superlearner library for sequential regression
lib.L = c("SL.glm","SL.earth","SL.ranger","SL.mean") # superlearner library for density ratio estimation via bayes rule for variables in L
lib.M = c("SL.glm","SL.earth","SL.ranger","SL.mean") # superlearner library for density ratio estimation via bayes rule for variables in M
lib.Y = c("SL.glm","SL.earth","SL.ranger","SL.mean") # superlearner library for outcome regression
lib.A = c("SL.glm","SL.earth","SL.ranger","SL.mean") # superlearner library for propensity score
formulaY="Y ~ ."
formulaA="A ~ ." # regression formula for outcome regression and propensity score if superlearner is not used
linkY_binary="logit"
linkA="logit" # link function for outcome regression and propensity score if superlearner is not used
n.iter=500
cvg.criteria=0.01
truncate_lower=0
truncate_upper=1



factorial <- function(n) {
  if (n == 0 || n == 1) {
    return(1)
  } else {
    return(n * factorial(n - 1))
  }
}

# Test the factorial function
factorial(5) 



des <- f.descendants(vertices=c("A","B","C","D"), di_edges=list(c('A','B'), c('B','C'), c('C','D'), c('A','C'), c('B','D')), nodes='A')
vertices=c("A","B","C","D"); di_edges=list(c('A','B'), c('B','C'), c('C','D'), c('A','C'), c('B','D')); nodes='A'


children <- f.children(vertices=c("A","B","C","D"), di_edges=list(c('A','B'), c('B','C'), c('C','D'), c('A','C'), c('B','D')), nodes='A')

m <- data.frame(fixed=rep(F,3), row.names=c("A","B","C"))

m[c("A","B"),"fixed"]


G <- list(vertices = data.frame(fixed=rep(FALSE, length(vertices)), row.names=vertices), di_edges=di_edges, bi_edges=bi_edges)


graph <- make.graph(vertices=c('T','M','L','Y'), bi_edges=list(c('M','Y'),c('L','Y')), di_edges=list(c('T','Y'), c('T','M'), c('T','L'), c('M','L'), c('L','Y')))

tmp <- is.np.saturated(graph)

graph <- make.graph(vertices=c('A','M','L','Y','X'), bi_edges=list(c('A','Y')), di_edges=list(c('X','A'), c('X','M'), c('X','L'),c('X','Y'), c('M','Y'), c('A','M'), c('A','L'), c('M','L'), c('L','Y')))
tmp <- is.np.saturated(make.graph(vertices=c('A','M','L','Y','X'), bi_edges=list(c('A','Y')), di_edges=list(c('X','A'), c('X','M'), c('X','L'),c('X','Y'), c('M','Y'), c('A','M'), c('A','L'), c('M','L'), c('L','Y'))))
tmp <- is.mb.shielded(make.graph(vertices=c('A','M','L','Y','X'), bi_edges=list(c('A','Y')), di_edges=list(c('X','A'), c('X','M'), c('X','L'),c('X','Y'), c('M','Y'), c('A','M'), c('A','L'), c('M','L'), c('L','Y'))))


tmp <- count_districts(make.graph(vertices=c('A','M','L','Y','X'), bi_edges=list(c('A','Y')), di_edges=list(c('X','A'), c('X','M'), c('X','L'),c('X','Y'), c('M','Y'), c('A','M'), c('A','L'), c('M','L'), c('L','Y'))))

graph <- make.graph(vertices=c('T','M','L','Y'), bi_edges=list(c('M','Y'),c('L','Y')), di_edges=list(c('T','M'),  c('T','L'),c('M','L'), c('L','Y')))
tmp <- is.np.saturated(make.graph(vertices=c('T','M','L','Y'), bi_edges=list(c('M','Y'),c('L','Y')), di_edges=list(c('T','M'),  c('T','L'),c('M','L'), c('L','Y'))))

