# Rerank the target list according to the reference list
# In our case, we re-rank the vector C, M, L according to topological order tau

rerank <- function(target,reference){
  
  # Find the indices of elements in M within V
  indices <- match(target, reference)
  
  # Reorder M according to the indices
  re_ranked_target <- target[order(indices)]
  
  return(re_ranked_target)
}

