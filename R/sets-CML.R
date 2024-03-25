# Extract the sets C, M, L from an ADMG givan a treatment

CML <- function(vertices, di_edges, bi_edges, treatment){
  
  # return topological ordering
  tau <- top_order(vertices, di_edges, bi_edges)
  
  # Get set C, M, L
  C <- tau[1:(which(tau==treatment)-1)] # everything comes before the treatment following topolofical order tau
  
  L <- intersect(district(vertices, di_edges, bi_edges, treatment),
                 tau[(which(tau==treatment)):length(tau)]) # variables within the district of treatment and comes after the treatment (including the treatment itself) following topolofical order tau
  
  M <- setdiff(vertices,c(C,L)) # everything else
  
  return(list(C=C, M=M, L=L))
}
