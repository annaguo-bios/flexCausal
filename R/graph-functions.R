###########################################
# make a graph object.
###########################################

make.graph <- function(vertices, bi_edges, di_edges) {
  
  graph <- list(vertices = vertices, fixed = data.frame(fixed=rep(FALSE, length(vertices)), row.names=vertices), di_edges=di_edges, bi_edges=bi_edges)
  
  return(graph)
}

###########################################
# put treatment to the end of the queue.
###########################################
# Called by the top_orderMAT() function
###########################################
treatment.queue <- function(adj.matrix, queue, treatment) {
  
  # find treatment corresponding to which index number
  treatment_index <- which(rownames(adj.matrix) == treatment)
  
  # Check if treatment is in the queue
  if (treatment_index %in% queue) {
    
    # Remove treatment from its current position
    queue <- queue[queue != treatment_index]
    
    # Append treatment to the end of the queue
    queue <- c(queue, treatment_index)
  }
  
  return(queue)
}

#######################################################################################################################################################
# get topological ordering of graph with adjacency matrix as input. Treatment is ranked as far back as possible by using the treatment.queue function.
########################################################################################################################################################
# Called by the top_order() function.
########################################################################################################################################################

# Function to perform topological sorting using Kahn's algorithm
f.top_orderMAT <- function(adj.matrix, treatment=NULL) {
  
  # Number of nodes
  n <- nrow(adj.matrix)
  
  # Initialize in-degree array
  in_degree <- rowSums(adj.matrix)
  
  # Initialize queue for nodes with in-degree 0
  queue <- which(in_degree == 0)
  
  if (!is.null(treatment)) {
    
    # reorder the queue such that treatment is always at the end
    queue <- treatment.queue(adj.matrix, queue, treatment)
    
  }
  
  
  # Initialize result list for topological order
  top_order <- vector("numeric", n)
  order_index <- 1
  
  while (length(queue) > 0) {
    
    # Dequeue a node
    node <- queue[1]
    queue <- queue[-1]
    
    # Add node to topological order
    top_order[order_index] <- node
    order_index <- order_index + 1
    
    # Update in-degree of adjacent nodes and enqueue nodes with in-degree 0
    for (adj_node in which(adj.matrix[, node] != 0)) {
      
      in_degree[adj_node] <- in_degree[adj_node] - 1
      
      if (in_degree[adj_node] == 0) {
        
        queue <- c(queue, adj_node)
        
        if (!is.null(treatment)) {
          
          # reorder the queue such that treatment is always at the end
          queue <- treatment.queue(adj.matrix, queue, treatment)
          
        }
        
        
      }
      
    }
  }
  
  # Check for cycles
  if (order_index <= n) {
    stop("The graph contains a cycle. Topological ordering not possible. Ordering is only available for directed acyclic graphs (DAGs).")
  }
  
  return(rownames(adj.matrix)[top_order]) 
}

##########################################
# get adjacency matrix from graph.
##########################################
# Called by the top_order() function
##########################################
f.adj_matrix <- function(graph){
  
  # extract vertices and di_edges from the graph
  vertices <- graph$vertices
  di_edges <- graph$di_edges
  
  # Create an empty adjacency matrix
  adj_matrix <- matrix(0, nrow = length(vertices), ncol = length(vertices), dimnames = list(vertices, vertices))
  
  # Update adjacency matrix based on directed edges
  for (edge in di_edges) {
    from <- edge[1]
    to <- edge[2]
    adj_matrix[to, from] <- 1
  }
  
  return(adj_matrix)
}


######################################################################################################################################
# get topological ordering of graphs from vertices and directed edges. Treatment variable is ranked as far back as possible.
######################################################################################################################################
f.top_order <- function(graph, treatment=NULL){
  
  f.top_orderMAT(f.adj_matrix(graph), treatment)

}



#################################################
# get parents of a node OR nodes in a graph
#################################################
# Called by the markov_pillow() function
#################################################

f.parents <- function(graph, nodes){
  
  vertices <- graph$vertices
  
  adj_matrix <- f.adj_matrix(graph)
  
  parents <- list()
  
  for (node in nodes) {
    
    parents[[node]] <- vertices[which(adj_matrix[node, ] != 0)]
    
  }
  
  return(unique(unlist(parents)))
  
}

#################################################
# get the children of a node OR nodes in a graph
#################################################
# Called by the is.p.fix() function
#################################################

f.children <- function(graph, nodes){
  
  vertices <- graph$vertices
  
  adj_matrix <- f.adj_matrix(graph)
  
  children <- list()
  
  for (node in nodes) {
    
    children[[node]] <- vertices[which(adj_matrix[, node] != 0)]
    
  }
  
  return(unique(unlist(children)))
  
}


#################################################
# get the descendants of a node OR nodes in a graph
#################################################
# Called by the reachable_closure() function
#################################################

f.descendants <- function(graph, nodes){

  descendants <- nodes # the set of descendants by default include the nodes themself
  
  # Recursively find descendants
  find_descendants <- function(nodes) {
    
    if (length(nodes) == 0) {

      return(NULL)

    }

    children <- f.children(graph, nodes)  # Get children of current nodes
    
    descendants <<- union(descendants, children)  # Add children to the list of descendants
    
    find_descendants(children)  # Recursively find descendants of children
  }
  
  # Start the recursive search
  find_descendants(nodes)
  
  return(descendants)
}

#################################################
# get district of a node in a graph
#################################################
# Called by the markov_pillow() function
#################################################

f.district <- function(graph, node){
  
  # extract bi_edges from the graph
  bi_edges <- graph$bi_edges
  
  connected_nodes <- c(node)
  
  for (edge in bi_edges) { # iterate over all bidirected edges
    
    if (node %in% edge) { # find out whether the given node is in this bidirected edge
      connected_nodes <- c(connected_nodes, edge[edge != node])
    }
    
  }
  
  return(unique(connected_nodes))
  

  
}

#################################################
# count the number of districts in a graph
#################################################
# Called by the is.np.saturated() function
#################################################

count.districts <- function(graph) {
  
  vertices <- graph$vertices
  di_edges <- graph$di_edges
  bi_edges <- graph$bi_edges
  fixed <- graph$fixed
  
  districts <- list()  # Initialize a list to store districts
  
  # Iterate over each vertex which is not fixed
  for (v in vertices[fixed[vertices,"fixed"]==F]) {
    
    # Add the current vertex itself to its district
    district.v <- c(v)
    
    # Find neighbors connected via bidirected edges
    siblings <- unlist(lapply(bi_edges, function(x) if (v %in% x) setdiff(x, v) else NULL))
    
    # Add neighbors to the district
    district.v <- c(district.v, siblings)
    
    # Store the district for the current vertex
    districts[[v]] <- sort(unique(district.v))
  }
  
  # Count the number of unique districts
  num_districts <- length(unique(districts))
  
  return(list(n.districts=num_districts, districts=districts))
}

#################################################
# get markov blanket of a node in a graph
#################################################
f.markov_blanket <- function(graph, node){
  
  # get the district of the node
  dist <- f.district(graph, node)
  
  # get the union of district and parents of the district
  union <- unique(c(dist, f.parents(graph, dist)))
  
  # get the markov pillow of the node, which is the subset of union that proceed the node in the topological ordering
  return(union)
  
}

#################################################
# get markov pillow of a node in a graph
#################################################
f.markov_pillow <- function(graph, node, treatment=NULL){

  # get topological ordering of the graph
  t.order <- f.top_order(graph, treatment)
  
  # get the markov blanket of the node
  union <- f.markov_blanket(graph, node)
  
  # get the markov pillow of the node, which is the subset of union that proceed the node in the topological ordering
  return(union[union %in% t.order[1:which(t.order == node)]])
  
}


#################################################
# whether the treatment is primal fixable
#################################################
is.p.fix <- function(graph, treatment){
  
  # get the children of the treatment
  ch <- f.children(graph, treatment)
  
  # get the bidirected edges of the graph
  bi_edges <- graph$bi_edges
  
  # whether there is no bidirected edge between the treatment and any of its children
  for (child in ch) { # loop over all the children of the treatment
    
    if (c(treatment, child) %in% bi_edges || c(child, treatment) %in% bi_edges) {
      return(FALSE)  # Bidirected edge found between node and child
    }
    
  }
  return(TRUE)  # No bidirected edge found between node and any child

  
}


#################################################
# perform the fix operation on a graph
#################################################
# Called by the reachable_closure() function
#################################################
fix <- function(graph, nodes) {

  # Iterate over nodes
  for (n in nodes) {
    
    # Mark node as fixed
    graph$fixed[n,"fixed"] <- TRUE
    
    # Delete incoming directed edges
    parents <- f.parents(graph, n)
    
    di_edges.remove <- lapply(parents, function(x) c(x, n))
    
    graph$di_edges <- setdiff(graph$di_edges, di_edges.remove)
    
    # Delete bidirected edges
    bi_edges.remove <- lapply(graph$bi_edges, function(x) if (n %in% x) x) # find all bidirected edges that contain the node

    graph$bi_edges <- setdiff(graph$bi_edges, bi_edges.remove)
    
  }
  
  return(graph)
}



##################################################
# define reachable closure
##################################################
f.reachable_closure <- function(graph, nodes) {
  
  # initialize set of vertices that can potentially be fixed
  remaining_vertices <- setdiff(graph$vertices, nodes)
  
  fixing_order <- character(0)
  
  fixed <- TRUE
  
  while (length(remaining_vertices) > 0 & fixed) {
    
    fixed <- FALSE
    
    for (v in remaining_vertices) {
      
      if (length(intersect(f.descendants(graph, v), f.district(graph, v) )) == 1) { # if there is no bi-directed edge between v and any of its descendants, then fix v
        
        graph <- fix(graph, v)
        
        remaining_vertices <- setdiff(remaining_vertices, v) # remove v from the remaining vertices
        
        fixing_order <- c(fixing_order, v)
        
        fixed <- TRUE
        
        break # breaks the for loop
        
      } # end of if statement
      
    } # end of for loop
    
  }
  
  reachable_closure <- setdiff(graph$vertices, fixing_order)
  
  return(list(reachable_closure = reachable_closure, fixing_order = fixing_order, graph = graph))
}


##################################################################
# whether the graph is nonparametrically saturated.
##################################################################
is.np.saturated <- function(graph) {
  
  # get the topological ordering of the graph
  top_order <- f.top_order(graph)

  # Iterate over all pairs of vertices
  for (pair in combn(graph$vertices, 2, simplify = FALSE)) {
    
    Vi <- pair[1]
    Vj <- pair[2]
    
    
    # order Vi and Vj
    if (which(top_order == Vi) > which(top_order == Vj)) {
      V1 <- Vj
      V2 <- Vi
    } else {
      V1 <- Vi
      V2 <- Vj
    }
    
    # Check if there is no dense inducing path between Vi and Vj
    
    # 1. V1 is not in the parents set of Di for all Di in D, where D is the district of V2 in conditional acyclic directed mixed graphs (CADMG) obtained by recursively fixing as many vertices as possible in V/V2
    # AND 2. The CADMG obtained by recursively fixing as many vertices as possible in V/{Vi, Vj} has more than one district.
    
    # (count.districts(f.reachable_closure(graph,c(Vi,Vj))[[3]])$n.districts > 1) 
    # !(V1 %in% f.district(f.reachable_closure(graph,c(V1,V2))[[3]], V2))
    
    if (!(V1 %in% f.parents(graph, f.district(f.reachable_closure(graph,V2)[[3]], V2))) &&  (count.districts(f.reachable_closure(graph,c(Vi,Vj))[[3]])$n.districts > 1)){
      
      print("The graph is not nonparametrically saturated.")
      print(c(V1,V2))
      return(FALSE)
      
    } # end of if statement
  
  } # end of for loop
  
  print("The graph is nonparametrically saturated.")
  
  return(TRUE)
}





#######################################################################################
# whether the graph is mb-shield, which is the graph contains only ordinary constrains.
########################################################################################
is.mb.shielded <- function(graph) {
  
  # Iterate over all pairs of vertices
  for (pair in combn(graph$vertices, 2, simplify = F)) {
    
    Vi <- pair[1]
    Vj <- pair[2]
    
    if (!( list(c(Vi,Vj)) %in% graph$bi_edges || list(c(Vj,Vi)) %in% graph$bi_edges || list(c(Vi,Vj)) %in% graph$di_edges || list(c(Vj,Vi)) %in% graph$di_edges) ) { # edge between Vi and Vj is absent in G
      
      if (Vi %in% f.markov_blanket(graph, Vj) || Vj %in% f.markov_blanket(graph, Vi)) { # Vi is in the Markov blanket of Vj or Vj is in the Markov blanket of Vi
        
        print("The graph is not mb-shield.")
        
        return(FALSE)
        
      } # end of if statement
      
    } # end of if statement

  } # end of for loop
  
  print("The graph is mb-shielded.")
  return(TRUE)
}
