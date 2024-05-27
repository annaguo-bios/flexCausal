#' @title Create graph object.
#' @description This function create a graph object that can be used in other functions in this package.
#' @param vertices A character vector of vertices in the graph.
#' @param bi_edges A list of vectors that record the bidirectional edges in the graph. For example, `bi_edges=list(c('A','B'))` means there is a bidirectional edge between vertex A and B.
#' @param di_edges A list of vectors that record the directed edges in the graph. For example, `di_edges=list(c('A','B'))` means there is a directed edge from vertex A to B.
#' @param multivariate.variables A list of variables that are multivariate in the causal graph.
#' For example, `multivariate.variables=list(M=c('M1,'M2'))` means M is bivariate and the corresponding columns in the dataframe are M1 and M2. Default is \code{NULL}.
#' @keywords graph ADMG
#' @return A graph object with the following components:
#' \describe{
#'       \item{\code{vertices}}{Equivalent to the input argument \code{vertices}.}
#'       \item{\code{fixed}}{A data frame with a column \code{fixed} that indicates whether the vertex is fixed or not. The vertices is not fixed initially.}
#'       \item{\code{bi_edges}}{Equivalent to the input argument \code{bi_edges}.}
#'       \item{\code{di_edges}}{Equivalent to the input argument \code{di_edges}.}
#'       \item{\code{multivariate.variables}}{A list of variables that are multivariate in the causal graph. For example, `multivariate.variables=list(M=c('M1,'M2'))` means M is bivariate and the corresponding columns in the dataframe are M1 and M2. Default is \code{NULL}.}}
#' @export
#' @examples
#' make.graph(vertices=c('A','M','L','Y','X'),
#' bi_edges=list(c('A','Y')),
#' di_edges=list(c('X','A'), c('X','M'), c('X','L'),
#' c('X','Y'), c('M','Y'), c('A','M'), c('A','L'), c('M','L'), c('L','Y')))
#'
make.graph <- function(vertices, bi_edges, di_edges, multivariate.variables=NULL) {

  graph <- list(vertices = vertices, fixed = data.frame(fixed=rep(FALSE, length(vertices)), row.names=vertices), di_edges=di_edges, bi_edges=bi_edges, multivariate.variables=multivariate.variables)

  return(graph)
}

###########################################
# put treatment to the end of the queue.
###########################################
# Called by the f.top_orderMAT() function
###########################################


#' Put treatment to the end of the queue.
#'
#' This function is called by the `f.top_orderMAT()` function to put the treatment to the end of the queue when employing Kahn's algorithm for topological ordering.
#' In deriving the topological ordering, there can be multiple vertices with the same in-degree. This function ensures that the treatment is ranked as far back as possible when there are ties.
#' @param adj.matrix An adjacency matrix of the graph.
#' @param queue A vector of integers that records the current queue of vertices.
#' @param treatment A character string indicating the treatment variable.
#' @keywords graph ADMG
#' @return A updated queue with the treatment variable at the end if it is in the queue.
#' @export
#'
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


##############################################################################
# Replace element in a specific position of vector A with another vector B
##############################################################################
# Called by the replace.vector() function
##############################################################################

#' Replace element at a specific position.
#'
#' This function replace element at a specific position of a vector called vec with another vector called value.
#' @param vec A vector.
#' @param pos An integer indicating the position to replace element in vec with value.
#' @param value A vector to insert into vec.
#' @keywords vector replace element
#' @return A vector with value inserted at position pos.
#' @export
#' @examples
#' insert(c('A','B','C'), 2, 'D')
insert <- function(vec, pos, value) {

  if (pos == 1) { # if the position is the first element of the vector

    return(c(value, vec[-1]))

  } else if (pos == length(vec)) { # if the position is the last element of the vector

    return(c(vec[-length(vec)], value))

  } else { # if the position is in the middle of the vector
    return(c(vec[1:(pos-1)], value, vec[(pos+1):length(vec)]))
  }

}


##############################################################################
# Replace element in a specific position of vector A with vector under
##############################################################################

#' Replace elements at a vector with other vectors.
#'
#' This function replace the elements in vec with the vectors in multivariate.variables under the same name.
#' @param vec A vector.
#' @param multivariate.variables A list of vectors to replace the elements in vec.
#' @keywords vector replace element
#' @return A vector with elements replaced by the vectors in multivariate.variables.
#' @export
#' @examples
#' replace.vector(c('A','B','C'), list(A=c('D','E')))
replace.vector <- function(vec, multivariate.variables=NULL) {

  if (!is.null(multivariate.variables)) { # if multivariate variables are provided, replace the variable names in the order vector with it's individual component

    for (var_name in names(multivariate.variables)) {

      replacements <- multivariate.variables[[var_name]]

      if (var_name %in% vec) {

        vec <- insert(vec, which(vec == var_name), replacements)

      } # end of if condition

    } # end of for loop

  }

  return(vec)

}

#######################################################################################################################################################
# get topological ordering of graph with adjacency matrix as input. Treatment is ranked as far back as possible by using the treatment.queue function.
########################################################################################################################################################
# Called by the top_order() function.
########################################################################################################################################################

#' Topological ordering of a graph from adjacency matrix.
#'
#' Function to perform topological sorting using Kahn's algorithm. The treatment variable is ranked with the largest order possible if there are ties.
#' @param adj.matrix An adjacency matrix of the graph.
#' @param treatment A character string indicating the treatment variable.
#' @keywords graph ADMG ordering
#' @return A vector of vertices ranked with rank from small to large.
#' @export
#'
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

  order <- rownames(adj.matrix)[top_order]

  return(order)
}

##########################################
# get adjacency matrix from graph.
##########################################
# Called by the f.top_order() function
##########################################

#' Get adjacency matrix from graph.
#'
#' Function to perform topological sorting using Kahn's algorithm. The treatment variable is ranked with the largest order possible if there are ties.
#' @param graph A graph object generated by the `make.graph()` function.
#' @keywords graph ADMG ordering
#' @return An adjacency matrix of the graph.
#' @export
#' @examples
#' graph <- make.graph(vertices=c('A','M','L','Y','X'),
#' bi_edges=list(c('A','Y')),
#' di_edges=list(c('X','A'), c('X','M'), c('X','L'),
#' c('X','Y'), c('M','Y'), c('A','M'), c('A','L'), c('M','L'), c('L','Y')))
#' f.adj_matrix(graph)
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

#' Get the topological ordering of a graph from graph object.
#'
#' Function to extract the adjacency matrix from a graph object.
#' @param graph A graph object generated by the `make.graph()` function.
#' @param treatment A character string indicating the treatment variable.
#' If NULL, this function will rank vertices according to their input order in the vertices vector when there are ties.
#' @keywords graph ADMG adjancency.matrix
#' @return A vector of vertices ranked with rank from small to large.
#' @export
#' @examples
#' graph <- make.graph(vertices=c('A','M','L','Y','X'),
#' bi_edges=list(c('A','Y')),
#' di_edges=list(c('X','A'), c('X','M'), c('X','L'),
#' c('X','Y'), c('M','Y'), c('A','M'), c('A','L'), c('M','L'), c('L','Y')))
#' f.top_order(graph)
#'
f.top_order <- function(graph, treatment=NULL){

  f.top_orderMAT(f.adj_matrix(graph), treatment)

}



#################################################
# get parents of a node OR nodes in a graph
#################################################
# Called by the markov_pillow() function
#################################################

#' Get the parents of a node OR nodes in a graph.
#'
#' Function to extract the parents of a node OR nodes in a graph object.
#' Parents of a node are the nodes that have directed edges pointing to the node.
#' @param graph A graph object generated by the `make.graph()` function.
#' @param nodes A character vector of nodes for which to extract parents.
#' @keywords graph ADMG parents
#' @return A vector of vertices contains parents set of the given nodes.
#' @export
#' @examples
#' graph <- make.graph(vertices=c('A','M','L','Y','X'),
#' bi_edges=list(c('A','Y')),
#' di_edges=list(c('X','A'), c('X','M'), c('X','L'),
#' c('X','Y'), c('M','Y'), c('A','M'), c('A','L'), c('M','L'), c('L','Y')))
#' f.parents(graph, c('Y','L'))
#'
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

#' Get the children of a node OR nodes in a graph.
#'
#' Function to extract the children of a node OR nodes in a graph object.
#' Children of a node are the nodes that have edges from the given node.
#' @param graph A graph object generated by the `make.graph()` function.
#' @param nodes A character vector of nodes for which to extract children.
#' @keywords graph ADMG children
#' @return A vector of vertices contains children set of the given nodes.
#' @export
#' @examples
#' graph <- make.graph(vertices=c('A','M','L','Y','X'),
#' bi_edges=list(c('A','Y')),
#' di_edges=list(c('X','A'), c('X','M'), c('X','L'),
#' c('X','Y'), c('M','Y'), c('A','M'), c('A','L'), c('M','L'), c('L','Y')))
#' f.children(graph, c('A'))
#'
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

#' Get the descendants of a node OR nodes in a graph.
#'
#' Function to extract the descendants of a node OR nodes in a graph object.
#' Descendants of a node Vi are set Vj such that there is a directed path Vi->...->Vj. Descendants set including Vi itself by convention.
#' @param graph A graph object generated by the `make.graph()` function.
#' @param nodes A character vector of nodes for which to extract children.
#' @keywords graph ADMG descendants
#' @return A vector of vertices contains descendants set of the given nodes.
#' @export
#' @examples
#' graph <- make.graph(vertices=c('A','M','L','Y','X'),
#' bi_edges=list(c('A','Y')),
#' di_edges=list(c('X','A'), c('X','M'), c('X','L'),
#' c('X','Y'), c('M','Y'), c('A','M'), c('A','L'), c('M','L'), c('L','Y')))
#' f.descendants(graph, c('A'))
#'
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
# Called by the f.markov_pillow() function
#################################################

#' Get the district of a vertex in a graph.
#'
#' Function to extract the name of vertices that is in the district of a given vertex in a graph object.
#' District of a unfixed vertex Vi is the set of vertices that are connected to Vi by bidirected edges, including Vi itself by convention.
#' @param graph A graph object generated by the `make.graph()` function.
#' @param node A character string of a vertex for which to extract district.
#' @keywords graph ADMG district
#' @return A vector of vertices that is in the district of the given vertex.
#' @export
#' @examples
#' graph <- make.graph(vertices=c('A','M','L','Y','X'),
#' bi_edges=list(c('A','Y')),
#' di_edges=list(c('X','A'), c('X','M'), c('X','L'),
#' c('X','Y'), c('M','Y'), c('A','M'), c('A','L'), c('M','L'), c('L','Y')))
#' f.district(graph, c('A'))
#'
f.district <- function(graph, node){

  # extract bi_edges from the graph
  bi_edges <- graph$bi_edges

  connected_nodes <- c(node)

  # Recursively find descendants
  find_district <- function(nodes) {

    if (length(nodes) == 0) {

      return(NULL)

    }

    new_connected_nodes <- c() # Initialize a new set of connected nodes

    for (node in nodes){

      for (edge in bi_edges) { # iterate over all bidirected edges

        if (node %in% edge) { # find out whether the given node is in this bidirected edge

          if (!(edge[edge!=node] %in% connected_nodes)){ # the new node is not yet in the connected nodes set

            new_connected_nodes <- c(new_connected_nodes, edge[edge != node])

          }



        }

      }

    }

    new_connected_nodes <- unique(new_connected_nodes)

    connected_nodes <<- union(connected_nodes, new_connected_nodes) # update connected_nodes in the global environment

    find_district(new_connected_nodes)  # Recursively find district of the new nodes
  }

  # Start the recursive search
  find_district(node)

  return(connected_nodes)
}

#################################################
# count the number of districts in a graph
#################################################
# Called by the is.np.saturated() function
#################################################

#' Get the number of districts in a graph.
#'
#' Function to count the number of districts in a graph object.
#' Number of districts in a graph is the number of bidirected connected components in the graph.
#' A vertex that is not bidirected connected to any other vertex is also considered as a district by convention.
#' For example, a graph like A B <-> C has two districts: `{{A, B}}` and `{{C}`.
#' @param graph A graph object generated by the `make.graph()` function.
#' @keywords graph ADMG district
#' @return A list with the following two components:
#' \describe{
#'       \item{districts}{A list of vectors, each vector contains the vertices in a district.}
#'       \item{n.districts}{The number of districts in the graph.}}
#' @export
#'
cnt.districts <- function(graph) {

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

#' Get the Markov blanket of a vertex in a graph.
#'
#' Function to get the Markov blanket of a vertex in a graph object.
#' Markov blanket of a vertex Vi is the union of vertices that is in the district of Vi and their parents set.
#' Mb= union(district(Vi), parents(district(Vi))). THIS IS PSUEDO CODE.
#' @param graph A graph object generated by the `make.graph()` function.
#' @param node A character string of a vertex for which to extract Markov blanket.
#' @keywords graph ADMG markov blanket
#' @return A vector of vertices that is in the Markov blanket of the given vertex.
#' @export
#' @examples
#' graph <- make.graph(vertices=c('A','M','L','Y','X'),
#' bi_edges=list(c('A','Y')),
#' di_edges=list(c('X','A'), c('X','M'), c('X','L'),
#' c('X','Y'), c('M','Y'), c('A','M'), c('A','L'), c('M','L'), c('L','Y')))
#' f.markov_blanket(graph, 'A')
#'
f.markov_blanket <- function(graph, node){

  # get the district of the node
  dist <- f.district(graph, node)

  # get the union of district and parents of the district
  union <- unique(c(dist, f.parents(graph, dist)))

  # get the markov pillow of the node, which is the subset of union that proceed the node in the topological ordering
  return(setdiff(union, node))

}

#################################################
# get markov pillow of a node in a graph
#################################################

#' Get the Markov pillow of a vertex in a graph.
#'
#' Function to get the Markov pillow of a vertex in a graph object.
#' Markov pillow of a vertex Vi is the subset of the Markov blanket of Vi that proceed Vi in the topological ordering of the graph.
#' Mp= `{{Vj in union(district(Vi), parents(district(Vi))): Vj proceed Vi}}`. THIS IS PSUEDO CODE.
#' @param graph A graph object generated by the `make.graph()` function.
#' @param node A character string of a vertex for which to extract Markov pillow.
#' @param treatment A character string specifying the treatment variable in the graph.
#' @keywords graph ADMG markov pillow
#' @return A vector of vertices that is in the Markov pillow of the given vertex.
#' @export
#' @examples
#' graph <- make.graph(vertices=c('A','M','L','Y','X'),
#' bi_edges=list(c('A','Y')),
#' di_edges=list(c('X','A'), c('X','M'), c('X','L'),
#' c('X','Y'), c('M','Y'), c('A','M'), c('A','L'), c('M','L'), c('L','Y')))
#' f.markov_pillow(graph, 'A')
#'
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

#' Primal fixability of a treatment variable in a graph.
#'
#' Function to check if a treatment variable is primal fixable in a graph object.
#' If the treatment is primal fixable, then the average causal effect of the treatment on any choice of the outcome in the given graph is always identified.
#' @param graph A graph object generated by the `make.graph()` function.
#' @param treatment A character string specifying the treatment variable in the graph.
#' @keywords graph ADMG primal fixable
#' @return A logical value indicating whether the treatment is primal fixable.
#' @export
#' @examples
#' graph <- make.graph(vertices=c('A','M','L','Y','X'),
#' bi_edges=list(c('A','Y')),
#' di_edges=list(c('X','A'), c('X','M'), c('X','L'),
#' c('X','Y'), c('M','Y'), c('A','M'), c('A','L'), c('M','L'), c('L','Y')))
#' is.p.fix(graph, 'A')
#'
is.p.fix <- function(graph, treatment){

  # get the children of the treatment
  ch <- f.children(graph, treatment)

  # get the bidirected edges of the graph
  bi_edges <- graph$bi_edges

  # whether there is no bidirected edge between the treatment and any of its children
  for (child in ch) { # loop over all the children of the treatment

    if (list(c(treatment, child)) %in% bi_edges || list(c(child, treatment)) %in% bi_edges) {

      print("The treatment is not primal fixable in the provided graph.")

      return(FALSE)  # Bidirected edge found between node and child
    }

  }

  print("The treatment is primal fixable in the provided graph.")
  return(TRUE)  # No bidirected edge found between node and any child


}




#################################################
# whether the treatment is fixable
#################################################

#' Fixability of a treatment variable in a graph.
#'
#' Function to check if a treatment variable is fixable in a graph object.
#' If the treatment is fixable, then the average causal effect of the treatment on any choice of the outcome in the given graph is always identified via backdoor adjustment.
#' @param graph A graph object generated by the `make.graph()` function.
#' @param treatment A character string specifying the treatment variable in the graph.
#' @keywords graph ADMG primal fixable
#' @return A logical value indicating whether the treatment is primal fixable.
#' @export
#' @examples
#' graph <- make.graph(vertices=c('A','M','L','Y','X'),
#' bi_edges=list(c('A','Y')),
#' di_edges=list(c('X','A'), c('X','M'), c('X','L'),
#' c('X','Y'), c('M','Y'), c('A','M'), c('A','L'), c('M','L'), c('L','Y')))
#' is.p.fix(graph, 'A')
#'
is.fix <- function(graph, treatment){

  # get the children of the treatment
  de <- f.descendants(graph, treatment)

  # get the bidirected edges of the graph
  bi_edges <- graph$bi_edges

  # whether there is no bidirected edge between the treatment and any of its children
  for (d in de) { # loop over all the children of the treatment

    if (list(c(treatment, d)) %in% bi_edges || list(c(d, treatment)) %in% bi_edges) {

      print("The treatment is not fixable in the provided graph.")

      return(FALSE)  # Bidirected edge found between node and descendants
    }

  }

  print("The treatment is primal fixable in the provided graph.")
  return(TRUE)  # No bidirected edge found between node and any descendant


}


#################################################
# perform the fix operation on a graph
#################################################
# Called by the reachable_closure() function
#################################################

#' Perform fixing operation on the given vertices in a graph.
#'
#' Function to perform the fixing operation on the given vertices in a graph object.
#' Fixing a vertex in a graph means deleting all the incoming directed edges to the vertex and deleting all the bidirected edges concerning this vertex. The vertex is also marked as fixed.
#' This function is called by the `f.reachable_closure()` function.
#' @param graph A graph object generated by the `make.graph()` function.
#' @param nodes A chacater vector of vertices to fix in the graph.
#' @keywords graph ADMG fix
#' @return A graph object with the given vertices fixed.
#' @export
#' @examples
#' graph <- make.graph(vertices=c('A','M','L','Y','X'),
#' bi_edges=list(c('A','Y')),
#' di_edges=list(c('X','A'), c('X','M'), c('X','L'),
#' c('X','Y'), c('M','Y'), c('A','M'), c('A','L'), c('M','L'), c('L','Y')))
#' fix(graph, 'A')
#'
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

#' Reachable closure of a set of vertices in a graph.
#'
#' Function to return the reachable closure of a set of vertices in a graph object.
#' First obtain a Conditional ADMG (CADMG) via recursively fixing as many vertices as possible in the set of all vertices (V) excluding the set of vertices specified by the `nodes` parameter (S), i.e. V \ S.
#' The reachable closure is the subset of V \ S, where each vertex is not fixable even upon fixing other vertices.
#' @param graph A graph object generated by the `make.graph()` function.
#' @param nodes A character vector of vertices.
#' @keywords graph ADMG reachable closure
#' @return A list containing the following components:
#' \describe{
#'  \item{reachable_closure}{A character vector containing the reachable closure of the given vertices.}
#'  \item{fixing_order}{A character vector of vertices telling the order in which the vertices were fixed.}
#'  \item{graph}{The CADMG obtained via recursively fixing as many vertices as possible in the set of all vertices (V) excluding the set of vertices specified by the `nodes` parameter (S), i.e. V \ S. }}
#' @export
#' @examples
#' graph <- make.graph(vertices=c('A','M','L','Y','X'),
#' bi_edges=list(c('A','Y')),
#' di_edges=list(c('X','A'), c('X','M'), c('X','L'),
#' c('X','Y'), c('M','Y'), c('A','M'), c('A','L'), c('M','L'), c('L','Y')))
#' f.reachable_closure(graph, 'A')
#'
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

#' Check if a graph is nonparametrically saturated.
#'
#' Function to check if a graph is nonparametrically saturated.
#' A graph being nonparametrically saturated means that the graph implies NO equality constraints on the observed data distribution
#' @param graph A graph object generated by the `make.graph()` function.
#' @keywords graph ADMG nonparametrically saturated
#' @return A logical value indicating whether the graph is nonparametrically saturated.
#' @export
#' @importFrom utils combn
#' @examples
#' graph <- make.graph(vertices=c('A','M','L','Y','X'),
#' bi_edges=list(c('A','Y')),
#' di_edges=list(c('X','A'), c('X','M'), c('X','L'),
#' c('X','Y'), c('M','Y'), c('A','M'), c('A','L'), c('M','L'), c('L','Y')))
#' is.np.saturated(graph)
#'
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

    # (cnt.districts(f.reachable_closure(graph,c(Vi,Vj))[[3]])$n.districts > 1)
    # !(V1 %in% f.district(f.reachable_closure(graph,c(V1,V2))[[3]], V2))

    if (!(V1 %in% f.parents(graph, f.district(f.reachable_closure(graph,V2)[[3]], V2))) &&  (cnt.districts(f.reachable_closure(graph,c(Vi,Vj))[[3]])$n.districts > 1)){

      message("The graph is not nonparametrically saturated.")
      # print(c(V1,V2))
      return(FALSE)

    } # end of if statement

  } # end of for loop

  message("The graph is nonparametrically saturated.")

  return(TRUE)
}





#######################################################################################
# whether the graph is mb-shield, which is the graph contains only ordinary constrains.
########################################################################################

#' Check if a graph is mb-shielded.
#'
#' Function to check if a graph is mb-shielded.
#' A graph being mb-shielded means that the graph only implies ordinary equality constraints on the observed data distribution.
#' @param graph A graph object generated by the `make.graph()` function.
#' @keywords graph ADMG mb-shielded
#' @return A logical value indicating whether the graph is mb-shielded.
#' @export
#' @examples
#' graph <- make.graph(vertices=c('A','M','L','Y','X'),
#' bi_edges=list(c('A','Y')),
#' di_edges=list(c('X','A'), c('X','M'), c('X','L'),
#' c('X','Y'), c('M','Y'), c('A','M'), c('A','L'), c('M','L'), c('L','Y')))
#' is.mb.shielded(graph)
#' @importFrom utils combn
#'
is.mb.shielded <- function(graph) {

  # Iterate over all pairs of vertices
  for (pair in combn(graph$vertices, 2, simplify = F)) {

    Vi <- pair[1]
    Vj <- pair[2]

    if (!( list(c(Vi,Vj)) %in% graph$bi_edges || list(c(Vj,Vi)) %in% graph$bi_edges || list(c(Vi,Vj)) %in% graph$di_edges || list(c(Vj,Vi)) %in% graph$di_edges) ) { # edge between Vi and Vj is absent in G

      if (Vi %in% f.markov_blanket(graph, Vj) || Vj %in% f.markov_blanket(graph, Vi)) { # Vi is in the Markov blanket of Vj or Vj is in the Markov blanket of Vi

        message("The graph is not mb-shield.")

        return(FALSE)

      } # end of if statement

    } # end of if statement

  } # end of for loop

  message("The graph is mb-shielded.")
  return(TRUE)
}


#######################################################################################
# Function to get sets C, M, L
########################################################################################

#' A function for getting the sets C, M, L from an ADMG given a treatment
#'
#' Set C is the pre-treatment variable
#' Set L is the post-treatment variables that are within the district of the treatment, including the treatment itself.
#' Set M is the rest of the variables
#' @param graph A graph object generated via function `make.graph()`.
#' @param treatment A character string indicating the treatment variable.
#' @keywords CML
#' @return A list of sets C, M, L
#' @export
#' @examples
#' graph <- make.graph(vertices=c('A','M','L','Y','X'),
#' bi_edges=list(c('A','Y')),
#' di_edges=list(c('X','A'), c('X','M'), c('X','L'),
#' c('X','Y'), c('M','Y'), c('A','M'), c('A','L'), c('M','L'), c('L','Y')))
#' CML(graph, treatment='A')
#'
CML <- function(graph, treatment){

  vertices <- graph$vertices

  # return topological ordering, rank treatment variable as far back as possible
  tau <- f.top_order(graph, treatment)

  # Get set C, M, L
  C <- tau[1:(which(tau==treatment)-1)] # everything comes before the treatment following topological order tau

  L <- intersect(f.district(graph, treatment),
                 tau[(which(tau==treatment)):length(tau)]) # variables within the district of treatment and comes after the treatment (including the treatment itself) following topological order tau

  M <- setdiff(vertices,c(C,L)) # everything else

  return(list(C=C, M=M, L=L))
}


#######################################################################################
# Function for re-rank a target vector based a reference vector
########################################################################################

#' Rerank the target vector according to the variable order in reference vector
#'
#' In our case, we re-rank the vector C, M, L according to topological order tau
#' @param target A vector to be re-ranked
#' @param reference A vector that provides the order for re-ranking
#' @keywords rerank ranking
#' @return A re-ranked vector
#' @export
#' @examples
#' rerank(target=c('A','C','B'), reference=c('A','B','C','D'))
#'
rerank <- function(target,reference){

  # Find the indices of elements in M within V
  indices <- match(target, reference)

  # Reorder M according to the indices
  re_ranked_target <- target[order(indices)]

  return(re_ranked_target)
}
