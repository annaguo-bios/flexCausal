% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graph-functions.R
\name{treatment.queue}
\alias{treatment.queue}
\title{Put treatment to the end of the queue.}
\usage{
treatment.queue(adj.matrix, queue, treatment)
}
\arguments{
\item{adj.matrix}{An adjacency matrix of the graph.}

\item{queue}{A vector of integers that records the current queue of vertices.}

\item{treatment}{A character string indicating the treatment variable.}
}
\value{
A updated queue with the treatment variable at the end if it is in the queue.
}
\description{
This function is called by the \code{f.top_orderMAT()} function to put the treatment to the end of the queue when employing Kahn's algorithm for topological ordering.
In deriving the topological ordering, there can be multiple vertices with the same in-degree. This function ensures that the treatment is ranked as far back as possible when there are ties.
}
\keyword{ADMG}
\keyword{graph}
