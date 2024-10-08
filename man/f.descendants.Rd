% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graph-functions.R
\name{f.descendants}
\alias{f.descendants}
\title{Get the descendants of a node OR nodes in a graph.}
\usage{
f.descendants(graph, nodes)
}
\arguments{
\item{graph}{A graph object generated by the \code{make.graph()} function.}

\item{nodes}{A character vector of nodes for which to extract children.}
}
\value{
A vector of vertices contains descendants set of the given nodes.
}
\description{
Function to extract the descendants of a node OR nodes in a graph object.
Descendants of a node Vi are set Vj such that there is a directed path Vi->...->Vj. Descendants set including Vi itself by convention.
}
\examples{
graph <- make.graph(vertices=c('A','M','L','Y','X'),
bi_edges=list(c('A','Y')),
di_edges=list(c('X','A'), c('X','M'), c('X','L'),
c('X','Y'), c('M','Y'), c('A','M'), c('A','L'), c('M','L'), c('L','Y')))
f.descendants(graph, c('A'))

}
\keyword{ADMG}
\keyword{descendants}
\keyword{graph}
