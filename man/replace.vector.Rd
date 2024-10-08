% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graph-functions.R
\name{replace.vector}
\alias{replace.vector}
\title{Replace elements at a vector with other vectors.}
\usage{
replace.vector(vec, multivariate.variables = NULL)
}
\arguments{
\item{vec}{A vector.}

\item{multivariate.variables}{A list of vectors to replace the elements in vec.}
}
\value{
A vector with elements replaced by the vectors in multivariate.variables.
}
\description{
This function replace the elements in vec with the vectors in multivariate.variables under the same name.
}
\examples{
replace.vector(c('A','B','C'), list(A=c('D','E')))
}
\keyword{element}
\keyword{replace}
\keyword{vector}
