#' This is data generated from the simulation model in Figure 4(a) of the paper "Targeted Learning of Average Causal Effects in Graphical Models with Unmeasured Variables: Extensions of the Front-Door Model".
#' The model can also be found in Figure (a) of the github repository: \url{https://github.com/annaguo-bios/ADMGtmle/tree/main}
#'
#'
#' @format A data frame with 2000 rows and 6 variables: X, A, M, L, Y.
#' \describe{
#'  \item{X}{a variable that follows uniform distribution at the interval of \eqn{[0,1]}}
#'  \item{A}{a binary treatment}
#'  \item{M}{a variable that contains M.1 and M.2, generated from a bivariate normal distribution}
#'  \item{L}{a variable that follows normal distribution}
#'  \item{Y}{a variable that follows normal distribution}
#'  }

#' @name data_fig_4a
#' @docType data
#' @references \url{https://github.com/annaguo-bios/ADMGtmle/tree/main}
#' @keywords data
"data_fig_4a"


#' This is data generated from the simulation model in Figure 4(b) of the paper "Targeted Learning of Average Causal Effects in Graphical Models with Unmeasured Variables: Extensions of the Front-Door Model".
#' The model can also be found in Figure (b) of the github repository: \url{https://github.com/annaguo-bios/ADMGtmle/tree/main}
#'
#'
#' @format A data frame with 2000 rows and 6 variables: X, A, M, L, Y.
#' \describe{
#'  \item{X}{a variable that follows uniform distribution at the interval of \eqn{[0,1]}}
#'  \item{A}{a binary treatment}
#'  \item{M}{a variable that contains M.1 and M.2, generated from a bivariate normal distribution}
#'  \item{L}{a variable that follows normal distribution}
#'  \item{Y}{a variable that follows normal distribution}
#'  }

#' @name data_fig_4b
#' @docType data
#' @references \url{https://github.com/annaguo-bios/ADMGtmle/tree/main}
#' @keywords data
"data_fig_4b"



#' This is data generated from a classic backdoor model, where the pre-treatment variable X cause A and Y, and the treatment A cause Y.
#'
#'
#' @format A data frame with 2000 rows and 3 variables: X, A, Y.
#' \describe{
#'  \item{X}{a variable that follows uniform distribution at the interval of \eqn{[0,1]}}
#'  \item{A}{a binary treatment}
#'  \item{Y}{a variable that follows normal distribution}
#'  }

#' @name data_backdoor
#' @docType data
#' @references \url{https://github.com/annaguo-bios/ADMGtmle/tree/main}
#' @keywords data
"data_backdoor"
