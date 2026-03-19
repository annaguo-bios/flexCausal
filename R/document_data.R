#' Data generated from the simulation model in Figure 4(a) of the paper \url{https://arxiv.org/abs/2409.03962}.
#' The model can also be found in Figure (a) of the github repository: \url{https://github.com/annaguo-bios/flexCausal}
#'
#'
#' @format A data frame with 2000 rows and 7 variables: X, A, U, M=(M.1,M.2), L, Y.
#' \describe{
#'  \item{X}{a variable that follows uniform distribution at the interval of \eqn{[0,1]}}
#'  \item{A}{a binary treatment}
#'  \item{U}{an unmeasured confounder following a normal distribution}
#'  \item{M.1}{first component of the bivariate mediator M, following a normal distribution}
#'  \item{M.2}{second component of the bivariate mediator M, following a normal distribution}
#'  \item{L}{a variable that follows normal distribution}
#'  \item{Y}{a variable that follows normal distribution}
#'  }
#' @name data_example_a
#' @docType data
#' @references \url{https://github.com/annaguo-bios/flexCausal}
#' @keywords data
"data_example_a"


#' Data generated from the simulation model in Figure 4(b) of the paper \url{https://arxiv.org/abs/2409.03962}.
#' The model can also be found in Figure (b) of the github repository: \url{https://github.com/annaguo-bios/flexCausal}
#'
#'
#' @format A data frame with 2000 rows and 6 variables: X, A, U1, U2, M=(M.1,M.2), L, Y.
#' \describe{
#'  \item{X}{a variable that follows uniform distribution at the interval of \eqn{[0,1]}}
#'  \item{A}{a binary treatment}
#'  \item{U1}{first unmeasured confounder following a normal distribution}
#'  \item{U2}{second unmeasured confounder following a normal distribution}
#'  \item{M.1}{first component of the bivariate mediator M, following a normal distribution}
#'  \item{M.2}{second component of the bivariate mediator M, following a normal distribution}
#'  \item{L}{a variable that follows normal distribution}
#'  \item{Y}{a variable that follows normal distribution}
#'  }
#' @name data_example_b
#' @docType data
#' @references \url{https://github.com/annaguo-bios/flexCausal}
#' @keywords data
"data_example_b"



#' Data generated from a classic backdoor model, where the pre-treatment variable X cause A and Y, and the treatment A cause Y.
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
#' @references \url{https://github.com/annaguo-bios/flexCausal}
#' @keywords data
"data_backdoor"


#' Data generated from a classic front-door model.
#'
#' Data generated from a front-door model where a pre-treatment variable X and
#' an unmeasured confounder U affect both the treatment A and outcome Y, and
#' the treatment A affects Y through a binary mediator M.
#'
#' @format A data frame with 2000 rows and 5 variables:
#' \describe{
#'  \item{X}{a pre-treatment variable following a uniform distribution on \eqn{[0, 1]}}
#'  \item{U}{an unmeasured confounder following a normal distribution}
#'  \item{A}{a binary treatment variable (0/1)}
#'  \item{M}{a binary mediator variable (0/1)}
#'  \item{Y}{a continuous outcome variable following a normal distribution}
#' }
#'
#' @name data_frontdoor
#' @docType data
#' @references \url{https://github.com/annaguo-bios/flexCausal}
#' @keywords data
"data_frontdoor"
