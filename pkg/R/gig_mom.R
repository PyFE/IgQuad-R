#' Analytic moments of GIG distribution.
#' Eq. (9)
#'
#' @param order moment order
#' @param gamma GIG parameter in the paper
#' @param delta GIG parameter in the paper
#' @param p GIG parameter in the paper
#'
gig.mom <- function(order=2, gamma=1, delta=1, p=-0.5){
  return( (delta/gamma)^order * besselK(gamma*delta, p+order)/besselK(gamma*delta, p) )
}
