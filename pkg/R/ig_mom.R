#' IG moment
#'
#' @param order moment order
#' @param mu IG parameter
#' @param lambda IG parameter
#'
ig.mom <- function(order=2, mu=1, lambda=1){
  delta <- sqrt(lambda)
  gamma <- sqrt(lambda)/mu
  return( gig.mom(order, gamma, delta, -0.5) )
}
