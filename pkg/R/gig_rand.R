#' GIG random variate
#'
#' @param n number of RVs
#' @param gamma GIG parameter in the paper
#' @param delta GIG parameter in the paper
#' @param p GIG parameter in the paper
#' @param n.quad number of quadrature points
#' @param correct normalize weights if TRUE
#' @param antithetic
#' @import stats
#'
gig.rand <- function(n=99, gamma=1, delta=1, p=-0.5, n.quad=50, correct=T, antithetic=F) {

  gig.quad <- gig.quad(gamma=gamma, delta=delta, p=p, n.quad=n.quad, correct=correct)
  x_gig <- gig.quad$nodes
  w_gig <- c(-999, cumsum(gig.quad$weights))
  w_gig[length(w_gig)] <- 999
  if(antithetic) {
    u <- runif(n/2)
    u <- c(u, 1.0-u)
  } else {
    u <- runif(n)
  }
  x <- x_gig[.bincode(u, w_gig)]
  return( x )
}
