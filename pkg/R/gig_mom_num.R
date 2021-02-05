#' Numerical moment of GIG distribution.
#'
#' @param order moment order
#' @param gamma GIG parameter in the paper
#' @param delta GIG parameter in the paper
#' @param p GIG parameter in the paper
#' @param n.quad # of quadrature points
#' @param correct normalize weights if TRUE
#'
gig.mom.num <- function(order=2, gamma=1, delta=1, p=-0.5, n.quad=8, correct=T){
  quad <- gig.quad(gamma=gamma, delta=delta, p=p, n.quad=n.quad, correct=correct)
  x.gig <- quad$nodes
  w.gig <- quad$weights
  x.gig.power <- outer(x.gig, order, "^")

  return( colSums( x.gig.power*w.gig ) )
}
