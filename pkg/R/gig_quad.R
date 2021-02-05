#' GIG quadrature
#' Corollary 1 in the paper
#'
#' @param gamma GIG parameter in the paper
#' @param delta GIG parameter in the paper
#' @param p GIG parameter in the paper
#' @param n.quad number of quadrature points
#' @param correct normalize weights if TRUE
#'
gig.quad <- function(gamma=1, delta=1, p=-0.5, n.quad=8, correct=T) {
  scale <- delta/gamma
  ig.quad <- ig.quad(mu=1, lambda=gamma*delta, n.quad=n.quad)
  x_ig <- scale * ig.quad$nodes
  w_ig <- ig.quad$weights

  if (p == -0.5) {
    ratio <- 1
  } else {
    ratio =  sqrt(pi/2) / scale^p/delta * exp(-gamma*delta)/besselK(gamma*delta, p)
  }

  w_gig = w_ig * x_ig^(p+0.5) * ratio

  if(correct) {
    w_gig = w_gig / sum(w_gig)
  }

  return( list(nodes=x_ig, weights=w_gig) )
}
