#' Moment generating function of GIG distribution
#' Numerically evaluated at t
#' Eq. (11)
#'
#' @param t dummy variable value
#' @param gamma GIG parameter in the paper
#' @param delta GIG parameter in the paper
#' @param p GIG parameter in the paper
#' @param n.quad number of quadrature points
#' @param correct normalize weights if TRUE
#'
gig.mgf.num <- function(t, gamma=1, delta=1, p=-0.5, n.quad=50, correct=F) {

  gig.quad <- gig.quad(delta=delta, gamma=gamma, p=-abs(p), n.quad=n.quad, correct=correct)
  x_gig <- gig.quad$nodes
  w_gig <- gig.quad$weights

  mgf = colSums(exp(x_gig %o% t)*w_gig)
  tmp = (1 - 2*t/gamma^2)^min(-p,0)

  return( tmp*mgf )
}
