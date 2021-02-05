#' Inverse Gaussian Quadrature
#' Theorem 1 in the paper
#'
#' @param mu IG parameter
#' @param lambda IG parameter
#' @param n.quad # of quadrature points
#'
ig.quad <- function(mu=1, lambda=1, n.quad=8){

  quad <- statmod::gauss.quad.prob(n.quad, "normal")
  z <- quad$nodes
  w <- quad$weights

  fac <- 0.5*mu/lambda

  y_hat <- fac*z*z

  x_p <- 1 + y_hat + sqrt(2*fac)*z*sqrt(1 + 0.5*y_hat)
  x_ig <- x_p*mu
  w_ig <- 2*w/(1+x_p)

  return( list(nodes=x_ig, weights=w_ig) )
}
