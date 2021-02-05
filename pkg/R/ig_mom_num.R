#' IG moment numerically evaluated with IG quadrature
#'
#' @param order moment order
#' @param mu IG parameter
#' @param lambda IG parameter
#' @param n.quad number of quadrature points
#'
ig.mom.num <- function(order=2, mu=1, lambda=1, n.quad=8){
  quad <- ig.quad(mu, lambda, n.quad)
  x.ig <- quad$nodes
  w.ig <- quad$weights
  x.ig.power <- outer(x.ig, order, "^")

  return( colSums( x.ig.power*w.ig ) )
}
