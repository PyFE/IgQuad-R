#' Double quadrature of the GH distribution
#'
#' @param mu GH parameter in the paper
#' @param alpha GH parameter in the paper
#' @param beta GH parameter in the paper
#' @param delta GH parameter in the paper
#' @param p GH parameter in the paper
#' @param param (mu, delta, alpha, beta, p) instead of above
#' @param n.quad # of quadrature points
#' @param correct normalize weights if TRUE
#'
gh.quad <- function(x, mu=0, alpha=1, beta=0, delta=1, p=-0.5,
                    param=c(mu, delta, alpha, beta, p),
                    n.quad=c(8,4), correct=T) {

  mu <- param[1]
  delta <- param[2]
  gamma <- sqrt(param[3]^2 - param[4]^2)
  beta <- param[4]
  p <- param[5]

  gig.quad <- gig.quad(delta=delta, gamma=gamma, p=p, n.quad=n.quad[1], correct=correct)
  x_gig <- gig.quad$nodes
  w_gig <- gig.quad$weights

  z.quad <- statmod::gauss.quad.prob(n.quad[2], "normal")
  z <- z.quad$nodes
  h <- z.quad$weights

  y_gh <- mu + beta*x_gig + sqrt(x_gig) %o% z
  w_gh <- w_gig %o% h

  return( list(nodes=y_gh, weights=w_gh) )
}
