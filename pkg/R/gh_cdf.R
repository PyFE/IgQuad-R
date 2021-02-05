#' CDF of GH distribution
#' Eq. (3)
#'
#' @param mu GH parameter in the paper
#' @param alpha GH parameter in the paper
#' @param beta GH parameter in the paper
#' @param delta GH parameter in the paper
#' @param p GH parameter in the paper
#' @param param (mu, delta, alpha, beta, p) instead of above
#' @param n.quad number of quadrature points
#' @param correct normalize weights if TRUE
#' @param lower.tail
#' @import stats
#'
gh.cdf <- function(x, mu=0, alpha=1, beta=0, delta=1, p=-0.5,
                   param=c(mu, delta, alpha, beta, p),
                   n.quad=8, correct=T, lower.tail=T) {

  mu <- param[1]
  delta <- param[2]
  gamma <- sqrt(param[3]^2 - param[4]^2)
  beta <- param[4]
  p <- param[5]

  gig.quad <- gig.quad(delta=delta, gamma=gamma, p=p, n.quad=n.quad, correct=correct)
  x_gig <- gig.quad$nodes
  w_gig <- gig.quad$weights

  p = colSums(pnorm( outer(-mu-beta*x_gig, x, "+")/sqrt(x_gig), lower.tail=lower.tail) * w_gig)

  return( p )
}
