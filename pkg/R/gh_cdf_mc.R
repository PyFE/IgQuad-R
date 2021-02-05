#' CDF of GH distribution at x from simulation
#'
#' @param mu GH parameter in the paper
#' @param alpha GH parameter in the paper
#' @param beta GH parameter in the paper
#' @param delta GH parameter in the paper
#' @param p GH parameter in the paper
#' @param param (mu, delta, alpha, beta, p) instead of above
#' @param n.quad # of quadrature points
#' @param correct normalize weights if TRUE
#' @param GIGrvg If TRUE, use GIGrvg::rgig instead
#' @param antithetic
#'
gh.cdf.mc <- function(x, mu=0, alpha=1, beta=0, delta=1, p=-0.5,
                      param=c(mu, delta, alpha, beta, p),
                      n.quad=50, n=1e6, correct=T, GIGrvg=F,
                      antithetic=F) {
  y <- ghyp.rand(param=param, n=n, p=p, n.quad=n.quad, correct=correct,
                 GIGrvg=GIGrvg, antithetic=antithetic)
  p = colSums( outer(y, x, '<') )/n
  return( p )
}
