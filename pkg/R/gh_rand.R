#' GH random variate
#' Eq. (5)
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
gh.rand <- function(n=99, mu=0, alpha=1, beta=0, delta=1, p=-0.5,
                    param=c(mu, delta, alpha, beta, p),
                    n.quad=50, correct=T, GIGrvg=F, antithetic=T) {
  mu <- param[1]
  delta <- param[2]
  gamma <- sqrt(param[3]^2 - param[4]^2)
  beta <- param[4]
  p <- param[5]

  n_nig <- ifelse(antithetic, n/2, n)

  if(GIGrvg){
    x <- GIGrvg::rgig(n_nig, lambda=p, psi=gamma^2, chi=delta^2)
  } else {
    x <- gig.rand(n_nig, gamma=gamma, delta=delta, p=p,
                  n.quad=n.quad, correct=correct, antithetic=F)
  }

  if(antithetic){
    z <- rnorm(n/2)
    y <- mu + beta*x + sqrt(x)*c(z, -z)
  } else {
    z <- rnorm(n)
    y <- mu + beta*x + sqrt(x)*z
  }

  return( y )
}
