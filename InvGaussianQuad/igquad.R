library(statmod)
library(GIGrvg)

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

#' GIG quadrature
#' Corollary 1 in the paper
#' 
#' @param gamma GIG parameter in the paper
#' @param delta GIG parameter in the paper
#' @param p GIG parameter in the paper
#' @param n.quad # of quadrature points
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

#' IG moment numerically evaluated with IG quadrature
#'
#' @param order moment order
#' @param mu IG parameter
#' @param lambda IG parameter
#' @param n.quad # of quadrature points
#' 
ig.mom.num <- function(order=2, mu=1, lambda=1, n.quad=8){
  quad <- ig.quad(mu, lambda, n.quad)
  x.ig <- quad$nodes
  w.ig <- quad$weights
  x.ig.power <- outer(x.ig, order, "^")
  
  return( colSums( x.ig.power*w.ig ) )
}

#' IG moment
#'
#' @param order moment order
#' @param mu IG parameter
#' @param lambda IG parameter
#' 
ig.mom <- function(order=2, mu=1, lambda=1){
  delta <- sqrt(lambda)
  gamma <- sqrt(lambda)/mu
  return( gig.mom(order, gamma, delta, -0.5) )
}

#' Analytic moments of GIG distribution.
#' Eq. (9)
#' 
#' @param order moment order
#' @param gamma GIG parameter in the paper
#' @param delta GIG parameter in the paper
#' @param p GIG parameter in the paper
#' 
gig.mom <- function(order=2, gamma=1, delta=1, p=-0.5){
  return( (delta/gamma)^order * besselK(gamma*delta, p+order)/besselK(gamma*delta, p) )
}

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

#' GIG random variate
#' 
#' @param n # of RVs
#' @param gamma GIG parameter in the paper
#' @param delta GIG parameter in the paper
#' @param p GIG parameter in the paper
#' @param n.quad # of quadrature points
#' @param correct normalize weights if TRUE
#' @param antithetic
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

#' Moment generating function of GIG distribution
#' Eq. (10) 
#'
#' @param t dummy variable value 
#' @param gamma GIG parameter in the paper
#' @param delta GIG parameter in the paper
#' @param p GIG parameter in the paper
#'
gig.mgf <- function(t, gamma=1, delta=1, p=-0.5) {
  tmp = sqrt(gamma^2 - 2*t)
  mgf = (gamma/tmp)^p * besselK(delta*tmp, p) / besselK(gamma*delta, p)
  return(mgf)
}

#' Moment generating function of GIG distribution
#' Numerically evaluated at t
#' Eq. (11) 
#' 
#' @param t dummy variable value 
#' @param gamma GIG parameter in the paper
#' @param delta GIG parameter in the paper
#' @param p GIG parameter in the paper
#' @param n.quad # of quadrature points
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
ghyp.quad <- function(x, mu=0, alpha=1, beta=0, delta=1, p=-0.5,
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

#' CDF of GH distribution
#' Eq. (3)
#' 
#' @param mu GH parameter in the paper
#' @param alpha GH parameter in the paper
#' @param beta GH parameter in the paper
#' @param delta GH parameter in the paper
#' @param p GH parameter in the paper
#' @param param (mu, delta, alpha, beta, p) instead of above
#' @param n.quad # of quadrature points
#' @param correct normalize weights if TRUE
#' @param lower.tail
#' 
ghyp.cdf <- function(x, mu=0, alpha=1, beta=0, delta=1, p=-0.5, 
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
ghyp.rand <- function(n=99, mu=0, alpha=1, beta=0, delta=1, p=-0.5, 
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
ghyp.cdf.mc <- function(x, mu=0, alpha=1, beta=0, delta=1, p=-0.5, 
                        param=c(mu, delta, alpha, beta, p),
                        n.quad=50, n=1e6, correct=T, GIGrvg=F,
                        antithetic=F) {
  y <- ghyp.rand(param=param, n=n, p=p, n.quad=n.quad, correct=correct, 
                 GIGrvg=GIGrvg, antithetic=antithetic) 
  p = colSums( outer(y, x, '<') )/n
  return( p )
}