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
