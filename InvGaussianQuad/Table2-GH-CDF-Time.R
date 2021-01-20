library(GeneralizedHyperbolic)
source("igquad.R")

mu <- c(0, 0.00029, 0.000666, 0.000048)
alpha <- c(1, 138.78464, 214.4, 9.0)
beta <- c(0, -4.90461, -6.17, 2.73)
delta <- c(1, 0.00646, 0.0022, 0.0161)
p <- c(-0.5, -0.5, 0.8357, -1.6630)

params <- cbind(mu, delta, alpha, beta, p)

#### Check computation time with similar accuracy
####
param <- params[2,]  # <-- Change the set #
n.grid = 100
q = (1:n.grid-1)/n.grid
n.run = 1000

x <- GeneralizedHyperbolic::qghyp(q, param=param, intTol=1e-8)
system.time( p <- GeneralizedHyperbolic::pghyp(x, param=param, intTol=1e-14, subdivisions=500) )
system.time( for(k in 1:n.run) p_int <- GeneralizedHyperbolic::pghyp(x, param=param, intTol=2e-3) )
system.time( for(k in 1:n.run) p_quad <- ghyp.cdf(x, param=param, n.quad=50) )

# Error
print( c(max(abs(p_int - p)), max(abs(p_quad - p))) )

