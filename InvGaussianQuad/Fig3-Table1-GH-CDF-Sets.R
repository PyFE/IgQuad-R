library(GeneralizedHyperbolic)
source("igquad.R")

### Simple test of quadrature
q1 <- gig.quad(n.quad=15, p=-10, delta=1)
q1$nodes*rev(q2$nodes)-1
sum(exp(0.5*q1$nodes)*q1$weights)

### Four Parameter Sets
# c(mu, delta, alpha, beta, lambda))
# Standard (NIG), EURUSD (NIG), NYSE Composite(PhD), BMW(GH), 
mu <- c(0, 0.00029, 0.000666, 0.000048)
alpha <- c(1, 138.78464, 214.4, 9.0)
beta <- c(0, -4.90461, -6.17, 2.73)
delta <- c(1, 0.00646, 0.0022, 0.0161)
p <- c(-0.5, -0.5, 0.8357, -1.6630)
params <- cbind(mu, delta, alpha, beta, p)

gamma <- sqrt(alpha^2-beta^2)
sigma <- sqrt(gamma*delta)

mom = matrix(0, ncol=4, nrow=length(beta))
mom_quad = matrix(0, ncol=4, nrow=length(beta))

for(n in 1:length(beta)){
  y_quad <- ghyp.quad(param=params[n,], n.quad=c(20,3))
  for(k in 1:4) {
    m.type <- if (k==1) "raw" else "central"
    mom[n,k] <- GeneralizedHyperbolic::ghypMom(k, param=params[n,], momType=m.type)
    mom_quad[n,k] = sum(y_quad$nodes^k * y_quad$weights)
  }
}
mom_quad/mom
write.table(t(params), "clipboard", sep="\t", row.names=FALSE, col.names=FALSE)
write.table(t(mom), "clipboard", sep="\t", row.names=FALSE, col.names=FALSE)


#### Compute CDF of NIG and compare against integration from `GeneralizedHyperbolic` package
####
n.quad.arr <- seq(20,100, by=5)
p_err = matrix(0, ncol=length(n.quad.arr), nrow=nrow(params))

for(n in 1:nrow(params)){
  param <- params[n,]
  # , method="integrate"
  x <- GeneralizedHyperbolic::qghyp((1:99)/100, param=param, intTol=1e-8)
  p_int <- GeneralizedHyperbolic::pghyp(x, param=param, intTol=1e-14, subdivisions=500)

    for (k in 1:length(n.quad.arr)) {
      p_quad <- ghyp.cdf(x, param=param, n.quad=n.quad.arr[k], correct=T)
      p_err[n,k] <- max(abs(p_quad - p_int))
      #p_err[n,k] <- sqrt(mean((p_quad - p_int)^2))
    }
}

### Plot
par(lwd=2, mar=c(3,3,0.5,0.5)+0.1, mgp=c(2,0.6,0), cex=1.5) 
plot(NA, xlim=c(min(n.quad.arr),max(n.quad.arr)), ylim=c(-12.5,-3.5), 
     ylab=latex2exp::TeX("Log_{10}(|error|)"), xlab=latex2exp::TeX("n (quadrature size)"))
grid(lty=6)

points(n.quad.arr, log10(p_err[1,]), col="blue", pch=19, cex=1) # circle
points(n.quad.arr, log10(p_err[2,]), col="red", pch=17, cex=1) # triangle
points(n.quad.arr, log10(p_err[3,]), col="darkgreen", pch=3, cex=1) # cross
points(n.quad.arr, log10(p_err[4,]), col="black", pch=4, cex=1) # x
legend(103, -3.2, c("Set 1", "Set 2", "Set 3", "Set 4"), cex=0.8, ncol=2, xjust=1, y.intersp=0.75,x.intersp=0.75,
       pch=c(19, 17, 3, 4), col=c("blue", "red", "darkgreen", "black"))
