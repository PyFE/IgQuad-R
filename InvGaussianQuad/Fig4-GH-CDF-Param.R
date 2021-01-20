source("igquad.R")

### Error for varying |Beta|
mu <- 0
n.quad.arr = c(60, 80, 100)

beta.arr <- seq(0,5, by=0.3)
param.arr <- beta.arr

sigma = 1
alpha = sqrt(sigma^2 + beta.arr^2)
p = -0.5

#### Exact GH CDF is from `GeneralizedHyperbolic` package
p_err = matrix(0, nrow=length(param.arr), ncol=length(n.quad.arr))

for(n in 1:length(param.arr)){
  param <- c(mu, sigma, alpha[n], beta.arr[n], p)
  
  # , method="integrate"
  x <- GeneralizedHyperbolic::qghyp((1:99)/100, param=param, intTol=1e-8)
  p_int <- GeneralizedHyperbolic::pghyp(x, param=param, intTol=1e-14, subdivisions=500)
  
  for (k in 1:length(n.quad.arr)) {
    p_quad <- ghyp.cdf(x, param=param, n.quad=n.quad.arr[k], correct=T)
    p_err[n,k] <- max(abs(p_quad - p_int))
  }
}

### Plot
par(lwd=2, mar=c(3,3,0.5,0.5)+0.1, mgp=c(2,0.6,0), cex=1.5) 
plot(NA, xlim=c(min(beta.arr),max(beta.arr)), ylim=c(-15.4,-2.6), 
     ylab=latex2exp::TeX("Log_{10}(|error|)"), xlab=latex2exp::TeX("$|\\tilde{\\beta}|$"))
grid(lty=6)

points(param.arr, log10(p_err[,1]), col="blue", pch=19, cex=1) # circle
points(param.arr, log10(p_err[,2]), col="red", pch=17, cex=1) # triangle
points(param.arr, log10(p_err[,3]), col="darkgreen", pch=3, cex=1) # cross

legend(4.8, -11.5, c("n= 60", "n= 80", "n=100"), cex=0.6, ncol=1, xjust=1, y.intersp=0.75,x.intersp=0.75,
       pch=c(19, 17, 3), col=c("blue", "red", "darkgreen"))


### Error for varying Sigma
mu <- 0
n.quad.arr = c(60, 80, 100)

beta = 0.0

sigma.arr <- seq(0.05, 1, by=0.05)
param.arr <- sigma.arr

alpha = sqrt(sigma.arr^2 + beta^2)
p = -0.5

#### Exact GH CDF is from `GeneralizedHyperbolic` package
p_err = matrix(0, nrow=length(param.arr), ncol=length(n.quad.arr))

for(n in 1:length(param.arr)){
  param <- c(mu, sigma.arr[n], alpha[n], beta, p)
  
  # , method="integrate"
  x <- GeneralizedHyperbolic::qghyp((1:99)/100, param=param, intTol=1e-8)
  p_int <- GeneralizedHyperbolic::pghyp(x, param=param, intTol=1e-14, subdivisions=500)
  
  for (k in 1:length(n.quad.arr)) {
    p_quad <- ghyp.cdf(x, param=param, n.quad=n.quad.arr[k], correct=T)
    p_err[n,k] <- max(abs(p_quad - p_int))
  }
}

### Plot
par(lwd=2, mar=c(3,3,0.5,0.5)+0.1, mgp=c(2,0.6,0), cex=1.5) 
plot(NA, xlim=c(min(param.arr),max(param.arr)), ylim=c(-12.5,-1.5), 
     ylab=latex2exp::TeX("Log_{10}(|error|)"), xlab=latex2exp::TeX("$\\sigma$"))
grid(lty=6)

points(param.arr, log10(p_err[,1]), col="blue", pch=19, cex=1) # circle
points(param.arr, log10(p_err[,2]), col="red", pch=17, cex=1) # triangle
points(param.arr, log10(p_err[,3]), col="darkgreen", pch=3, cex=1) # cross
legend(1.0, -1.5, c("n= 60", "n= 80", "n=100"), cex=0.6, ncol=1, xjust=1, y.intersp=0.75,x.intersp=0.75,
       pch=c(19, 17, 3), col=c("blue", "red", "darkgreen"))


### Error for varying p
mu <- 0
n.quad.arr = c(60, 80, 100)

beta = 0.0
sigma <- 1
p.arr <- seq(-15, 15, by=1)
param.arr <- p.arr
alpha = sqrt(sigma^2 + beta^2)

#### Exact GH CDF is from `GeneralizedHyperbolic` package
p_err = matrix(0, nrow=length(param.arr), ncol=length(n.quad.arr))

for(n in 1:length(param.arr)){
  param <- c(mu, sigma, alpha, beta, p.arr[n])
  
  # , method="integrate"
  x <- GeneralizedHyperbolic::qghyp((1:99)/100, param=param, intTol=1e-8)
  p_int <- GeneralizedHyperbolic::pghyp(x, param=param, intTol=1e-14, subdivisions=500)
  
  for (k in 1:length(n.quad.arr)) {
    p_quad <- ghyp.cdf(x, param=param, n.quad=n.quad.arr[k], correct=T)
    p_err[n,k] <- max(abs(p_quad - p_int))
  }
}

### Plot
par(lwd=2, mar=c(3,3,0.5,0.5)+0.1, mgp=c(2,0.6,0), cex=1.5) 
plot(NA, xlim=c(min(param.arr),max(param.arr)), ylim=c(-16.5,-9.5), 
     ylab=latex2exp::TeX("Log_{10}(|error|)"), xlab=latex2exp::TeX("p"))
grid(lty=6)

points(param.arr, log10(p_err[,1]), col="blue", pch=19, cex=1) # circle
points(param.arr, log10(p_err[,2]), col="red", pch=17, cex=1) # triangle
points(param.arr, log10(p_err[,3]), col="darkgreen", pch=3, cex=1) # cross
legend(-7, -9.5, c("n= 60", "n= 80", "n=100"), cex=0.6, ncol=1, xjust=1, y.intersp=0.75,x.intersp=0.75,
       pch=c(19, 17, 3), col=c("blue", "red", "darkgreen"))