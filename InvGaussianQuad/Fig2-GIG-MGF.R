library(latex2exp)
source("igquad.R")

#### Compute MGF and compare against analytic value for varing gamma=delta values

p <- -0.5 # -0.5, 1 or 90
n.quad.arr <- seq(21,101, by=5)  # this is to avoid a kink at n=40
sig.arr = c(0.5, 0.75, 1, 2)
p_err = matrix(0, ncol=length(n.quad.arr), nrow=length(sig.arr))

for (n in 1:length(sig.arr)){
  sigma <- sig.arr[n]
  for (k in 1:length(n.quad.arr)) {
    t=0.8 * 0.5*sigma^2
    mgf = gig.mgf(t, gamma=sigma, delta=sigma, p=p)
    mgf_num = gig.mgf.num(t, gamma=sigma, delta=sigma, p=p, n.quad=n.quad.arr[k], correct=T)
    p_err[n,k] <- mgf_num / mgf - 1
  }
}

### Plot
par(lwd=2, mar=c(3,3,0.5,0.5)+0.1, mgp=c(2,0.6,0), cex=1.5) 
plot(NA, xlim=c(min(n.quad.arr),max(n.quad.arr)), ylim=c(-16.5,-3.5), #ylim=c(-14.5,1), # 
     ylab=latex2exp::TeX("Log_{10}(|relative error|)"), xlab=latex2exp::TeX("n (quadrature size)"))
grid(lty=6)

points(n.quad.arr, log10(abs(p_err[1,])), col="blue", pch=19, cex=1) # circle
points(n.quad.arr, log10(abs(p_err[2,])), col="red", pch=17, cex=1) # triangle
points(n.quad.arr, log10(abs(p_err[3,])), col="darkgreen", pch=3, cex=1) # cross
points(n.quad.arr, log10(abs(p_err[4,])), col="black", pch=4, cex=1) # x
legend(110, -3.2+4, c(latex2exp::TeX("$\\sigma$=0.5"), latex2exp::TeX("$\\sigma$=0.75"), 
                    latex2exp::TeX("$\\sigma$=1"), latex2exp::TeX("$\\sigma$=2")), 
       cex=0.7, ncol=1, xjust=1, y.intersp=0.75, x.intersp=0.75,
       pch=c(19, 17, 3, 4), col=c("blue", "red", "darkgreen", "black"))