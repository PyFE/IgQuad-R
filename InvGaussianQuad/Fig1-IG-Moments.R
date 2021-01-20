source('igquad.R')
### Check momnents of IG against analytic 

n.quad = 10
sigma <- 1
order <- -n.quad:(n.quad+1)
mom.num = gig.mom.num(order=order, gamma=sigma, delta=sigma, n.quad=n.quad)
mom.ana = gig.mom(order=order, gamma=sigma, delta=sigma)
abs(mom.num/mom.ana - 1) < 1e-13

n.quad=10
order1 <- seq(0,n.quad+2,by=0.01)-0.00001
order2 <- seq(0,n.quad+2,by=1)+0.00001
order <- sort(c(order1, order2))
mom.num = gig.mom.num(order=order, gamma=sigma, delta=sigma, n.quad=n.quad)
mom.ana = gig.mom(order=order, gamma=sigma, delta=sigma)
abs(mom.num/mom.ana - 1) < 1e-13

# plot
par(lwd=2, mar=c(3,3,0.5,0.5)+0.1, mgp=c(2,0.6,0), cex=1.5) 
plot(NA, xlim=c(min(order),max(order)), ylim=c(-7,-1.9), 
     #plot(NA, xlim=c(min(order),max(order)), ylim=c(-10.5,-5.6), 
     ylab=latex2exp::TeX("Log_{10}(|relative error|)"), xlab="r")
#title("...")
grid(lty=6)

ind_pos = (mom.num/mom.ana>1)
log_err = log10(abs(mom.num/mom.ana-1))
log_err_pos = log_err
log_err_pos[!ind_pos] <- NA
log_err_neg = log_err
log_err_neg[ind_pos] <- NA

lines(order, log_err_pos, col="blue" )
lines(order, log_err_neg, col="red", lty="dashed" )
