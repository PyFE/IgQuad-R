library(GeneralizedHyperbolic)
source("igquad.R")

#### Compute CDF using MC
####

n.run = 1000
q_arr = c(0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99)
n_q_arr = length(q_arr)
p_err_mc <- matrix(0, ncol=n_q_arr, nrow=nrow(params))
p_err_sd <- matrix(0, ncol=n_q_arr, nrow=nrow(params))

for(n in 1:nrow(params)){
  param <- params[n,]
  
  # , intTol=1e-12, subdivisions=500, method='integrate'
  x <- GeneralizedHyperbolic::qghyp(q_arr, param=param, intTol=1e-8)
  print(x)
  p_int <- GeneralizedHyperbolic::pghyp(x, param=param, intTol=1e-14, subdivisions=500)
  
  set.seed(755)
  p_mc = matrix(0, nrow=n.run, ncol=length(x))
  for(k in 1:n.run) { 
    p_mc[k,] <- ghyp.cdf.mc(x, param=param, n.quad=50, n=1e6, GIGrvg=T,
                            antithetic=F) 
  }
  ##for(k in 1:100) { p_mc[k,] <- ghyp.cdf.mchyb(x, param=param, n.quad=25, n=1e5) }
  p_err_mc[n,] <- colMeans(p_mc) - p_int
  p_err_sd[n,] <- apply(p_mc, 2, sd)
}

write.table(t(p_err_mc), "clipboard", sep="\t", row.names=FALSE, col.names=FALSE)
write.table(t(p_err_sd), "clipboard", sep="\t", row.names=FALSE, col.names=FALSE)
