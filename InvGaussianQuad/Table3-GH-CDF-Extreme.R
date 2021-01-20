source("igquad.R")

#### Compute CDF on fixed points 
####

#q_arr = c(0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99)
lower.tail = F
q_arr = c(1e-9, 1e-6, 0.001)
if(!lower.tail) {q_arr = 1-rev(q_arr)}
n_q_arr = length(q_arr)
p_err_quad <- matrix(0, ncol=n_q_arr, nrow=nrow(params))

for(n in 1:nrow(params)){
  param <- params[n,]
  
  # , intTol=1e-12, subdivisions=500, method='integrate'
  x <- GeneralizedHyperbolic::qghyp(q_arr, param=param, intTol=1e-14, method='integrate')
  print(x)
  p_int <- GeneralizedHyperbolic::pghyp(x, param=param, 
                                        intTol=1e-14, subdivisions=500, lower.tail=lower.tail)
  p_quad <- ghyp.cdf(x, param=param, n.quad=50, lower.tail=lower.tail)
  p_err_quad[n,] <- (p_quad - p_int)*ifelse(lower.tail,1,-1)
}

write.table(t(p_err_quad), "clipboard", sep="\t", row.names=FALSE, col.names=FALSE)
