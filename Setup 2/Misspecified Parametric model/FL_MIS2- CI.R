### BOOTSTRAP FUNCTIONS (INCREASE) ###
library(boot)
tau = log(3)
p_c1   = function(data, indices) {
  p.c     = ifelse( mean( data[indices] > tau, na.rm = T ) < 0.001, # prob of success in the control grp
                    0.001,
                    mean( data[indices] > tau, na.rm = T ))
  var.pc  = var( data[indices] )                                       # variance of sampled control values
  mean.pc = mean( data[indices] )                                      # mean of sampled control values
  return(c(p.c, var.pc, mean.pc))
}

p_t1 = function(data, indices) {
  p.t     = ifelse( mean( data[indices] > tau, na.rm = T ) > 0.999, # prob of success in the treatment grp
                    0.999,
                    mean( data[indices] > tau, na.rm = T  ) )
  var.pt  = var( data[indices] )                                       # variance of sampled treatment values
  mean.pt = mean( data[indices] )                                      # mean of sampled treatment values
  return(c(p.t, var.pt, mean.pt))
}

### Initializing 
n <- 50
I <- 10
treat <- list(mode="vector",length=I)
control <- list(mode="vector",length=I)
p_c.boot <- list(mode="vector",length=1000)
p_t.boot <- list(mode="vector",length=1000)
p_c <- numeric(I)
nnt.f.bs <- vector("list", I)
nnt.f <- vector("list", I)
ci.bs <- vector("list", I)
yc_bar <- numeric(I)
yt_bar <- numeric(I)
s_ml.bs <- list(mode="vector",length=I)
d.bs <- list(mode="vector",length=I)
diff2 <- numeric(I)
c2 <- numeric(I)
Bias <- numeric(I)



set.seed(123)
for (i in 1:I){
  treat[[i]] <- rgamma(n, shape = 1.5, rate = 1)
  control[[i]] <- rgamma(n, shape = 1.5, rate = 2)
  
  yc_bar[i]  = mean( control[[i]] )
  yt_bar[i]  = mean( treat[[i]] )
  
  n_c <- length(control[[i]])
  n_t <- length(treat[[i]])
  

# 95% quantile BS confidence interval
p_c.boot[[i]] = boot(data = control[[i]], statistic = p_c1, R = 1000)
p_t.boot[[i]] = boot(data = treat[[i]],   statistic = p_t1, R = 1000)

s_ml.bs[[i]]    = ( 1 / ( n_c + n_t )  * ( (n_c - 1) * p_c.boot[[i]]$t[ ,2] + (n_t - 1) * p_t.boot[[i]]$t[ ,2] ) ) ^ ( 1 / 2 )

d.bs[[i]]       = ( p_t.boot[[i]]$t[ ,3] - p_c.boot[[i]]$t[ ,3] ) / s_ml.bs[[i]]



nnt.f.bs[[i]]         = ifelse( 1 / ( pnorm(  d.bs[[i]] + qnorm(p_c.boot[[i]]$t[ ,1]) )  - p_c.boot[[i]]$t[ ,1] ) > 0,
                          1 / ( pnorm(  d.bs[[i]] + qnorm(p_c.boot[[i]]$t[ ,1]) )  - p_c.boot[[i]]$t[ ,1] ),
                          Inf )

# BS CI

ci.bs[[i]]            =  c( max( quantile(nnt.f.bs[[i]], .025), 1), quantile(nnt.f.bs[[i]], .975) )

## Coverage & Median Length

diff2[i] <- ci.bs[[i]][2] - ci.bs[[i]][1]

c2[i] <- ifelse(ci.bs[[i]][1] <= 3.22 & ci.bs[[i]][2] >= 3.22, 1, 0)


}
median(diff2)
mean(c2)
