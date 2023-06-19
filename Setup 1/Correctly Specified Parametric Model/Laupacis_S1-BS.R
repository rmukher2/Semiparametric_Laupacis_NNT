### BOOTSTRAP FUNCTIONS (INCREASE) ###
library(boot)
tau = 100
p_c1   = function(data, indices) {
  p.c     = ifelse( mean( data[indices] > 100, na.rm = T ) < 0.001, # prob of success in the control grp
                    0.001,
                    mean( data[indices] > 100, na.rm = T ))
  var.pc  = var( data[indices] )                                       # variance of sampled control values
  mean.pc = mean( data[indices] )                                      # mean of sampled control values
  return(c(p.c, var.pc, mean.pc))
}

p_t1 = function(data, indices) {
  p.t     = ifelse( mean( data[indices] > 100, na.rm = T ) > 0.999, # prob of success in the treatment grp
                    0.999,
                    mean( data[indices] > 100, na.rm = T  ) )
  var.pt  = var( data[indices] )                                       # variance of sampled treatment values
  mean.pt = mean( data[indices] )                                      # mean of sampled treatment values
  return(c(p.t, var.pt, mean.pt))
}

# 95% quantile BS confidence interval
n <- 50
I <- 10
treat <- list(mode="vector",length=I)
control <- list(mode="vector",length=I)
p_c.boot <- list(mode="vector",length=1000)
p_t.boot <- list(mode="vector",length=1000)
nntl.bs <- vector("list", I)
ci_bs <- vector("list", I)
bs_diff <- numeric(I)
coverage <- numeric(I)

set.seed(123)
for (i in 1:I){
treat[[i]] <- rnorm(n,mean = 110, sd = 10)
control[[i]] <- rnorm(n,mean = 100,sd = 10)
p_c.boot[[i]] = boot(data = control[[i]], statistic = p_c1, R = 1000)
p_t.boot[[i]] = boot(data = treat[[i]],   statistic = p_t1, R = 1000)

nntl.bs[[i]]     = ifelse( 1 / ( p_t.boot[[i]]$t[ ,1] - p_c.boot[[i]]$t[ ,1] ) > 0,
                         1 / ( p_t.boot[[i]]$t[ ,1] - p_c.boot[[i]]$t[ ,1] ),
                         Inf )
##Bootstrap CI
ci_bs[[i]]        = c( max(quantile(nntl.bs[[i]], .025), 1),  quantile(nntl.bs[[i]], .975) )
##Length of CI
bs_diff[i] <- ci_bs[[i]][2] - ci_bs[[i]][1]
## Coverage
coverage[i] <- ifelse(ci_bs[[i]][1] <= 2.93 & ci_bs[[i]][2] >= 2.93, 1, 0)
}
median(bs_diff)
mean(coverage)
