tau = 100
### Initializing 
n <- 50
I <- 10
treat <- list(mode="vector",length=I)
control <- list(mode="vector",length=I)
p_c <- numeric(I)
nnt.f <- vector("list", I)
yc_bar <- numeric(I)
yt_bar <- numeric(I)
s_ml <- numeric(I)
d <- numeric(I)
Bias <- numeric(I)



set.seed(123)
for (i in 1:I){
  treat[[i]] <- rexp(n, rate = 0.0065)
  control[[i]] <- rexp(n, rate = 0.02)
  
  yc_bar[i]  = mean( control[[i]] )
  yt_bar[i]  = mean( treat[[i]] )
  
  n_c <- length(control[[i]])
  n_t <- length(treat[[i]])
  
  s_ml[i]    = ( 1 / ( n_c + n_t)  * ( (n_c - 1) * var( control[[i]] ) + (n_t - 1) * var( treat[[i]] ) ) ) ^ ( 1/2 )
  
  d[i]       = ( mean(treat[[i]]) - mean(control[[i]]) ) / s_ml[i]
# probability control
p_c[i]     =  ifelse( mean( control[[i]] > 100, na.rm = T ) < 0.001,
                   0.001,
                   mean( control[[i]] > 100, na.rm = T) )


# point est.

nnt.f[i]            = ifelse( 1 / ( pnorm(  d[i]   + qnorm(p_c[i]) )  - p_c[i] ) > 0,
                           1 / ( pnorm(  d[i]   + qnorm(p_c[i]) )  - p_c[i] ),
                           Inf )

}

sum(is.infinite(unlist(nnt.f)))
nnt.f1 <- unlist(nnt.f[is.finite(unlist(nnt.f))])
Bias <- nnt.f1 - 2.58
mean(Bias)
est <- abs(mean(Bias))
sd_est <- sd(nnt.f1)
sd_est
RMSE <- sqrt(mean((nnt.f1 - 2.58)^2))
RMSE
bias_percent <- (est/sd_est)*100
bias_percent
