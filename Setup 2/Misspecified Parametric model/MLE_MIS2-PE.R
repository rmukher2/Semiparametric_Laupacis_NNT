tau = log(3)

n <- 50
I <- 10
treat <- list(mode="vector",length=I)
control <- list(mode="vector",length=I)
nnt.v <- vector("list", I)
yc_bar <- numeric(I)
yt_bar <- numeric(I)
s_ml <- numeric(I)
d <- numeric(I)
s_t <- list("vector", I)
s_c <- list("vector", I)
s_t1 <- numeric(I)
s_c1 <- numeric(I)
Bias <- numeric(I)
nntl_fin <- numeric(I)


set.seed(123)
for (i in 1:I){
  treat[[i]] <- rgamma(n,shape = 1.5, rate = 1)
  control[[i]] <- rgamma(n,shape = 1.5, rate = 2)

### Initial values
  yc_bar[i]  = mean( control[[i]] )
  yt_bar[i]  = mean( treat[[i]] )
  
  n_c <- length(control[[i]])
  n_t <- length(treat[[i]])

s_t[i]     = ( ( n_t - 1 ) * var( treat[[i]] ) /  n_t ) ^ ( 1 / 2 )
s_c[i]     = ( ( n_c - 1 ) * var( control[[i]] ) / n_c ) ^ ( 1 / 2 )
s_t1[i]    = unlist(s_t[i])
s_c1[i]    = unlist(s_c[i])

## Estimator

nnt.v[[i]]      = ifelse( ( exp( - yt_bar[i] ^ (-1) * tau ) - exp( - yc_bar[i] ^ (-1) * tau ) ) ^ (-1) > 1,
                     ( exp( - yt_bar[i] ^ (-1) * tau ) - exp( - yc_bar[i] ^ (-1) * tau ) ) ^ (-1),
                     Inf )


}

sum(is.infinite(unlist(nnt.v)))

nntl_fin <- unlist(nnt.v[is.finite(unlist(nnt.v))])
Bias <- nntl_fin -3.22

est <- abs(mean(Bias))
sd_est <- sd(nntl_fin)
sd_est
RMSE <- sqrt(mean((nntl_fin - 3.22)^2))
RMSE
bias_percent <- (est/sd_est)*100
bias_percent
mean(Bias)

