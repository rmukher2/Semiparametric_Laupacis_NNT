tau = 100
n <- 50
I = 10
treat <- list(mode="vector",length=I)
control <- list(mode="vector",length=I)
p_t <- numeric(I)
p_c <- numeric(I)
tau = 100
nntl <- numeric(I)

Bias <- numeric(I)

set.seed(123)
for (i in 1:I){
  treat[[i]] <- rexp(n, rate = 0.0065)
  control[[i]] <- rexp(n, rate = 0.02)
  p_t[i]<- mean( treat[[i]]   > tau, na.rm = T )
  p_c[i]<- mean( control[[i]]   > tau, na.rm = T ) 
  
  
  ## point est
  nntl[i]        = ifelse( 1 / ( p_t[i] - p_c[i] ) > 0,
                           1 / ( p_t[i] - p_c[i] ),
                           Inf )
  

  n_t <- length(treat[[i]])
  n_c <- length(control[[i]])
  
}
sum(is.infinite(unlist(nntl)))
nntl1 <- unlist(nntl[is.finite(unlist(nntl))])
Bias <- nntl1 - 2.58
mean(Bias)


est <- abs(mean(Bias))
sd_est <- sd(nntl1)
sd_est
RMSE <- sqrt(mean((nntl1 - 2.58)^2))
RMSE
bias_percent <- (est/sd_est)*100
bias_percent


