n <- 50
I = 10
treat <- list(mode="vector",length=I)
control <- list(mode="vector",length=I)
p_t <- numeric(I)
p_c <- numeric(I)
tau = 100
nntl <- numeric(I)
diff <- numeric(I)

sd.wald <- numeric(I)
sd.delta <- numeric(I)

cov_wald <- numeric(I)
cov_dl <- numeric(I)
ci_w <- vector("list", I)
ci_d <- vector("list", I)
ci_bs <- vector("list", I)
ci_w1 <- vector("list", I)
ddd <- vector("list", I)
Bias <- numeric(I)
ntl.bs <- numeric(I)


wald_diff <- numeric(I)
delta_diff <- numeric(I)

set.seed(123)
for (i in 1:I){
treat[[i]] <- rnorm(n,mean = 110, sd = 10)
control[[i]] <- rnorm(n,mean = 100, sd = 10)
n_t <- length(treat[[i]])
n_c <- length(control[[i]])
p_t[i]<- mean( treat[[i]]   > tau, na.rm = T )
p_c[i]<- mean( control[[i]]   > tau, na.rm = T ) 


## point est
nntl[i]        = ifelse( 1 / ( p_t[i] - p_c[i] ) > 0,
                         1 / ( p_t[i] - p_c[i] ),
                         Inf )

}

for (i in 1:I){

Bias[i] <- nntl[i] - 2.93

n_t <- length(treat[[i]])
n_c <- length(control[[i]])

##Wald's CI

sd.wald[i]       = sqrt( p_t[i] * ( 1 - p_t[i] ) / n_t + p_c[i] * ( 1 - p_c[i] ) / n_c )


ci_w[[i]]          =    c( max( 1 / (  p_t[i] - p_c[i] + qnorm(.975) * sd.wald[i] ), 1),
                      ifelse( 1 / (  p_t[i] - p_c[i] - qnorm(.975) * sd.wald[i] ) > 0,
                              1 / (  p_t[i] - p_c[i] - qnorm(.975) * sd.wald[i] ),
                              Inf ))

ci_w1[[i]]          =   c(max( nntl[i] - qnorm(.975) * sd.wald[i], 1), nntl[i] + qnorm(.975) * sd.wald[i])  

# DELTA's CI
sd.delta[i]      = ( 1 / (p_t[i] - p_c[i]) ^ 2 ) * sqrt(   p_t[i] * (1 - p_t[i]) / n_t
                                                           + p_c[i] * (1 - p_c[i]) / n_c )


ci_d[[i]]          =  c(max( nntl[i] - qnorm(.975) * sd.delta[i], 1), nntl[i] + qnorm(.975) * sd.delta[i]) 


wald_diff[i] <- ci_w[[i]][2] - ci_w[[i]][1]
delta_diff[i] <- ci_d[[i]][2] - ci_d[[i]][1]


ddd[[i]] <- c(ci_w[[i]][1] - (1 / (  p_t[i] - p_c[i] )), ci_w[[i]][2] - (1 / (  p_t[i] - p_c[i] )))



cov_wald[i] <- ifelse(ci_w[[i]][1] <= 2.93 & ci_w[[i]][2] >= 2.93, 1, 0)
cov_dl[i] <- ifelse(ci_d[[i]][1] <= 2.93 & ci_d[[i]][2] >= 2.93, 1, 0)

}
sum(is.infinite(unlist(nntl)))
nntl1 <- unlist(nntl[is.finite(unlist(nntl))])
Bias <- nntl1 - 2.93
mean(Bias)
#median(nntl1)

est <- abs(mean(Bias))
sd_est <- sd(nntl1)
sd_est
RMSE <- sqrt(mean((nntl1 - 2.93)^2))
RMSE
bias_percent <- (est/sd_est)*100
bias_percent
median(wald_diff)
median(delta_diff)


median(delta_diff[delta_diff!="NaN"])
mean(cov_wald)
mean(cov_dl)

