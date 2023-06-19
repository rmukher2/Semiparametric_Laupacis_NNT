n <- 50
I = 10
treat <- list(mode="vector",length=I)
control <- list(mode="vector",length=I)
p_t <- numeric(I)
p_c <- numeric(I)
tau = log(3)
nntl <- list(mode="vector",length=I)
nntl_fin <- list(mode="vector",length=I)
diff <- numeric(I)

sd.wald <- numeric(I)
sd.delta <- numeric(I)

cov1 <- numeric(I)
cov2 <- numeric(I)
ci_w <- vector("list", I)
ci_d <- vector("list", I)
ci_bs <- vector("list", I)
Bias <- numeric(I)
ntl.bs <- numeric(I)


wald_diff <- numeric(I)
delta_diff <- numeric(I)

set.seed(123)
for (i in 1:I){
treat[[i]] <- rexp(n,0.4)
control[[i]] <- rexp(n,1)
p_t[i]<- mean( treat[[i]]   > tau, na.rm = T )
p_c[i]<- mean( control[[i]]   > tau, na.rm = T ) 

## point est
nntl[[i]]        = ifelse( 1 / ( p_t[[i]] - p_c[[i]] ) > 0,
                         1 / ( p_t[[i]] - p_c[[i]] ),
                         Inf )

}

nntl_fin <- unlist(nntl[is.finite(unlist(nntl))])

sum(is.infinite(unlist(nntl_fin)))

for (i in 1:I){
  
Bias <- nntl_fin - 3.21



n_t <- length(treat[[i]])
n_c <- length(control[[i]])

##Wald's CI

sd.wald[i]       = sqrt( p_t[i] * ( 1 - p_t[i] ) / n_t + p_c[i] * ( 1 - p_c[i] ) / n_c )


ci_w[[i]]          =    c( max( 1 / (  p_t[i] - p_c[i] + qnorm(.975) * sd.wald[i] ), 1),
                      ifelse( 1 / (  p_t[i] - p_c[i] - qnorm(.975) * sd.wald[i] ) > 0,
                              1 / (  p_t[i] - p_c[i] - qnorm(.975) * sd.wald[i] ),
                              Inf ))

# DELTA's CI
sd.delta[i]      = ( 1 / (p_t[i] - p_c[i]) ^ 2 ) * sqrt(   p_t[i] * (1 - p_t[i]) / n_t
                                                           + p_c[i] * (1 - p_c[i]) / n_c )
ci_d[[i]]          =  c(max( nntl[[i]] - qnorm(.975) * sd.delta[i], 1), nntl[[i]] + qnorm(.975) * sd.delta[i])


wald_diff[i] <- ci_w[[i]][2] - ci_w[[i]][1]
delta_diff[i] <- ci_d[[i]][2] - ci_d[[i]][1]

cov1[i] <- ifelse(ci_w[[i]][1] <= 3.21 & ci_w[[i]][2] >= 3.21, 1, 0)
cov2[i] <- ifelse(ci_d[[i]][1] <= 3.21 & ci_d[[i]][2] >= 3.21, 1, 0)

}

sum(is.infinite(unlist(nntl_fin)))


est <- abs(mean(Bias))
#nntl1 <- unlist(nntl)
sd_est <- sd(nntl_fin)
sd_est
RMSE <- sqrt(mean((nntl_fin - 3.21)^2))
RMSE
bias_percent <- (est/sd_est)*100
bias_percent
mean(Bias)
median(wald_diff)
sum(is.na(delta_diff))
median(delta_diff[delta_diff != 'NaN'])
mean(cov1)
mean(cov2)




