n1 <- 50
n0 <- 50
n <- n1 + n0
rho <- n1/n0
pi <- n1/n
I <- 10
tau <- log(3)

treat <- list(mode="vector",length=I)
control <- list(mode="vector",length=I)
max.len <- list(mode="vector",length=I)
t <- list(mode="vector",length=I)
d <- list(mode="vector",length=I)
t_c <- list(mode="vector",length=I)
d_c <- list(mode="vector",length=I)
dat <- list(mode="vector",length=I)
mylogit <- list(mode="vector",length=I)
beta_c <- numeric(I)
alpha_c <- numeric(I)
pi_c <- list(mode="vector",length=I)
wt <- list(mode="vector",length=I)
p_c <- numeric(I)
p_t <- numeric(I)
NNT_D <- numeric(I)
m1 <-list(I)
A_0 <- numeric(I)
m2 <-list(I)
A_1 <- numeric(I)
m3 <-list(I)
A_2 <- numeric(I)
A <- list(I)
A_inv <- list(I)
A0_t <- list(I)
A1_t <- list(I)
ma1 <- list(I)
mul <- numeric(I)
f_t <- numeric(I)
g_t <- numeric(I)
C <- numeric(I)
NP_var <- numeric(I)
SP_var <- numeric(I)
sdsp <- numeric(I)
NNT_D <- numeric(I)
g_tilde <- numeric(I)
fhat_t <- numeric(I)
ghat_t <- numeric(I)
SP <- numeric(I)
mul1 <- numeric(I)
NP_v1 <- numeric(I)
ci_dr <- vector("list", I)
delta_diff <- numeric(I)
cov_dl <- numeric(I)
Bias <- numeric(I)

set.seed(123)
for (i in 1:I){
  treat[[i]] <- rgamma(n1,shape = 1.5, rate = 1)
  control[[i]] <- rgamma(n0,shape = 1.5, rate = 2)
  
  
  t[[i]]=rbind(cbind(treat[[i]],1),cbind(control[[i]],0))
  colnames(t[[i]])=c('t','d')
  dat[[i]]=as.data.frame(t[[i]])
  
  mylogit[[i]] <- glm(d ~ t, family=binomial(link="logit"), data = dat[[i]])
  beta_c[i] <- coef(mylogit[[i]])[2]
  alpha_c[i] <- ifelse((1-pi)/pi == 1 | (1-pi)/pi < 0 , coef(mylogit[[i]])[1], coef(mylogit[[i]])[1] + log((1-pi)/pi))
  
  pi_c[[i]] <- 1/(n0*(1 + rho*exp(alpha_c[i] + dat[[i]]$t*beta_c[i])))
  wt[[i]] <- exp(alpha_c[i] + dat[[i]]$t*beta_c[i])
  
  p_c[i] <- sum(ifelse(dat[[i]]$t <= tau, pi_c[[i]], 0))
  p_t[i] <- sum(ifelse(dat[[i]]$t <= tau, pi_c[[i]]*wt[[i]], 0))
  
  NNT_D[i] <- ifelse(1/(p_c[i] - p_t[i]) > 0,
                     1/(p_c[i] - p_t[i]),Inf)
  
  f_t[i]<- mean( treat[[i]]   <= log(3), na.rm = T )
  g_t[i]<- mean( control[[i]] <= log(3), na.rm = T ) 
  
  C[i] <- 1/(g_t[i] - f_t[i])^2
  
  
  A0_t[i] <- sum(ifelse(dat[[i]]$t <= log(3),  pi_c[[i]] *
                          (exp(alpha_c[i] + dat[[i]]$t*beta_c[i]))/(1+ rho *exp(alpha_c[i] + dat[[i]]$t*beta_c[i])),0))
  A1_t[i] <- sum(ifelse(dat[[i]]$t <= log(3),  pi_c[[i]] *
                          (exp(alpha_c[i] + dat[[i]]$t*beta_c[i])/(1+ rho *exp(alpha_c[i] + dat[[i]]$t*beta_c[i])))*dat[[i]]$t
                        ,0))
  
  A_0[i] <- sum(pi_c[[i]] * (exp(alpha_c[i] + dat[[i]]$t*beta_c[i]))/(1+ rho * exp(alpha_c[i] + dat[[i]]$t*beta_c[i])))
  A_1[i] <- sum(pi_c[[i]] *(exp(alpha_c[i] + dat[[i]]$t*beta_c[i])/(1+ rho * exp(alpha_c[i] + dat[[i]]$t*beta_c[i])))*dat[[i]]$t)
  A_2[i] <- sum(pi_c[[i]] *(exp(alpha_c[i] + dat[[i]]$t*beta_c[i])/(1+ rho *exp(alpha_c[i] + dat[[i]]$t*beta_c[i])))
                *t(dat[[i]]$t)*dat[[i]]$t)
  
  
  A[[i]] = matrix(c(A_0[i],A_1[i],t(A_1[i]),A_2[i]),nrow = 2,ncol = 2, byrow = TRUE)
  A_inv[[i]] <- solve(A[[i]])
  
  ma1[[i]] <- matrix(c(A0_t[[i]],A1_t[[i]]),nrow = 1, ncol = 2, byrow = T)
  mul[i] <- A0_t[[i]] - (ma1[[i]] %*% A_inv[[i]] %*%  t(ma1[[i]]))
  
  
  NP_var[i] <-  (  g_t[i] * (1 - g_t[i]) /n0 + f_t[i] * (1 - f_t[i])/n1 )
  
  SP_var[i] <-  (1 + 1/rho)^2 *(rho /  n0) * mul[i]
  
  sdsp[i] <- ( 1 / (p_c[i] - p_t[i]) ^ 2  ) *sqrt((NP_var[i] - SP_var[i]) )
  
  
  
  ci_dr[[i]]          =  c(max( NNT_D[i] - qnorm(.975) * sdsp[i], 1), NNT_D[i] + qnorm(.975) * sdsp[i]) 
  
  delta_diff[i] <- ci_dr[[i]][2] - ci_dr[[i]][1]
  cov_dl[i] <- ifelse(ci_dr[[i]][1] <= 3.22 & ci_dr[[i]][2] >= 3.22, 1, 0)
  
  
  
  Bias[i] <- NNT_D[i] - 3.22
}
est <- abs(mean(Bias))
sd_est <- sd(NNT_D)
sd_est
RMSE <- sqrt(mean((NNT_D - 3.22)^2))
RMSE
bias_percent <- (est/sd_est)*100
bias_percent
mean(Bias)

#median(NNT_D)
median(delta_diff, na.rm = T)
mean(cov_dl, na.rm = T)

median(delta_diff)
mean(cov_dl)










