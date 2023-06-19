logit.fn <- function(data, tau=100) {
  fit<-glm(d ~ t, data = data, family = "binomial")  
  n1 <- length(which(data$d==1))
  n0 <- length(which(data$d==0))
  rho <- n1/n0
  n <- n1 + n0
  pi <- n1/n
  a1 <- coef(fit)[1] + log((1-pi)/pi)
  b1 <- coef(fit)[2]
  pi_c <- 1/(n0*(1 + rho*exp(a1 + data$t*b1)))
  wt <- exp(a1 + data$t*b1)
  
  p_c <- sum(ifelse(data$t <= tau, pi_c, 0))
  p_t <- sum(ifelse(data$t<= tau, pi_c*wt, 0))
  A=  ifelse (1/(p_c-p_t) > 0 ,
              1/(p_c-p_t), Inf)
  #
  return(A)
}


simudata <- function(mn_t=110,mn_c=100,sd_t=10,sd_c=10,n0=50,n1=50){
  treat <- rnorm(n1,mean = mn_t, sd = sd_t)
  control <- rnorm(n0,mean = mn_c, sd = sd_c)
  
  
  t1=rbind(cbind(treat,1),cbind(control,0))
  colnames(t1)=c('t','d')
  # dat <- na.omit(data.frame(t,d))
  dat=as.data.frame(t1)
  return (dat)
}

logit.fn(simudata(n1= 50, n0 = 50))

btsimu<-function(I=10,J=10,n1=50,n0=50,mn_t=110,mn_c=100,sd_t=10,sd_c=10,seed1=123){
  #set.seed(seed1)
  res0=c()
  for (i in 1:I)
  {
    dat1=simudata(n1=n1,n0=n0)
    res1=c()
    for (j in 1:J){
      dat1bt=dat1[sample(1:(n1+n0), size = n1 + n0, replace = TRUE),]
      res1= c(res1, logit.fn(dat1bt))
    }
    ci1=c( max(quantile(res1, .025), 1),  quantile(res1, .975) )
    cid=c(ci1,ci1[2]-ci1[1])
    res0=rbind(res0,cid)
    
  }
  
  
  return(res0)
}

set.seed(123)
results <- btsimu(n1= 50, n0 = 50)
results
cov <- ifelse(results[,1] <= 2.93 & results[,2] >= 2.93,1,0)
median(results[,3])
mean(cov)

