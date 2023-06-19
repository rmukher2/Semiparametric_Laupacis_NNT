logit.fn <- function(data, tau=log(3)) {
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


simudata <- function(shape_t = 1.5, shape_c = 1.5, rate_t = 1, rate_c = 2, n0=50,n1=50){
  treat <- rgamma(n1, shape = shape_t, rate = rate_t)
  control <- rgamma(n0, shape = shape_c, rate = rate_c)
  
  
  
  t1=rbind(cbind(treat,1),cbind(control,0))
  colnames(t1)=c('t','d')
  # dat <- na.omit(data.frame(t,d))
  dat=as.data.frame(t1)
  return (dat)
}

logit.fn(simudata(n1= 50, n0 = 50))

btsimu<-function(I=10,J=10,shape_t = 1.5, shape_c = 1.5, rate_t = 1, rate_c = 2, n0=50,n1=50,seed1=123){
  
  res0=c()
  for (i in 1:I)
  {
    dat1=simudata(n1=n1,n0=n0)
    res1=c()
    for (j in 1:J){
      dat1bt=dat1[sample(1:(n1+n0), size = n1 + n0, replace = TRUE),]
      res1= c(res1, logit.fn(dat1bt))
    }
    ci1=c( max(quantile(res1, .025, na.rm = TRUE), 1),  quantile(res1, .975, na.rm = TRUE) )
    cid=c(ci1,ci1[2]-ci1[1])
    res0=rbind(res0,cid)
    
  }
  
  
  return(res0)
}

set.seed(123)
results <- btsimu(n1= 50, n0 = 50)
results
cov <- ifelse(results[,1] <= 3.22 & results[,2] >= 3.22,1,0)
median(results[,3])
mean(cov)

