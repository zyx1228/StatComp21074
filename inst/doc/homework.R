## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
Z<-function(n){
x<-seq(-7,7,0.01) 
truth<-dnorm(x,0,2*sqrt(2)) 
plot(density(rnorm(n,0,2)+rnorm(n,0,2)),main="density estimate of the normal addition model",ylim=c(0,0.25),lwd=2,lty=2) 
lines(x,truth,col="red",lwd=2) 
legend("topright",c("true","estimated"),col=c("red","black"),lwd=2,lty=c(1,2))
}
Z(1000)

## -----------------------------------------------------------------------------
library(xtable)
library(gt)
library(datasets)
library(knitr)

data(airquality)
table<-kable(airquality[1:7,])
table

gt_table<-gt(data=airquality[1:7,])
gt_table <-
	gt_table %>%
	tab_header(
   	 	title=md("**New York Air Quality Measurements**"),
    	subtitle=md("*The first seven days of May, 1973*") 
  	) %>%
 	tab_source_note(
 		source_note="Source: The data were obtained from the New York State Department of Conservation (ozone data) and the National Weather Service (meteorological data)."
	) %>%
	tab_source_note(
		source_note=md("Reference: Chambers, J. M., Cleveland, W. S., Kleiner, B. and Tukey, P. A. (1983) *Graphical Methods for Data Analysis*. Belmont, CA: Wadsworth.")
	) 

gt_table


## -----------------------------------------------------------------------------
n<-1e3;

rleigh<-function(sigma){
u<-runif(n)
x<-sqrt(-2*(sigma^2)*log(1-u))
y<-seq(0,100,0.1)
hist(x,probability = TRUE,main=paste ("sigma=",as.character(sigma),sep=""),border=NA)
lines(y,(y/(sigma^2))*exp(-y^2/(2*sigma^2)))
}

par(mfrow=c(2,3))
rleigh(4)
rleigh(5)
rleigh(10)
rleigh(15)
rleigh(20)

## -----------------------------------------------------------------------------
n<-1000
x1<-rnorm(n,0,1)
x2<-rnorm(n,3,1)
r1<-sample(c(1,0),n,replace=TRUE,prob=c(0.75,0.25))
z1<-r1*x1+(1-r1)*x2
r2<-sample(c(1,0),n,replace=TRUE,prob=c(0.6,0.4))
z2<-r2*x1+(1-r2)*x2
r3<-sample(c(1,0),n,replace=TRUE,prob=c(0.5,0.5))
z3<-r3*x1+(1-r3)*x2
r4<-sample(c(1,0),n,replace=TRUE,prob=c(0.4,0.6))
z4<-r4*x1+(1-r4)*x2
r5<-sample(c(1,0),n,replace=TRUE,prob=c(0.1,0.9))
z5<-r5*x1+(1-r5)*x2

par(mfrow=c(2,3))
hist(z1,prob=TRUE,main=expression(p==0.75));
hist(z2,prob=TRUE,main=expression(p==0.6));
hist(z3,prob=TRUE,main=expression(p==0.5));
hist(z4,prob=TRUE,main=expression(p==0.4));
hist(z5,prob=TRUE,main=expression(p==0.1));

## -----------------------------------------------------------------------------
n<-1e4;alpha<-4;beta<-3;lambda<-5;t<-10 
x<-numeric(n) 

for (i in 1:n){
  Nt<-rpois(1,lambda*t)
  y<-rgamma(Nt,alpha,beta)
  x[i]<-sum(y)
}

meanm<-mean(x)
varm<-var(x)
meanr<-lambda*t*(alpha/beta)
varr<-lambda*t*(alpha*(1+alpha)/beta^2)
meanm
varm
meanr
varr

## -----------------------------------------------------------------------------
set.seed(0)
#set parameters
x<-seq(0.1,0.9,0.1)
n<-1e3
alpha<-3
beta<-3
F1<-numeric(length(x))
B<-factorial(alpha+beta-1)/(factorial(alpha-1)*factorial(beta-1))

#function to caculate MC Integration
F<-function(t){
  y<-runif(n,0,t)
  m<-B*t*mean(y^(alpha-1)*(1-y)^(beta-1))
}
for (i in 1:length(x)){
  F1[i]<-F(x[i])
}

plot(x,F1,type='l',xlab='x',ylab='F(x)',col=1,lwd=2)
lines(x,pbeta(x,alpha,beta),lty=2,col=2)
exbeta<-c(expression(paste("Monte Carlo Estimate")),expression(paste("pbeta function in R")))
legend("bottomright",exbeta,lty=c(1,2),col=c(1,2),lwd=2)

## -----------------------------------------------------------------------------
set.seed(0)
n<-1e3
x<-2
sig<-4 #Here we set the upper limit x=2, parameter sigma=4. You can choose any value of these two parameters. 
u<-runif(n,0,x)
ua<-c(u[1:(n/2)])#use half of u to construct antithetic variable

mcm<-x*mean(u*exp(-u^2/2*sig^2)/sig^2)
mcv<-var(u*exp(-u^2/2*sig^2)/sig^2)/n
antimcm<-x*sum(ua*exp(-ua^2/2*sig^2)/sig^2+(x-ua)*exp(-(x-ua)^2/2*sig^2)/sig^2)/n
antimcv<-var(ua*exp(-ua^2/2*sig^2)/sig^2+(x-ua)*exp(-(x-ua)^2/2*sig^2)/sig^2)/(2*n)

cat("The percent reduction in variance is",antimcv/mcv)

## -----------------------------------------------------------------------------
x<-seq(1,10,0.1)
# we also tried x<-seq(1,100,1), but the g(x) close to 0 as x increase
w<-2
g<-x^2*exp(-x^2/2)/sqrt(2*pi)
f1<-sqrt(exp(1))*x^2*exp(-x/2)/10
f2<-sqrt(exp(1))*x*exp(-x^2/2)

#figure (a)
plot(x,g,type = "l",main = "(a)",ylab="",col=1,lwd = w)
lines(x,g/g,lty = 2,col=2,lwd = w)
lines(x,f1,lty = 3,col=3,lwd = w)
lines(x,f2,lty = 4,col=4,lwd = w)
legend("topright",legend=c("g","g/g","f1","f2"),lty=c(1:4),col=c(1:4),lwd = w,inset = 0.02)

#figure (b)
plot(x,g/g,type = "l",main = "(b)",lty=2,col=2,ylab="",lwd = w)
lines(x,g/f1,lty = 3,col=3,lwd = w)
lines(x,g/f2,lty = 4,col=4,lwd = w)
legend("topright",legend =c("g/g","g/f1","g/f2"),lty=c(2:4),col=c(2:4),lwd = w, inset = 0.02)

## -----------------------------------------------------------------------------
set.seed(0)
n<-1e3
u<-runif(n)
x<-sqrt(-2*log(1-u*(1-1/sqrt(exp(1)))))#inverse method
EX<-mean(x)*(1-1/sqrt(exp(1)))/sqrt(2*pi)
cat("The estimating integral is",1/2-EX)

## -----------------------------------------------------------------------------
set.seed(1)
n<-20
m<-1000 
s<-0 
alpha<-0.05
for (i in 1:m){
  x<-rchisq(n,df=2)
  y1<-mean(x)-qt(1-alpha/2, df=n-1)*sd(x)/sqrt(n) 
  y2<-mean(x)+qt(1-alpha/2, df=n-1)*sd(x)/sqrt(n)
  if (y1<2 & y2>2) {
    s<-s+1
  }
}
cat("The coverage probability of the t-interval is",s/m)

## -----------------------------------------------------------------------------
set.seed(1)
n<-20
m<-1000 
s<-0 
alpha<-0.05
mu0<-1
p<-numeric(m)
for (i in 1:m){
  x<-rchisq(n,df=1)
  y1<-mean(x)-qt(1-alpha/2, df=n-1)*sd(x)/sqrt(n) 
  y2<-mean(x)+qt(1-alpha/2, df=n-1)*sd(x)/sqrt(n)
  if (y1<1 & y2>1) {
    s<-s+1
  }
  ttest<-t.test(x,mu=mu0)
  p[i]<-ttest$p.value
}
p.hat<-mean(p<alpha)
se.hat<-round(sqrt(p.hat*(1-p.hat)/m),3)
cat("The coverage probability of the t-interval is",s/m,". The p-value is",p.hat,". And the standard error of the estimate is approximately",se.hat,".")

## -----------------------------------------------------------------------------
set.seed(1)
n<-20
m<-1000 
s<-0 
alpha<-0.05
mu0<-1
p<-numeric(m)
for (i in 1:m){
  x<-runif(n,0,2)
  y1<-mean(x)-qt(1-alpha/2, df=n-1)*sd(x)/sqrt(n) 
  y2<-mean(x)+qt(1-alpha/2, df=n-1)*sd(x)/sqrt(n)
  if (y1<1 & y2>1) {
    s<-s+1
  }
  ttest<-t.test(x,mu=mu0)
  p[i]<-ttest$p.value
}
p.hat<-mean(p<alpha)
se.hat<-round(sqrt(p.hat*(1-p.hat)/m),3)
cat("The coverage probability of the t-interval is",s/m,". The p-value is",p.hat,". And the standard error of the estimate is approximately",se.hat,".")

## -----------------------------------------------------------------------------
set.seed(1)
n<-20
m<-1000 
s<-0 
alpha<-0.05
mu0<-1
p<-numeric(m)
for (i in 1:m){
  x<-rexp(n,1)
  y1<-mean(x)-qt(1-alpha/2, df=n-1)*sd(x)/sqrt(n) 
  y2<-mean(x)+qt(1-alpha/2, df=n-1)*sd(x)/sqrt(n)
  if (y1<1 & y2>1) {
    s<-s+1
  }
  ttest<-t.test(x,mu=mu0)
  p[i]<-ttest$p.value
}
p.hat<-mean(p<alpha)
se.hat<-round(sqrt(p.hat*(1-p.hat)/m),3)
cat("The coverage probability of the t-interval is",s/m,". The p-value is",p.hat,". And the standard error of the estimate is approximately",se.hat,".")

## -----------------------------------------------------------------------------
library(MASS)
set.seed(123)
d <- 3 #dimension
m <- 100 #sample sizes of A
n <- 20  #sample sizes of X in each A
sktests <- numeric(m) 
cv <- qchisq(.95, d*(d+1)*(d+2)/6) #crit. values for each n

skn<-function(x){
    xbar <-  apply(x,1,mean)
    y <- matrix(nrow = d, ncol = n)
    s <- matrix(nrow = 1, ncol = n)
    for (a in 1:n){
        y[,a] <- x[,a] - xbar
    }
    cov <- cov(t(x))
    for (b in 1:n){
        s[b] <- 0
        for (c in 1:n){
            s[b] <- s[b]+(t(y[,b])%*%ginv(cov)%*%y[,c])^3
         }
    }
    return(sum(s)/n/6)
}

for (k in 1:m) {
    x<-matrix(rnorm(d*n), nrow = d, ncol = n)
   #test decision is 1 (reject) or 0
    sktests[k] <- as.integer(skn(x) >= cv)
}
p.reject <- mean(sktests) #proportion rejected
p.reject

## -----------------------------------------------------------------------------
result<-matrix(nrow = 1, ncol = 4)
result[1,] <- c(0.000,0.000,0.030,0.060)
dimnames(result)[[2]] <- c("n=10","n=20","n=50","n=100") 
knitr::kable(result)

## -----------------------------------------------------------------------------
library(MASS)
set.seed(1)
alpha <- .1
n <- 30
m <- 100
epsilon <- c(seq(0,1,0.05))
N <- length(epsilon)
pwr <- numeric(N)
cv <- qchisq(1-alpha, d*(d+1)*(d+2)/6)

for (j in 1:N) { #for each epsilon
  e <- epsilon[j]
  sktests <- numeric(m)
  for (i in 1:m) { #for each replicate
    sigma <- sample(c(1, 10), replace = TRUE, size = d*n, prob = c(1-e, e))
    x <- matrix(rnorm(d*n,0,sigma), nrow = d, ncol = n)
    sktests[i] <- as.integer(abs(skn(x)) >= cv)
  }
  pwr[j] <- mean(sktests)
}

#plot power vs epsilon
plot(epsilon, pwr, type = "b",xlab = bquote(epsilon), ylim = c(0,1))
abline(h = .1, lty = 3)
se <- sqrt(pwr * (1-pwr) / m) #add standard errors
lines(epsilon, pwr+se, lty = 3)
lines(epsilon, pwr-se, lty = 3)

## -----------------------------------------------------------------------------
library(boot); library(bootstrap); library(MASS); 
set.seed(12345)

b.pca <- function(x,i){
  lam <- eigen(cov(x[i,]))$values
  theta_hat <- lam[1] / sum(lam)
} 
x <- cbind(scor$mec, scor$vec, scor$alg, scor$ana, scor$sta)
obj <- boot(x,statistic = b.pca, R = 200)

round(c(original=obj$t0,bias.boot=mean(obj$t)-obj$t0,se.boot=sd(obj$t)),3)

## -----------------------------------------------------------------------------
n<-nrow(scor)
theta_hat <- b.pca(x,1:n)
theta_jack <- numeric(n)

for(i in 1:n){
theta_jack[i] <- b.pca(x,(1:n)[-i])
}

bias.jack <- (n-1)*(mean(theta_jack)-theta_hat)
se.jack <- sqrt((n-1)*mean((theta_jack-theta_hat)^2))

round(c(original=theta_hat,bias.jack=bias.jack,se.jack=se.jack),3)

## -----------------------------------------------------------------------------
ci <- matrix(NA,2,2)

CI <- boot.ci(obj, conf = 0.95, type = c('perc', 'bca'))
ci[1,]<-CI$percent[4:5];ci[2,]<-CI$bca[4:5]

result<-matrix(nrow = 1, ncol = 4)
result[1,] <- c(ci[1,1],ci[1,2],ci[2,1],ci[2,2])
dimnames(result)[[2]] <- c("perc lower","perc upper","BCa lower","BCa upper") 
knitr::kable(result)

## -----------------------------------------------------------------------------
#normal populations (skewness=0)
library(boot);library(moments)
set.seed(12345)
skew<-0;n<-20;m<-1e3;
ci.norm<-ci.basic<-ci.perc<-matrix(NA,m,2)

boot.sk <- function(x,i){
  skewness(x[i]) 
}

#bootstrap
for(i in 1:m){
U <- rnorm(n,0,1)
de <- boot(data=U , statistic=boot.sk, R = 200)
ci <- boot.ci(de , conf = 0.95 , type=c("norm","basic","perc"))
ci.norm[i,] <- ci$norm[2:3]
ci.basic[i,] <- ci$basic[4:5]
ci.perc[i,] <- ci$percent[4:5]
}

#coverage probabilities
cat('norm =',mean(ci.norm[,1]<=skew & ci.norm[,2]>=skew),
'basic =',mean(ci.basic[,1]<=skew & ci.basic[,2]>=skew),
'perc =',mean(ci.perc[,1]<=skew & ci.perc[,2]>=skew))

#miss on the left
cat('norm left=',mean(ci.norm[,1]>skew),
'basic left=',mean(ci.basic[,1]>skew),
'perc left=',mean(ci.perc[,1]>skew))

#miss on the right
cat('norm right=',mean(ci.norm[,2]<skew),
'basic right=',mean(ci.basic[,2]<skew),
'perc right=',mean(ci.perc[,2]<skew))

## -----------------------------------------------------------------------------
#chisq populations 
library(boot);library(moments)
set.seed(12345)
skew<-sqrt(8/5);n<-10;m<-1e3;
ci.norm<-ci.basic<-ci.perc<-matrix(NA,m,2)

boot.sk <- function(x,i){
  skewness(x[i]) 
}

#bootstrap 
for(i in 1:m){
U <- rchisq(n,5)
de <- boot(data=U , statistic=boot.sk, R = 200)
ci <- boot.ci(de , conf = 0.95 , type=c("norm","basic","perc"))
ci.norm[i,] <- ci$norm[2:3]
ci.basic[i,] <- ci$basic[4:5]
ci.perc[i,] <- ci$percent[4:5]
}

#coverage probabilities
cat('norm =',mean(ci.norm[,1]<=skew & ci.norm[,2]>=skew),
'basic =',mean(ci.basic[,1]<=skew & ci.basic[,2]>=skew),
'perc =',mean(ci.perc[,1]<=skew & ci.perc[,2]>=skew))

#miss on the left
cat('norm left=',mean(ci.norm[,1]>skew),
'basic left=',mean(ci.basic[,1]>skew),
'perc left=',mean(ci.perc[,1]>skew))

#miss on the right
cat('norm right=',mean(ci.norm[,2]<skew),
'basic right=',mean(ci.basic[,2]<skew),
'perc right=',mean(ci.perc[,2]<skew))

## -----------------------------------------------------------------------------
library(boot)
library(energy)
library(Ball)
library(RANN)

## -----------------------------------------------------------------------------
set.seed(12345); 
n1 <- n2 <- 50; 
K <- 1:100; 
m <- 999;
reps <- numeric(m); 
x <- rchisq(n1,2); 
y <- rgamma(n2,0.5,7);
plot(y~x,main="scatter diagram",xlab="X",ylab="y")

z <- c(x,y) #sample pool
s0 <- cor.test(x,y, method ="spearman",continuity=TRUE,conf.level=0.95)$statistic

for (i in 1:m) {
k <- sample(K, size = n1, replace = FALSE)
x1 <- z[k];y1 <- z[-k] #complement of x1
reps[i] <- cor.test(x1, y1, method ="spearman",continuity=TRUE,conf.level=0.95)$statistic
}
p <- mean(abs(c(s0, reps)) >= abs(s0))

round(c(p,cor.test(x,y,method ="spearman",continuity=TRUE,conf.level=0.95)$p.value),3)

## -----------------------------------------------------------------------------
set.seed(123)
m <- 1e2; k <- 3; p <- 2; 
n1 <- n2 <- 50; R <- 99; n <- n1+n2; N = c(n1,n2)

#KNN
Tn3 <- function(z, ix, sizes,k) {
  n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
  if(is.vector(z)) z <- data.frame(z,0);
  z <- z[ix, ];
  NN <- nn2(data=z, k=k+1) 
  block1 <- NN$nn.idx[1:n1,-1]
  block2 <- NN$nn.idx[(n1+1):n,-1]
  i1 <- sum(block1 <= n1); i2 <- sum(block2 > n1)
  (i1 + i2) / (k * n)
}
eqdist.nn <- function(z,sizes,k){
  boot.obj <- boot(data = z, statistic = Tn3, R = R, sim = "permutation", sizes = sizes, k = k)
  ts <- c(boot.obj$t0,boot.obj$t)
  p.value <- mean(ts >= ts[1])
  list(statistic=ts[1],p.value=p.value)
}

p.values <- matrix(NA,m,3)
for(i in 1:m){
  x <- matrix(rnorm(n1*p,0,1.5),ncol=p);
  y <- cbind(rnorm(n2),rnorm(n2));
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=x,y=y,num.permutations=99,seed=i*12)$p.value
}
alpha <- 0.1;
pow1 <- colMeans(p.values<alpha)
cat(" power of NN =",pow1[1],'\n',"power of energy =",pow1[2],'\n', "power of ball =", pow1[3])

## -----------------------------------------------------------------------------
set.seed(123)
p.values <- matrix(NA,m,3)
for(i in 1:m){
  x <- matrix(rnorm(n1*p,1,3),ncol=p);
  y <- cbind(rnorm(n2,0.5,4),rnorm(n2,0.5,4));
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=x,y=y,num.permutations=999,seed=i*12)$p.value
}
alpha <- 0.1;
pow2 <- colMeans(p.values<alpha)
cat(" power of NN =",pow2[1],'\n',"power of energy =",pow2[2],'\n', "power of ball =", pow2[3])

## -----------------------------------------------------------------------------
set.seed(123)
p.values <- matrix(NA,m,3)
for(i in 1:m){
  x <- matrix(rt(n1*p,1),ncol=p);
  mu <- sample(c(0, 1), replace = TRUE, size = n2*p, prob = c(0.3, 0.7))
  y <- matrix(rnorm(n2*p,mu,1), ncol = p)
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=x,y=y,num.permutations=99,seed=i*12)$p.value
}
alpha <- 0.1;
pow3 <- colMeans(p.values<alpha)
cat(" power of NN =",pow3[1],'\n',"power of energy =",pow3[2],'\n', "power of ball =", pow3[3])

## -----------------------------------------------------------------------------
set.seed(123)
n1<-10;n2<-100;n=n1+n2;N=c(n1,n2)
p.values <- matrix(NA,m,3)
for(i in 1:m){
  x <- matrix(rnorm(n1*p,1,1),ncol = p);
  y <- cbind(rnorm(n2),rnorm(n2))
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=x,y=y,num.permutations=99,seed=i*12)$p.value
}
alpha <- 0.1;
pow4 <- colMeans(p.values<alpha)
cat(" power of NN =",pow4[1],'\n',"power of energy =",pow4[2],'\n', "power of ball =", pow4[3])


## -----------------------------------------------------------------------------
power <- as.data.frame(cbind(pow1, pow2, pow3, pow4)) 
colnames(power) <- c("unequal var","unequal var&ex","Non-normal distributions","Unbalanced samples")
rownames(power) <- c("NN", "Energy", "Ball")
knitr::kable(power) 

## ----eval=FALSE---------------------------------------------------------------
#  #standard cauchy dist
#  f <- function(x, theta) {
#  return(1/(theta*pi*(1+(x/theta)^2)))
#  }
#  
#  #MCMC
#  set.seed(12)
#  m <- 5000
#  theta <- 1
#  sigma <- 2
#  burn <- 1000
#  x <- numeric(m)
#  
#  x[1] <- rnorm(1, mean = 0, sd = theta) #proposal dist
#  k <- 0
#  u <- runif(m)
#  
#  for (i in 2:m) {
#  xt <- x[i-1]
#  y <- rnorm(1, mean = xt, sd = sigma)
#  num <- f(y, theta) * dnorm(xt, y, sigma)
#  den <- f(xt, theta) * dnorm(y, xt, sigma)
#  if (u[i] <= num/den) x[i] <- y else {
#  x[i] <- xt
#  k <- k+1 #y is rejected
#  }
#  }
#  
#  cat('reject probability: ',round(k/m,3))
#  
#  #trace plot
#  index <- (burn+1):m
#  y <- x[index]
#  plot(index, y, type="l", main="", ylab="x")
#  
#  #QQ plot
#  a <- ppoints(10)
#  QR <- tan(pi*(a-1/2)) #deciles of standard cauchy
#  Q <- quantile(y, a)
#  qqplot(QR, Q, main="", xlab="standard cauchy Quantiles", ylab="Sample Quantiles")
#  abline(0,1,col='blue',lwd=2)

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(12)
#  #initialize constants and parameters
#  N <- 5000 #length of chain
#  burn <- 1000 #burn-in length
#  X <- matrix(0, N, 2) #the chain, a bivariate sample
#  
#  #conditional parameters
#  n <- 16
#  alpha <- 2
#  beta <- 4
#  
#  ###### generate the chain #####
#  X[1, ] <- c(0.5, 0.5) #initialize
#  for (i in 2:N) {
#  x2 <- X[i-1, 2]
#  X[i, 1] <- rbinom(1, n, x2)
#  x1 <- X[i, 1]
#  X[i, 2] <- rbeta(1, x1 + alpha, n - x1 + beta)
#  }
#  
#  b <- burn + 1
#  x <- X[b:N, ]
#  plot(x[,1],type='l',col=1,lwd=2,xlab='Index',ylab='Random numbers')
#  lines(x[,2],col=2,lwd=2)
#  legend('topright',c(expression(X[1]),expression(X[2])),col=1:2,lwd=2)
#  hist(x, breaks="scott", main="", xlab="", freq=FALSE)
#  plot(x, main="", cex=.5, xlab=bquote(X[1]),ylab=bquote(X[2]), ylim=range(x[,2]))
#  

## ----eval=FALSE---------------------------------------------------------------
#  #calculate G.R statistics
#  Gelman.Rubin <- function(psi) {
#          # psi[i,j] is the statistic psi(X[i,1:j])
#          # for chain in i-th row of X
#          psi <- as.matrix(psi)
#          n <- ncol(psi)
#          k <- nrow(psi)
#  
#          psi.means <- rowMeans(psi)     #row means
#          B <- n * var(psi.means)        #between variance est.
#          psi.w <- apply(psi, 1, "var")  #within variances
#          W <- mean(psi.w)               #within est.
#          v.hat <- W*(n-1)/n + (B/n)     #upper variance est.
#          r.hat <- v.hat / W             #G-R statistic
#          return(r.hat)
#  }
#  
#  #9.3 MC
#  f <- function(x, theta) {
#  return(1/(theta*pi*(1+(x/theta)^2)))
#  }
#  
#  cauchy.chain <- function(sigma, N, X1) {
#          x <- rep(0, N)
#          x[1] <- X1
#          u <- runif(N)
#          for (i in 2:N) {
#            xt <- x[i-1]
#            y <- rnorm(1, mean = xt, sd = sigma)
#            num <- f(y, theta = 1) * dnorm(xt, y, sigma)
#            den <- f(xt, theta = 1) * dnorm(y, xt, sigma)
#            if (u[i] <= num/den) x[i] <- y else {
#                x[i] <- xt
#            }
#          }
#          return(x)
#  }
#  
#  sigma <- 2      #parameter of proposal distribution
#  k <- 4          #number of chains to generate
#  n <- 18000      #length of chains
#  b <- 2000       #burn-in length
#  
#  #choose overdispersed initial values
#  x0 <- c(-8, -5, 5, 8)
#  
#  #generate the chains
#  set.seed(1)
#  X <- matrix(0, nrow=k, ncol=n)
#  for (i in 1:k)
#      X[i, ] <- cauchy.chain(sigma, n, x0[i])
#      #compute diagnostic statistics
#      psi <- t(apply(X, 1, cumsum))
#  for (i in 1:nrow(psi))
#      psi[i,] <- psi[i,] / (1:ncol(psi))
#  
#  #plot psi for the four chains
#  par(mfrow=c(2,2))
#  for (i in 1:k)
#      plot(psi[i, (b+1):n], type="l",xlab='Index', ylab=bquote(phi))
#  
#  #plot the sequence of R-hat statistics
#  rhat <- rep(0, n)
#  for (j in (b+1):n)
#  rhat[j] <- Gelman.Rubin(psi[,1:j])
#  
#  plot(rhat[(b+1):n], type="l", xlab="", ylab="R")
#  abline(h=1.2, lty=2)

## ----eval=FALSE---------------------------------------------------------------
#  #9.8 MC
#  bb.chain <- function(sigma, N, X1) {
#  #generates a Metropolis chain for Normal(0,1)
#  #with Normal(X[t], sigma) proposal distribution
#  #and starting value X1
#  x <- matrix(0, N, 2)
#  x[1,] <- X1
#  
#  for (j in 2:N) {
#  x2 <- x[j-1, 2]
#  x[j, 1] <- rbinom(1, size =16, prob = x2)
#  x1 <- x[j, 1]
#  x[j, 2] <- rbeta(1, shape1 = x1 + 2, shape2 = 20 - x1)
#  #n=16,alpna=2,beta=4
#  }
#  return(as.vector(t(x)))
#  }
#  
#  k <- 4
#  #number of chains to generate
#  N <- 15000
#  #length of chains
#  b <- 1000
#  #burn-in length
#  
#  #choose overdispersed initial values
#  x0 <- matrix(c(-3, -1, 1, 3, -1, 0, 0, 1), nrow = 4, ncol = 2)
#  
#  #generate the chains
#  set.seed(12345)
#  X <- matrix(0, nrow = 4, ncol = 2 * N)
#  for (i in 1:k)
#  X[i, ] <- bb.chain(sigma, N, x0[i,])
#  
#  #compute diagnostic statistics
#  psi <- t(apply(X, 1, cumsum))
#  for (i in 1:nrow(psi))
#  psi[i,] <- psi[i,] / (1:ncol(psi))
#  
#  #plot psi for the four chains
#  par(mfrow=c(2,2))
#  for (i in 1:k)
#  plot(psi[i, (b+1):2*N], type="l",xlab='Index', ylab=bquote(phi))
#  
#  #plot the sequence of R-hat statistics
#  rhat <- rep(0, 2*N)
#  for (l in (b+1):2*N)
#  rhat[l] <- Gelman.Rubin(psi[,1:l])
#  
#  plot(rhat[(b+1):2*N], type="l", xlab="", ylab="R")
#  abline(h=1.2, lty=2)

## -----------------------------------------------------------------------------
explg <- function(k,a,d){
  (-1)^k*dist(a)^(2*k+2)*exp(lgamma((d+1)/2)+lgamma(k+(3/2))-lgamma(k+1+(d/2)))/(exp(log(2*k+1)+log(2*k+2)+log(prod(1:k))+k*log(2)))
}

rk <- 3
a <- c(1, 2)
rd <- length(a)
z <- numeric(rd)
ra <- rbind(a, z)

cat("The third term is",round(explg(rk, ra, rd), 5))

## -----------------------------------------------------------------------------
rk <- 200
a <- c(1, 2)
rd <- length(a)
z <- numeric(rd)
ra <- rbind(a, z)

cat("The 200th term is",explg(rk, ra, rd))

## -----------------------------------------------------------------------------
a <- c(1, 2)
rd <- length(a)
z <- numeric(rd)
ra <- rbind(a, z)

sum <- 0
for (i in 1:400)
  sum <- sum + explg(i, ra, rd)
  
cat("The sum is",round(sum, 5))


## -----------------------------------------------------------------------------
f <- function(a,k){
  g1<-function(u) (1+(u^2/(k-1)))^(-k/2)
  g2<-function(u) (1+u^2/k)^(-(k+1)/2)
  g11 <- integrate(g1, lower = 0, upper = sqrt(a^2*(k-1)/(k-a^2)))$value
  g22 <- integrate(g2, lower=0, upper = sqrt(a^2*k/(k+1-a^2)))$value
  exp(lgamma(k/2)-lgamma((k-1)/2)-log(sqrt(k-1))) * g11 - exp(lgamma((k+1)/2)-lgamma(k/2)-log(sqrt(k))) * g22
}


r <- numeric(25)
i <- 1
for (k in c(4:17)){
  
  g <- function(a) f(a,k)  #求根必须是只关于a的函数
  
  res <- uniroot(g, lower = 0.1, upper = sqrt(k)-0.001)
  r[i] <- res$root
  i <- i + 1
}

i <- 15
for (k in c(18:25,100,500,1000)){
  
  g <- function(a) f(a,k) 
  
  res <- uniroot(g, lower = 0.1, upper = 3)
  r[i] <- res$root
  i <- i + 1
}

r 

## -----------------------------------------------------------------------------
####the observed data MLE
library(stats4)
y <- c(0.54, 0.48, 0.33, 0.43, 1.00, 1.00, 0.91, 1.00, 0.21, 0.85)

mlogL <- function(theta=1) {
        # minus log-likelihood
        return(-((length(y)-3)*log(theta)-theta*sum(y)))
    }

fit <- mle(mlogL)
as.numeric(c(fit@coef,sqrt(fit@vcov)))

## -----------------------------------------------------------------------------
y <- c(0.54, 0.48, 0.33, 0.43, 1.00, 1.00, 0.91, 1.00, 0.21, 0.85)
N <- 10000 #max. number of iterations
L <- 1 #initial est. for lambdas
tol <- .Machine$double.eps^0.5
L.old <- L + 1

repeat { 
  L <- 10/(sum(y)+3/L) #The iterative formula
  if((abs(L - L.old)/L.old) < tol) {
    break
  }
  L.old <- L
}
L

## -----------------------------------------------------------------------------
trims <- c(0, 0.1, 0.2, 0.5)
x <- rcauchy(100)

lapply(trims, function(trim) mean(x, trim = trim))
lapply(trims, mean, x = x)

## -----------------------------------------------------------------------------
####model in 204.3
formulas <- list(
mpg ~ disp,
mpg ~ I(1 / disp),
mpg ~ disp + wt,
mpg ~ I(1 / disp) + wt
)

m1 <- lapply(formulas, lm, data = mtcars)

## -----------------------------------------------------------------------------
###model in 204.4
bootstraps <- lapply(1:10, function(i) {
rows <- sample(1:nrow(mtcars), rep = TRUE)
mtcars[rows, ]
})

m2 <- lapply(bootstraps , lm, formula = mpg ~ disp )

## -----------------------------------------------------------------------------
###204.5
rsq <- function(mod) summary(mod)$r.squared
r1 <- lapply(m1, rsq)
r2 <- lapply(m2, rsq)
r1
r2

## -----------------------------------------------------------------------------
set.seed(12)
###the standard deviation of every column in a numeric data frame
n <- 20
a <- data.frame(rnorm(n,10,20), rexp(n,2), rt(n,6))
round(vapply(a, sd, FUN.VALUE=c(mean=0)), 4)

## -----------------------------------------------------------------------------
set.seed(12)
###the standard deviation of every numeric column in a mixed data frame
n <- 20
b <- data.frame(rnorm(n,10,20), letters[1:n], rep(c(1,2,3),c(4,4,12)))
round(vapply(b[vapply(b, is.numeric, FUN.VALUE=logical(1))], sd, FUN.VALUE=c(mean=0)), 4)
# the inner vapply used to extract numeric column

## -----------------------------------------------------------------------------
library(parallel)
cl.cores <- detectCores() #Check the number of cores currently available on computer
cl <- makeCluster(cl.cores) #Use the kernel parallelism just examined

boot_df <- function(x) x[sample(nrow(x), rep = T), ]
rsquared <- function(mod) summary(mod)$r.square
boot_lm <- function(i) {
dat <- boot_df(mtcars)
rsquared(lm(mpg ~ wt + disp, data = dat))
}

#Load custom functions into the parallel computing environment
clusterExport(cl, "boot_df")  
clusterExport(cl, "rsquared")
clusterExport(cl, "boot_lm")

n <- 1e3
a <- system.time(sapply(1:n, boot_lm))
b <- system.time(parSapply(cl = cl, 1:n, boot_lm))
a
b

## ----eval=FALSE---------------------------------------------------------------
#  #####Rcpp function for 9.8
#  #####myGibbsC.cpp
#  
#  #include <Rcpp.h>
#  using namespace Rcpp;
#  // [[Rcpp::export]]
#  NumericMatrix myGibbsC(int N, int b) {
#    NumericMatrix XC(N, 2);
#    NumericMatrix XCN(N-b, 2);
#    XC(0,0) = 0.5, XC(0,1) = 0.5;
#    for(int i = 1; i < N; i++) {
#      XC(i,0) = rbinom(1, 16, XC(i-1,1))[0];
#      XC(i,1) = rbeta(1, XC(i,0) + 2, 20 - XC(i,0))[0];
#    }
#    for(int j = 0; j < N-b;j++){
#      XCN(j,0) = XC(j+b,0);
#      XCN(j,1) = XC(j+b,1);
#    }
#    return(XCN);
#  }

## -----------------------------------------------------------------------------
library(Rcpp)
dir_cpp <- '../src/'
# Can create source file in Rstudio
sourceCpp(paste0(dir_cpp,'myGibbsC.cpp'))

#####Pure R function for 9.8
#####myGibbsR.R
myGibbsR <- function(N,b){
  XR <- matrix(nrow = N, ncol = 2)
  XRN <- matrix(nrow = N - b, ncol = 2)
  XR[1, ] <- c(0.5, 0.5) #initialize
  for (i in 2:N) {
    XR[i, 1] <- rbinom(1, 16, XR[i-1, 2])
    XR[i, 2] <- rbeta(1, XR[i, 1] + 2, 20 - XR[i, 1])
  }
  XRN <- XR[(b+1):N, ]
}

N <- 5000 #length of chain
b <- 1000 #burn-in length
set.seed(12)

#Compare the random numbers using the function “qqplot”
XR<-myGibbsR(N,b)
XC<-myGibbsC(N,b)
qqplot(XR[,1], XC[,1], xlab = deparse1(substitute(x1)),ylab = deparse1(substitute(X2)))
qqplot(XR[,2], XC[,2], xlab = deparse1(substitute(y1)),ylab = deparse1(substitute(y2)))


#Compare the computation time using function “microbenchmark”
library(microbenchmark)
ts <- microbenchmark(myGibbsR=myGibbsR(N,b),myGibbsC=myGibbsC(N,b))
summary(ts)[,c(1,3,5,6)]

