## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  #### Our function1
#  simudata <- function(beta0, Beta, gamma0, Time, sig = 1, N, alpha) {
#    W = W(N, alpha)
#    p=length(gamma0)
#    Z = Z(N, p)
#    G = Beta[1]*W + Beta[2]*diag(nrow(W))
#    mu = statMu(G, beta0 + Z %*% gamma0)
#    Gamma0 = statGamma(G, sigma = sig)
#    Y0 = mvrnorm(n = 1, mu = mu, Sigma = Gamma0)
#    Ysim = Y(Y0 = Y0, beta0 + Z %*% gamma0, Beta, W, Time = Time)
#    return(Ysim)
#  }

## ----eval=FALSE---------------------------------------------------------------
#  W <- function(N, alpha) {
#    Numfollowers = rpldis(N, 1, alpha)  # generate random numbers for each nodes following power-law(1, alpha)
#    A = sapply(Numfollowers, function(n) {
#      # for node i, randomly select di nodes to follow it
#      sel = rep(0, N)
#      sel[sample(1:N, min(n, N))] = 1
#      return(sel)
#    })
#    diag(A) = 0
#    for (i in 1:N){
#      if (sum(A[i,]) == 0)         # in case some row sums are zero
#        A[i,sample(1:N,3)] = 1   # for those node, randomly select 3 followees
#    }
#    W = A / rowSums(A) # row-normalized adjacency matrix
#    return(W)
#  }

## ----eval=FALSE---------------------------------------------------------------
#  Z <- function(N, p, s = 0.5) {
#    cov = diag(p)  # Set covariance matrix with the element of 0.5 ^ | i - j |
#    for (i in 1:(p-1)){
#      for (j in (i+1):p){
#        cov[i,j] = s^abs(i-j)
#      }
#    }
#    cov = t(cov)+cov
#    diag(cov)=1
#    Z = mvrnorm(N, mu = rep(0,p), Sigma = cov)
#    return(Z)
#  }

## ----eval=FALSE---------------------------------------------------------------
#  #using proposition1 and Case3 to obtain the mean and covariance matrix of stationary time series
#  statMu <- function(G, Beta0) {
#    mu = solve(diag(length(Beta0)) - G) %*% Beta0
#    return(mu)
#  }
#  statGamma <- function(G, sigma) {
#    Gamma0 = (diag(nrow(G)) + tcrossprod(G) + tcrossprod(G %*% G) + tcrossprod(G %*% G %*% G)) * sigma^2
#    # Tayler expansion in case3
#    return(Gamma0)
#  }

## ----eval=FALSE---------------------------------------------------------------
#  Y <- function(Y0, Beta0, Beta, W, sig = 1, Time) {
#    Ysim = matrix(0, nrow = length(Beta0), ncol = Time + 1)  # use Ysim to store the simulated data
#    Ysim[, 1] = as.vector(Y0)  # the first column of Ysim holds Y0
#    for (i in 1:Time) {
#      Ysim[, i + 1] = as.vector(Beta0 + Beta[1] * W %*% Ysim[, i] + Beta[2] *
#                                  Ysim[, i] + rnorm(length(Beta0), sd = sig))  # simulate
#    }
#    return(Ysim)
#  }

## ----eval=FALSE---------------------------------------------------------------
#  # The function requires the user to obtain simulated data
#  # and then input it as a parameter to the est()
#  est <- function(Ysim, N, alpha, p){
#    W = W(N, alpha)
#  
#    Z = Z(N, p)
#    WY = W %*% Ysim  # WY
#    Time = ncol(Ysim) - 1  # cause Ysim included Y0
#  
#    # Joining together the X
#    int = rep(1, nrow(Ysim) * Time)
#    self = as.vector(WY[, -ncol(Ysim)])
#    other = as.vector(Ysim[, -ncol(Ysim)])
#    realz = matrix(0, nrow = N*Time, ncol = p)
#    for (i in seq(1, N*Time,N)){
#      realz[i:(i+N-1),] = Z
#    }
#  
#    X = as.matrix(bind_cols(int, self, other, realz))
#    invXX = solve(crossprod(X))
#    Y = as.vector(Ysim[, -1])
#  
#    thetahat = invXX %*% t(X) %*% Y  # estimation equation (2.8)
#    return(list(theta = thetahat))
#  }

## -----------------------------------------------------------------------------
library(dplyr)
library(plyr)
library(Matrix)
library(MASS)
library(poweRlaw)
library(StatComp21074)

N = 50  #numbers of nodes
alpha = 1.2  #parameter of power-law
Time = 30   #The length of time
gamma0 = c(-0.5,0.3,0.8,0,0) #fixed parameter of Z
p = length(gamma0)
beta0 = rep(0.3,N) #Intercept term coefficient
Beta = c(0.1,0.5)  #Regression coefficient
Ysim = simudata(beta0, Beta, gamma0, Time, sig = 1, N, alpha)
est(Ysim,N,alpha,p)

