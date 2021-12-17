#' @title generate the desired simulation data
#' @description Use network matrix, node-specific random vector and others to generate the simulation data
#' @param N the number of nodes
#' @param alpha the parameter of power-law distribution
#' @param beta0 Initial value of intercept item
#' @param Beta Initial value of network effect
#' @param gamma0 Initial value of nodal effect
#' @param Time The observation time
#' @param sig Variance of the error term
#' @return a simulated dataset of N*(p+1)
#' @examples
#' \dontrun{
#' N=50  
#' alpha=1.2  
#' Time=30   
#' gamma0=c(-0.5,0.3,0.8,0,0) 
#' beta0=rep(0.3,N) 
#' Beta=c(0.1,0.5)  
#' Ysim=simudata(beta0, Beta, gamma0, Time, sig = 1, N, alpha)
#' }
#' @import plyr
#' @import Matrix
#' @importFrom MASS mvrnorm
#' @importFrom poweRlaw rpldis
#' @importFrom stats rnorm
#' @export
simudata <- function(beta0, Beta, gamma0, Time, sig = 1, N, alpha) {
  W = W(N, alpha)
  p = length(gamma0)
  Z = Z(N, p)
  G = Beta[1]*W + Beta[2]*diag(nrow(W))
  mu = statMu(G, beta0 + Z %*% gamma0)
  Gamma0 = statGamma(G, sigma = sig)
  Y0 = mvrnorm(n = 1, mu = mu, Sigma = Gamma0)
  Ysim = Y(Y0 = Y0, beta0 + Z %*% gamma0, Beta, W, Time = Time)
  return(Ysim)
}




#' @title calculate estimators
#' @description Use the simulation data to calculate estimators
#' @param N the number of nodes
#' @param alpha the parameter of power-law distribution
#' @param Ysim the simulation data
#' @param p the dimension of node-specific random vector Z
#' @return an estimator of coefficients
#' @examples
#' \dontrun{
#' N=50  
#' alpha=1.2  
#' Time=30   
#' gamma0=c(-0.5,0.3,0.8,0,0) 
#' beta0=rep(0.3,N) 
#' Beta=c(0.1,0.5)  
#' Ysim=simudata(beta0, Beta, gamma0, Time, sig = 1, N, alpha)
#' est(Ysim,N,alpha,p)
#' }
#' @import plyr
#' @import Matrix
#' @importFrom MASS mvrnorm
#' @importFrom dplyr bind_cols
#' @importFrom poweRlaw rpldis
#' @export
est <- function(Ysim, N, alpha, p){
  W = W(N, alpha)
  Z = Z(N, p)
  WY = W %*% Ysim  
  Time = ncol(Ysim) - 1  
  
  int = rep(1, nrow(Ysim) * Time)
  self = as.vector(WY[, -ncol(Ysim)])
  other = as.vector(Ysim[, -ncol(Ysim)])
  realz = matrix(0, nrow = N*Time, ncol = p)
  for (i in seq(1, N*Time,N)){
    realz[i:(i+N-1),] = Z
  }
  
  X = as.matrix(bind_cols(int, self, other, realz))
  invXX = solve(crossprod(X)) 
  Y = as.vector(Ysim[, -1])  
  
  thetahat = invXX %*% t(X) %*% Y 
  return(list(theta = thetahat))
}


W <- function(N, alpha) {
  Numfollowers = rpldis(N, 1, alpha)  
  A = sapply(Numfollowers, function(n) {
    
    sel = rep(0, N)
    sel[sample(1:N, min(n, N))] = 1
    return(sel)
  })
  diag(A) = 0
  for (i in 1:N){
    if (sum(A[i,]) == 0)         
      A[i,sample(1:N,3)] = 1   
  }
  W = A / rowSums(A) 
  return(W)
} 


Z <- function(N, p, s = 0.5) {
  cov = diag(p)   
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      cov[i,j] = s^abs(i-j)
    }
  }
  cov = t(cov)+cov
  diag(cov)=1
  Z = mvrnorm(N, mu = rep(0,p), Sigma = cov) 
  return(Z)
}


statMu <- function(G, Beta0) {
  mu = solve(diag(length(Beta0)) - G) %*% Beta0  
  return(mu)
}
statGamma <- function(G, sigma) {
  Gamma0 = (diag(nrow(G)) + tcrossprod(G) + tcrossprod(G %*% G) + tcrossprod(G %*% G %*% G)) * sigma^2
  return(Gamma0)
} 


Y <- function(Y0, Beta0, Beta, W, sig = 1, Time) {
  Ysim = matrix(0, nrow = length(Beta0), ncol = Time + 1)  
  Ysim[, 1] = as.vector(Y0) 
  for (i in 1:Time) {
    Ysim[, i + 1] = as.vector(Beta0 + Beta[1] * W %*% Ysim[, i] + Beta[2] * 
                                Ysim[, i] + rnorm(length(Beta0), sd = sig))  
  }
  return(Ysim)
}
