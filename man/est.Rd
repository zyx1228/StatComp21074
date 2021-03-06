% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simudata.R
\name{est}
\alias{est}
\title{calculate estimators}
\usage{
est(Ysim, N, alpha, p)
}
\arguments{
\item{Ysim}{the simulation data}

\item{N}{the number of nodes}

\item{alpha}{the parameter of power-law distribution}

\item{p}{the dimension of node-specific random vector Z}
}
\value{
an estimator of coefficients
}
\description{
Use the simulation data to calculate estimators
}
\examples{
\dontrun{
N=50  #numbers of nodes
alpha=1.2  #parameter of power-law
Time=30   #The length of time
gamma0=c(-0.5,0.3,0.8,0,0) #fixed parameter of Z
beta0=rep(0.3,N) #Intercept term coefficient
Beta=c(0.1,0.5)  #Regression coefficient
Ysim=simudata(beta0, Beta, gamma0, Time, sig = 1, N, alpha)
est(Ysim,N,alpha,p)
}
}
