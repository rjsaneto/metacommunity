% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meanLocal.R
\name{meanLocal}
\alias{meanLocal}
\title{Calculate the relative species abundance of the local communities.}
\usage{
meanLocal(df, meanRSA = TRUE)
}
\arguments{
\item{df}{A dataframe containing the results of the simulations.}

\item{meanRSA}{If set to TRUE, the RSA of the simulations will be used for the mean; if set to FALSE, the RSA will be calculated from the mean abundances of all simulations. Generally, the first method may result in values different from 1 when the RSA is summed, whereas in the second method, this does not occur.}
}
\value{
The meanLocal function returns an object containing the species rank, RSA (Relative Species Abundance), RSA standard deviation, and the simulation's sigma value.
}
\description{
Function to extract the arithmetic mean of the RSA from all simulations performed by the modelMetacomm function.
}
\examples{
mymodel <- modelMetacomm(J=100,theta = 10,sigma=0.03,Jl=10,nloc=3,mig=0.25,loccomm="diversity",abonly=F,nSim=2)
results <- meanLocal(mymodel, meanRSA = TRUE) #The final RSA will be obtained by the means of the RSAs of all simulations.
results2 <- meanLocal(mymodel, meanRSA = FALSE) #The final RSA will be obtained by the  mean of the abundances from all simulations.
}
