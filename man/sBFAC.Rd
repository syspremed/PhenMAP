\name{sBFAC}
\alias{sBFAC}
\title{
The sparse Bayesian FAC model.
}
\description{
A statistical tool for clustering expression data, finding associated phenotypes and molecular features.
}
\usage{
sBFAC(Y, X, fmin = 1, fmax = 2, gmin = 1, gmax = 2, w0 = 1, g0 = nrow(Y), mcmc = 20000)
}
\arguments{
  \item{Y}{
Expression data in samples by features format.
}
  \item{X}{
Design matrix from model.matrix function which defines how covariates are included into sBFAC.
}
  \item{fmin}{
Minimum number of latent variables to be considered. Default is 1.
}
  \item{fmax}{
Maximum number of latent variables to be considered. Default is 2. fmax should be greater or equal to fmin.
}
  \item{gmin}{
Minimum number of subtypes to be considered. Default is 1.
}
  \item{gmax}{
Maximum number of subtypes to be considered. Default is 2.
}
  \item{w0}{
Controls sparseness in features, minimum value is 0.01 for non-informative priors.
}
  \item{g0}{
Decreasing this regularize regression coefficients, i.e., Very large value of g0 provides non-informative prior == Classical estimates.
}
  \item{mcmc}{
Number of MCMC samples. Default is 20000 but can be increased to improve convergence.
}
}
\details{
 The function first select the optimal number of latent variables and fit sBFAC model via MCMC to identify subtypes and their associated
 phenotypes and molecular features. 
}
\value{
\item{u.store }{This is an array of q x N x mcmc containing scores for the q latent variables on N obaservations over mcmc samples.}
\item{w.store }{This is an array of p x q x mcmc containing loadings for the p genes/features on q latent variables over mcmc samples.}
\item{b.store }{This is an array of q x L x mcmc containing regression coefficients of L covariates on q latent variables over mcmc samples.}
\item{s.store }{This is an matrix of p x mcmc containing residual variablity of each feature.}
\item{s2.store }{This is an matrix of q x mcmc containing residual variablity of each latent variable.}
\item{MBClustering }{Model based clustering of q latent variables.}
\item{KMeansClustering }{KMeans clustering of q latent variables.}
}
\references{
Nyamundanda et al (2016). A next generation tool enables identification of functional cancer subtypes with associated biological phenotypes.
}
\author{
Gift Nyamundanda, Katie Eason, Pawan Poudel and Anguraj Sadanandam.
}

\examples{
## help(sBFAC)
## Expression and phenotype data are required for this function 
data(Y)
data(Cv)
## Create design matrix for the phenotye data
X<-model.matrix(~Etoposide+Fascaplysin+Bortezomib+Geldanamycin,Cv)  
## Run the function
sBFAC(Y,X,fmin=1,fmax=3,gmin=1,gmax=7,w0=1,g0=1000,mcmc=10000)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
