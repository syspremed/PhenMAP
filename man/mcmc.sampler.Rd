\name{mcmc.sampler}
\alias{mcmc.sampler}
\title{
Function for fitting sBFAC model via MCMC methods.
}
\description{
Function for fitting sBFAC model via MCMC methods.
}
\usage{
mcmc.sampler(facfit, Y, X, w0, g0, mcmc, thinstep)
}
\arguments{
  \item{facfit}{
MAP estimates of sBFAC used to initialize MCMC algorithm.}
  \item{Y}{
Expression data in samples by features format
}
  \item{X}{
Design matrix from model.matrix function which defines how covariates are included into sBFAC.
}
  \item{w0}{
Controls sparseness in features, minimum value should be 0.01 for non-informative priors.
}
  \item{g0}{
Decreasing this regularize regression coefficients, i.e., g0=1000 provides non-informative prior == Classical estimates
}
  \item{mcmc}{
Number of MCMC samples. Default is 20000 but can be increased to improve convergence.
}
  \item{thinstep}{
Number of iterations to thin to remove autocorrelation between samples.
}
}
\details{
Function for fitting sBFAC model via MCMC methods.
}
\value{
\item{scores.store }{This is an array of qxNxmcmc containing scores for the q latent variables on N obaservations over mcmc samples.}
\item{loadings.store }{This is an array of pxqxmcmc containing loadings for the p genes/features on q latent variables over mcmc samples.}
\item{beta.store }{This is an array of qxLxmcmc containing regression coefficients of L covariates on q latent variables over mcmc samples.}
\item{sigma.store }{This is an matrix of pxmcmc containing residual variablity of each feature.}
\item{s2.store }{This is an matrix of qxmcmc containing residual variablity of each latent variable.}
}
\references{
Nyamundanda et al (2015). A next generation tool enables identification of functional cancer subtypes with associated biological phenotypes.
}
\author{
Gift Nyamundanda, Katie Eason, Pawan Poudel and Anguraj Sadanandam.
}
\examples{
## help(mcmc.sampler)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line


