\name{fac}
\alias{fac}
\title{
Selecting the optimal number of latent variables.
}
\description{
Select the optimal number of latent variables by varying the number of latent variables from fmin to fmax.
}
\usage{
fac(Y, X, minq = 1, maxq = 2, epsilon = 0.01, plot.BIC = FALSE, printout = TRUE)
}
\arguments{
  \item{Y}{
Expression data in samples by features format.
}
  \item{X}{
Design matrix from model.matrix function which defines how covariates are included into sBFAC.
}
  \item{minq}{
Minimum number of latent variables to be considered. Default is 1.
}
  \item{maxq}{
Maximum number of latent variables to be considered. Default is 2. fmax should be greater or equal to fmin.
}
  \item{epsilon}{
Assessing convergence.
}
  \item{plot.BIC}{
Plot BIC values for selecting the optimal number of latent variables.
}
  \item{printout}{
Print when the model have converged.
}
}
\details{
 Select the optimal number of latent variables by varying the number of latent variables from fmin to fmax.
}
\value{
\item{qopt }{The optimal number of latent variables selected by BIC.}
\item{sig }{Residual variability associated with each feature.}
\item{sig }{Residual variability associated with each latent variable.}
\item{factors }{Scores of observations on each latent variable.}
\item{loadings }{Loadings of features on each latent variable.}
\item{coefficients }{Regression coefficient for each covariate.}
}
\references{
Nyamundanda et al (2015). A next generation tool enables identification of functional cancer subtypes with associated biological phenotypes.
}
\author{
Gift Nyamundanda, Katie Eason, Pawan Poudel and Anguraj Sadanandam.
}
\examples{
## help(fac) 
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
