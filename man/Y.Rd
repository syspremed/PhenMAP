\name{Gene Expression data}
\alias{Expression data}
\docType{data}
\title{
Gene expression data for thirty six breast cancer cell lines.
}
\description{
Gene expression data for thirty six breast cancer cell lines.
}
\usage{data("Y")}
\format{
  The format is:
 num [1:36, 1:994] 7.71 6.68 7.56 8.24 5.79 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:36] "AU565" "BT474" "BT483" "BT549" ...
  ..$ : chr [1:994] "TSPAN6" "CFH" "TMEM176A" "ICA1" ...
}
\details{
 Gene expression data for thirty six breast cancer cell lines.
}
\source{
A next generation tool enables identification of functional cancer subtypes with associated biological phenotypes.
}
\references{
Nyamundanda et al (2015). A next generation tool enables identification of functional cancer subtypes with associated biological phenotypes.
}
\examples{
data(Y)
## maybe str(Y) ; plot(Y) ...
}
\keyword{datasets}
