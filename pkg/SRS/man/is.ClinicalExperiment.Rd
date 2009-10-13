\name{is.ClinicalExperiment}
\alias{is.ClinicalExperiment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Checks if an object a ClinicalExperiment object }
\description{
  Checks if an object a ClinicalExperiment object
}

\usage{
is.ClinicalExperiment(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ any object }
}
\details{
 }
\value{
  Returns true or false
}
\references{ }
\author{ Balasubramanian Narasimhan }
\note{ }

\seealso{ \code{\link{ClinicalExperiment}} }
\examples{
expt <- ClinicalExperiment(number.of.factors = 3,
                           number.of.factor.levels = c(2, 2, 3),
                           number.of.treatments = 3)
is.ClinicalExperiment(expt)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ design }
\keyword{ distribution }% __ONLY ONE__ keyword per line
