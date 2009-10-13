\name{lastRandomization}
\alias{lastRandomization}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ A function to return details of the last randomized subject }
\description{
  A function to return details of the last randomized subject
}
\usage{
lastRandomization(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ An object of class \code{PocockSimonRandomizer} }
}
\details{
}
\value{
  A one row data frame with named variables.
}
\references{ }
\author{ Balasubramanian Narasimhan }
\note{ } 
\seealso{ }
\examples{
expt <- ClinicalExperiment(number.of.factors = 3,
                           number.of.factor.levels = c(2, 2, 3),
                           number.of.treatments = 3)
r.obj <- new("PocockSimonRandomizer", expt, as.integer(12345))
r.obj <- randomize(r.obj, "Subject-001", c("1","2","2"))
lastRandomization(r.obj)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ design }
\keyword{ distribution }% __ONLY ONE__ keyword per line
