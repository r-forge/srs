\name{ClinicalExperiment}
\alias{ClinicalExperiment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ClinicalExperiment constructor function}
\description{
  This is a constructor function for the \code{ClinicalExperiment}
}
\usage{
ClinicalExperiment(number.of.factors = 2, factor.names,
number.of.factor.levels, factor.level.names, number.of.treatments=2, treatment.names)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{number.of.factors}{the number of factors, defaults to 2 }
  \item{factor.names}{ the names of the factors; default names are
    constructed if not specified }
  \item{number.of.factor.levels}{ the number of levels of each factor;
    defaults to 2 for each}
  \item{factor.level.names}{ the names of the factor levels; default names are
    constructed if not specified }
  \item{number.of.treatments}{ the number of treatments; defaults to 2 }
  \item{treatment.names}{ the names of the treatments; defaults are
    constructed if not specified}
}
\details{
  This function merely constructs a \code{ClinicalExperiment} object for
  use later on with a randomizer.
}
\value{
  An object of class \code{ClinicalExperiment}
}
\references{ }
\author{ Balasubramanian Narasimhan}
\note{ }

\seealso{  }
\examples{
expt <- ClinicalExperiment(number.of.factors = 3,
                           number.of.factor.levels = c(2, 2, 3),
                           number.of.treatments = 3)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ design }
\keyword{ distribution }% __ONLY ONE__ keyword per line
