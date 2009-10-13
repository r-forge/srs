\name{randomize}
\alias{randomize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ A function to randomize a new subject }
\description{
  A function to randomize a new subject to a treatment given factor variables
}
\usage{
randomize(object, subject.id, factor.values)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ An object of class \code{PocockSimonRandomizer}}
  \item{subject.id}{ A string identifying the subject; no check is made
    on uniqueness }
  \item{factor.values}{ A vector of factor values for the subject }
}
\details{
}
\value{
  Returns an object of \code{\linkS4class{PocockSimonRandomizer}} class
  that contains the orginal object with updated slots
}
\references{ }
\author{ Balasubramanian Narasimhan }
\note{ }

\seealso{ \code{\link{lastRandomization}}}
\examples{

expt <- ClinicalExperiment(number.of.factors = 3,
                           number.of.factor.levels = c(2, 2, 3),
                           number.of.treatments = 3)
r.obj <- new("PocockSimonRandomizer", expt, as.integer(12345))
r.obj <- randomize(r.obj, "Subject-001", c("1","2","2"))
## Try once again with same subject and you should get an error
## r.obj <- randomize(r.obj, "Subject-001", c("1","2","2"))
## Or try with a bad factor value and you should get an error
## r.obj <- randomize(r.obj, "Subject-002", c("1","2","5"))
lastRandomization(r.obj)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ design }
\keyword{ distribution }% __ONLY ONE__ keyword per line
