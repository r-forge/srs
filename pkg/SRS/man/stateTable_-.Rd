\name{stateTable<-}
\alias{stateTable<-}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ A function to change the state of a Randomizer }
\description{
  This function allows you to set the slot named \code{stateTable} of an
  object of class \code{PocockSimonRandomizer} to a value.  This allows
  you to start the randomization from arbitrary points in time.
}
\usage{
stateTable(x)<- value
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ An object of class \code{PocockSimonRandomizer} }
  \item{value}{ a matrix of marginal counts of treatment allocations by
    the levels of each factor. It is important to ensure that the
    dimensions match with the experimental setup.  At the moment, no
    checks are done, although that might change in the future.}
}
\details{
}
\value{
}
\references{ }
\author{ Balasubramanian Narasimhan } 
\note{ }
\seealso{ }
\examples{
state.matrix <- matrix (c(9,8,8,9,8,4,5,
                       10,7,6,11,8,5,4,
                       9,7,7,9,8,3,5), nrow=3, byrow=TRUE)

rownames(state.matrix) <- c("Tr1", "Tr2", "Tr3")
colnames(state.matrix) <- c("F1:1", "F1:2", "F2:1", "F2:2", "F3:1", "F3:2", "F3:3")

expt <- ClinicalExperiment(number.of.factors = 3,
                           number.of.factor.levels = c(2, 2, 3),
                           number.of.treatments = 3)

y <- new("PocockSimonRandomizer", expt, as.integer(12345))

y@stateTable <- state.matrix

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ design }
\keyword{ distribution }% __ONLY ONE__ keyword per line
