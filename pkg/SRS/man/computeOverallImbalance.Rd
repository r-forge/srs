\name{computeOverallImbalance}
\alias{computeOverallImbalance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ A function to compute overall imbalances for each treatment }
\description{
  This function computes the overall treatement imbalances for each treatment 
}
\usage{
computeOverallImbalance(object, imbalances)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ An object of class \code{PocockSimonRandomizer}}
  \item{imbalances}{ \code{imbalances} a matrix of imbalance measures,
    usually the result of the function call to
    \code{\link{computeImbalances}}}
}

\details{
  
}
\value{
  A vector of imbalances resulting from assigning each of the
  treatments, in turn.
}
\references{ Sequential Treatment Assigment with Balancing for Prognostic
  Factors in the Controlled Clinical Trial, by S.~J.~Pocock and
  R.~Simon, Biometrics, 31, 103-115}
\author{ Balasubramanian Narasimhan }
\note{ }
\seealso{ }
\examples{
expt <- ClinicalExperiment(number.of.factors = 3,
                           number.of.factor.levels = c(2, 2, 3),
                           number.of.treatments = 3)
r.obj <- new("PocockSimonRandomizer", expt, as.integer(12345))
computeOverallImbalance(r.obj, computeImbalances(r.obj, c("1","2","2")))
##
## Another example
##
ex.matrix <- matrix (c(9,8,8,9,8,4,5,
                       10,7,6,11,8,5,4,
                       9,7,7,9,8,3,5), nrow=3, byrow=TRUE)

rownames(ex.matrix) <- c("Tr1", "Tr2", "Tr3")
colnames(ex.matrix) <- c("F1:1", "F1:2", "F2:1", "F2:2", "F3:1", "F3:2", "F3:3")
stateTable(r.obj) <- ex.matrix
computeOverallImbalance(r.obj, computeImbalances(r.obj, c("1","2","2")))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ design }
\keyword{ distribution }% __ONLY ONE__ keyword per line
