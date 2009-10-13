\name{computeImbalances}
\alias{computeImbalances}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ A function to compute imbalances resulting from assigning each
  of the treatments in turn}
\description{
  This function computes the imbalances resulting from assigning each
  of the treatments in turn for a prospective subject coming in with the
  factors specified
}
\usage{
computeImbalances(object, factor.values)
}

\arguments{
  %- maybe also 'usage' for other objects documented here.  \arguments{
  \item{object}{An object of class \code{PocockSimonRandomizer} }
  \item{factor.values}{A vector of character values specifying the level of each
    factor }
}

\details{ }

\value{ The result is an m by n matrix of values, where m is the number
  of factors and n is the number of treatments. The entries are what the
  counts of subjects assigned to each treatment would be for the levels of
  the factors specified. Suppose the factor levels are first level of
  factor 1, \code{F1:1}, third level of factor 2, \code{F2:3}, and second
  level of factor 3, \code{F3:2} and there are four treatments \code{Tr1},
  \code{Tr2}, \code{Tr3}, and \code{Tr4}.  Then the 3 by 4 matrix has row
  names \code{F1-1, F2-3, F3-2} and column names \code{Tr1, Tr2, Tr3,
    Tr4}. The first row would indicate the range of imbalance of subject
  counts (assuming default imbalance function of range) for each of the
  four treatments.
}
\references{ Sequential Treatment Assigment with Balancing for Prognostic
  Factors in the Controlled Clinical Trial, by S.~J.~Pocock and
  R.~Simon, Biometrics, 31, 103-115}
\author{ Balasubramanian Narasimhan }
\note{ }

\seealso{ \code{\link{computeOverallImbalance}} }
\examples{
expt <- ClinicalExperiment(number.of.factors = 3,
                           number.of.factor.levels = c(2, 2, 3),
                           number.of.treatments = 3)
r.obj <- new("PocockSimonRandomizer", expt, as.integer(12345))
computeImbalances(r.obj, c("1","2","2"))
##
## Another example
##
ex.matrix <- matrix (c(9,8,8,9,8,4,5,
                       10,7,6,11,8,5,4,
                       9,7,7,9,8,3,5), nrow=3, byrow=TRUE)

rownames(ex.matrix) <- c("Tr1", "Tr2", "Tr3")
colnames(ex.matrix) <- c("F1:1", "F1:2", "F2:1", "F2:2", "F3:1", "F3:2", "F3:3")
stateTable(r.obj) <- ex.matrix
computeImbalances(r.obj, c("1","2","2"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ design }
\keyword{ distribution }% __ONLY ONE__ keyword per line
