\name{PRS-package}
\alias{PRS-package}
\alias{PRS}
\docType{package}
\title{
Randomized patients based on Pocock-Simon minimization method
}
\description{
This package allows one to randomize subjects based on prognostic
factors using the Pocock Simon (1975) minimization method
}
\details{
\tabular{ll}{
Package: \tab PRS\cr
Type: \tab Package\cr
Version: \tab 0.01\cr
Date: \tab 2007-08-20\cr
License: \tab GPL (version 2 or later)\cr
}
The package is used as follows. First create a \code{ClinicalExperiment}
object using the function \code{ClinicalExperiment} class and then a
randomizer using a \code{PocockSimonRandomizer} class. The function
\code{randomize} will randomize subjects to treatments.  Read the
vignette for a complete example.
}
\author{
Balasubramanian Narasimhan

Maintainer: Balasubramanian Narasimhan<naras@stanford.edu>
}
\references{
  Sequential Treatment Assigment with Balancing for Prognostic
  Factors in the Controlled Clinical Trial, by S.~J.~Pocock and
  R.~Simon, Biometrics, 31, 103-115
}
\keyword{ package }
\keyword{ design }
\keyword{ distribution }
\seealso{
}
\examples{
##
## Please refer to the vignette
##
}
