\name{ClinicalExperiment-class}
\docType{class}
\alias{ClinicalExperiment-class}

\title{Class "ClinicalExperiment"}
\description{A class to encapsulate a clinical experiment used in
  clinical trials}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("ClinicalExperiment", ...)}.
  The arguments include the number of factors in the experiment, the
  names of the factors, the number of levels of each factor, the number
  of treatments and the treatment names. Many, if not all, of the
  arguments are optional, in which case a two-treatment, two-factor
  experiment with two levels of each factor is assumed.
}
\section{Slots}{
  \describe{
    \item{\code{number.of.factors}:}{Object of class \code{"integer"} An
      integer specifying the number of factors, defaults to 2}
    \item{\code{factor.names}:}{Object of class \code{"character"} a
      character vector specifying the character names, defaults to
      \code{F1} and \code{F2}}
    \item{\code{number.of.factor.levels}:}{Object of class
      \code{"integer"} A vector specifying the number of levels of each
      factor, defaults to \code{c(2,2)}}
    \item{\code{factor.level.names}:}{Object of class \code{"list"} a
      list of character vectors specifying the factor level names, defaults to
      \code{1}, \code{2}, etc.}
    \item{\code{number.of.treatments}:}{Object of class \code{"integer"}
      The number of treatments, defaults to 2}
    \item{\code{treatment.names}:}{Object of class \code{"character"}
      The names of the treatments; defaults to \code{c("Tr1", "Tr2")}}
  }
}
\section{Methods}{
  No methods defined with class "ClinicalExperiment" in the signature
  for now. These will be added in future versions.
}
\references{ }
\author{ Balasubramanian Narasimhan }
\note{  }

\seealso{
  \code{\link{ClinicalExperiment}} for a constructor function and
  \code{\linkS4class{PocockSimonRandomizer}} for an example of how this
  class is used. 
}
\examples{
showClass("ClinicalExperiment")

##
## Construct a Clinical Experiment with 3 factors, with levels 2, 2, 3
## respectively, and three treatments with default names for all.
##
expt <- ClinicalExperiment(number.of.factors = 3,
                           number.of.factor.levels = c(2, 2, 3),
                           number.of.treatments = 3)
}
\keyword{classes}
\keyword{design}
