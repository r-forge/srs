\name{PocockSimonRandomizer-class}
\docType{class}
\alias{PocockSimonRandomizer-class}
\alias{computeImbalances,PocockSimonRandomizer,integer-method}
\alias{computeOverallImbalance,PocockSimonRandomizer,matrix-method}
\alias{initialize,PocockSimonRandomizer-method}
\alias{lastRandomization,PocockSimonRandomizer-method}
\alias{randomize,PocockSimonRandomizer,ANY,integer-method}
\alias{stateTable<-,PocockSimonRandomizer-method}
\alias{tr.assignments<-,PocockSimonRandomizer-method}

\title{Class "PocockSimonRandomizer"}

\description{ An encapsulation of the minimization randomization
  detailed by Pocock and Simon (1975) }

\section{Objects from the Class}{
  Objects can be created by calls of the form
  \code{new("PocockSimonRandomizer", expt, seed, stateTable, tr.ratios,
    d.func, g.func, p.func)}.  Arguments include a ClinicalExperiment
  object, a random number seed, (optionally) a state table of marginal
  counts if the randomization is to begin at a particular point in time,
  the ratio of treatment counts, a function to use for computing
  imbalance for each treatment, a function to compute the overall
  imbalance, and a function that provides the probability allocation.
}

\section{Slots}{
  \describe{
    \item{\code{expt}:}{Object of class \code{"ClinicalExperiment"}
      Defines the clinical experiment context for this randomizer}
    \item{\code{seed}:}{Object of class \code{"integer"} An integer used
    for seeding the random number generator for reproducibility}
    \item{\code{stateTable}:}{Object of class \code{"matrix"} A matrix
      of counts indicating the marginal distribution of each factor
      level per treatment; see example below }
    \item{\code{tr.assignments}:}{Object of class \code{"data.frame"}
      The treatment assignments so far}
    \item{\code{tr.ratios}:}{Object of class \code{"integer"} The ratio
      of treatment counts, example 2:1 in a two treatment experiment}
    \item{\code{d.func}:}{Object of class \code{"function"} A function
      that computes the imbalance; see example below}
    \item{\code{g.func}:}{Object of class \code{"function"} A function
      that computes the overall imbalance to be minimized; see example below}
    \item{\code{p.func}:}{Object of class \code{"function"} A function
      that computes the probability vector of treatment assignments; see
    example below}
  }
}

\section{Methods}{
  \describe{
    \item{computeImbalances}{\code{signature(object =
	"PocockSimonRandomizer", factor.values = "list")}: Given a
      set of factor values associated with a subject, compute imbalances
      that occur if each of the treatments in turn is assigned to the subject } 
    \item{computeOverallImbalance}{\code{signature(object =
	"PocockSimonRandomizer", imbalances = "matrix")}: Given a vector
      of imbalances resulting from assigning each of the treatments in
      turn, compute the overall imbalance} 
    \item{initialize}{\code{signature(.Object =
	"PocockSimonRandomizer")}: Create an instance of this object }
    \item{lastRandomization}{\code{signature(object =
	"PocockSimonRandomizer")}: Return the details of the last
      randomization that was done}
    \item{randomize}{\code{signature(object = "PocockSimonRandomizer",
	subject.id = "character", factor.values = "character")}: Given a subject
      id and a set of factor values, randomize the subject to one of the treatments}
    \item{stateTable<-}{\code{signature(x = "PocockSimonRandomizer")}:
      Set the \code{stateTable} slot}
    \item{tr.assignments<-}{\code{signature(x =
	"PocockSimonRandomizer")}: Set the \code{}tr.assignments} slot}
}

\references{ This implementation is based directly on the paper
  \emph{Sequential Treatment Assigment with Balancing for Prognostic
    Factors in the Controlled Clinical Trial}, by S.~J.~Pocock and
  R.~Simon, Biometrics, 31, 103-115
}

\author{ Balasubramanian Narasimhan }

\note{  }

\seealso{
  \code{\link{ClinicalExperiment}},
  \code{\linkS4class{ClinicalExperiment}} class
}

\examples{
showClass("PocockSimonRandomizer")
##
## Create a simple PocockSimonRandomizer class
##
expt <- ClinicalExperiment(number.of.factors = 2,
                           number.of.factor.levels = c(2, 2))
randomizer <- new("PocockSimonRandomizer", expt, as.integer(12345))
randomizer <- randomize(randomizer, "Subject 1", c("1", "2"))
randomizer <- randomize(randomizer, "Subject 2", c("2", "1"))
randomizer

}


\keyword{classes}
