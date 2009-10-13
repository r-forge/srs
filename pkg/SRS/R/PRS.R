##
## Define Clinical Experiment
##
setClass("ClinicalExperiment",
         representation(number.of.factors="integer",
                        factor.names="character",
                        factor.level.names="list",
                        number.of.factor.levels="integer",
                        number.of.treatments="integer",
                        treatment.names="character"),
         prototype(number.of.factors=as.integer(2),
                   factor.names=c("F1", "F2"),
                   factor.level.names=list(),
                   number.of.factor.levels=as.integer(c(2,2)),
                   number.of.treatments=as.integer(2),
                   treatment.names=c("Tr1", "Tr2")))

##
## Validate a clinical experiment object
##
.validClinicalExperiment <- function(object) {
    retval <- NULL
    if (length(object@number.of.factors) != 1) {
        retval <- c(retval, "number.of.factors is not an integer")
    }
    
    if (object@number.of.factors <= 0) {
        retval <- c(retval, "number.of.factors is not positive")
    }
    
    if (length(object@factor.names) != object@number.of.factors) {
        retval <- c(retval, "factor.names length mismatch with number.of.factors")
    }
    
    if (any(sapply(object@factor.names, function(x) nchar(x) == 0))) {
        retval <- c(retval, "some factor.names are empty")
    }
    
    if (length(object@factor.level.names) != object@number.of.factors) {
        retval <- c(retval, "factor.level.names length mismatch with number.of.factors")
    }
    
    if (length(object@number.of.factor.levels) != object@number.of.factors) {
        retval <- c(retval, "number.of.factor.levels length mismatch with number.of.factors")
    }
    
    for (i in 1:object@number.of.factors) {
      if (length(object@factor.level.names[[i]]) != object@number.of.factor.levels[i]) {
        retval <- c(retval, "factor.level.names mismatch with number.of.factors.levels")        
      }
    }

    if (length(object@number.of.treatments) != 1) {
        retval <- c(retval, "number.of.treatments is not a single integer")
    }
    
    if (object@number.of.treatments < 2) {
        retval <- c(retval, "number.of.treatments must be at least 2")
    }
    
    if (length(object@treatment.names) != object@number.of.treatments) {
        retval <- c(retval, "treatment names length mismatch with number.of.treatment")
    }
    
    if (is.null(retval)) {
        return(TRUE)
    } else {
        return(retval)
    }
}

setValidity("ClinicalExperiment", .validClinicalExperiment)


##
## Constructor function for ClinicalExperiement
##
"ClinicalExperiment" <- function(number.of.factors=2,
                                 factor.names,
                                 number.of.factor.levels,
                                 factor.level.names,                                 
                                 number.of.treatments=2,
                                 treatment.names) {
  
  if (length(number.of.factors) != 1) {
    stop("number.of.factors is not a single integer")
  } else if (number.of.factors < 2) {
    stop("number.of.factors must be at least 2!")
  }
  
  if (missing(factor.names)) {
    factor.names <- sapply(1:number.of.factors,
                           function(x) paste("F", x, sep=""))
  }
  
  if (missing(number.of.factor.levels)) {
    number.of.factor.levels <- rep(2, number.of.factors)
  }
  
  if (missing(treatment.names)) {
    treatment.names <- sapply(1:number.of.treatments,
                              function(x) paste("Tr", x, sep=""))
  }

  if (missing(factor.level.names)) {
    factor.level.names <- lapply(number.of.factor.levels,
                                 function(x) as.character(1:x))
  }

  new("ClinicalExperiment",
      number.of.factors = as.integer(number.of.factors),
      factor.names = factor.names,
      factor.level.names = factor.level.names,
      number.of.factor.levels = as.integer(number.of.factor.levels),
      number.of.treatments = as.integer(number.of.treatments),
      treatment.names = treatment.names)
}
        
is.ClinicalExperiment <- function(x) is(x, "ClinicalExperiment")

    
##
## PocockSimonRandomizer class
##
setClass("PocockSimonRandomizer",
         representation(expt = "ClinicalExperiment",
                        seed = "integer",
                        stateTable = "matrix",
                        tr.assignments = "data.frame",
                        tr.ratios = "numeric",
                        d.func = "function",
                        g.func = "function",
                        p.func = "function")
         )

##
## Helper default functions. Internal
##

##
## Default measure of imbalance is range, as in paper.
## Result is vector of length = number of treatments
##
.default.d.func <- function(x) {
  diff(range(x))
}

##
## Default measure of overall imbalance. Just sum of imbalances.
##
.default.g.func <- function(x) {
  sum(x)
}

##
## Default probablity assignments based on overall imbalance
##
.default.p.func <- function(overallImbalance) {
  number.of.treatments <- length(overallImbalance)
  p.star <- 2/3
  k <- which(overallImbalance == min(overallImbalance))
  if (length(k) > 1) {
    k <- sample(k, 1)
  }
  p.vec <- rep((1-p.star)/(number.of.treatments-1), number.of.treatments)
  p.vec[k] <- p.star
  p.vec
}

##
## Another alternative probability function based on overall imbalance
## Deterministically assign the treatment that will minimize overall imbalance
##
.alt.p.func.best <- function(overallImbalance) {
    number.of.treatments <- length(overallImbalance)
    k <- which(overallImbalance == min(overallImbalance))
    if (length(k) > 1) {
        k <- sample(k, 1)
    }
    p.vec <- rep(0, number.of.treatments)
    p.vec[k] <- 1
    p.vec
}

##
## Another alternative probability function based on overall imbalance
## Load the die so that the treatment that will minimize overall imbalance
## is heavily favored and the remaining ones are all equiprobable
## .alt.p.func.best (FAVORED.PROB = 1) and .default.p.func (FAVORED.PROB = 0.75)
## are special cases of this.
##
.alt.p.func.prob <- function(overallImbalance) {
    FAVORED.PROB <- 0.75
    number.of.treatments <- length(overallImbalance)
    k <- which(overallImbalance == min(overallImbalance))
    if (length(k) > 1) {
        k <- sample(k, 1)
    }
    p.vec <- rep((1-FAVORED.PROB)/(number.of.treatments-1), number.of.treatments)
    p.vec[k] <- FAVORED.PROB
    p.vec
}

##
## Another alternative probability function based on overall imbalance
## Load the die proportionately according to imbalance.
## Have to ensure overallImbalance is never negative!
##
.alt.p.func.prop <- function(overallImbalance) {
    p.vec <- overallImbalance / sum(overallImbalance)
    p.vec
}



setMethod("initialize", "PocockSimonRandomizer",
          function (.Object, expt, seed, stateTable, tr.ratios, d.func=.default.d.func,
                    g.func=.default.g.func, p.func=.default.p.func) {
            ##print("I am called")
            if (missing(seed)) seed <- 12345
            if (missing(stateTable)) {
              m <- expt@number.of.treatments
              n <- sum(expt@number.of.factor.levels)
              stateTable <- matrix(0, nrow=m, ncol=n)
              rownames(stateTable) <- expt@treatment.names
              colnames(stateTable) <- unlist(sapply(1:expt@number.of.factors,
                                                    function(i) 
                                                    sapply(1:expt@number.of.factor.levels[i],
                                                           function(j) paste(expt@factor.names[i],
                                                                             expt@factor.level.names[[i]][j], sep=":"))))
            }
            if (missing(tr.ratios)) {
              m <- expt@number.of.treatments
              tr.ratios <- rep(1, m)
            }
            
            .Object@expt <- expt
            .Object@seed <- as.integer(seed)
            .Object@stateTable <- stateTable
            .Object@tr.ratios <- tr.ratios / sum (tr.ratios)
            .Object@d.func <- d.func
            .Object@g.func <- g.func
            .Object@p.func <- p.func
            set.seed(seed)
            .Object
          })

##
## Validator for PocockSimonRandomizer
##
.validPocockSimonRandomizer <- function(object) {
    retval <- NULL
    expt <- object@expt
    m <- expt@number.of.treatments
    n <- sum(expt@number.of.factor.levels)
    
    if (nrow(object@stateTable) != m) {
        retval <- c(retval, "State matrix has the wrong number of rows")
    }
    
    if (ncol(object@stateTable) != n) {
        retval <- c(retval, "State matrix has the wrong number of columns")
    }

    treatment.names <- expt@treatment.names
    
    if (!all(rownames(object@stateTable) == treatment.names)) {
        retval <- c(retval, "State matrix has the wrong row names")
    }
        
    expected.colnames <- unlist(sapply(1:expt@number.of.factors,
                                       function(i) 
                                       sapply(1:expt@number.of.factor.levels[i],
                                              function(j) paste(expt@factor.names[i],
                                                                expt@factor.level.names[[i]][j], sep=":"))))
    
    if (!all(colnames(object@stateTable) == expected.colnames)) {
        retval <- c(retval, "State matrix has the wrong column names")
    }
        

    if (any(object@tr.ratios) <= 0) {
        retval <- c(retval, "Treatment ratios should be positive")
    }
    
    if (is.null(retval)) {
        return(TRUE)
    } else {
        return(retval)
    }
}

setValidity("PocockSimonRandomizer", .validPocockSimonRandomizer)

## ## Constructor function Don't need this because I have the initializer!
## "PocockSimonRandomizer" <- function(expt, seed,
##                                stateTable) {
##     my.expt <- as(expt, "ClinicalExperiment")
##     new("PocockSimonRandomizer",
##         my.expt,
##         as.integer(seed),
##         as.matrix(stateTable))
## }

##
## The State table contains all the relevant marginal information.
## For testing purposes, it is convenient to be able to start the
## Randomizer at particular states. Hence this generic method
##
if (!isGeneric("stateTable<-")) {
    if (is.function("stateTable<-")) {
        setGeneric("stateTable<-", "stateTable<-")
    } else {
    setGeneric("stateTable<-",
               function(x, value) standardGeneric("stateTable<-"))
    }

}

setReplaceMethod("stateTable", "PocockSimonRandomizer", function(x, value) {
    expt <- x@expt
    treatment.names <- expt@treatment.names
    ##print("In statetable")
    ##print(treatment.names)
    if (!all(rownames(value) == treatment.names)) {
        stop("Matrix has the wrong row names")
    }
    ##print("check passed")
    
    expected.colnames <- unlist(sapply(1:expt@number.of.factors,
                                       function(i) 
                                       sapply(1:expt@number.of.factor.levels[i],
                                              function(j) paste(expt@factor.names[i],
                                                                expt@factor.level.names[[i]][j], sep=":"))))
    
    if (!all(colnames(value) == expected.colnames)) {
        stop("Matrix has the wrong column names")
    }
    
    x@stateTable <- value
    x
})

##
## This sets the list of treatment assignments so far
##
if (!isGeneric("tr.assignments<-")) {
    if (is.function("tr.assignments<-")) {
        setGeneric("tr.assignments<-", "tr.assignments<-")
    } else {
        setGeneric("tr.assignments<-",
                   function(x, value) standardGeneric("tr.assignments<-"))
    }
}

setReplaceMethod("tr.assignments", "PocockSimonRandomizer", function(x, value) {
  ##
  ## Check for appropriate values
  ##
  if (!is.data.frame(value)) {
    stop("Expecting data frame for treatment assignments!")
  }
  expt <- x@expt
  if (ncol(value) != (expt@number.of.factors + 1)) {
    stop("Wrong dimension for treatment assignment data frame!")
  }
  if (any(names(value) != c(expt@factor.names, "Treatment"))) {
    stop("Wrong variable names for treatment assignments!")
  }

  if (nrow(value) > 0) {
    ##
    ## Now check for the values of the factor levels
    ##
    for (i in 1:expt@number.of.factors) {
      if (! all(value[, i] %in% expt@factor.level.names[[i]])) {
        stop(paste("Wrong value for variable", names(value)[i], "in treatment assignments!"))
      }
    }
    if (! all(value[, ncol(value)] %in% expt@treatment.names)) {
      stop(paste("Bad treatement names in treatment assignments!"))
    }
  }             
  x@tr.assignments <- value
  x
})

##
## We need a way of getting at the imbalance for each possible treatment, so we need a
## generic function that computes imbalances
##

if (!isGeneric("computeImbalances")) {
    if (is.function("computeImbalances")) {
        setGeneric("computeImbalances", computeImbalances)
    } else {
        setGeneric("computeImbalances",
                   function(object, factor.values) standardGeneric("computeImbalances"))
    }

}

setMethod("computeImbalances",
          signature(object="PocockSimonRandomizer", factor.values="character"),
          function (object, factor.values) {
              ##   factor.values <- as.integer(factor.values)
              if (missing(factor.values)) {
                  stop("Need factor values")
              }
              expt <- object@expt
              number.of.factors <- expt@number.of.factors
              if (length(factor.values) != number.of.factors) {
                  stop("Not correct number of factors")
              }
              number.of.factor.levels <- expt@number.of.factor.levels
              factor.level.names <- expt@factor.level.names
              factor.values.kosher <- sapply(1:number.of.factors,
                                             function(x) factor.values[x] %in% factor.level.names[[x]])
              if (! all(factor.values.kosher)) {
                stop("Incorrect factor values provided")
              }
              
              treatment.names <- expt@treatment.names
              number.of.treatments <- expt@number.of.treatments
              d.func <- object@d.func
              factor.names <- expt@factor.names
              state.matrix <- object@stateTable
              tr.ratios <- object@tr.ratios
              named.factors <- paste(factor.names, factor.values, sep=":")
              f.mat <- state.matrix[, named.factors]
              ##print("our mat")
              ##print(f.mat)
              ##print(treatment.names)

              result <- sapply(treatment.names, function(x) {
                ##print("fmat")
                ##print(f.mat)
                new.mat <- f.mat
                new.mat[x, ] <- new.mat[x, ] + 1
                ##print(new.mat)                  
                exp.mat <- apply(new.mat, 2, function(x) tr.ratios * sum(x))
                ## The D function in the paper
                ##print(exp.mat)
                apply(new.mat - exp.mat, 2, d.func)
              })
              colnames(result) <- treatment.names
              ##print("Imbalances")
              ##print(result)
              result
          })



if (!isGeneric("computeOverallImbalance")) {
    if (is.function("computeOverallImbalance")) {
        setGeneric("computeOverallImbalance", computeOverallImbalance)
    } else {
        setGeneric("computeOverallImbalance",
                   function(object, imbalances) standardGeneric("computeOverallImbalance"))
    }

}

##
## The G function in the paper
##
setMethod("computeOverallImbalance",
          signature(object="PocockSimonRandomizer", imbalances="matrix"),
          function (object, imbalances) {
              g.func <- object@g.func
              apply(imbalances, 2, g.func)
          })


if (!isGeneric("randomize")) {
    if (is.function("randomize")) {
        setGeneric("randomize", randomize)
    } else {
        setGeneric("randomize",
                   function(object, subject.id, factor.values) standardGeneric("randomize"))
    }              
    
}

setMethod("randomize",
          signature(object="PocockSimonRandomizer", subject.id="character", factor.values="character"),
          function (object, subject.id, factor.values) {
              if (missing(subject.id)) {
                  stop("Need subject id for randomization")
              }
              
              if (length(subject.id) > 1) {
                stop("Need a single subject id for randomization")
              }
              
              if (subject.id %in% rownames(object@tr.assignments)) {
                stop(paste("Subject ID", subject.id, "already randomized!"))
              }
                

              ##   factor.values <- as.integer(factor.values)
              if (missing(factor.values)) {
                  stop("Need factor values")
              }
              expt <- object@expt
              number.of.factors <- expt@number.of.factors
              if (length(factor.values) != number.of.factors) {
                stop("Not correct number of factors")
              }
              factor.level.names <- expt@factor.level.names
              factor.values.kosher <- sapply(1:number.of.factors,
                                             function(x) factor.values[x] %in% factor.level.names[[x]])
              if (! all(factor.values.kosher)) {
                stop("Incorrect factor values provided")
              }
              
              number.of.treatments <- expt@number.of.treatments
              treatment.names <- expt@treatment.names
              factor.names <- expt@factor.names
              named.factors <- paste(factor.names, factor.values, sep=":")
              imbalances <- computeImbalances(object, factor.values)
              overallImbalance <- computeOverallImbalance(object, imbalances)
              ##print("Overall imbalance")
              ##print(overallImbalance)
              p.func <- object@p.func
              tr.ratios <- object@tr.ratios
              
              p.vec <- p.func(overallImbalance)
              ##print("pvec")
              ##print(p.vec)
              tr.index <- sample(number.of.treatments, 1, prob=p.vec)
              tr.name <- treatment.names[tr.index]
              ##print(paste("Treatment is", tr.name))
              
              ## Update state
              state.matrix <- object@stateTable
              state.matrix[tr.name, named.factors] <- state.matrix[tr.name, named.factors] + 1
              stateTable(object) <- state.matrix
              
              ## Update assignment table
              current.assignment <- data.frame(as.list(factor.values), tr.name, stringsAsFactors=FALSE)
              rownames(current.assignment) <- subject.id
              colnames(current.assignment) <- c(factor.names, "Treatment")
              if (nrow(object@tr.assignments) == 0) {
                  assignments <- current.assignment
              } else {
                  assignments <- rbind(object@tr.assignments, current.assignment)
              }
              tr.assignments(object) <- assignments
              object
          })


if (!isGeneric("lastRandomization")) {
  if (is.function("lastRandomization")) {
      setGeneric("lastRandomization", lastRandomization)
  } else {
      setGeneric("lastRandomization",
                 function(object) standardGeneric("lastRandomization"))
  }

}

setMethod("lastRandomization",
          signature(object="PocockSimonRandomizer"),
          function (object) {
              n <- nrow(object@tr.assignments)
              if (n < 1) {
                  stop("No assignment yet!")
              } else {
                  object@tr.assignments[n, ]
              }
          })
