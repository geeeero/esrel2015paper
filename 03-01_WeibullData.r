# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
#
# S4 implementation of generalized iLUCK models
# -- Classes for inference on data from a Weibull distribution 
#    with fixed shape parameter k: W(k,lambda)  --
# (relation of parametrization with pweibull: shape a <-> k, scale b <-> lambda^(1/k))
# (furthermore, E[T|lambda] = lambda^(1/k) * gamma(1 + 1/k).)
#
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #


# ---------------------------------------------------------------------------- #
# Class for data: ScaledNormalData
# ---------------------------------------------------------------------------- #

# data class definition
setClass ("WeibullData",
          contains = "LuckModelData")
# TODO: add k as slot for WeibullData

# constructor function handling as input
#    an object of class LuckModelData (or some subclass)
# or an observed data vector
# or a mean and sample size to simulate data according to.
# shape parameter k defaults to 2.
WeibullData <- function (arg1 = NULL, arg2 = NULL, data = NULL,
                         mean = NULL, n = NULL, sim = FALSE, k = 2){
 if (all(is.null(c(arg1, arg2, data, mean, n)))) {
  stop("No arguments given for WeibullData()!")
 } else {
  # if arg1 is a LuckModelData object, turn it into a WeibullData object directly
  if (is(arg1, "LuckModelData")) {
   return (new("WeibullData", tauN = tauN(arg1), rawData = rawData(arg1)))
  } else {
  # see if a "ready" data vector is given either named or as a single unnamed argument
   if (is.vector(data) | (is.vector(arg1) & is.null(c(arg2, data, mean, n)))) {
    if (is.vector(data))
     arg1 <- data
    return (.WeibullDataVector(data = arg1, k = k))
   } else { # no data vector is given
    if (!sim){
      stop("If no data is given, set sim = TRUE to simulate data.")
    } else {
     # check if some of the named arguments "mean" and "n" are given.
     # if so, write them in arg1 (mean) and arg2 (n)
     # if not, arg1 and arg2 give already the mean and sample size
     if (!is.null(mean) & is.samplesize(n)) { # mean and n are given
      arg1 <- mean  
      arg2 <- n
     }
     if (is.null(mean) & is.samplesize(n)) { # unnamed argument in arg1 is mean
      arg2 <- n     
     }
     if (!is.null(mean) & is.null(n)) {  # unnamed argument in arg1 is n
      arg2 <- arg1 
      arg1 <- mean
     }
     # check (again) if arg1 and arg2 are ok
     if (!is.vector(arg1))
      stop("mean must be a single value or, if used to simulate an inhomogeneous sample, a vector.")
     if (!is.samplesize(arg2))
      stop("n must be a proper sample size, i.e., a positive integer.")
     # take mean and n and simulate data according to them
     return (.WeibullDataSimulate(mean = arg1, n = arg2, k = k))
    }
   } 
  }
 }
}

# constructor function if an observed data vector is given:
.WeibullDataVector <- function (data, k){
 # create tauN with LuckModelData constructor function
 .tau <- sum(data^k)
 .dataobj <- LuckModelData(tau = .tau, n = length(data))
 .tauN <- tauN(.dataobj)
 # create the object with rawdata = data
 new("WeibullData", tauN = .tauN, rawData = matrix(data, ncol = 1))
}

# constructor function if data is simulated:
.WeibullDataSimulate <- function (mean, n, k){
 # simulate data from a Weibull distribution
 .data <- rweibull (n, shape = k, scale = mean / gamma(1 + 1/k))
 # use then the data constructor function
 .WeibullDataVector(data = .data, k = k)
}


# accessor methods are inherited from LuckModelData 

# replacement method for slot rawData
# (other replacement methods are inherited from LuckModelData)
if(!isGeneric("rawData<-"))
 setGeneric("rawData<-", function(object, ...) standardGeneric("rawData<-"))
setReplaceMethod("rawData", "WeibullData", function(object, value){
 # create entirely new object with the new data
 WeibullData(data = value)
})


# show method for class WeibullData (no plot method)
setMethod("show", "WeibullData", function(object){
 if (is.null(tauN(object))) {
  cat("Default data object with no data specified.\n")
 } else {
  .rawData <- rawData(object)
  .tau <- tau(object)
  .n <- n(object)
  if (is.null(.rawData)) {
   cat("WeibullData object containing sufficient statistic tau =", .tau, "and sample size", .n, ".\n")
  } else { # object has rawData
   cat("WeibullData object containing", .n, "observations: \n",
       as.vector(.rawData), "\n")
  }
 }
})


#