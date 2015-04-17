# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
#
# S4 implementation of generalized iLUCK models
# -- Classes for inference on data from a Weibull distribution 
#    with fixed shape parameter k W(k,lambda)  --
# (relation of parametrization with pweibull: a <-> k, b <-> lambda^(1/k))
#
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

# provides dinvgamma(x, shape, scale), where
# shape = alpha = n0 + 1 and 
# scale = beta = n0*y0
library(actuar)

# ---------------------------------------------------------------------------- #
# Class for model: WeibullLuckModel
# ---------------------------------------------------------------------------- #

# class definition: extends class LuckModel
setClass ("WeibullLuckModel",
          contains = "LuckModel")

# ---------------------------------------------------------------------------- #
# Constructor function for WeibullLuckModel
# ---------------------------------------------------------------------------- #

# arguments may either be the arguments for LuckModel()
# or an object of class LuckModel()
WeibullLuckModel <- function (arg1 = NULL, n0 = NULL, y0 = NULL, data = new("WeibullData")){
 if (all(is.null(c(arg1, n0, y0)))) {
  stop("No arguments given for WeibullLuckModel()!")
 } else {
  # .data will stay a default object if no data argument was given.
  .data <- WeibullData(data)
  # if arg1 is not a LuckModel, see if n0 and y0 are given.
  # If so, read out n0 and y0.
  if (!is(arg1, "LuckModel")){
   if (any(is.null(c(n0, y0))))
    stop("To define a WeibullLuckModel both arguments n0 and y0 must be given!")
   .n0 <- n0
   .y0 <- y0
  } else { # if arg1 is a LuckModel, read out its slots
   .n0 <- n0(arg1)
   .y0 <- y0(arg1)
   if (!is.null(tauN(data(arg1)))) { # LuckModel contains data    
    if (!is.null(tauN(.data))) # data was given also in the constructor call
     stop("Object to be converted contains already some data.\n Replace this data later if it should be changed.")
    # overwrite default data object with data slot from arg1
    .data <- WeibullData(data(arg1))
   }
  }
  # build a LuckModel object from the slots for checking inputs
  .object <- LuckModel(n0 = .n0, y0 = .y0, data = .data)
  # check if y0 is one-dimensional, n0 > 1 and y0 > 0
  if (dim(y0(.object))[1] != 1)
   stop("For inference on Weibull data with given shape parameter k, y0 must be one-dimensional!")
  if (y0(.object)[1] <= 0)
   stop("For inference on Weibull data with given shape parameter k, y0 must be positive!")
  if (n0(.object)[1] <= 1)
   stop("For inference on Weibull data with given shape parameter k, n0 must larger than 1!")
  # now create the WeibullLuckModel object from the LuckModel object
  new("WeibullLuckModel", n0 = n0(.object), y0 = y0(.object), data = .data)
 }
}


# accessor and replacement methods are inherited from LuckModel


# ---------------------------------------------------------------------------- #
# show (= print in S3) method for WeibullLuckModel objects                                                             
# ---------------------------------------------------------------------------- #

# helper function to produce Weibull specific text 
.pItextWei <- function(object) {
 .oneNoneY = paste("corresponding to an inverse gamma prior\n with mean", #
                   y0(object)[1], #                            
                   "and variance", #                             
                   (y0(object)[1])^2/(1 - 1/n0(object)[1]))
 .oneNtwoY = paste("corresponding to a set of inverse gamma priors\n with means in [", #
                   y0(object)[1,1], ";", y0(object)[1,2], #                    
                   "] and variances in [", #                                         
                   (y0(object)[1])^2/(1 - 1/n0(object)[1]), ";", #
                   (y0(object)[2])^2/(1 - 1/n0(object)[1]), "]")    
 .twoNoneY = paste("corresponding to a set of inverse gamma priors\n with mean", #
                   y0(object)[1,1], #                                    
                   "and variances in [", #                               
                   (y0(object)[1])^2/(1 - 1/n0(object)[2]), ";", #
                   (y0(object)[1])^2/(1 - 1/n0(object)[1]), "]")    
 .twoNtwoY = paste("corresponding to a set of inverse gamma priors\n with means in [", #                                   
                   y0(object)[1,1], ";", y0(object)[1,2], #                    
                   "] and variances in [", #                                   
                   (y0(object)[1])^2/(1 - 1/n0(object)[2]), ";", #
                   (y0(object)[2])^2/(1 - 1/n0(object)[1]), "]")    
 .genilucktree(object = object, oneNoneY = .oneNoneY, #
                                oneNtwoY = .oneNtwoY, #
                                twoNoneY = .twoNoneY, #
                                twoNtwoY = .twoNtwoY)
}


# show method uses helper function .showLuckModels (see show method for LuckModel)
setMethod("show", "WeibullLuckModel", function(object){
 .showLuckModels(object = object, #
                 forInference = "for inference from weibull data with fixed shape parameter k\n", #
                 parameterInterpretation = .pItextWei(object))
})


# ---------------------------------------------------------------------------- #
# singleHdi: function returning symmetric credibility intervals for WeibullLuckModel
# ---------------------------------------------------------------------------- #
# shape = alpha = n0 + 1 and 
# scale = beta = n0*y0

if(!isGeneric("singleHdi"))
 setGeneric("singleHdi", function(object, ...) standardGeneric("singleHdi"))

# argument object is used only for method dispatching
setMethod("singleHdi", "WeibullLuckModel", function(object, n, y, gamma) {
 # TODO: proper highest density interval instead of symmetric interval?
 # using the quantiles for symmetric credibility interval
 .lqua <- (1-gamma)/2    # lower quantile
 .uqua <- gamma + .lqua  # upper quantile  
 .lhd <- qinvgamma (.lqua, shape = n + 1, scale = n*y) # lower border
 .uhd <- qinvgamma (.uqua, shape = n + 1, scale = n*y) # upper border
 c(.lhd,.uhd) # return the credibility interval
})


# ---------------------------------------------------------------------------- #
# singleCdf: function returning cdf values for WeibullLuckModel
# ---------------------------------------------------------------------------- #

if(!isGeneric("singleCdf"))
 setGeneric("singleCdf", function(object, ...) standardGeneric("singleCdf"))

# argument object is used only for method dispatching
setMethod("singleCdf", "WeibullLuckModel", function(object, n, y, x) {
 pinvgamma (x, shape = n + 1, scale = n*y)
})


#