###############################################################################
#                       Code for ISIPTA'15 system poster                      #
###############################################################################

# function to calculate P(C_t = l | n^(0), y^(0), data),
# the posterior predictive probability that l components function at time t
# in the system observed until t_now (notion "of type k" dropped here)
# this implements (19)
# n0y0    pair c(n0,y0) of prior parameters
# t       time t for which to calculate P(C_t), t > t_now
# l       number of functioning components, \in {0, 1, ..., n-e}
# n       number of components in the system
# e       number of failed components at time t_now
# tnow    time until the system is observed
# fts     vector of length e giving the observed failure times
# kappa   fixed weibull shape parameter
postpredC <- function(n0y0, t, l, n, e, tnow, fts, kappa){
  if (length(fts) != e)
    stop("fts must be a vector of e failure times")
  if (l > n-e)
    stop("l can be at most n-e")
  nn <- n0y0[1] + e
  nnyn <- n0y0[1]*n0y0[2] + (n-e)*tnow + sum(fts^kappa)
  j <- seq(0, n-e-l)
  choose(n-e, l) * sum( (-1)^j * choose(n-e-l, j) * (nnyn/(nnyn + (l+j)*(t^kappa - tnow^kappa)))^(nn + 1))
}








#