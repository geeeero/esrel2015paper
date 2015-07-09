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
# tnow    time until the system is observed
# fts     vector of length e giving the observed failure times
# kappa   fixed weibull shape parameter
postpredC <- function(n0y0, kappa, n, fts, tnow, t, l){
  if (t < tnow)
    stop("t must be larger than tnow")
  e <- length(fts)
  if (n < e)
    stop("there can be at most n failure times!")
  if (l > n-e)
    stop("l can be at most n-e")
  nn <- n0y0[1] + e
  nnyn <- n0y0[1]*n0y0[2] + (n-e)*tnow + sum(fts^kappa)
  j <- seq(0, n-e-l)
  choose(n-e, l) * sum( (-1)^j * choose(n-e-l, j) * (nnyn/(nnyn + (l+j)*(t^kappa - tnow^kappa)))^(nn + 1))
}

kappa <- 2
n0y0_1 <- c(2, failuretolambda(9,kappa))
              # n0 = 2, y0 corresponding to mean failure time 9
n_1 <- 4      # four components of type 1
fts_1 <- c(1,2) # two failures at times 1 and 2, so l = 0,1,2
tnow <- 3     # observed until tnow=3, so t > 3 
postpredC(n0y0_1, kappa, n_1, fts_1, tnow, t = 4, l = 2)
  
postpredCpmf <- function(n0y0, kappa, n, fts, tnow, t){
  l <- seq(0, n-length(fts))
  res <- numeric(length(l))
  for (i in l) res[i+1] <- postpredC(n0y0, kappa, n, fts, tnow, t, i)
  res <- array(res)
  dimnames(res)[[1]] <- l
  res
}

postpredCpmf(n0y0_1, kappa, n_1, fts_1, tnow, t = 4)

postpredCpmf <- function(n0y0, kappa, n, fts, tnow, t){
  l <- seq(0, n-length(fts))
  res <- numeric(length(l))
  for (i in l) res[i+1] <- postpredC(n0y0, kappa, n, fts, tnow, t, i)
  res <- array(res)
  dimnames(res)[[1]] <- l
  res
}

postpredCcmf <- function(n0y0, kappa, n, fts, tnow, t){
  pmf <- postpredCpmf(n0y0, kappa, n, fts, tnow, t)
  cmf <- cumsum(pmf)
  cmf
}

Ccmfplot <- function(cmf, add = FALSE, ylim = c(0,1), xlab = "l", ylab = "F(C = l)",...){
  if(add)
    lines(as.numeric(names(cmf)), cmf, type="s", ...)
  else
    plot(as.numeric(names(cmf)), cmf, type="s", ylim = ylim, xlab = xlab, ylab = ylab, ...)
}

asdf <- postpredCcmf(n0y0_1, kappa, n_1, fts <- c(7,8), tnow=8, t = 9)
asdf2 <- postpredCcmf(n0y0_1, kappa, n_1, fts <- c(7,8), tnow=8, t = 9)
Ccmfplot(asdf)
Ccmfplot(asdf2, add=TRUE, lty=2)

#