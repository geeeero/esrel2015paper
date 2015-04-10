# --------------------------------------------------------------------------- #
# ESREL 2015 paper - code for Fig 1
# --------------------------------------------------------------------------- #

library(luck)

# translate lambda parameter of weibull to expected failure time
# lambda - scale parameter of weibull
# k      - shape parameter of weibull
lambdatofailure <- function(lambda, k = 2){
  lambda^(1/k) * gamma(1 + 1/k)
}

# translate failure time to lambda parameter of weibull
# ft - expected failure time
# k  - shape parameter of weibull
failuretolambda <- function(ft, k){
  (ft/gamma(1 + 1/k))^k
}

# example for Fig 1: assume mean ft = 9
k <- 2
y0 <- failuretolambda(9,k)
n0 <- 2
t1 <- 1
t2 <- 2

n2 <- n0 + 2
y2 <- (n0*y0 + t1^k + t2^k)/n2

eft2 <- lambdatofailure(y2,k)

#