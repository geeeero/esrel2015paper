# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
#                                                 
# S4 implementation of generalized iLUCK models   
# -- Tests for Class WeibullLuckModel incl. show method, unionHdi and cdfplot --             
#                                                 
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

library(luck)
library(testthat)
source("/home/gero/work/QdC+n/R/r-forge/luck/pkg/luck/R/00-02_LuckModel.r")
source("03-02_Weibull.r")

# constructor function
bsp1 <- LuckModel(n0 = c(2,5), y0 = c(103, 154))
weibsp1 <- WeibullLuckModel(bsp1)
n0(weibsp1)
y0(weibsp1)

bsp1e <- LuckModel(n0 = c(1,5), y0 = c(1, 2))
expect_error(WeibullLuckModel(bsp1e))
bsp1e <- LuckModel(n0 = c(2,5), y0 = c(-1, 2))
expect_error(WeibullLuckModel(bsp1e))
bsp2 <- WeibullLuckModel(n0=2, y0=2)
bsp2
bsp3 <- WeibullLuckModel(n0=c(2,3), y0=2)
bsp3
bsp4 <- WeibullLuckModel(n0=2, y0=c(2,5))
bsp4
bsp5 <- WeibullLuckModel(n0=c(2,3), y0=c(2,5))
bsp5
bsp6 <- WeibullLuckModel(n0=c(2,5), y0=c(1,10), data = WeibullData(1:2))
bsp6
n0(bsp6)
y0(bsp6)
data(bsp6)
tauN(data(bsp6))

# singleHdi 
singleHdi(weibsp1, n=2, y=103, gamma=0.95)
singleHdi(weibsp1, n=5, y=103, gamma=0.95)
singleHdi(weibsp1, n=2, y=154, gamma=0.95)
singleHdi(weibsp1, n=5, y=154, gamma=0.95)

# unionHdi
unionHdi(weibsp1)
data(weibsp1) <- 1:2
unionHdi(weibsp1, posterior=TRUE)$borders

# singleCdf
singleCdf(weibsp1, n = 2, y = 103, x = 50)

# cdfplot
cdfplot(weibsp1)
cdfplot(weibsp1, xvec = 20)
cdfplot(weibsp1, xvec = seq(0, 200, length.out = 300))
cdfplot(weibsp1, xvec = seq(100, 250, length.out = 300))

cdfplot(weibsp1, xvec = seq(0, 200, length.out = 300), control = controlList(posterior = TRUE))
cdfplot(weibsp1, xvec = seq(0, 200, length.out = 300), control = controlList(posterior = TRUE))
#