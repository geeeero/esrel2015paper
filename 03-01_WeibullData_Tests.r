# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
#                                                 
# S4 implementation of generalized iLUCK models   
# -- Tests for Class WeibullData incl. show method --             
#                                                 
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

library(luck)
library(testthat)
source("03-01_WeibullData.r")

weidata1 <- WeibullData(1:10)
weidata1 <- WeibullData(data=1:10)
weidata2 <- WeibullData(data=1:10, k = 3)

set.seed(8552)
weidata3 <- WeibullData(mean = 10, n = 5, sim = TRUE)
weidata4 <- WeibullData(mean = 10, n = 5, sim = TRUE, k = 3)

data1 <- LuckModelData(5, 2)
weidata5 <- WeibullData(data1)

weidata1
weidata2
weidata3
weidata4
weidata5
#