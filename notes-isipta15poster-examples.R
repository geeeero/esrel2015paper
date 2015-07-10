###############################################################################
#                     Examples for ISIPTA'15 system poster                    #
###############################################################################

kappa <- 2
n0y0_1 <- c(2, failuretolambda(9,kappa))
# n0 = 2, y0 corresponding to mean failure time 9
n_1 <- 4      # four components of type 1
fts_1 <- c(1,2) # two failures at times 1 and 2, so l = 0,1,2
tnow <- 3     # observed until tnow=3, so t > 3 
postpredC(n0y0_1, kappa, n_1, fts_1, tnow, t = 4, l = 2)

asdf <- postpredCcmf(n0y0_1, kappa, n_1, fts <- c(7,8), tnow=8, t = 9)
asdf2 <- postpredCcmf(n0y0_1, kappa, n_1, fts <- c(7,8), tnow=8, t = 9)
Ccmfplot(asdf)
Ccmfplot(asdf2, add=TRUE, lty=2)

fc1 <- LuckModel(n0 = c(2,5), y0 = c(failuretolambda(9,kappa), failuretolambda(11,kappa)))
#n_1 <- 4      # four components of type 1
#fts_npdc <- c(10,11) 
#tnow <- 11
fourCornersCcmf(fc1, kappa, n=4, fts=c(10,11), tnow=11, t=12) # F_lower = bl, F_upper = tr
fourCornersCcmf(fc1, kappa, n=4, fts=c(10,11), tnow=11, t=15) # F_lower = bl-br!!, F_upper = tr
fourCornersCcmf(fc1, kappa, n=4, fts=c(10,11), tnow=11, t=20) # F_lower = br, F_upper = tr-tl!!!

fourCornersCcmf(fc1, kappa, n=4, fts=c(1,2), tnow=2, t=3)  # F_lower = bl, F_upper = tr
fourCornersCcmf(fc1, kappa, n=4, fts=c(1,2), tnow=2, t=5)  # F_lower = bl, F_upper = tr
fourCornersCcmf(fc1, kappa, n=4, fts=c(1,2), tnow=2, t=10) # F_lower = bl, F_upper = tr
fourCornersCcmf(fc1, kappa, n=4, fts=c(1,2), tnow=2, t=20) # F_lower = bl, F_upper = tr

fourCornersCcmf(fc1, kappa, n=4, fts=c(1,2), tnow=8, t=10) # F_lower = bl, F_upper = tr

library(ReliabilityTheory)

sys1 <- graph.formula(s -- 1:2 -- 3:4:5 -- t)
V(sys1)$compType <- NA # This just creates the attribute compType
V(sys1)$compType[match(c("1","3"), V(sys1)$name)] <- "Type 1"
V(sys1)$compType[match(c("2","4"), V(sys1)$name)] <- "Type 2"
V(sys1)$compType[match(c("5"), V(sys1)$name)] <- "Type 3"
V(sys1)$compType[match(c("s","t"), V(sys1)$name)] <- NA
sys1sign <- computeSystemSurvivalSignature(sys1)

n0y0_1 <- c(2, failuretolambda(9,kappa))
n0y0_2 <- c(5, failuretolambda(4,kappa))
n0y0_3 <- c(3, failuretolambda(10,kappa))


sysrel(n0y0 = list(n0y0_1,n0y0_2,n0y0_3), survsign = sys1sign, kappa=rep(2,3),
       fts = list(c(7), c(3,4), NULL), tnow = 7, t = 8)
sysrel(n0y0 = list(n0y0_1,n0y0_2,n0y0_3), survsign = sys1sign, kappa=rep(2,3),
       fts = list(c(7), c(3,4), NULL), tnow = 7, t = 8, table = T)

tvec <- seq(7,15,by=0.1)
rvec <- sapply(tvec, FUN = sysrel, n0y0 = list(n0y0_1,n0y0_2,n0y0_3),
               survsign = sys1sign, kappa=rep(2,3), fts = list(c(7), c(3,4), NULL), tnow = 7)
plot(tvec, rvec, type="l", ylim=c(0,1), xlab="t", ylab=expression(R[sys](t)))
# starts with 0.5 as one of type 1 and one of type 3 still function,
# and the working type 1 component could be either of the two in the system.

rvec <- sapply(tvec, FUN = sysrel, n0y0 = list(n0y0_1,n0y0_2,n0y0_3),
               survsign = sys1sign, kappa=rep(2,3), fts = list(c(7), c(3), NULL), tnow = 7)
plot(tvec, rvec, type="l", ylim=c(0,1), xlab="t", ylab=expression(R[sys](t)))

rvec <- sapply(tvec, FUN = sysrel, n0y0 = list(n0y0_1,n0y0_2,n0y0_3),
               survsign = sys1sign, kappa=rep(2,3), fts = list(NULL, c(4), c(7)), tnow = 7)
plot(tvec, rvec, type="l", ylim=c(0,1), xlab="t", ylab=expression(R[sys](t)))

# more complex system:
fig3 <- graph.formula(s -- 1:4 -- 2:5 -- 3:6 -- t, s -- 7:8, 8 -- 9, 7:9 -- t)
plot(fig3)
V(fig3)$compType <- NA # This just creates the attribute compType
V(fig3)$compType[match(c("1"), V(fig3)$name)] <- "Type 1"
V(fig3)$compType[match(c("2","3","4","7"), V(fig3)$name)] <- "Type 2"
V(fig3)$compType[match(c("5","6","8","9"), V(fig3)$name)] <- "Type 3"
V(fig3)$compType[match(c("s","t"), V(fig3)$name)] <- NA
sys3sign <- computeSystemSurvivalSignature(fig3)

tvec <- seq(7,20,by=0.1)
rvec <- sapply(tvec, FUN = sysrel, n0y0 = list(n0y0_1,n0y0_2,n0y0_3),
               survsign = sys3sign, kappa=rep(2,3), fts = list(NULL, c(4,5), c(7)), tnow = 7)
plot(tvec, rvec, type="l", ylim=c(0,1), xlab="t", ylab=expression(R[sys](t)))






#