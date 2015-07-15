###############################################################################
#                     Examples for ISIPTA'15 system poster                    #
###############################################################################

source("notes-isipta15poster.R")

kappa <- 2
n0y0_1 <- c(2, failuretolambda(9,kappa))
# n0 = 2, y0 corresponding to mean failure time 9
n_1 <- 4      # four components of type 1
fts_1 <- c(1,2) # two failures at times 1 and 2, so l = 0,1,2
tnow <- 3     # observed until tnow=3, so t > 3 
postpredC(n0y0_1, kappa = 2, n_1, fts_1, tnow, t = 4, l = 2)

asdf <- postpredCcmf(n0y0_1, 2, n_1, fts = c(7,8), tnow=8, t = 9)
asdf2 <- postpredCcmf(n0y0_1, 2, n_1, fts = c(7,8), tnow=8, t = 9)
Ccmfplot(asdf)
Ccmfplot(asdf2, add=TRUE, lty=2)

# -------------------------------------------------------

fcc <- LuckModel(n0 = c(1,10), y0 = c(failuretolambda(9,2), failuretolambda(11,2)))
applyfun <- function(n0, y0, kappa, n, fts, tnow, t)
  postpredCcmf(n0y0 = c(n0, y0), kappa = kappa, n = n, fts = fts, tnow = tnow, t = t)

n0vec <- seq(n0(fcc)[1], n0(fcc)[2], length.out=101)
ccmfs <- sapply(n0vec, applyfun, y0 = y0(fcc)[1], kappa = 2, n = 4, fts = c(10), tnow = 11, t = 15)
ccmfs <- t(rbind(n0vec, ccmfs))
matplot(ccmfs[,1],  ccmfs[,2:5], type="l", ylim=c(0,1), lty = 1, col = as.numeric(dimnames(ccmfs)[[2]][-1]) + 1)

fourCornersCcmf(fcc, kappa = 2, n=4, fts=c(10), tnow=11, t=15)

# -------------------------------------------------------

fc1 <- LuckModel(n0 = c(2,5), y0 = c(failuretolambda(9,2), failuretolambda(11,2)))
#n_1 <- 4      # four components of type 1
#fts_npdc <- c(10,11) 
#tnow <- 11
pdf("cmf1.pdf")
fourCornersCcmf(fc1, kappa = 2, n=4, fts=c(10,11), tnow=11, t=12) # F_lower = bl, F_upper = tr
dev.off()
pdf("cmf2.pdf")
fourCornersCcmf(fc1, kappa = 2, n=4, fts=c(10,11), tnow=11, t=15) # F_lower = bl-br!!, F_upper = tr
dev.off()
pdf("cmf3.pdf")
fourCornersCcmf(fc1, kappa = 2, n=4, fts=c(10,11), tnow=11, t=20) # F_lower = br, F_upper = tr-tl!!!
dev.off()

fourCornersCcmf(fc1, kappa = 2, n=4, fts=c(1,2), tnow=2, t=3)  # F_lower = bl, F_upper = tr
fourCornersCcmf(fc1, kappa = 2, n=4, fts=c(1,2), tnow=2, t=5)  # F_lower = bl, F_upper = tr
fourCornersCcmf(fc1, kappa = 2, n=4, fts=c(1,2), tnow=2, t=10) # F_lower = bl, F_upper = tr
fourCornersCcmf(fc1, kappa = 2, n=4, fts=c(1,2), tnow=2, t=20) # F_lower = bl, F_upper = tr

fourCornersCcmf(fc1, kappa, n=4, fts=c(1,2), tnow=8, t=10) # F_lower = bl, F_upper = tr

sys1 <- graph.formula(s -- 1:2 -- 3:4:5 -- t)
V(sys1)$compType <- NA # This just creates the attribute compType
V(sys1)$compType[match(c("1","3"), V(sys1)$name)] <- "Type 1"
V(sys1)$compType[match(c("2","4"), V(sys1)$name)] <- "Type 2"
V(sys1)$compType[match(c("5"), V(sys1)$name)] <- "Type 3"
V(sys1)$compType[match(c("s","t"), V(sys1)$name)] <- NA
sys1sign <- computeSystemSurvivalSignature(sys1)

n0y0_1 <- c(2, failuretolambda(9, 2))
n0y0_2 <- c(5, failuretolambda(4, 2))
n0y0_3 <- c(3, failuretolambda(10,2))


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

# simpler system again
fc1 <- LuckModel(n0 = c(2,5), y0 = c(failuretolambda(9, 2), failuretolambda(11, 2)))
fc2 <- LuckModel(n0 = c(4,8), y0 = c(failuretolambda(4, 2), failuretolambda(6, 2)))
fc3 <- LuckModel(n0 = c(1,3), y0 = c(failuretolambda(9, 2), failuretolambda(11, 2)))
fclist <- list(fc1,fc2,fc3)

tvec <- seq(7,15,by=0.1)
r1 <- fourKcornersSysrel(luckobjlist = fclist, survsign = sys1sign, kappa = rep(2,3),
                         fts = list(c(7), c(3), NULL), tnow = 7, tvec = tvec)
fourKcornersSysrelPlot(tvec = tvec, rframe = r1$lower, legend = TRUE, add = FALSE, col = c(1,2,4), lty = c(1,1,1,2,2,2,3,3))
fourKcornersSysrelPlot(tvec = tvec, rframe = r1$lower, legend = TRUE, add = FALSE, col = 1:8, lwd = 2)
fourKcornersSysrelPlot(tvec = tvec, rframe = r1$upper, legend = FALSE, add = TRUE, col = 1:8, lwd = 2, lty = 2)

tvec <- seq(7,15,by=0.1)
r3 <- fourKcornersSysrel(luckobjlist = fclist, survsign = sys3sign, kappa = rep(2,3),
                         fts = list(NULL, c(3,4), 7), tnow = 7, tvec = tvec)
fourKcornersSysrelPlot(tvec = tvec, rframe = r3$lower, legend = TRUE, add = FALSE, col = 1:8, lwd = 2)
fourKcornersSysrelPlot(tvec = tvec, rframe = r3$upper, legend = FALSE, add = TRUE, col = 1:8, lwd = 2, lty = 2)
# there might be a switch in the lower here, too.

r1luck10 <- sysrelLuck(luckobjlist = fclist, survsign = sys1sign, kappa = rep(2,3),
                       fts= list(c(7), c(3), NULL), tnow = 7, t = 10)


# --------------------------------------------------------------------------
# --------------------------------------------------------------------------

par(mar = c(4, 4, 3.5, 0.2) + 0.1)
# example system with full signature
tvec <- seq(7,15,by=0.1)
r1 <- fourKcornersSysrel(luckobjlist = fclist, survsign = sys1sign, kappa = rep(2,3),
                         fts = list(c(7), c(3), NULL), tnow = 7, tvec = tvec)
pdf("fourKcornerSys1.pdf", width=5, height=4)
par(mar = c(4, 4, 3.5, 0.2) + 0.1)
fourKcornersSysrelPlot(tvec = tvec, rframe = r1$lower, legend = TRUE, add = FALSE, col = 1:8, lwd = 2)
fourKcornersSysrelPlot(tvec = tvec, rframe = r1$upper, legend = FALSE, add = TRUE, col = 1:8, lwd = 2, lty = 2)
#lines(rep(10,2), r1luck10$rel, lwd=3, lend=2)
dev.off()
pdf("rsyspbox1.pdf", width=5, height=4)
par(mar = c(4, 4, 3.5, 0.2) + 0.1)
pbox1 <- sysrelPbox(luckobjlist = fclist, survsign = sys1sign, kappa = rep(2,3),
                    fts = list(c(7), c(3), NULL), tnow = 7, tvec = tvec, returnres = T)
fourKcornersSysrelPlot(tvec = tvec, rframe = r1$lower, add = TRUE)
fourKcornersSysrelPlot(tvec = tvec, rframe = r1$upper, add = TRUE)
dev.off()
#look in r1 to see whether the switch between 111 and 121 in the lower bound really happens
checkswitchl <- r1$lower[-c(1:3),c(1,3)]
checkswitchl <- cbind(checkswitchl, checkswitchl[,1] > checkswitchl[,2], tvec)
table(checkswitchl[,3])
#look in r1 to see whether there is also a switch between 212 and 222 in the upper bound
checkswitchu <- r1$upper[-c(1:3),c(6,8)]
checkswitchu <- cbind(checkswitchu, checkswitchu[,1] > checkswitchu[,2], tvec)
table(checkswitchu[,3]) # no switch


# example system with reduced signature
sys1f <- graph.formula(s -- 1 -- 2:3 -- t)
V(sys1f)$compType <- NA # This just creates the attribute compType
V(sys1f)$compType[match(c("1"), V(sys1f)$name)] <- "Type 1"
V(sys1f)$compType[match(c("2"), V(sys1f)$name)] <- "Type 2"
V(sys1f)$compType[match(c("3"), V(sys1f)$name)] <- "Type 3"
V(sys1f)$compType[match(c("s","t"), V(sys1f)$name)] <- NA
sys1fsign <- computeSystemSurvivalSignature(sys1f)

tvec <- seq(7,15,by=0.1)
r1f <- fourKcornersSysrel(luckobjlist = fclist, survsign = sys1fsign, kappa = rep(2,3),
                          fts = list(c(7), c(3), NULL), tnow = 7, tvec = tvec, nk = c(2,2,1))
pdf("fourKcornerSys1failed.pdf", width=5, height=4)
par(mar = c(4, 4, 3.5, 0.2) + 0.1)
fourKcornersSysrelPlot(tvec = tvec, rframe = r1f$lower, legend = TRUE, add = FALSE, col = 1:8, lwd = 2)
fourKcornersSysrelPlot(tvec = tvec, rframe = r1f$upper, legend = FALSE, add = TRUE, col = 1:8, lwd = 2, lty = 2)
#lines(rep(10,2), r1luck10$rel, lwd=3, lend=2)
dev.off()
pdf("rsyspbox1failed.pdf", width=5, height=4)
par(mar = c(4, 4, 3.5, 0.2) + 0.1)
pbox1 <- sysrelPbox(luckobjlist = fclist, survsign = sys1fsign, kappa = rep(2,3),
                    fts = list(c(7), c(3), NULL), tnow = 7, tvec = tvec, returnres = T, nk = c(2,2,1))
fourKcornersSysrelPlot(tvec = tvec, rframe = r1f$lower, add = TRUE)
fourKcornersSysrelPlot(tvec = tvec, rframe = r1f$upper, add = TRUE)
dev.off()
#look in r1f to verify the switch between 111 and 121 in the lower bound
checkswitchlf <- r1f$lower[-c(1:3),c(1,3)]
checkswitchlf <- cbind(checkswitchlf, checkswitchlf[,1] > checkswitchlf[,2], tvec)
table(checkswitchlf[,3])

#-------------------------------------------------------------------------------------

# check in example system sys1 that when only the single component of type 3 fails,
# the calculations with the full and the reduced signature lead to the same results

# determine the reduced signature when type 3 comoponent fails
sys13f <- graph.formula(s -- 1:2 -- 3:4 -- t)
V(sys13f)$compType <- NA # This just creates the attribute compType
V(sys13f)$compType[match(c("1","3"), V(sys13f)$name)] <- "Type 1"
V(sys13f)$compType[match(c("2","4"), V(sys13f)$name)] <- "Type 2"
#V(sys13f)$compType[match(c("5"), V(sys1)$name)] <- "Type 3"
V(sys13f)$compType[match(c("s","t"), V(sys13f)$name)] <- NA
sys13fsign <- computeSystemSurvivalSignature(sys13f)
tvec <- seq(7,20,by=0.1)
# with full signature when type 3 component fails
r13 <- fourKcornersSysrel(luckobjlist = list(fc1,fc2,fc3), survsign = sys1sign, kappa = rep(2,3),
                          fts = list(NULL, NULL, c(7)), tnow = 7, tvec = tvec, nk = c(2,2,1))
# with reduced signature when type 3 component fails
r13f <- fourKcornersSysrel(luckobjlist = list(fc1,fc2), survsign = sys13fsign, kappa = rep(2,3),
                           fts = list(NULL, NULL), tnow = 7, tvec = tvec, nk = c(2,2))
# plot
pdf("rsys3failedtest.pdf", width=5, height=4)
par(mfrow=c(1,2), mar = c(4, 4, 3.5, 0.2) + 0.1)
pbox13 <- sysrelPbox(luckobjlist = list(fc1,fc2,fc3), survsign = sys1sign, kappa = rep(2,3),
                     fts = list(NULL, NULL, c(7)), tnow = 7, tvec = tvec, nk = c(2,2,1), returnres = T)
fourKcornersSysrelPlot(tvec = tvec, rframe = r13$lower, add = TRUE)
fourKcornersSysrelPlot(tvec = tvec, rframe = r13$upper, add = TRUE)
pbox13f <- sysrelPbox(luckobjlist = list(fc1,fc2), survsign = sys13fsign, kappa = rep(2,3),
                      fts = list(NULL, NULL), tnow = 7, tvec = tvec, nk = c(2,2), returnres = T)
fourKcornersSysrelPlot(tvec = tvec, rframe = r13f$lower, add = TRUE)
fourKcornersSysrelPlot(tvec = tvec, rframe = r13f$upper, add = TRUE)
dev.off()
table(pbox13[,1] == pbox13f[,1]) ; table(pbox13[,2] == pbox13f[,2]) # test succesfull
# additional check whether the min of the 'corner' is the same as the optimized
corner13lower <- r13$lower[-(1:3),]
corner13lower <- apply(corner13lower, 1, min)
table(corner13lower == pbox13[,1]) # due to numeric?
all.equal(corner13lower, pbox13[,1])
# indeed, pbox13 has some non-extreme values for n0 at the end
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------


# system with four components to see whether that problem with the Ccmfs carries over
fc4 <- LuckModel(n0 = c(2,5), y0 = c(failuretolambda(9, 2), failuretolambda(11, 2)))
sys4 <- graph.formula(s -- 2 -- 3:4 -- t, s -- 1 -- t)
V(sys4)$compType <- NA # This just creates the attribute compType
V(sys4)$compType[match(c("1","2","3","4"), V(sys4)$name)] <- "Type 1"
V(sys4)$compType[match(c("s","t"), V(sys4)$name)] <- NA
sys4sign <- computeSystemSurvivalSignature(sys4)

pdf("cmfssys4.pdf", width=12, height=4)
par(mfrow=c(1,3))
fourCornersCcmf(fc4, kappa = 2, n=4, fts=c(10,11), tnow=11, t=12)
fourCornersCcmf(fc4, kappa = 2, n=4, fts=c(10,11), tnow=11, t=15)
fourCornersCcmf(fc4, kappa = 2, n=4, fts=c(10,11), tnow=11, t=20)
par(mfrow=c(1,1))
dev.off()


tvec <- seq(11, 27, by=0.1)
r4 <- fourKcornersSysrel(luckobjlist = list(fc4), survsign = sys4sign, kappa = 2,
                         fts = list(c(10, 11)), tnow = 11, tvec = tvec)
pdf("rsyspbox2.pdf", width=5, height=4)
pbox4 <- sysrelPbox(luckobjlist = list(fc4), survsign = sys4sign, kappa = 2,
                    fts = list(c(10,11)), tnow = 11, tvec = tvec, returnres = T)
fourKcornersSysrelPlot(tvec = tvec, rframe = r4$lower, add = TRUE)
fourKcornersSysrelPlot(tvec = tvec, rframe = r4$upper, add = TRUE)
dev.off()

tveczoom <- seq(14, 16, by=0.02)
r4zoom <- fourKcornersSysrel(luckobjlist = list(fc4), survsign = sys4sign, kappa = 2,
                             fts = list(c(10, 11)), tnow = 11, tvec = tveczoom)
pdf("rsyspbox2zoom.pdf", width=5, height=4)
pbox4zoom <- sysrelPbox(luckobjlist = list(fc4), survsign = sys4sign, kappa = 2,
                        fts = list(c(10,11)), tnow = 11, tvec = tveczoom, returnres = T, ylim = c(0.1,0.3))
fourKcornersSysrelPlot(tvec = tveczoom, rframe = r4zoom$lower, add = TRUE)
fourKcornersSysrelPlot(tvec = tveczoom, rframe = r4zoom$upper, add = TRUE)
dev.off()

###############################################################################
# the three poster examples
###############################################################################
posterprior1 <- LuckModel(n0 = c(2,10), y0 = c(failuretolambda(9,2), failuretolambda(11,2)))
posterprior2 <- LuckModel(n0 = c(8,16), y0 = c(failuretolambda(4,2), failuretolambda(5,2)))
posterprior3 <- LuckModel(n0 = c(1,5), y0 = c(failuretolambda(9,2), failuretolambda(11,2)))
posterprior <- list(posterprior1, posterprior2, posterprior3)

# prior system reliability
posterfts0 <- list(NULL, NULL, NULL)
# full system signature (for prior system reliability)
postersys <- graph.formula(s -- 1:2 -- 3:4:5 -- t)
V(postersys)$compType <- NA # This just creates the attribute compType
V(postersys)$compType[match(c("1","3"), V(postersys)$name)] <- "Type 1"
V(postersys)$compType[match(c("2","4"), V(postersys)$name)] <- "Type 2"
V(postersys)$compType[match(c("5"), V(postersys)$name)] <- "Type 3"
V(postersys)$compType[match(c("s","t"), V(postersys)$name)] <- NA
postersyssign <- computeSystemSurvivalSignature(postersys)

postertvec0 <- seq(0, 10, by=0.02)
postercorners0 <- fourKcornersSysrel(luckobjlist = posterprior, survsign = postersyssign, nk = c(2,2,1),
                                     kappa = rep(2,3), fts = posterfts0, tnow = 0, tvec = postertvec0)
pdf("poster0.pdf", width=5, height=4)
posterpbox0 <- sysrelPbox(luckobjlist = posterprior, survsign = postersyssign,
                          kappa = rep(2,3), fts = posterfts0, tnow = 0, tvec = postertvec0,
                          returnres = T, polygonFillCol = rgb(255,154,0, max =255))
fourKcornersSysrelPlot(tvec = postertvec0, rframe = postercorners0$lower, add = TRUE)
fourKcornersSysrelPlot(tvec = postertvec0, rframe = postercorners0$upper, add = TRUE)
dev.off()


# example 1: failure times as expected: no failures type 1, failures type 2 c(4,5), no failure type 3
posterfts1 <- list(NULL, NULL)
# reduced signature
postersys1 <- graph.formula(s -- 1 -- 2:3 -- t)
V(postersys1)$compType <- NA # This just creates the attribute compType
V(postersys1)$compType[match(c("1", "2"), V(postersys1)$name)] <- "Type 1"
V(postersys1)$compType[match(c("3"), V(postersys1)$name)] <- "Type 2"
#V(postersys1)$compType[match(c("3"), V(postersys1)$name)] <- "Type 3"
V(postersys1)$compType[match(c("s","t"), V(postersys1)$name)] <- NA
postersys1sign <- computeSystemSurvivalSignature(postersys1)

postertvec1 <- seq(5, 15, by=0.02)
postercorners1 <- fourKcornersSysrel(luckobjlist = list(posterprior1, posterprior3), survsign = postersys1sign,
                                     kappa = rep(2,2), fts = posterfts1, tnow = 5, tvec = postertvec1, nk = c(2,1))
pdf("./poster-system/poster1.pdf", width=5, height=4)
par(mar=c(3,3,0,0)+0.1)
posterpbox1 <- sysrelPbox(luckobjlist = list(posterprior1, posterprior3), survsign = postersys1sign,
                          kappa = rep(2,2), fts = posterfts1, tnow = 5, tvec = postertvec1, nk = c(2,1),
                          returnres = T, polygonFillCol = rgb(255,154,0, max =255))
fourKcornersSysrelPlot(tvec = postertvec1, rframe = postercorners1$lower, add = TRUE)
fourKcornersSysrelPlot(tvec = postertvec1, rframe = postercorners1$upper, add = TRUE)
dev.off()

# example 2: early failure times: one failure of type 1 at 3, one failure of type 2 at 1, no failure type 3
posterfts2 <- list(c(3), c(1), NULL)
# reduced signature
postersys2 <- graph.formula(s -- 1 -- 2:3 -- t)
V(postersys2)$compType <- NA # This just creates the attribute compType
V(postersys2)$compType[match(c("1"), V(postersys2)$name)] <- "Type 1"
V(postersys2)$compType[match(c("2"), V(postersys2)$name)] <- "Type 2"
V(postersys2)$compType[match(c("3"), V(postersys2)$name)] <- "Type 3"
V(postersys2)$compType[match(c("s","t"), V(postersys2)$name)] <- NA
postersys2sign <- computeSystemSurvivalSignature(postersys2)

postertvec2 <- seq(3, 13, by=0.02)
postercorners2 <- fourKcornersSysrel(luckobjlist = posterprior, survsign = postersys2sign,
                                     kappa = rep(2,3), fts = posterfts2, tnow = 3, tvec = postertvec2, nk = c(2,2,1))
pdf("./poster-system/poster2.pdf", width=5, height=4)
par(mar=c(3,3,0,0)+0.1)
posterpbox2 <- sysrelPbox(luckobjlist = posterprior, survsign = postersys2sign,
                          kappa = rep(2,3), fts = posterfts2, tnow = 3, tvec = postertvec2, nk = c(2,2,1),
                          returnres = T, polygonFillCol = rgb(255,154,0, max =255))
fourKcornersSysrelPlot(tvec = postertvec2, rframe = postercorners2$lower, add = TRUE)
fourKcornersSysrelPlot(tvec = postertvec2, rframe = postercorners2$upper, add = TRUE)
dev.off()

# example 3: late failure times: no failure of type 1, one failure of type 2 at 10, no failure type 3
posterfts3 <- list(NULL, c(10), NULL)
# reduced signature
postersys3 <- graph.formula(s -- 1:2 -- 3:4 -- t)
V(postersys3)$compType <- NA # This just creates the attribute compType
V(postersys3)$compType[match(c("1", "3"), V(postersys3)$name)] <- "Type 1"
V(postersys3)$compType[match(c("2"), V(postersys3)$name)] <- "Type 2"
V(postersys3)$compType[match(c("4"), V(postersys3)$name)] <- "Type 3"
V(postersys3)$compType[match(c("s","t"), V(postersys3)$name)] <- NA
postersys3sign <- computeSystemSurvivalSignature(postersys3)

postertvec3 <- seq(10, 20, by=0.02)
postercorners3 <- fourKcornersSysrel(luckobjlist = posterprior, survsign = postersys3sign,
                                     kappa = rep(2,3), fts = posterfts3, tnow = 10, tvec = postertvec3, nk = c(2,2,1))
pdf("./poster-system/poster3.pdf", width=5, height=4)
par(mar=c(3,3,0,0)+0.1)
posterpbox3 <- sysrelPbox(luckobjlist = posterprior, survsign = postersys3sign,
                          kappa = rep(2,3), fts = posterfts3, tnow = 10, tvec = postertvec3, nk = c(2,2,1),
                          returnres = T, polygonFillCol = rgb(255,154,0, max =255))
fourKcornersSysrelPlot(tvec = postertvec3, rframe = postercorners3$lower, add = TRUE)
fourKcornersSysrelPlot(tvec = postertvec3, rframe = postercorners3$upper, add = TRUE)
dev.off()

#