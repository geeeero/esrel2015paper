# --------------------------------------------------------------------------- #
# ESREL 2015 paper - code for talk
# --------------------------------------------------------------------------- #

source("figure1.R")

tuecyan <- function(alpha=1)
  rgb(0.000,0.635,0.871,alpha=alpha)
tuered  <- function(alpha=1)
  rgb(0.839,0.000,0.290,alpha=alpha)
tueyellow <- function(alpha=1)
  rgb(1,154/255,0, alpha)
tuegreen <- function(alpha=1)
  rgb(0,172/255,130/255, alpha)

xseq <- seq(0.1, 200, length.out = 300)

pdf("talk-esrel/pdc1.pdf",width=5,height=3)
par(mar=c(3,3,0,0)+0.1)
plot(xseq, cdinvgamma(xseq, n0, y0), type="l", ylim = c(0,0.024), xlab="", ylab="", col = tuered(), lwd = 2)
mtext(expression(lambda), 1, 2)
mtext(expression(f(lambda)), 2, 2)
points(y0, 0, pch = 17, col = tuered())
#legend("topright", legend = c("prior", "posterior"), lty= 1, col = c(tuered, tuecyan), inset = 0.02)
dev.off()

pdf("talk-esrel/pdc2.pdf",width=5,height=3)
par(mar=c(3,3,0,0)+0.1)
plot(xseq, cdinvgamma(xseq, n0, y0), type="l", ylim = c(0,0.024), xlab="", ylab="", col = tuered(), lwd = 2)
mtext(expression(lambda), 1, 2)
mtext(expression(f(lambda)), 2, 2)
points(y0, 0, pch = 17, col = tuered())
points(2.5, 0, pch = 17)
#legend("topright", legend = c("prior", "posterior"), lty= 1, col = c(tuered, tuecyan), inset = 0.02)
dev.off()

pdf("talk-esrel/pdc3.pdf",width=5,height=3)
par(mar=c(3,3,0,0)+0.1)
plot(xseq, cdinvgamma(xseq, n0, y0), type="l", ylim = c(0,0.024), xlab="", ylab="", col = tuered(), lwd = 2)
mtext(expression(lambda), 1, 2)
mtext(expression(f(lambda)), 2, 2)
points(y0, 0, pch = 17, col = tuered())
points(2.5, 0, pch = 17)
lines(xseq, cdinvgamma(xseq, n2, y2), col = tuecyan(), lwd = 2)
points(y2, 0, pch = 17, col = tuecyan())
#legend("topright", legend = c("prior", "posterior"), lty= 1, col = c(tuered, tuecyan), inset = 0.02)
dev.off()


# almost same posterior for no-pdc example
y0a <- failuretolambda(7,k); y0a
invgammasd(n2, y0a)
t1a <- 6
t2a <- 7
(t1a^k + t2a^k)/2

y2a <- (n0*y0a + t1a^k + t2a^k)/n2; y2a
lambdatofailure(y2a)
invgammasd(n2, y2a)


pdf("talk-esrel/nopdc1.pdf",width=5,height=3)
par(mar=c(3,3,0,0)+0.1)
plot(xseq, cdinvgamma(xseq, n0, y0a), type="l", ylim = c(0,0.024), xlab="", ylab="", col = tuered(), lwd = 2)
mtext(expression(lambda), 1, 2)
mtext(expression(f(lambda)), 2, 2)
points(y0a, 0, pch = 17, col = tuered())
#legend("topright", legend = c("prior", "posterior"), lty= 1, col = c(tuered, tuecyan), inset = 0.02)
dev.off()

pdf("talk-esrel/nopdc2.pdf",width=5,height=3)
par(mar=c(3,3,0,0)+0.1)
plot(xseq, cdinvgamma(xseq, n0, y0a), type="l", ylim = c(0,0.024), xlab="", ylab="", col = tuered(), lwd = 2)
mtext(expression(lambda), 1, 2)
mtext(expression(f(lambda)), 2, 2)
points(y0a, 0, pch = 17, col = tuered())
points(42.5, 0, pch = 17)
#legend("topright", legend = c("prior", "posterior"), lty= 1, col = c(tuered, tuecyan), inset = 0.02)
dev.off()

pdf("talk-esrel/nopdc3.pdf",width=5,height=3)
par(mar=c(3,3,0,0)+0.1)
plot(xseq, cdinvgamma(xseq, n0, y0a), type="l", ylim = c(0,0.024), xlab="", ylab="", col = tuered(), lwd = 2)
mtext(expression(lambda), 1, 2)
mtext(expression(f(lambda)), 2, 2)
points(y0a, 0, pch = 17, col = tuered())
points(42.5, 0, pch = 17)
lines(xseq, cdinvgamma(xseq, n2, y2a), col = tuecyan(), lwd = 2)
points(y2a, 0, pch = 17, col = tuecyan())
#legend("topright", legend = c("prior", "posterior"), lty= 1, col = c(tuered, tuecyan), inset = 0.02)
dev.off()


# --------------------------------------------------------------------------- #
# code for pdc - no pdc with luck
# --------------------------------------------------------------------------- #

library(luck)
library(actuar)
source("00-06_cdfplotLuckModel.r")
source("03-01_WeibullData.r")
source("03-02_Weibull.r")

luckrect <- function(luck){
  n <- updateLuckN(n0(luck), n(data(luck)))
  yl <- min(updateLuckY(n0(luck), y0(luck)[1], tau(data(luck)), n(data(luck))))
  yu <- max(updateLuckY(n0(luck), y0(luck)[2], tau(data(luck)), n(data(luck))))
  list(n = n, y = c(yl, yu))
}

varlambda <- function(luck, posterior = FALSE){
  optimfu <- function(.n0y0, luck){
    if(posterior){
      ny <- c(updateLuckN(.n0y0[1], n(data(luck))),
              updateLuckY(.n0y0[1], .n0y0[2], tau(data(luck)), n(data(luck))))
    } else {
      ny <- .n0y0
    }
    ny[2]^2/(1 - 1/ny[1])
  }
  lres <- wrapOptim(par = c(mean(n0(luck)), mean(y0(luck))), fn = optimfu, control = list(fnscale = 1),
                    lower = c(n0(luck)[1], y0(luck)[1]), upper = c(n0(luck)[2], y0(luck)[2]), luck = luck)
  ures <- wrapOptim(par = c(mean(n0(luck)), mean(y0(luck))), fn = optimfu, control = list(fnscale = -1),
                    lower = c(n0(luck)[1], y0(luck)[1]), upper = c(n0(luck)[2], y0(luck)[2]), luck = luck)
  c(lres$value, ures$value)
}

luck12 <- WeibullLuckModel(n0 = c(2,5), y0 = c(failuretolambda(8,2), failuretolambda(10,2)),
                             data = WeibullData(1:2))
luck67 <- WeibullLuckModel(n0 = c(2,5), y0 = c(failuretolambda(6,2), failuretolambda(8,2)),
                             data = WeibullData(6:7))
xseq <- seq(0, 200, length.out = 201)

y0(luck67)
sqrt(varlambda(luck67))
luckrect(luck67)
lambdatofailure(luckrect(luck67)$y)
sqrt(varlambda(luck67, posterior=T))

y0(luck12)
sqrt(varlambda(luck12))
luckrect(luck12)
lambdatofailure(luckrect(luck12)$y)
sqrt(varlambda(luck12, posterior=T))

pdf("talk-esrel/lucknopdc1.pdf",width=5,height=3)
par(mar=c(3,3,0,0)+0.1)
cdfplot(luck67, xvec = xseq, vertdist = TRUE,
        control = controlList(posterior = FALSE, borderCol = tuered(), polygonCol = tuered(0.4)))
mtext(expression(lambda), 1, 2)
mtext(expression(F(lambda)), 2, 2)
lines(y0(luck67), rep(0.01,2), lwd = 2, lend = 2, col = tuered())
dev.off()

pdf("talk-esrel/lucknopdc2.pdf",width=5,height=3)
par(mar=c(3,3,0,0)+0.1)
cdfplot(luck67, xvec = xseq, vertdist = TRUE,
        control = controlList(posterior = FALSE, borderCol = tuered(), polygonCol = tuered(0.4)))
lines(y0(luck67), rep(0.01,2), lwd = 2, lend = 2, col = tuered())
points(42.5, 0, pch = 17)
mtext(expression(lambda), 1, 2)
mtext(expression(F(lambda)), 2, 2)
dev.off()

pdf("talk-esrel/lucknopdc3.pdf",width=5,height=3)
par(mar=c(3,3,0,0)+0.1)
cdfplot(luck67, xvec = xseq, vertdist = TRUE,
        control = controlList(posterior = FALSE, borderCol = tuered(), polygonCol = tuered(0.4)))
lines(y0(luck67), rep(0.01,2), lwd = 2, lend = 2, col = tuered())
points(42.5, 0, pch = 17)
cdfplot(luck67, xvec = xseq, add = TRUE, vertdist = TRUE,
        control = controlList(posterior = TRUE, borderCol = tuecyan(), polygonCol = tuecyan(0.8)))
lines(luckrect(luck67)$y, rep(-0.01,2), lwd = 2, lend = 2, col = tuecyan())
mtext(expression(lambda), 1, 2)
mtext(expression(F(lambda)), 2, 2)
dev.off()

pdf("talk-esrel/luckpdc1.pdf",width=5,height=3)
par(mar=c(3,3,0,0)+0.1)
cdfplot(luck12, xvec = xseq, vertdist = TRUE,
        control = controlList(posterior = FALSE, borderCol = tuered(), polygonCol = tuered(0.4)))
mtext(expression(lambda), 1, 2)
mtext(expression(F(lambda)), 2, 2)
lines(y0(luck12), rep(0.01,2), lwd = 2, lend = 2, col = tuered())
dev.off()

pdf("talk-esrel/luckpdc2.pdf",width=5,height=3)
par(mar=c(3,3,0,0)+0.1)
cdfplot(luck12, xvec = xseq, vertdist = TRUE,
        control = controlList(posterior = FALSE, borderCol = tuered(), polygonCol = tuered(0.4)))
lines(y0(luck12), rep(0.01,2), lwd = 2, lend = 2, col = tuered())
points(2.5, 0, pch = 17)
mtext(expression(lambda), 1, 2)
mtext(expression(F(lambda)), 2, 2)
dev.off()

pdf("talk-esrel/luckpdc3.pdf",width=5,height=3)
par(mar=c(3,3,0,0)+0.1)
cdfplot(luck12, xvec = xseq, vertdist = TRUE,
        control = controlList(posterior = FALSE, borderCol = tuered(), polygonCol = tuered(0.4)))
lines(y0(luck12), rep(0.01,2), lwd = 2, lend = 2, col = tuered())
points(2.5, 0, pch = 17)
cdfplot(luck12, xvec = xseq, add = TRUE, vertdist = TRUE,
        control = controlList(posterior = TRUE, borderCol = tuecyan(), polygonCol = tuecyan(0.8)))
lines(luckrect(luck12)$y, rep(-0.01,2), lwd = 2, lend = 2, col = tuecyan())
mtext(expression(lambda), 1, 2)
mtext(expression(F(lambda)), 2, 2)
dev.off()

pdf("talk-esrel/pdc-nopdc.pdf",width=5,height=3)
par(mar=c(3,3,0,0)+0.1)
cdfplot(luck12, xvec = xseq, vertdist = TRUE,
        control = controlList(posterior = TRUE, borderCol = tuecyan(), polygonCol = tuecyan(0.8)))
lines(luckrect(luck12)$y, rep(-0.01,2), lwd = 2, lend = 2, col = tuecyan())
cdfplot(luck67, xvec = xseq, add = TRUE, vertdist = TRUE,
        control = controlList(posterior = TRUE, borderCol = tueyellow(), polygonCol = tueyellow(0.8)))
lines(luckrect(luck67)$y, rep( 0.01,2), lwd = 2, lend = 2, col = tueyellow())
mtext(expression(lambda), 1, 2)
mtext(expression(F(lambda)), 2, 2)
dev.off()


# --------------------------------------------------------------------------- #
# code for parallel system graphs
# --------------------------------------------------------------------------- #

source("notes-isipta15poster.R")

esrelpriorluck <- LuckModel(n0 = c(2,5), y0 = c(failuretolambda(9,2), failuretolambda(11,2)))
esrelprior <- list(esrelpriorluck)
esrelfts <- list(NULL)

# signature full parallel system
esrelsys <- graph.formula(s -- 1:3 -- t)
V(esrelsys)$compType <- NA # This just creates the attribute compType
V(esrelsys)$compType[match(c("1", "2", "3"), V(esrelsys)$name)] <- "Type 1"
V(esrelsys)$compType[match(c("s","t"), V(esrelsys)$name)] <- NA
esrelsyssign <- computeSystemSurvivalSignature(esrelsys)

# esrel paper example 1: failure times as expected:
esrel1fts <- list(c(10,11))
# reduced signature
esrel1sys <- graph.formula(s -- 1 -- t)
V(esrel1sys)$compType <- NA # This just creates the attribute compType
V(esrel1sys)$compType[match(c("1"), V(esrel1sys)$name)] <- "Type 1"
V(esrel1sys)$compType[match(c("s","t"), V(esrel1sys)$name)] <- NA
esrel1syssign <- computeSystemSurvivalSignature(esrel1sys)

esreltvec <- seq(0, 70, length.out=101)
esrel1tvec <- seq(11, 70, length.out=101)

esrel1prior <- sysrelPbox(luckobjlist = esrelprior, survsign = esrelsyssign, nk = 3,
                          kappa = 2, fts = esrelfts, tnow = 0, prior = TRUE, tvec = esreltvec,
                          returnres = TRUE, polygonFillCol = "gray")
esrel1corners <- fourKcornersSysrel(luckobjlist = esrelprior, survsign = esrel1syssign,
                                    kappa = 2, fts = esrel1fts, tnow = 11, tvec = esrel1tvec, nk = 1)
#posterpbox1 <- sysrelPbox(luckobjlist = list(posterprior1, posterprior3), survsign = postersys1sign,
#                          kappa = rep(2,2), fts = posterfts1, tnow = 5, tvec = postertvec1, nk = c(2,1),
#                          returnres = T, polygonFillCol = rgb(255,154,0, max =255))

#pdf("./poster-system/poster1prior.pdf", width=5, height=4)
#par(mar=c(3,3,0,0)+0.1)
#plotSysrelPbox(tvec = postertvec1, res = priorpbox1, add = FALSE, ylim = c(0,1), xlab = "t",
#               ylab = expression(R[sys](t)), polygonBorderCol = rgb(255,154,0,160, max =255),
#               polygonFillCol = rgb(255,154,0, 80, max =255))
#plotSysrelPbox(tvec = postertvec1, res = posterpbox1, add = TRUE,
#               polygonBorderCol = NA,
#               polygonFillCol = rgb(16,16,115, 80, max =255))
#fourKcornersSysrelPlot(tvec = postertvec1, rframe = postercorners1$lower, add = TRUE, col = rgb(16,16,115, max =255))
#fourKcornersSysrelPlot(tvec = postertvec1, rframe = postercorners1$upper, add = TRUE, col = rgb(16,16,115, max =255))
#dev.off()

parallelRsys <- function(n0y0, t, fts, tnow, l, m, k = 2){
  sumover <- 0:(l-m)
  signs <- (-1)^(sumover)
  bincoef <- choose(l-m, sumover)
  enumerator <- n0y0[1]*n0y0[2] + (l-m)*(tnow^k) + sum(fts)
  thing <- (enumerator/(enumerator + sumover*(t^k - tnow^k)))^(n0y0[1] + m + 1)
  1 - sum(signs * bincoef * thing)
}

#parallelRsys(n0y0 = c(2,103), t = 10, fts = 0, tnow = 0, l = 3, m = 0, k = 2)

parallelRsysLuck <- function(luck, t, fts, tnow, l, m, k = 2, returnasvec = FALSE){
  parl <- c(n0(luck)[1], y0(luck)[1])
  paru <- c(n0(luck)[2], y0(luck)[2])
  optl <- optim(par = (parl+paru)/2, fn = parallelRsys, method = "L-BFGS-B", lower = parl, upper = paru,
                t = t, fts = fts, tnow = tnow, l = l, m = m, k = k)
  optu <- optim(par = (parl+paru)/2, fn = parallelRsys, method = "L-BFGS-B", lower = parl, upper = paru,
                t = t, fts = fts, tnow = tnow, l = l, m = m, k = k,
                control = list(fnscale = -1))
  if(returnasvec)
    return(c(optl$value, optu$value, optl$par, optu$par))
  else
    return(list(rel = c(optl$value, optu$value), lpar = optl$par, upar = optu$par))
}

parallelRsysPbox <- function(luck, tvec, fts, tnow, l, m, k = 2,
                             add = FALSE, ylim = c(0,1), xlab = "t", ylab = expression(R[sys](t)),
                             polygonBorderCol = NA, polygonFillCol = "grey",
                             returnres = FALSE, noplot = FALSE, ...){
  res <- sapply(tvec, FUN = parallelRsysLuck, luck = luck, fts = fts, tnow = tnow, l = l, m = m, returnasvec = TRUE)
  res <- t(res)
  colnames(res) <- c("lower", "upper", "ln0", "ly0", "un0", "uy0")
  if(noplot){
    add = TRUE
    returnres = TRUE
  }
  if(!add){
    plot(c(min(tvec), max(tvec)), c(0,0), type = "n", ylim = ylim, xlab = "", ylab = "", ...)
    mtext(xlab, 1, 2)
    mtext(ylab, 2, 2)
  }
  tvecrev <- numeric(length(tvec))
  lower <- res[,1]
  upper <- numeric(length(tvec))
  for (t in 1:length(tvec)){
    tvecrev[length(tvec)-t+1] <- tvec[t]
    upper[length(tvec)-t+1] <- res[t,2]
  }
  if(!noplot)
    polygon(c(tvec, tvecrev), c(lower, upper), border = polygonBorderCol, col = polygonFillCol)
  if(returnres)
    return(res)
}

plotparallelRsysPbox <- function(tvec, res, add = FALSE, ylim = c(0,1), xlab = "t", ylab = expression(R[sys](t)),
                                 polygonBorderCol = NA, polygonFillCol = "grey", ...){
  if(!add){
    plot(c(min(tvec), max(tvec)), c(0,0), type = "n", ylim = ylim, xlab = "", ylab = "", ...)
    mtext(xlab, 1, 2)
    mtext(ylab, 2, 2)
  }
  tvecrev <- numeric(length(tvec))
  lower <- res[,1]
  upper <- numeric(length(tvec))
  for (t in 1:length(tvec)){
    tvecrev[length(tvec)-t+1] <- tvec[t]
    upper[length(tvec)-t+1] <- res[t,2]
  }
  polygon(c(tvec, tvecrev), c(lower, upper), border = polygonBorderCol, col = polygonFillCol)
}



esrelpriorluck <- LuckModel(n0 = c(2,5), y0 = c(failuretolambda(9,2), failuretolambda(11,2)))
priortvec <- seq(0, 70, by = 0.1)
post1tvec <- seq(11, 70, by = 0.1)
post2tvec <- seq(2, 70, by = 0.1)
post3tvec <- seq(15, 70, by = 0.1)

esrelpriorpbox <- parallelRsysPbox(luck = esrelpriorluck, tvec = priortvec, fts = 0, tnow = 0, l = 3, m = 0,
                                   polygonBorderCol = NA, polygonFillCol = "grey",
                                   returnres = TRUE, noplot = FALSE)
  
esrelpost1pbox <- parallelRsysPbox(luck = esrelpriorluck, tvec = post1tvec, fts = c(10, 11), tnow = 11,
                                   l = 3, m = 2, polygonBorderCol = NA, polygonFillCol = "red",
                                   returnres = TRUE, noplot = FALSE, add = TRUE)

esrelpost2pbox <- parallelRsysPbox(luck = esrelpriorluck, tvec = post2tvec, fts = c(1, 2), tnow = 2,
                                   l = 3, m = 2, polygonBorderCol = NA, polygonFillCol = "red",
                                   returnres = TRUE, noplot = FALSE, add = TRUE)

esrelpost3pbox <- parallelRsysPbox(luck = esrelpriorluck, tvec = post3tvec, fts = 0, tnow = 15,
                                   l = 3, m = 0, polygonBorderCol = NA, polygonFillCol = "red",
                                   returnres = TRUE, noplot = FALSE, add = TRUE)

plotparallelRsysPbox(tvec = priortvec, res = esrelpriorpbox,
                     polygonBorderCol = tuered(), polygonFillCol = tuered(0.4))
plotparallelRsysPbox(tvec = post1tvec, res = esrelpost1pbox, 
                     polygonBorderCol = tuecyan(), polygonFillCol = tuecyan(0.8), add = TRUE)

pdf("talk-esrel/sysex1-1.pdf",width=5,height=3)
par(mar=c(3,3,0,0)+0.1)
plotparallelRsysPbox(tvec = priortvec, res = esrelpriorpbox,
                     polygonBorderCol = tuered(), polygonFillCol = tuered(0.4))
dev.off()

pdf("talk-esrel/sysex1-2.pdf",width=5,height=3)
par(mar=c(3,3,0,0)+0.1)
plotparallelRsysPbox(tvec = priortvec, res = esrelpriorpbox,
                     polygonBorderCol = tuered(), polygonFillCol = tuered(0.4))
plotparallelRsysPbox(tvec = post1tvec, res = esrelpost1pbox, 
                     polygonBorderCol = tuecyan(), polygonFillCol = tuecyan(0.8), add = TRUE)
dev.off()

pdf("talk-esrel/sysex1-3.pdf",width=5,height=3)
par(mar=c(3,3,0,0)+0.1)
plotparallelRsysPbox(tvec = priortvec, res = esrelpriorpbox,
                     polygonBorderCol = tuered(), polygonFillCol = tuered(0.4))
plotparallelRsysPbox(tvec = post1tvec, res = esrelpost1pbox, 
                     polygonBorderCol = tuecyan(), polygonFillCol = tuecyan(0.8), add = TRUE)
plotparallelRsysPbox(tvec = post2tvec, res = esrelpost2pbox, 
                     polygonBorderCol = tueyellow(), polygonFillCol = tueyellow(0.8), add = TRUE)
dev.off()

pdf("talk-esrel/sysex1-4.pdf",width=5,height=3)
par(mar=c(3,3,0,0)+0.1)
plotparallelRsysPbox(tvec = priortvec, res = esrelpriorpbox,
                     polygonBorderCol = tuered(), polygonFillCol = tuered(0.4))
plotparallelRsysPbox(tvec = post1tvec, res = esrelpost1pbox, 
                     polygonBorderCol = tuecyan(), polygonFillCol = tuecyan(0.8), add = TRUE)
plotparallelRsysPbox(tvec = post2tvec, res = esrelpost2pbox, 
                     polygonBorderCol = tueyellow(), polygonFillCol = tueyellow(0.8), add = TRUE)
plotparallelRsysPbox(tvec = post3tvec, res = esrelpost3pbox, 
                     polygonBorderCol = tuegreen(), polygonFillCol = tuegreen(0.6), add = TRUE)
dev.off()
#