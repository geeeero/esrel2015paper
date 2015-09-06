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



luck12 <- WeibullLuckModel(n0 = c(2,5), y0 = c(failuretolambda(8,2), failuretolambda(10,2)),
                             data = WeibullData(1:2))
luck67 <- WeibullLuckModel(n0 = c(2,5), y0 = c(failuretolambda(6,2), failuretolambda(8,2)),
                             data = WeibullData(6:7))
xseq <- seq(0, 200, length.out = 201)

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


#