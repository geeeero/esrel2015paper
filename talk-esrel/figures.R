# --------------------------------------------------------------------------- #
# ESREL 2015 paper - code for talk
# --------------------------------------------------------------------------- #

source("figure1.R")

tuecyan <- rgb(0.000,0.635,0.871)
tuered  <- rgb(0.839,0.000,0.290)

xseq <- seq(0.1, 200, length.out = 300)

pdf("talk-esrel/pdc1.pdf",width=5,height=3)
par(mar=c(3,3,0,0)+0.1)
plot(xseq, cdinvgamma(xseq, n0, y0), type="l", ylim = c(0,0.024), xlab="", ylab="", col = tuered, lwd = 2)
mtext(expression(lambda), 1, 2)
mtext(expression(f(lambda)), 2, 2)
points(y0, 0, pch = 17, col = tuered)
#legend("topright", legend = c("prior", "posterior"), lty= 1, col = c(tuered, tuecyan), inset = 0.02)
dev.off()

pdf("talk-esrel/pdc2.pdf",width=5,height=3)
par(mar=c(3,3,0,0)+0.1)
plot(xseq, cdinvgamma(xseq, n0, y0), type="l", ylim = c(0,0.024), xlab="", ylab="", col = tuered, lwd = 2)
mtext(expression(lambda), 1, 2)
mtext(expression(f(lambda)), 2, 2)
points(y0, 0, pch = 17, col = tuered)
points(2.5, 0, pch = 17)
#legend("topright", legend = c("prior", "posterior"), lty= 1, col = c(tuered, tuecyan), inset = 0.02)
dev.off()

pdf("talk-esrel/pdc3.pdf",width=5,height=3)
par(mar=c(3,3,0,0)+0.1)
plot(xseq, cdinvgamma(xseq, n0, y0), type="l", ylim = c(0,0.024), xlab="", ylab="", col = tuered, lwd = 2)
mtext(expression(lambda), 1, 2)
mtext(expression(f(lambda)), 2, 2)
points(y0, 0, pch = 17, col = tuered)
points(2.5, 0, pch = 17)
lines(xseq, cdinvgamma(xseq, n2, y2), col = tuecyan, lwd = 2)
points(y2, 0, pch = 17, col = tuecyan)
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
plot(xseq, cdinvgamma(xseq, n0, y0a), type="l", ylim = c(0,0.024), xlab="", ylab="", col = tuered, lwd = 2)
mtext(expression(lambda), 1, 2)
mtext(expression(f(lambda)), 2, 2)
points(y0a, 0, pch = 17, col = tuered)
#legend("topright", legend = c("prior", "posterior"), lty= 1, col = c(tuered, tuecyan), inset = 0.02)
dev.off()

pdf("talk-esrel/nopdc2.pdf",width=5,height=3)
par(mar=c(3,3,0,0)+0.1)
plot(xseq, cdinvgamma(xseq, n0, y0a), type="l", ylim = c(0,0.024), xlab="", ylab="", col = tuered, lwd = 2)
mtext(expression(lambda), 1, 2)
mtext(expression(f(lambda)), 2, 2)
points(y0a, 0, pch = 17, col = tuered)
points(42.5, 0, pch = 17)
#legend("topright", legend = c("prior", "posterior"), lty= 1, col = c(tuered, tuecyan), inset = 0.02)
dev.off()

pdf("talk-esrel/nopdc3.pdf",width=5,height=3)
par(mar=c(3,3,0,0)+0.1)
plot(xseq, cdinvgamma(xseq, n0, y0a), type="l", ylim = c(0,0.024), xlab="", ylab="", col = tuered, lwd = 2)
mtext(expression(lambda), 1, 2)
mtext(expression(f(lambda)), 2, 2)
points(y0a, 0, pch = 17, col = tuered)
points(42.5, 0, pch = 17)
lines(xseq, cdinvgamma(xseq, n2, y2a), col = tuecyan, lwd = 2)
points(y2a, 0, pch = 17, col = tuecyan)
#legend("topright", legend = c("prior", "posterior"), lty= 1, col = c(tuered, tuecyan), inset = 0.02)
dev.off()

#