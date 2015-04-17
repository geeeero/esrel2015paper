# --------------------------------------------------------------------------- #
# ESREL 2015 paper - code for Figs 2 & 3
# --------------------------------------------------------------------------- #

source("figure1.R")

forfig23 <- list(n0 = c(2,5), y0 = c(failuretolambda(9,2),failuretolambda(11,2)), k = 2)
datafig2 <- c(1,2)
datafig3 <- c(10,11)

xseq <- seq(0.1, 300, length.out = 300)
xseq <- seq(120, 140, length.out = 300) # for lower y0
xseq <- seq(180, 200, length.out = 300) # for upper y0

plot( xseq, cpinvgamma(xseq, forfig23$n0[1], forfig23$y0[1]), type="l", xlab="", ylab="")
lines(xseq, cpinvgamma(xseq, forfig23$n0[2], forfig23$y0[1]))
#plot( xseq, cpinvgamma(xseq, forfig23$n0[1], forfig23$y0[2]), type="l", xlab="", ylab="")
lines(xseq, cpinvgamma(xseq, forfig23$n0[1], forfig23$y0[2]))
lines(xseq, cpinvgamma(xseq, forfig23$n0[2], forfig23$y0[2]))

lines(xseq, cpinvgamma(xseq, 2.5, forfig23$y0[1]), lty = 2)
lines(xseq, cpinvgamma(xseq, 3.0, forfig23$y0[1]), lty = 2)
lines(xseq, cpinvgamma(xseq, 3.5, forfig23$y0[1]), lty = 2)
lines(xseq, cpinvgamma(xseq, 4.0, forfig23$y0[1]), lty = 2)
lines(xseq, cpinvgamma(xseq, 4.5, forfig23$y0[1]), lty = 2)

lines(xseq, cpinvgamma(xseq, 2.5, forfig23$y0[2]), lty = 2)
lines(xseq, cpinvgamma(xseq, 3.0, forfig23$y0[2]), lty = 2)
lines(xseq, cpinvgamma(xseq, 3.5, forfig23$y0[2]), lty = 2)
lines(xseq, cpinvgamma(xseq, 4.0, forfig23$y0[2]), lty = 2)
lines(xseq, cpinvgamma(xseq, 4.5, forfig23$y0[2]), lty = 2)

# The problem:
xseq <- seq(165, 215, length.out = 300) # for upper y0
plot( xseq, cpinvgamma(xseq, forfig23$n0[1], forfig23$y0[2]), type="l", xlab="", ylab="")
lines(xseq, cpinvgamma(xseq, 2.5, forfig23$y0[2]), col = 2)
lines(xseq, cpinvgamma(xseq, 3.0, forfig23$y0[2]), col = 3)
lines(xseq, cpinvgamma(xseq, 3.5, forfig23$y0[2]), col = 4)
lines(xseq, cpinvgamma(xseq, 4.0, forfig23$y0[2]), col = 5)
lines(xseq, cpinvgamma(xseq, 4.5, forfig23$y0[2]), col = 6)
lines(xseq, cpinvgamma(xseq, forfig23$n0[2], forfig23$y0[2]), col = 7)
legend("topleft", legend = seq(2,5,by=0.5), lty=1, col=1:7)

# ny = c(n0, y0)
pigforoptim <- function(ny, t, ...)
  cpinvgamma(t, ny[1], ny[2], ...)




# cpinvgamma(x, n0, y0)


setEPS()
postscript("fig2.eps",width=5,height=3)
par(mar=c(3,3,1,1)+0.1)

dev.off()

#