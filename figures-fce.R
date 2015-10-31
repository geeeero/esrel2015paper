###############################################################################
#                                plots for FCE talk                           #
###############################################################################

source("notes-isipta15poster.R")

tuecyan <- function(alpha=1)
  rgb(0.000,0.635,0.871,alpha=alpha)
tuered  <- function(alpha=1)
  rgb(0.839,0.000,0.290,alpha=alpha)
tueyellow <- function(alpha=1)
  rgb(1,154/255,0, alpha)
tuegreen <- function(alpha=1)
  rgb(0,172/255,130/255, alpha)


fceprior1 <- LuckModel(n0 = c(2,10), y0 = c(failuretolambda(9,2), failuretolambda(11,2)))
fceprior2 <- LuckModel(n0 = c(8,16), y0 = c(failuretolambda(4,2), failuretolambda(5,2)))
fceprior3 <- LuckModel(n0 = c(1,5), y0 = c(failuretolambda(9,2), failuretolambda(11,2)))
fceprior <- list(fceprior1, fceprior2, fceprior3)

# prior system reliability at time 0
fcetvec0 <- seq(0, 30, length.out=301)
fcefts0 <- list(NULL, NULL, NULL)
# full system signature (for prior system reliability)
fcesys0 <- graph.formula(s -- 1:2:3 -- 4:5:6 -- t)
V(fcesys0)$compType <- NA # This just creates the attribute compType
V(fcesys0)$compType[match(c("1","4"), V(fcesys0)$name)] <- "Type 1"
V(fcesys0)$compType[match(c("2","3","5"), V(fcesys0)$name)] <- "Type 2"
V(fcesys0)$compType[match(c("6"), V(fcesys0)$name)] <- "Type 3"
V(fcesys0)$compType[match(c("s","t"), V(fcesys0)$name)] <- NA
fcesyssign0 <- computeSystemSurvivalSignature(fcesys0)
fcecorners0 <- fourKcornersSysrel(luckobjlist = fceprior, survsign = fcesyssign0, nk = c(2,2,1),
                                  kappa = rep(2,3), fts = fcefts0, tnow = 0, prior = TRUE, tvec = fcetvec0)
fcepbox0 <- sysrelPbox(luckobjlist = fceprior, survsign = fcesyssign, nk = c(2,2,1),
                       kappa = rep(2,3), fts = fcefts0, tnow = 0, prior = TRUE, tvec =fcetvec0, noplot = TRUE)

# posterior 1: failure times as expected: no failures type 1, failures type 2 c(4,5), no failure type 3
fcetvec1 <- seq(5, 30, length.out=251)
fcefts1 <- list(NULL, c(4,5), NULL)
# reduced signature
fcesys1 <- graph.formula(s -- 1:2 -- 3:4 -- t)
V(fcesys1)$compType <- NA # This just creates the attribute compType
V(fcesys1)$compType[match(c("1", "3"), V(fcesys1)$name)] <- "Type 1"
V(fcesys1)$compType[match(c("2"), V(fcesys1)$name)] <- "Type 2"
V(fcesys1)$compType[match(c("4"), V(fcesys1)$name)] <- "Type 3"
V(fcesys1)$compType[match(c("s","t"), V(fcesys1)$name)] <- NA
fcesyssign1 <- computeSystemSurvivalSignature(fcesys1)
fcecorners1 <- fourKcornersSysrel(luckobjlist = fceprior, survsign = fcesyssign1, nk = c(2,1,1),
                                  kappa = rep(2,3), fts = fcefts1, tnow = 5, prior = TRUE, tvec = fcetvec1)
fcepbox1 <- sysrelPbox(luckobjlist = fceprior, survsign = fcesyssign1, nk = c(2,1,1),
                       kappa = rep(2,3), fts = fcefts1, tnow = 5, prior = TRUE, tvec = fcetvec1, noplot = TRUE)


# ----------- plots ----------------------------------

pdf("./talk-esrel/fcesysex-0.pdf", width=5,height=3)
par(mar=c(3,3,0,0)+0.1)
plotSysrelPbox(tvec = fcetvec0, res = fcepbox0, add = FALSE, ylim = c(0,1), xlab = "t",
               ylab = expression(R[sys](t)), polygonBorderCol = tuered(), polygonFillCol = tuered(0.4))
fourKcornersSysrelPlot(tvec = fcetvec0, rframe = fcecorners0$lower, col = tuered(), add = TRUE)
fourKcornersSysrelPlot(tvec = fcetvec0, rframe = fcecorners0$upper, col = tuered(), add = TRUE)
dev.off()

pdf("./talk-esrel/fcesysex-1.pdf", width=5,height=3)
par(mar=c(3,3,0,0)+0.1)
plotSysrelPbox(tvec = fcetvec0, res = fcepbox0, add = FALSE, ylim = c(0,1), xlab = "t",
               ylab = expression(R[sys](t)), polygonBorderCol = tuered(), polygonFillCol = tuered(0.4))
fourKcornersSysrelPlot(tvec = fcetvec0, rframe = fcecorners0$lower, col = tuered(), add = TRUE)
fourKcornersSysrelPlot(tvec = fcetvec0, rframe = fcecorners0$upper, col = tuered(), add = TRUE)
plotSysrelPbox(tvec = fcetvec1, res = fcepbox1, add = TRUE, polygonBorderCol = tuecyan(), polygonFillCol = tuecyan(0.8))
fourKcornersSysrelPlot(tvec = fcetvec1, rframe = fcecorners1$lower, col = tuecyan(), add = TRUE)
fourKcornersSysrelPlot(tvec = fcetvec1, rframe = fcecorners1$upper, col = tuecyan(), add = TRUE)
dev.off()

#######

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