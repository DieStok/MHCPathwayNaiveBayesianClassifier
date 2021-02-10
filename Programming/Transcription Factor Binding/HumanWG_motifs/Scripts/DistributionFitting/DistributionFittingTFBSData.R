##Fit Distributions to the TFBS data
pacman::p_load(fitdistrplus, DiscreteWeibull, VarianceGamma, SpatialExtremes, gPdtest, spd, ismev, car)

#SpatialExtremes = generalised pareto (among others)

# data("groundbeef", package = "fitdistrplus")
# head(groundbeef)
# nrow(groundbeef)
# str(groundbeef)
# plotdist(groundbeef$serving, histo = TRUE, demp = TRUE)
# ?plotdist
# descdist(groundbeef$serving, boot = 10000)
# 
# fitWeibull <- fitdist(groundbeef$serving, distr = "weibull")
# fitLogNormal <- fitdist(groundbeef$serving, distr = "lnorm")
# fitGamma     <- fitdist(groundbeef$serving, distr = "gamma")
# fitNormal    <- fitdist(groundbeef$serving, distr = "norm")
# fitUniform   <- fitdist(groundbeef$serving, distr = "unif")
# fitExp       <- fitdist(groundbeef$serving, distr = "exp")

toPlot <- list(fitWeibull, fitLogNormal, fitGamma, fitNormal, fitUniform, fitExp)
plot(fitWeibull)
summary(fitWeibull)

par(mfrow = c(2,2))
plotLegend <- c("Weibull", "lognormal", "gamma", "norm", "uniform", "exponential")
fitdistrplus::denscomp(toPlot, legendtext = plotLegend)
fitdistrplus::qqcomp(toPlot, legendtext = plotLegend)
fitdistrplus::cdfcomp(toPlot, legendtext = plotLegend)
fitdistrplus::ppcomp(toPlot, legendtext = plotLegend)

##qqplot: emphasises lack of fit at distribution tails
##ppplot: emphasis lack of fit at distribution center



##now on my data
plotCrossValHistogramAndDensityPlot$plottingData
positiveSetValues <- totalCrossValScoresAllGenes %>% filter(EnsemblId %in% dataStructureTFBS$positiveSetTable$EnsemblId)
negativeSetValues <- totalCrossValScoresAllGenes %>% filter(EnsemblId %ni% dataStructureTFBS$positiveSetTable$EnsemblId)

positiveSetActualScores <- positiveSetValues$totalScore
negativeSetActualScores <- negativeSetValues$totalScore

#this is wrong, because my data is not, in fact, discrete
#fitdist(positiveSetActualScores, "dweibull", start = list(q = 2, beta = 4), method = "mle", discrete = TRUE, lower = c(2,2))

#weibull and gamma
fitVarianceGammaTFBSPositive <- fitdist(positiveSetActualScores, distr = "vg", method = "mle", start = list(vgC = 2.5, sigma = 10, theta = 1, nu = 1 ), lower = c(-Inf, 0, -Inf, 0), control=list(trace=1, REPORT=1))
fitNormalDist                <- fitdist(positiveSetActualScores, distr = "norm")
plot(fitNormalDist)
fitWeibullTFBSPositive <- fitdist(positiveSetActualScores, distr = "weibull")
#not working. So, shift the data?

min(positiveSetActualScores)
min(negativeSetActualScores)
#shift by the largest negative value
positiveSetActualScoresShifted <- positiveSetActualScores + (1 + -(min(min(positiveSetActualScores), min(negativeSetActualScores))))
min(positiveSetActualScoresShifted)
negativeSetActualScoresShifted <- negativeSetActualScores + (1 + -(min(min(positiveSetActualScores), min(negativeSetActualScores))))
min(negativeSetActualScoresShifted)

#get beta distribution compatible values
positiveSetActualScoresShiftedBeta <- positiveSetActualScoresShifted/max(positiveSetActualScoresShifted)

plotdist(positiveSetActualScoresShifted)

plotdist(negativeSetActualScoresShifted)
descdist(positiveSetActualScoresShifted, boot = 1000)
descdist(negativeSetActualScoresShifted, boot = 1000)
fitWBPos <- fitdist(positiveSetActualScoresShifted, distr = "weibull")
fitLNPos <- fitdist(positiveSetActualScoresShifted, distr = "lnorm")
fitGammaPos <- fitdist(positiveSetActualScoresShifted, distr = "gamma")
fitNormPos     <- fitdist(positiveSetActualScoresShifted, distr = "norm")
#unfortunately, fitdist does not get this
fitgpdPos      <- fitdist(positiveSetActualScoresShifted, distr = "gpd", start = list(shape = -1, scale = 39))
#so let us try then, the spd package, which uses a kernel function in the middle (normal dist) and pareto in the tails
fitgpdPos <- spdfit(data = as.matrix(positiveSetActualScoresShifted), kernelfit = "normal", type = "mle", information = "expected")
plot(fitgpdPos)
#strange error messages, see what data(rain) does
data(rain)
hey <- gpd.fit(rain, 10)
hey2 <- gpd.fit(rain, 50)
gpd.diag.result = gpd.diag(hey)
gpd.diag.result2 <- gpd.diag(hey2)
?gpd.fit
#no idea what this all means, too complex

#fitBeta     <- fitdist(positiveSetActualScoresShiftedBeta, distr = "beta", lower = c(0,0))

fitWBNeg <- fitdist(negativeSetActualScoresShifted, distr = "weibull")
fitLNNeg <- fitdist(negativeSetActualScoresShifted, distr = "lnorm")
fitGammaNeg <- fitdist(negativeSetActualScoresShifted, distr = "gamma")
fitNormNeg     <- fitdist(negativeSetActualScoresShifted, distr = "norm")
fitgpdNeg      <- fitdist(negativeSetActualScoresShifted, distr = "gpd", start = list(loc = 0, scale = 1, shape = 0))


#fitUniform  <- fitdist(positiveSetActualScoresShifted, distr = "unif")
#fitExp      <- fitdist(positiveSetActualScoresShifted, distr = "exp")

par(mfrow = c(2,2))
toPlotPositive = list(fitWBPos, fitLNPos, fitGammaPos, fitNormPos, fitgpdPos)
plotLegend     = c("Weibull", "LogNormal", "Gamma", "Norm", "gdp")
fitdistrplus::denscomp(toPlotPositive, legendtext = plotLegend)
fitdistrplus::qqcomp(toPlotPositive, legendtext = plotLegend)
fitdistrplus::cdfcomp(toPlotPositive, legendtext = plotLegend)
fitdistrplus::ppcomp(toPlotPositive, legendtext = plotLegend)

positiveSetActualScores <- sort (positiveSetActualScores)
kernelPositive <- ksmooth(x = seq(1,length(positiveSetActualScores)) ,y = positiveSetActualScores, bandwidth = 2)
plot(positiveSetActualScores)
lines(kernelPositive)
points(positiveSetActualScores)
plot(kernelPositive)
plot(kernelPositive$x, kernelPositive$y)
pizzaBakker   <- lm(kernelPositive$y ~ poly(kernelPositive$y, 5))
?predict



par(mfrow = c(2,2))
toPlotNegative = list(fitWBNeg, fitLNNeg, fitGammaNeg, fitNormNeg, fitgpdNeg )
plotLegend     = c("Weibull", "LogNormal", "Gamma", "Norm", "gdp")
fitdistrplus::denscomp(toPlotNegative, legendtext = plotLegend)
fitdistrplus::qqcomp(toPlotNegative, legendtext = plotLegend)
fitdistrplus::cdfcomp(toPlotNegative, legendtext = plotLegend)
fitdistrplus::ppcomp(toPlotNegative, legendtext = plotLegend)

maakDfMetCriteria <- function(lijst) {
  
  banaan = data.frame()
  for(i in lijst) {
    
    banaan = rbind(banaan, c(summary(i)[6], summary(i)[7], summary(i)[8]))
  }
  banaan
}
freek <- maakDfMetCriteria(toPlotNegative)
freekPositief <- maakDfMetCriteria(toPlotPositive)
summary(fitWBNeg)
summary(fitLNNeg)
summary(fitGammaNeg)
summary(fitNormNeg)
#find initial parameters
plotdist(positiveSetActualScores, distr = "vg", para = list(vgC = 2.5, sigma = 10, theta = 1, nu = 2 ))

#
fitNorm <- fitdistrplus::fitdist(negativeSetActualScores, distr = 'norm')
plot(fitNorm)
