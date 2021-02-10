positiveSetValues <- totalCrossValScoresAllGenes %>% filter(EnsemblId %in% dataStructureTFBS$positiveSetTable$EnsemblId)
negativeSetValues <- totalCrossValScoresAllGenes %>% filter(EnsemblId %ni% dataStructureTFBS$positiveSetTable$EnsemblId)

positiveSetActualScores <- positiveSetValues$totalScore
negativeSetActualScores <- negativeSetValues$totalScore

positiveSetActualScores
negativeSetActualScores
densityPositive <- density(positiveSetActualScores, bw = 2.5)
densityNegative <- density(negativeSetActualScores, bw = 2.5)
plot(densityPositive)
plot(densityNegative)

ksmoothpositive <- ksmooth(densityPositive$x, densityPositive$y)


?ksmooth
lo <- loess(densityPositive$y ~ densityPositive$x, family = "gaussian")
dfPositive <- data.frame(x = densityPositive$x, y = densityPositive$y)
polyP <- lm(y ~ poly(x, 20), data = dfPositive)

dfNegative <- data.frame(x=densityNegative$x, y=densityNegative$y)
polyN <- lm(y ~ poly(x, 20), data = dfNegative)
?loess
plot(densityPositive$x, densityPositive$y)
lines(predict(lo))
plot(predict(lo))
plot(predict(poly))
?loess


plot(densityPositive$y ~ densityPositive$x)
lines(predict(polyP) ~ densityPositive$x)

plot(densityNegative$y ~ densityNegative$x)
lines(predict(polyN) ~ densityNegative$x)


polyFunction <- function(x, polyModel) {
  
  k <- polyModel$coefficients
  print(k)
  print(k[21])
  print(x^20 * k[21])
  
  func <- k[1] + x * k[2] + x^2 * k[3] + x^3 * k[4] + x^4 * k[5] + x^5 * k[6] +
    x^6 * k[7] + x^7 * k[8] + x^8 * k[9] + x^9 * k[10] + x^10 * k[11] +
    x^11 * k[12] + x^12 * k[13] + x^13 * k[14] + x^14 * k[15] +
    x^15 * k[16] + x^16 * k[17] + x^17 * k[18] + x^18 * k[19] +
    x^19 * k[20] + x^20 * k[21]
  print(func)
  
  func
  
}

freek <- polyFunction(3.5484762252, polyP)
freek
positivePolyScoresAllGenes <- predict(polyP, data.frame(x = totalCrossValScoresAllGenes$totalScore))
negativePolyScoresAllGenes <- predict(polyN, data.frame(x = totalCrossValScoresAllGenes$totalScore))
totalLog2ScoresPoly        <- log2(positivePolyScoresAllGenes/negativePolyScoresAllGenes)
totalLog2ScoresPoly
dfTotalLog2Scores <- data.frame(EnsemblId = totalCrossValScoresAllGenes$EnsemblId,
                                totalScore = totalCrossValScoresAllGenes$totalScore,
                                Log2Scores = totalLog2ScoresPoly)

dfTotalLog2Scores <- dfTotalLog2Scores %>% arrange(totalScore)
dfTotalLog2Scores[dfTotalLog2Scores$Log2Scores > 10,]


tail(dfTotalLog2Scores)
#totalPlot
plot(densityNegative$y ~ densityNegative$x, col = "red", cex = 0.2)
lines(predict(polyN) ~ densityNegative$x, col = "darkred")
points(densityPositive$y ~ densityPositive$x, col = "green", cex = 0.2)
lines(predict(polyP) ~ densityPositive$x, col = "darkgreen")
plot(Log2Scores ~ totalScore, data = dfTotalLog2Scores, col = "black")

dfTotalLog2Scores[dfTotalLog2Scores$totalScore >= 3.95 & dfTotalLog2Scores$totalScore <= 4.05,]
dfTotalLog2Scores[dfTotalLog2Scores$totalScore >= 5.95 & dfTotalLog2Scores$totalScore <= 6.05,]
dfTotalLog2Scores[dfTotalLog2Scores$totalScore >= 7.95 & dfTotalLog2Scores$totalScore <= 8.05,]

#find the local maximum that I will use as the max and min values
head(dfTotalLog2Scores, 30)
#minimal value = -3.363383, from ENSG00000167822, totalScore = -10.83280
tail(dfTotalLog2Scores, 50)
#maximal value = 2.406856, from ENSG00000205659, totalScore = 13.22579

#edit the scores
dfTotalLog2ScoresBounded <- dfTotalLog2Scores
minimalBoundGene <- dfTotalLog2ScoresBounded[dfTotalLog2ScoresBounded$EnsemblId == "ENSG00000167822",]
maximalBoundGene <- dfTotalLog2ScoresBounded[dfTotalLog2ScoresBounded$EnsemblId == "ENSG00000205659",]
dfTotalLog2ScoresBounded[dfTotalLog2ScoresBounded$totalScore < minimalBoundGene$totalScore, ]$Log2Scores <- minimalBoundGene$Log2Scores
dfTotalLog2ScoresBounded[dfTotalLog2ScoresBounded$totalScore > maximalBoundGene$totalScore, ]$Log2Scores <- maximalBoundGene$Log2Scores

plot(densityNegative$y ~ densityNegative$x, col = "red", cex = 0.2, main = "score density P and N sets")
lines(predict(polyN) ~ densityNegative$x, col = "darkred")
points(densityPositive$y ~ densityPositive$x, col = "green", cex = 0.2)
lines(predict(polyP) ~ densityPositive$x, col = "darkgreen")
plot(Log2Scores ~ totalScore, data = dfTotalLog2Scores, col = "black", main = "unbounded Log2Scores")
plot(Log2Scores ~ totalScore, data = dfTotalLog2ScoresBounded, col = "black", main = "bounded Log2Scores")

orderedLogScores <- dfTotalLog2ScoresBounded %>% arrange(EnsemblId)
head(orderedLogScores)
colnames(orderedLogScores)[3] <- "polynomial-generated bounded log 2 scores"
head(orderedLogScores)
write_csv(orderedLogScores, path = "~/Documents/Project/Programming/Transcription Factor Binding/HumanWG_motifs/Data/singleLineFinalDf and Unified table/FinalLog2Scores/polynomialLog2Scores_10_11_2017.csv")
