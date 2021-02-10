##
##
## Started on 13th of February. 
##
##



library("pacman")
pacman::p_load(ggplot2, dplyr, stringr, purrr, readr, LaF, ff, data.table, R.utils, base, utils, devtools, microbenchmark, reshape2, magrittr, tidyr)
options(stringsAsFactors = FALSE)

celLineData <- fread("/home/dieter/Documents/Project/Programming/Human tissue enrichment/Data/rna_celline.tsv")
head(celLineData)
wideCelLineData <- celLineData %>% spread(Sample, Value)
head(wideCelLineData)

#encoded in the .ods in the data directory
ImmuneLines <- c("AN3-CA", "Daudi", "HAP1", "HDLM-2", "HL-60", "HMC-1", "K-562",
                    "Karpas-707", "MOLT-4", "NB-4", "REH", "RPMI-8226", "THP-1",
                    "U-266/70", "U-266/84", "U-698", "U-937")
ImmuneLinesData <- wideCelLineData %>% select(ImmuneLines)
OtherLinesData  <- wideCelLineData %>% select(-one_of(ImmuneLines)) %>% select(-one_of(c("Gene", "Gene name", "Unit")))
head(ImmuneLinesData); head(OtherLinesData)


tTestDataFrame <- as_tibble(data.frame("Ensembl-ID" = wideCelLineData$Gene, "t-statistic" = 0, "df" = 0, "p-value" = 0,
                                       "lower95%CI" = 0, "upper95%CI" = 0, "sampleEstDifference" = rep(0, nrow(wideCelLineData))))
head(tTestDataFrame); names(tTestDataFrame)

#ugly for-loop to do the t-tests
##still need to check for equality of variances and stuff
for(row in seq_len(nrow(wideCelLineData))) {
  
  tTest                 <- t.test(ImmuneLinesData[row,], OtherLinesData[row,])
  confIntLow            <- as.numeric(tTest$conf.int[[1]])
  confIntHigh           <- as.numeric(tTest$conf.int[[2]])
  tTestDF               <- as.numeric(tTest$parameter[[1]])
  tTestStat             <- as.numeric(tTest$statistic[[1]])
  tTsig                 <- as.numeric(tTest$p.value)
  sampleEstDifference   <- as.numeric(abs(tTest$estimate[[1]]-tTest$estimate[[2]]))
  ensemblID             <- wideCelLineData[row,]$Gene
  tTestDataFrame[row,]  <- c(ensemblID, tTestStat, tTestDF, tTsig, confIntLow, confIntHigh, sampleEstDifference) 
  
}
head(tTestDataFrame)
tTestDataFrame$t.statistic <- as.numeric(tTestDataFrame$t.statistic)
tTestDataFrame$df <- as.numeric(tTestDataFrame$df); tTestDataFrame$p.value <- as.numeric(tTestDataFrame$p.value)
tTestDataFrame$lower95.CI <- as.numeric(tTestDataFrame$lower95.CI); tTestDataFrame$upper95.CI <- as.numeric(tTestDataFrame$upper95.CI)
tTestDataFrame$sampleEstDifference <- as.numeric(tTestDataFrame$sampleEstDifference)
head(tTestDataFrame)
tTestDataFrame$p.value.BH = p.adjust(tTestDataFrame$p.value, method = "BH")
head(tTestDataFrame)

#zijn genen die juist lager tot expressie komen in immuuncellen per definitie oninteressant? niet-tot expressie komen kan ook belangrijk
#zijn voor goede functionliteit, non?

higherExpressedImmuneTissues <- tTestDataFrame %>% filter(p.value.BH <= 0.05, lower95.CI > 0)
nrow(higherExpressedImmuneTissues)
head(higherExpressedImmuneTissues, 30)
tail(higherExpressedImmuneTissues, 30)

####
####    Enrichment Stuff
####

positiveSet <- fread("~/Documents/Project/Programming/DataForAll/PositiveList_2.csv")
positiveSetIdCol <- "Ensembl_accession"
allEnsemblGenes <- fread("~/Documents/Project/Programming/DataForAll/EnsemblBioMartIDs_GRCh38.p10.txt", sep = ",", header = T)
allEnsemblGenesIDCol <- "Gene stable ID"
allEnsemblGenesNonPositiveTissueData <- allEnsemblGenes %>% filter(`Gene stable ID` %ni% positiveSet$Ensembl_accession) 
tissueIDCol <- "Ensembl.ID"
tissueValueCol <- "sampleEstDifference"
higherExpressedImmuneTissues <- as.data.table(higherExpressedImmuneTissues)

nrich_nalyse_yesno <- function(positiveSet, psIDCol = "Ensembl_accession", dataSet, dsValueCol, dsIDCol, allEnsemblGenes, aegIDCol) {
  
  
  pctCoverageDataSet          <- sum(allEnsemblGenes[, get(aegIDCol)] %in% dataSet[, get(dsIDCol)])/length(allEnsemblGenes[, get(aegIDCol)])*100
  print(paste0("In total, ", pctCoverageDataSet, "% of the genes are covered in this dataset."))
  
  
  presentPSGenes              <- positiveSet[, get(psIDCol)][positiveSet[, get(psIDCol)] %in% dataSet[, get(dsIDCol)]]
  presentDfPSGenes            <- positiveSet[positiveSet[, get(psIDCol)] %in% dataSet[, get(dsIDCol)]]
  
  absentPSGenes               <- positiveSet[, get(psIDCol)][positiveSet[, get(psIDCol)] %ni% dataSet[, get(dsIDCol)]]
  absentDfPSGenes             <- positiveSet[positiveSet[, get(psIDCol)] %ni% dataSet[, get(dsIDCol)]]
  amountPresentPSGenes        <- nrow(presentDfPSGenes)
  totalPSGenes                <- length(positiveSet[, get(psIDCol)])
  totalWGGenes                <- length(allEnsemblGenes[, get(aegIDCol)])
  pctPresentPSGenes           <- amountPresentPSGenes/totalPSGenes*100
  cat(paste0(pctPresentPSGenes, "% of the positive set genes are present in the dataset, i.e. ",
             amountPresentPSGenes, " out of ", totalPSGenes, " positive set genes."))
  
  
  presentDataSet             <- allEnsemblGenes[, get(aegIDCol)][allEnsemblGenes[, get(aegIDCol)] %in% dataSet[, get(dsIDCol)]]
  presentDfDataSet           <- allEnsemblGenes[allEnsemblGenes[,get(aegIDCol)] %in% presentDataSet, ]
  absentDataSet              <- allEnsemblGenes[, get(aegIDCol)][allEnsemblGenes[, get(aegIDCol)] %ni% dataSet[, get(dsIDCol)]]
  absentDfDataSet            <- allEnsemblGenes[allEnsemblGenes[,get(aegIDCol)] %ni% presentDataSet, ]
  
  
  allEnsemblGenesDf <- data.table(ENSEMBL_ID = allEnsemblGenes[, get(aegIDCol)],
                                  Value = "no")
  head(allEnsemblGenesDf)
  tail(allEnsemblGenesDf)
  allEnsemblGenesDf[ENSEMBL_ID %in% presentDataSet, ][, "Value"] <- "yes"
  
  #analysis
  
  bins <- c("no", "yes")
  
  
  
  fisherDataFrame <- data.frame(bins = character(),
                                pValue=numeric(),
                                sampleEstimate=numeric(),
                                lowerConfidenceInterval=numeric(),
                                upperConfidenceInterval=numeric())
  
  
  
  
  PSyes         <- filter(allEnsemblGenesDf, Value == "yes")
  PSyes         <- PSyes[ PSyes$ENSEMBL_ID %in% positiveSet[, get(psIDCol)],]
  
  PSno          <- filter(allEnsemblGenesDf, Value == "no")
  PSno          <- PSno[ PSno$ENSEMBL_ID %in% positiveSet[, get(psIDCol)],]
  
  WGyes         <- filter(allEnsemblGenesDf, Value == "yes")
  WGyes         <- WGyes[ WGyes$ENSEMBL_ID %ni% positiveSet[, get(psIDCol)],]
  
  WGno          <- filter(allEnsemblGenesDf, Value == "no")
  WGno          <- WGno[ WGno$ENSEMBL_ID %ni% positiveSet[, get(psIDCol)],]
  
  
  PSGenes       <- dataSet[dataSet[, get(dsIDCol)] %in% presentPSGenes, ]
  WGGenes       <- dataSet[dataSet[, get(dsIDCol)] %ni% presentPSGenes, ]
  
  
  listMatrices <- list()
  fisherMatrix <- matrix(c(nrow(PSyes),
                           totalPSGenes - nrow(PSyes),
                           nrow(WGyes),
                           totalWGGenes - nrow(WGyes)),
                         nrow = 2,
                         dimnames = list(Bin = c("yes", "no"),
                                         Set = c("Positive set", "Whole genome")))
  fisherTestResult <- fisher.test(fisherMatrix)
  
  pValue = fisherTestResult["p.value"]
  lowerConfInt <- unlist(fisherTestResult["conf.int"])[1]
  upperConfInt <- unlist(fisherTestResult["conf.int"])[2]
  sampleEst <- fisherTestResult["estimate"]
  
  listMatrices <- append(listMatrices,fisherMatrix)
  addVector <- c("viral PPI",
                 pValue,
                 sampleEst,
                 lowerConfInt,
                 upperConfInt)
  names(addVector) <- c("bins", "pValue", "sampleEstimate",
                        "lowerConfidenceInterval", "upperConfidenceInterval")
  #add to the Df the results of the fisherTest
  fisherDataFrame <- rbind(fisherDataFrame, addVector, make.row.names = FALSE)
  fisherDataFrame
}



tissueEnrichmentYesNo <- nrich_nalyse_yesno(positiveSet, positiveSetIdCol, higherExpressedImmuneTissues, 
                                            tissueValueCol, tissueIDCol, allEnsemblGenes, allEnsemblGenesIDCol)
tissueEnrichmentYesNo


###now for the actual bins

range(higherExpressedImmuneTissues$sampleEstDifference)
##perhaps I should select for at least 2-fold difference? we shall see

PositiveSetTissue <- higherExpressedImmuneTissues[higherExpressedImmuneTissues$Ensembl.ID %in% positiveSet$Ensembl_accession,]$sampleEstDifference
PositiveSetTissue
OtherTissue       <- higherExpressedImmuneTissues[higherExpressedImmuneTissues$Ensembl.ID %ni% positiveSet$Ensembl_accession,]$sampleEstDifference
OtherTissue
sum(OtherTissue <10)
sum(OtherTissue >= 10)
#ZScoresNotIndata(hence, 0)



dataForPlottingTissue <- data.frame(higherExpressedImmuneTissues$sampleEstDifference) %>%
  mutate(set = if_else(higherExpressedImmuneTissues$Ensembl.ID %in% positiveSet$Ensembl_accession, "PS", "NS")) %>%
  rename(sampleEstDifference = higherExpressedImmuneTissues.sampleEstDifference) %>%
  mutate(EnsemblID = higherExpressedImmuneTissues$Ensembl.ID)  %>%
  arrange(EnsemblID)
head(dataForPlottingTissue, 50); tail(dataForPlottingTissue)

range(dataForPlottingTissue$sampleEstDifference)
breaksTissue <- c(-1, 10, 1500)
breaksTissueLog10 <- c(1, 10, 10000) %>% log10()

totalPlotTissueBins <- dataForPlottingTissue %>%
  ggplot(aes(x = sampleEstDifference, fill = set, y = ..density..)) + geom_histogram(breaks = breaksTissue, colour = "black") 
totalPlotTissueBins



#get at the density data for calculating enrichments
dataTotalHistTissue <- ggplot_build(totalPlotTissueBins)
sum(dataTotalHistTissue$data[[1]]$density)

getScoresAdapted <- function(plotBuild) {
  
  logScores <- numeric()
  for(i in seq(2, nrow(plotBuild$data[[1]]), by = 2)) {
    
    logScores = c(logScores, log2((plotBuild$data[[1]][(i-1),]$count/nrow(positiveSet))/(plotBuild$data[[1]][i,]$count/nrow(allEnsemblGenesNonPositiveTissueData))))
    
  }
  newPlot <- plotBuild$plot + annotate("text", x = unique(plotBuild$data[[1]]$x), y = 0.25, label = round(logScores, digits = 3), angle = 90)
  
  list(plot = newPlot, scores = logScores)
  
}
tissuePlot <- getScoresAdapted(dataTotalHistTissue)
tissuePlot$plot 
tissuePlot$scores

kaas <- ggplot_build(tissuePlot$plot)
kaas$data[[1]]$xmin <- kaas$data[[1]]$xmin + 1
kaas$plot 

#log10
totalPlotTissueBinsLog10 <- dataForPlottingTissue %>% mutate(sampleEstDifference = log10(sampleEstDifference + 1)) %>%
  ggplot(aes(x = sampleEstDifference)) + geom_histogram(aes(fill = set, y = ..density..),  breaks = breaksTissueLog10, colour = "black") 
tissueLog10PlotContinued <- ggplot_build(totalPlotTissueBinsLog10) 
tissueLog10PlotContinued$plot + scale_x_continuous(labels = c(0, 10, 100, 1000, 10000), expand = c(0,0), name = "sample estimate of difference in TPM")  +
  scale_y_continuous(name = "density", expand = c(0,0), limits = c(0,1)) + annotate("text", x = c(0.5,2.5), y = 0.75, label = round(tissuePlot$scores, digits = 3), angle = 90) +
  theme_bw() + geom_segment(aes(x=3.153731, xend=3.153731, y=0.1, yend=0.017), 
                            arrow = arrow(length = unit(0.15, "cm"), type = "closed")) 




finalTissueLog2DataAllGenes <- dataForPlottingTissue %>% group_by(gr = cut(sampleEstDifference, breaksTissue)) %>%
  mutate(log2Scores = tissuePlot$scores[gr])
finalTissueLog2DataAllGenes
unique(finalTissueLog2DataAllGenes$log2Scores)
unique(finalTissueLog2DataAllGenes$gr)


#what positive set genes pop up?
dataForPlottingTissue %>% filter(EnsemblID %in% positiveSet$Ensembl_accession, sampleEstDifference > 0)
dataForPlottingTissue %>% filter( sampleEstDifference > 0)
dataForPlottingTissue %>% filter(EnsemblID %in% positiveSet$Ensembl_accession, sampleEstDifference >= 10)
dataForPlottingTissue %>% filter(sampleEstDifference >= 100) %>% arrange(sampleEstDifference)

##Ik weet niet hoezeer ik te spreken ben over de resultaten. Ze berusten allemaal op zeer weinig positieve set-genen.
#alles boven de 100 rust op 1 gen. Dat is eigenlijk onzin, want bij een leave-one-out analysis krijg je dus allemaal Inf scores.
#Ik denk dat ik derhalve maar een grote bin maak van alles boven de 10. Zie plot hieronder en data hierboven.
oldLog10Plot #<- tissueLog10PlotContinued




