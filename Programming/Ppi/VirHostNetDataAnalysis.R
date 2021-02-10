#AnalyseEnrichmentVirHostData
#written 19 juli 2017 DIeter Stoker

#packages
library("pacman")
pacman::p_load(ggplot2, dplyr, stringr, purrr, readr, LaF, ff, data.table, R.utils, base, utils, devtools, microbenchmark)

finalVirHostData <- fread("~/Documents/Project/Programming/Ppi/Data/virhostnetuniqueForAnalysis.csv", header = T, sep = ",")
head(finalVirHostData)
str(finalVirHostData)

#is this present in the positiveset?
positiveSet <- fread("~/Documents/Project/Programming/Transcription Factor Binding/PositiveList_2.csv", header = T, sep = ",")
head(positiveSet)
sum(positiveSet$Ensembl_accession %in% finalVirHostData$`Gene stable ID`)
length(positiveSet$Ensembl_accession)
#seems to be. 33/62 of the positive set, while only 5000 out of 30000 genes are here. Let's do this.

`%ni%` <- Negate(`%in%`)




#Testing
allEnsemblGenes <- fread("~/Documents/Project/Programming/DataForAll/EnsemblBioMartIDs_GRCh38.p10.txt", sep = ",", header = T)
head(allEnsemblGenes)
aegIDCol <- "Gene stable ID"

positiveSet <- fread("~/Documents/Project/Programming/Transcription Factor Binding/PositiveList_2.csv", sep = ",", header = T)
head(positiveSet)
psIDCol <- "Ensembl_accession"

dataSet <- finalVirHostData
head(dataSet)
dsIDCol <- "ensembl_gene_id"
dsValueCol <- "measuredPPI"




#nrich_nalyse_yesno <- function(positiveSet, psIDCol = "Ensembl_ID", dataSet, dsValueCol, dsIDCol, allEnsemblGenes, aegIDCol) {


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
#note: what is this LRG_10 stuff? perhaps I should remove that?
allEnsemblGenesDf[ENSEMBL_ID %in% presentDataSet, ][, "Value"] <- "yes"

#analysis

bins <- c("no", "yes")



fisherDataFrame <- data.frame(bins = character(),
                              pValue=numeric(),
                              sampleEstimate=numeric(),
                              lowerConfidenceInterval=numeric(),
                              upperConfidenceInterval=numeric())


#for(i in 1:length(bins)) {
  
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
  
  #binPSGenes    <- subset(dataSet[dataSet[, get(dsIDCol)] %in% presentPSGenes, ], dsValueCol == bin[i])
  #binWGGenes    <- subset(dataSet[dataSet[, get(dsIDCol)] %ni% presentPSGenes, ], dsValueCol == bin[i])
  
  
  
  
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
#adjust for multiple testing with Benjamini-Hochberg methodology (utilises FDR)
#fisherDataFrame$qValue <- p.adjust(fisherDataFrame$pValue, method="BH")
return(fisherDataFrame)
#}