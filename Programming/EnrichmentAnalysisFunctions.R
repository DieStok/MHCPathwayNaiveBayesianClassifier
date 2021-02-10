##Script for the analysis of enrichments of immune genes in a certain group##
##make 2 functions nrich_nalyse_discrete and nrich_nalyse_continuous, whereby the continuous one makes
##bins itself and also allows for user set bins. Discrete treats every unique case of the values for that 
##data type as a bin (i.e. yes/no, but also some, few, many).

##Should:
#1. take a data frame of immune genes with other data on those immune genes
#2. take a data frame of PPI or PS or domains and gene IDs (with column names for the Values and IDs column supplied by the user)
#3. take a list with all the Ensembl IDs considered in this study (i.e. the output of BioMart when asked for all Ensembl Ids)
#Get which amount of set 1 is in set 2
#also which amount of 3 is in 2 (how well does this dataset cover all the genes?)
#Do the enrichment
#do Benjamini-Hochberg-correction
#allow the user to see which positive set genes and which genes from the whole genome data had and did not have the feature.



#operators
`%ni%` <- Negate(`%in%`)

#packages
library("pacman")
pacman::p_load(ggplot2, dplyr, stringr, purrr, readr, LaF, ff, data.table, R.utils, base, utils, devtools, microbenchmark)


#Testing
allEnsemblGenes <- fread("~/Documents/Project/Programming/Positive selection/Data/AllHumanGenesPositiveSelection_Version78Ensembl.txt", sep = ",", header = T)
head(allEnsemblGenes)
aegIDCol <- "Ensembl Gene ID"

positiveSet <- fread("~/Documents/Project/Programming/Transcription Factor Binding/PositiveList_2.csv", sep = ",", header = T)
head(positiveSet)
psIDCol <- "Ensembl_accession"

dataSet <- fread("~/Documents/Project/Programming/Positive selection/Data/TableS4__331_PSG__info_statistics.txt", sep = "\t", header = T)
head(dataSet)
dsIDCol <- "Ensembl.Gene.ID"
dsValueCol <- "PSR.total"




#nrich_nalyse_yesno <- function(positiveSet, psIDCol = "Ensembl_ID", dataSet, dsValueCol, dsIDCol, allEnsemblGenes, aegIDCol) {
  
  
  pctCoverageDataSet          <- sum(allEnsemblGenes[, get(aegIDCol)] %in% dataSet[, get(dsIDCol)])/length(allEnsemblGenes[, get(aegIDCol)])*100
  print(paste0("In total, ", pctCoverageDataSet, "% of the genes are covered in this dataset."))
  
  
  presentPSGenes              <- positiveSet[, get(psIDCol)][positiveSet[, get(psIDCol)] %in% dataSet[, get(dsIDCol)]]
  presentDfPSGenes            <- positiveSet[positiveSet[, get(psIDCol)] %in% dataSet[, get(dsIDCol)]]
  
  absentPSGenes               <- positiveSet[, get(psIDCol)][positiveSet[, get(psIDCol)] %ni% dataSet[, get(dsIDCol)]]
  absentDfPSGenes             <- positiveSet[positiveSet[, get(psIDCol)] %ni% dataSet[, get(dsIDCol)]]
  amountPresentPSGenes        <- nrow(presentDfPSGenes)
  totalPSGenes                <- length(positiveSet[, get(psIDCol)])
  pctPresentPSGenes           <- amountPresentPSGenes/totalPSGenes*100
  cat(paste0(pctPresentPSGenes, "% of the positive set genes are present in the dataset, i.e. ",
               amountPresentPSGenes, " out of ", totalPSGenes, " positive set genes."))
  
  
  presentDataSet             <- allEnsemblGenes[, get(aegIDCol)][allEnsemblGenes[, get(aegIDCol)] %in% dataSet[, get(dsIDCol)]]
  presentDfDataSet           <- allEnsemblGenes[allEnsemblGenes[,get(aegIDCol)] %in% presentDataSet, ]
  absentDataSet              <- allEnsemblGenes[, get(aegIDCol)][allEnsemblGenes[, get(aegIDCol)] %ni% dataSet[, get(dsIDCol)]]
  absentDfDataSet            <- allEnsemblGenes[allEnsemblGenes[,get(aegIDCol)] %ni% presentDataSet, ]
  
  
  allEnsemblGenesDf <- data.table(ENSEMBL_ID = allEnsemblGenes$aegIDCol,
                                  Value = "")
  allEnsemblGenesDf[ENSEMBL_ID == presentPSGenes, ]["Value"] = "yes"
  
  #analysis
  
  bins <- unique(dataSet[, get(dsValueCol)])
  
  
  
  fisherDataFrame <- data.frame(bins = character(),
                                pValue=numeric(),
                                sampleEstimate=numeric(),
                                lowerConfidenceInterval=numeric(),
                                upperConfidenceInterval=numeric())
  
  
  for(i in 1:length(bins)) {
    
    PSGenes       <- dataSet[dataSet[, get(dsIDCol)] %in% presentPSGenes, ]
    WGGenes       <- dataSet[dataSet[, get(dsIDCol)] %ni% presentPSGenes, ]
    
    binPSGenes    <- subset(dataSet[dataSet[, get(dsIDCol)] %in% presentPSGenes, ], dsValueCol == bin[i])
    binWGGenes    <- subset(dataSet[dataSet[, get(dsIDCol)] %ni% presentPSGenes, ], dsValueCol == bin[i])
    
    
    
     
    listMatrices <- list()
    fisherMatrix <- matrix(c(nrow(binPSGenes),
                             PSGenes - nrow(binPSGenes),
                             binWGGenes,
                             nrow(binWGGenes) - binWGGenes),
                           nrow = 2,
                           dimnames = list(Bin = c(bin[i], "Other bins"),
                                           Set = c("Positive set", "Whole genome")))
    fisherTestResult <- fisher.test(fisherMatrix)
    pValue = fisherTestResult["p.value"]
    lowerConfInt <- unlist(fisherTestResult["conf.int"])[1]
    upperConfInt <- unlist(fisherTestResult["conf.int"])[2]
    sampleEst <- fisherTestResult["estimate"]
    
    listMatrices <- append(listMatrices,fisherMatrix)
    addVector <- c(as.character(bin[i]),
                   pValue,
                   sampleEst,
                   lowerConfInt,
                   upperConfInt)
    names(addVector) <- c("bins", "pValue", "sampleEstimate",
                          "lowerConfidenceInterval", "upperConfidenceInterval")
    #add to the Df the results of the fisherTest
    fisherDataFrame <- rbind(fisherDataFrame, addVector, make.row.names = FALSE)
    
  }
  #adjust for multiple testing with Benjamini-Hochberg methodology (utilises FDR)
  fisherDataFrame$qValue <- p.adjust(fisherDataFrame$pValue, method="BH")
  return(fisherDataFrame)
#}




  
  
  
  








computeFishersTest = function(wholeGenomeTFBSCountData, positiveSetTFBSCountData, run = FALSE) {
  
  #don't compute again unless explicitly told to. Otherwise, just open the .rds
  
    
    #Get only those genes from the whole genome that have TFBS that are also found in the positive MHC set
    relevantWholeGenome <- wholeGenomeTFBSCountData[wholeGenomeTFBSCountData$TFBS %in% positiveSetTFBSCountData$TFBS,] 
    totalCountsWholeGenome <- sum(wholeGenomeTFBSCountData$Counts)
    totalCountsMHCGenes <- sum(positiveSetTFBSCountData$Counts)
    
    #listMatrices is for debug purposes, to check whether correct contingency tables are created
    listMatrices <- list()
    fisherDataFrame <- data.frame(TFBS = character(),
                                  pValue=numeric(),
                                  sampleEstimate=numeric(),
                                  lowerConfidenceInterval=numeric(),
                                  upperConfidenceInterval=numeric())
    
    #cycle over all the rows of 
    for(bindingSites in 1:nrow(positiveSetTFBSCountData)) {
      
      fisherMatrix <- matrix(c(positiveSetTFBSCountData$Counts[bindingSites],
                               totalCountsMHCGenes - positiveSetTFBSCountData$Counts[bindingSites],
                               relevantWholeGenome$Counts[bindingSites],
                               totalCountsWholeGenome - relevantWholeGenome$Counts[bindingSites]),
                             nrow = 2,
                             dimnames = list(TFBS = c(relevantWholeGenome$TFBS[bindingSites], "TotalTFBS"),
                                             Set = c("Positive set", "Whole genome")))
      fisherTestResult <- fisher.test(fisherMatrix)
      pValue = fisherTestResult["p.value"]
      lowerConfInt <- unlist(fisherTestResult["conf.int"])[1]
      upperConfInt <- unlist(fisherTestResult["conf.int"])[2]
      sampleEst <- fisherTestResult["estimate"]
      
      listMatrices <- append(listMatrices,fisherMatrix)
      addVector <- c(as.character(relevantWholeGenome$TFBS[bindingSites]),
                     pValue,
                     sampleEst,
                     lowerConfInt,
                     upperConfInt)
      names(addVector) <- c("TFBS", "pValue", "sampleEstimate",
                            "lowerConfidenceInterval", "upperConfidenceInterval")
      #add to the Df the results of the fisherTest
      fisherDataFrame <- rbind(fisherDataFrame, addVector, make.row.names = FALSE)
      
    }
    #adjust for multiple testing with Benjamini-Hochberg methodology (utilises FDR)
    fisherDataFrame$qValue <- p.adjust(fisherDataFrame$pValue, method="BH")
    saveRDS(fisherDataFrame, file = paste(getwd(), "/fisherDataFrame.rds", sep = ""))
    return(fisherDataFrame)
  }
  
  
}