######################################################################################################
# Script started 4-08-2017. Author: Dieter Stoker. Function: Enrichment in the intersects and unions
# of the viral PPi data
######################################################################################################

library("pacman")
pacman::p_load(ggplot2, dplyr, stringr, purrr, readr, LaF, ff, data.table, R.utils, base, utils, devtools, microbenchmark, reshape2, rlang)



##Change the original code into a function that allows multiple positive sets
##NOTE: this was unnecessary. The different positive sets were only used for calculations here, not in the final Bayesian classifier.


generateViralPPICSVForBayesian <- function(positiveSet = fread("~/Documents/Project/Programming/DataForAll/PositiveList_2.csv"),
                                           positiveSetIDCol = "Ensembl_accession", 
                                           PositiveSetName = "MHCIAndII",
                                           allEnsemblGenes = fread("~/Documents/Project/Programming/DataForAll/EnsemblBioMartIDs_GRCh38.p10.txt", sep = ",", header = T),
                                           allEnsemblGenesIDCol = "Gene stable ID",
                                           IntersectData = fread("~/Documents/Project/Programming/Ppi/Data/IntersectionPPiData.csv"),
                                           IntersectDataIDCol = "ensembl_gene_id",
                                           UnionData = fread("~/Documents/Project/Programming/Ppi/Data/UnionPPiData.csv"),
                                           UnionDataIDCol = "ensembl_gene_id",
                                           valueCol = "measuredPPI") {
  
  IntersectDataForBayesian <- IntersectData %>% rename(IntersectMeasuredPPI = measuredPPI ) %>% select(-uniprot_gn) %>% arrange(!! rlang::sym(IntersectDataIDCol))
  UnionDataForBayesian <- UnionData %>% rename(UnionMeasuredPPI = measuredPPI) %>% select(-uniprot_gn) %>% arrange(!! rlang::sym(UnionDataIDCol))
  head(UnionDataForBayesian)
  head(IntersectDataForBayesian)
  UnionDataForBayesian %>% filter(presentInPhisto == "yes", presentInHPIDB == "yes", presentInVirHostNet == "yes")
  UnionDataForBayesian %>% filter(UnionMeasuredPPI == "no")
  
  totalDataFrame <- data.frame(ensembl_ID = as.data.frame(allEnsemblGenes)[,allEnsemblGenesIDCol], UnionMeasuredPPI = "", IntersectMeasuredPPI = "",
                               presentInPhisto = "", presentInHPIDB = "", presentInVirHostNet = "") %>% arrange(ensembl_ID)
  totalDataFrame$IntersectMeasuredPPI <- ifelse( totalDataFrame$ensembl_ID %in% IntersectDataForBayesian[,IntersectDataIDCol], IntersectDataForBayesian$IntersectMeasuredPPI, "no")
  totalDataFrame$UnionMeasuredPPI     <- ifelse(totalDataFrame$ensembl_ID %in% UnionDataForBayesian[,UnionDataIDCol], UnionDataForBayesian$UnionMeasuredPPI, "no")
  totalDataFrame <- totalDataFrame %>% left_join(UnionDataForBayesian, by = c("ensembl_ID" ="ensembl_gene_id")) %>%
    select(-presentInPhisto.x, -presentInHPIDB.x, -presentInVirHostNet.x, -UnionMeasuredPPI.y) %>%
    rename(presentInPhisto = presentInPhisto.y, presentInVirHostNet = presentInVirHostNet.y, presentInHPIDB = presentInHPIDB.y,
           UnionMeasuredPPI = UnionMeasuredPPI.x)

  totalDataFrame <- totalDataFrame %>% as.data.frame()
  #diagnostics
  #totalDataFrame[1,]
  #class(totalDataFrame[1,4])
  #totalDataFrame[1,4]
  #nrow(totalDataFrame)
  
  for(rows in 1:nrow(totalDataFrame)) {
    
    for (cols in 1:ncol(totalDataFrame)) {
      
      if(is.na(totalDataFrame[rows,cols])) {
        totalDataFrame[rows, cols] = "no"
      }
      
      
    }
  }
  nrow(totalDataFrame)
  head(totalDataFrame)
  #remove rows that are there double
  totalDataFrame <- totalDataFrame[!duplicated(totalDataFrame$ensembl_ID),]
  #check for positive set coverage
  positiveSet$Ensembl_accession %in% totalDataFrame$ensembl_ID
  nrow(totalDataFrame)
  write_csv(totalDataFrame, path = paste0("~/Documents/Project/Programming/Ppi/Data/ViralPPIForBayesianCSV", PositiveSetName, ".csv"))
  totalDataFrame
  
  
}




# positiveSet <- fread("~/Documents/Project/Programming/DataForAll/PositiveList_2.csv")
# positiveSetIdCol <- "Ensembl_accession"
# allEnsemblGenes <- fread("~/Documents/Project/Programming/DataForAll/EnsemblBioMartIDs_GRCh38.p10.txt", sep = ",", header = T)
# allEnsemblGenesIDCol <- "Gene stable ID"
# IntersectData <- fread("~/Documents/Project/Programming/Ppi/Data/IntersectionPPiData.csv")
# IntersectIDCol <- "ensembl_gene_id"
# head(IntersectData)
# UnionData <- fread("~/Documents/Project/Programming/Ppi/Data/UnionPPiData.csv")
# UnionIDCol <- "ensembl_gene_id"
# dataValueCol <- "measuredPPi"
# head(UnionData)
# tail(UnionData)

#################################################################
##
##
##            Output file for the .csv for the Bayesian
##
##
##
#################################################################


#Note: actual enrichment is calculated in scriptForRunningTheBayesian. Thus, the positive sets do not really matter here.
#Shame I figured that out after writing code as if they do. All code for different PS is therefore commented out.
MHCIAndIIData            <- generateViralPPICSVForBayesian()
#MHCIPositiveSetLocation  <- "/home/dieter/Documents/Project/Programming/DataForAll/PositiveList_MHCIOnly.csv"
#MHCIIPositiveSetLocation <- "/home/dieter/Documents/Project/Programming/DataForAll/PositiveList_MHCIIOnly.csv"
#MHCIData                 <- generateViralPPICSVForBayesian(positiveSet = fread(MHCIPositiveSetLocation),
#                                                PositiveSetName = "MHCIOnly")
#MHCIIData                <- generateViralPPICSVForBayesian(positiveSet = fread(MHCIIPositiveSetLocation),
#                                                PositiveSetName = "MHCIIOnly")

#head(MHCIData); head(MHCIIData); head(MHCIAndIIData)
#str(MHCIData); str(MHCIIData)
#sum(MHCIData$ensembl_ID %in% MHCIIData$ensembl_ID)
###############################################
###
### Code below is not used in the final classifier.
### It was used to check for enrichment without the x-gold cross-validation used
### in the final classifier. The code is therefore deprecated and there for reference only.
###
###############################################



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

IntersectEnrichment<- nrich_nalyse_yesno(positiveSet, positiveSetIdCol, IntersectData, dataValueCol, IntersectIDCol,
                                         allEnsemblGenes, allEnsemblGenesIDCol )
IntersectEnrichment


UnionEnrichment<- nrich_nalyse_yesno(positiveSet, positiveSetIdCol, UnionData, dataValueCol, UnionIDCol,
                                         allEnsemblGenes, allEnsemblGenesIDCol )
UnionEnrichment



##Include both in final score file
FinalViralScoreDf <- allEnsemblGenes %>% mutate(Log2ViralUnionPPiScore = if_else(allEnsemblGenes$`Gene stable ID` %in% UnionData$ensembl_gene_id, log2(UnionEnrichment$sampleEstimate), 0))
FinalViralScoreDf <- FinalViralScoreDf %>% mutate(Log2ViralIntersectPPiScore = if_else(FinalViralScoreDf$`Gene stable ID` %in% IntersectData$ensembl_gene_id, log2(IntersectEnrichment$sampleEstimate), 0))
head(FinalViralScoreDf); tail(FinalViralScoreDf)
FinalViralScoreDf <- FinalViralScoreDf %>% arrange(`Gene stable ID`)
head(FinalViralScoreDf); tail(FinalViralScoreDf)
#save
write_csv(FinalViralScoreDf, path = "~/Documents/Project/Programming/Ppi/Data/ViralPpiEnrichmentScoresWholeGenome.csv")



###Code Graveyard###

# IntersectDataForBayesian <- IntersectData %>% rename(IntersectMeasuredPPI = measuredPPI ) %>% select(-uniprot_gn) %>% arrange(ensembl_gene_id)
# UnionDataForBayesian <- UnionData %>% rename(UnionMeasuredPPI = measuredPPI) %>% select(-uniprot_gn) %>% arrange(ensembl_gene_id)
# head(UnionDataForBayesian)
# head(IntersectDataForBayesian)
# UnionDataForBayesian %>% filter(presentInPhisto == "yes", presentInHPIDB == "yes", presentInVirHostNet == "yes")
# UnionDataForBayesian %>% filter(UnionMeasuredPPI == "no")
# 
# totalDataFrame <- data.frame(ensembl_ID = allEnsemblGenes$`Gene stable ID`, UnionMeasuredPPI = "", IntersectMeasuredPPI = "",
#                              presentInPhisto = "", presentInHPIDB = "", presentInVirHostNet = "") %>% arrange(ensembl_ID)
# totalDataFrame$IntersectMeasuredPPI <- ifelse( totalDataFrame$ensembl_ID %in% IntersectDataForBayesian$ensembl_gene_id , IntersectDataForBayesian$IntersectMeasuredPPI, "no")
# totalDataFrame$UnionMeasuredPPI     <- ifelse(totalDataFrame$ensembl_ID %in% UnionDataForBayesian$ensembl_gene_id, UnionDataForBayesian$UnionMeasuredPPI, "no")
# totalDataFrame <- totalDataFrame %>% left_join(UnionDataForBayesian, by = c("ensembl_ID" ="ensembl_gene_id")) %>%
#   select(-presentInPhisto.x, -presentInHPIDB.x, -presentInVirHostNet.x, -UnionMeasuredPPI.y) %>%
#   rename(presentInPhisto = presentInPhisto.y, presentInVirHostNet = presentInVirHostNet.y, presentInHPIDB = presentInHPIDB.y,
#          UnionMeasuredPPI = UnionMeasuredPPI.x)
# totalDataFrame[1,]
# totalDataFrame <- totalDataFrame %>% as.data.frame()
# totalDataFrame[1,]
# class(totalDataFrame[1,4])
# totalDataFrame[1,4]
# nrow(totalDataFrame)
# 
# for(rows in 1:nrow(totalDataFrame)) {
#   
#   for (cols in 1:ncol(totalDataFrame)) {
#     
#     if(is.na(totalDataFrame[rows,cols])) {
#       totalDataFrame[rows, cols] = "no"
#     }
#     
#     
#   }
# }
# nrow(totalDataFrame)
# head(totalDataFrame)
# #remove rows that are there double
# totalDataFrame <- totalDataFrame[!duplicated(totalDataFrame$ensembl_ID),]
# #check for positive set coverage
# positiveSet$Ensembl_accession %in% totalDataFrame$ensembl_ID
# nrow(totalDataFrame)
# write_csv(totalDataFrame, path = "~/Documents/Project/Programming/Ppi/Data/ViralPPIForBayesianCSV.csv")
# ################################


