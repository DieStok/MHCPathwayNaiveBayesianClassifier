#####8 FEBRUARY 2018
###Script to get ensemblIDs for Z scores and incorporate final scores
###Z-scores = Neefjes lab screen for surface expression perturbation of MHC II


#Purpose: Staining was performed with CLIP and L243. CLIP measured MHCII-CLIP concentrations (Signalling a problem with peptide loading)
#L243 measures MHC II-peptide at the cell surface. As a metric for maximum potential of MHC II perturbation in any form by a gene, 
#I wish to select, per gene, the highest score of the two. These Z-scores will then be subjected to enrichment analysis in the Bayesian
#(those genes with more of an effect on MHC II surface expression might, themselves, be involved in regulating it)

###########################################
###########################################
##                                       ##
##      Loading packages and data        ##
##                                       ##
###########################################
###########################################
library("pacman")
pacman::p_load(ggplot2, dplyr, stringr, purrr, readr, LaF, ff, data.table, R.utils, base, utils, devtools, microbenchmark, reshape2, magrittr, biomaRt)
options(stringsAsFactors = FALSE)

#What does DO mean?? --> answer after e-mail correspondence: HLA-DO, MHC II molecule
#siDM = positive control
#siControl = negative control (nonsense siRNA)
#wt = wild type (i.e. nothing altered, used to see what is significantly different)
#note, for some genes, triplicates only have a value for one of the three triplicates. These were included.
#However, be aware that medians for those genes are thus not informative.


##########
#L243Data#
##########
L243Data <- fread("/home/dieter/Documents/Project/Programming/MHC2loadingandSurfaceExpressionNeefjes/Data/L243_rescreen_Zscores.csv") %>% arrange(GeneID) %>%
  filter(GeneID != 'wt', GeneID != 'DO', GeneID != "siControl", Median != "#VALUE!", Significance != 0)
L243Data_ModelRef <- L243Data[startsWith(L243Data$GeneID, "XM_"), ]
L243Data_KnownRef <- L243Data[startsWith(L243Data$GeneID, "NM_"), ]


#is this data unique?
duplicated(L243Data_KnownRef$GeneID)
#no...
L243Data_KnownRef[10:25,]
#take median of 6 values for CORT
medianCort <- median(c(as.numeric(L243Data_KnownRef[19, 4:6, drop = TRUE]), as.numeric(L243Data_KnownRef[20, 4:6, drop = TRUE])))
#assign to CORT, remove other instance
L243Data_KnownRef[19,"Median"] <- medianCort
L243Data_KnownRef <- L243Data_KnownRef[-20,]
L243Data_KnownRef[110:125,]
#Remove TOLLIP with 2 NA's
L243Data_KnownRef <- L243Data_KnownRef[-118,]
duplicated(L243Data_KnownRef$GeneID)

#that is done. Do the same for CLIP

print("model:"); head(L243Data_ModelRef); tail(L243Data_ModelRef); print("known:")
head(L243Data_KnownRef); tail(L243Data_KnownRef); 
nrow(L243Data)
nrow(L243Data_KnownRef) + nrow(L243Data_ModelRef)
#inconsistency is correct,2 missing because they were duplicate

##########
#CLIPData#
##########
CLIPData <- fread("/home/dieter/Documents/Project/Programming/MHC2loadingandSurfaceExpressionNeefjes/Data/CLIP_rescreen_Zscores.csv") %>%
  arrange(GeneID) %>%
  filter(GeneID != 'wt', GeneID != 'DO', GeneID != "siControl", Median != "#VALUE!", Significance != 0)
CLIPData_ModelRef <- CLIPData[startsWith(CLIPData$GeneID, "XM_"), ]
CLIPData_KnownRef <- CLIPData[startsWith(CLIPData$GeneID, "NM_"), ]



duplicated(CLIPData_KnownRef$GeneID)
#no...
CLIPData_KnownRef[1:10,]
CLIPData_KnownRef[47:52,]
#CORT and TOLLIP, just like in L243. Take Median of 6 values for CORT, remove one for TOLLIP
medianCortCLIP <- median(c(as.numeric(CLIPData_KnownRef[7, 4:6, drop = TRUE]), as.numeric(CLIPData_KnownRef[8, 4:6, drop = TRUE])))
CLIPData_KnownRef[7,"Median"] <- medianCortCLIP

#remove unneeded rows
CLIPData_KnownRef <- CLIPData_KnownRef[-c(8, 48),]
duplicated(CLIPData_KnownRef$GeneID)
#done
nrow(CLIPData)
nrow(CLIPData_KnownRef) + nrow(CLIPData_ModelRef)

nrow(CLIPData_ModelRef) + nrow(L243Data_ModelRef)
nrow(CLIPData_KnownRef) + nrow(L243Data_KnownRef)






#vignette("biomaRt")
listMarts()
ensembl = useMart("ENSEMBL_MART_ENSEMBL", "hsapiens_gene_ensembl", host = "apr2018.archive.ensembl.org")
#Check welke datasets er zijn
listDatasets(ensembl)
#check de filters die je kan toepassen
listFilters(ensembl)

modelAttrQuery   <- c("ensembl_gene_id", "refseq_mrna_predicted", "description")
knownAttrQuery   <- c("ensembl_gene_id", "refseq_mrna", "description")
modelFilterQuery <- c("refseq_mrna_predicted", "biotype")
knownFilterQuery <- c("refseq_mrna", "biotype")
valueL243Model <- list(L243Data_ModelRef$GeneID, "protein_coding")
valueL243Known <- list(L243Data_KnownRef$GeneID, "protein_coding")
valueCLIPModel <- list(CLIPData_ModelRef$GeneID, "protein_coding")
valueCLIPKnown <- list(CLIPData_KnownRef$GeneID, "protein_coding")

dataZscoreL243Model = getBM(modelAttrQuery,modelFilterQuery,valueL243Model,ensembl)
dataZscoreL243Known = getBM(knownAttrQuery,knownFilterQuery,valueL243Known,ensembl)
dataZscoreCLIPModel = getBM(modelAttrQuery,modelFilterQuery,valueCLIPModel,ensembl)
dataZscoreCLIPKnown = getBM(knownAttrQuery,knownFilterQuery,valueCLIPKnown,ensembl)
#check results
dataZscoreL243Model; dataZscoreL243Known ; dataZscoreCLIPModel; dataZscoreCLIPKnown
# --> Modelled are all not protein_coding. Will disregard. I am interested in protein-coding genes only in this analysis

#is this a unique mapping for CLIP?
length(unique(dataZscoreCLIPKnown$ensembl_gene_id)); length(unique(dataZscoreCLIPKnown$refseq_mrna))
#nope. So, get the relevant data.
duplicated(dataZscoreCLIPKnown$ensembl_gene_id)
duplicated(dataZscoreCLIPKnown$refseq_mrna)
dataZscoreCLIPKnown[49:52,]

#CLIP without HLA-DMA score (siDM) added in
dataZscoreCLIPKnown <- dataZscoreCLIPKnown %>% as_tibble() %>% inner_join(CLIPData_KnownRef, by = c("refseq_mrna" = "GeneID"))
length(unique(dataZscoreCLIPKnown$ensembl_gene_id))
#add in siDM
siDMStuffCLIP <- CLIPData %>% filter(GeneID == "siDM") %>% group_by(GeneID) %>%
  summarize(medianOfAllMedians = median(as.numeric(Median), na.rm = TRUE))
siDMStuffCLIP$GeneID = "ENSG00000204257"
siDMStuffCLIP$refseq_mrna = "NM_006120"
colnames(siDMStuffCLIP) <- c("ensembl_gene_id", "MedianCLIP", "refseq_mrna")
siDMStuffCLIP
dataZscoreCLIPKnown[nrow(dataZscoreCLIPKnown) +1,] <- c(siDMStuffCLIP$ensembl_gene_id, siDMStuffCLIP$refseq_mrna, "HLA-DMA",
                                                        NA, NA, NA, NA, NA, NA, "HLA-DMA", siDMStuffCLIP$MedianCLIP, 1, "CLIP")
tail(dataZscoreCLIPKnown) %>% as.data.frame()
#siDM EnsemblID = ENSG00000204257, refseq = NM_006120


#is this a unique mapping for L243?
duplicated(dataZscoreL243Known$ensembl_gene_id)
duplicated(dataZscoreL243Known$refseq_mrna)
#Nope, but it is just many transcripts mapping to different genes. That is fine
dataZscoreL243Known <- dataZscoreL243Known %>% as_tibble() %>% inner_join(L243Data_KnownRef, by = c("refseq_mrna" = "GeneID"))
length(unique(dataZscoreL243Known$ensembl_gene_id))

#add in HLA-DMA score here as well
siDMStuffL243 <- L243Data %>% filter(GeneID == "siDM") %>% group_by(GeneID) %>%
  summarize(medianOfAllMedians = median(as.numeric(Median), na.rm = TRUE))
#Okay, so none of these are significant. For comparison, will add in a line with a median of 0.
dataZscoreL243Known[nrow(dataZscoreL243Known) +1,] <- c(siDMStuffCLIP$ensembl_gene_id, siDMStuffCLIP$refseq_mrna, "HLA-DMA",
                                                        NA, NA, NA, NA, NA, NA, "HLA-DMA", 0, 0, "L243")


##############################################################
##
##    Make a table for determining which Z score is highest
##
##############################################################

dataZscoreL243Known$Median <- as.numeric(dataZscoreL243Known$Median)
names(dataZscoreL243Known)[11] <- "MedianL243"

dataZscoreCLIPKnown$Median <- as.numeric(dataZscoreCLIPKnown$Median)
names(dataZscoreCLIPKnown)[11] <- "MedianCLIP"

combinedTable <- dataZscoreCLIPKnown %>% inner_join(dataZscoreL243Known, by = c("refseq_mrna", "ensembl_gene_id"))
combinedTable
nrow(combinedTable)




#nu wil ik het absolute maximum van de twee
#en dan de echte waarde van degene die het hoogst is
combinedTable$HighScore <- ifelse(abs(combinedTable$MedianCLIP) > abs(combinedTable$MedianL243), combinedTable$MedianCLIP, combinedTable$MedianL243)
combinedTable$HighScoreAntibody <- ifelse(abs(combinedTable$MedianCLIP) > abs(combinedTable$MedianL243), combinedTable$Antibody.x, combinedTable$Antibody.y)
combinedTable
head(combinedTable)
combinedTable$HighScore
combinedTable$HighScoreAntibody
combinedTable$MedianCLIP

combinedTableForIncorporation <- combinedTable %>% dplyr::select(ensembl_gene_id, refseq_mrna, description = description.x, HighScore, Antibody = HighScoreAntibody) 
combinedTableForIncorporation

#get the entries which only affect CLIP or only affect L243 
L243Only <- dataZscoreL243Known[dataZscoreL243Known$refseq_mrna %ni% dataZscoreCLIPKnown$refseq_mrna,]
head(L243Only)
nrow(L243Only)
CLIPOnly <- dataZscoreCLIPKnown[dataZscoreCLIPKnown$refseq_mrna %ni% dataZscoreL243Known$refseq_mrna,]
nrow(CLIPOnly)

CLIPAndL243 <- L243Only %>% full_join(CLIPOnly, by = c("ensembl_gene_id", "refseq_mrna", "description", "plate", "position", "well", "Z-score R1: normalized_r1_ch1", "Z-score R2: normalized_r2_ch1", "Z-score R3: normalized_r3_ch1", "GeneSymbol", "Significance", "Antibody", "MedianL243" = "MedianCLIP"))
names(CLIPAndL243)[11] <- "HighScore"


finalDfZscores <- CLIPAndL243 %>% dplyr::select(ensembl_gene_id, refseq_mrna, description, HighScore, Antibody) %>%
  full_join(combinedTableForIncorporation) %>% arrange(ensembl_gene_id)
finalDfZscores
head(finalDfZscores)
tail(as.data.frame(finalDfZscores))
range(finalDfZscores$HighScore)
nrow(finalDfZscores)
####
####     the below is without filtering for immune expression and deduplexing, and is also used for plotting below
####

ZscoresRescreenOnly <- finalDfZscores

#######################################
#
#     Now check whether in immune cells and true after deconvolution
#
#######################################

deconvolutionData <- fread("~/Documents/Project/Programming/MHC2loadingandSurfaceExpressionNeefjes/Data/DeconvolutionScreenResults.csv")
names(deconvolutionData)[3] <- "L243EffectDuplexAmount"; names(deconvolutionData)[5] <- "CLIPEffectDuplexAmount"
deconvolutionData <- deconvolutionData %>% as_tibble()

finalFinalDfZscores <- finalDfZscores %>% inner_join(deconvolutionData, by = c("refseq_mrna" = "Accession #")) %>% arrange(ensembl_gene_id)
head(as.data.frame(finalFinalDfZscores))
nrow(finalFinalDfZscores)
nrow(finalDfZscores)
#Okay, so now I have a table of vetted hits in protein_coding genes that are expressed in immune cells! 143 candidates with a score.
#need to put hla_DMA in again because it was not in the deconvolutionscreen
finalFinalDfZscores[nrow(finalFinalDfZscores)+1, ] <- c(siDMStuffCLIP$ensembl_gene_id, siDMStuffCLIP$refseq_mrna, "HLA-DMA",
                                                        siDMStuffCLIP$MedianCLIP, "CLIP", "-", NA, "CLIP up", NA, NA, NA)
#144 candidates!
write_csv(finalFinalDfZscores, path = "~/Documents/Project/Programming/MHC2loadingandSurfaceExpressionNeefjes/Data/MaximalZscoreData.csv")
finalFinalDfZscores

#on to the murky quagmires of overrepresentation:

#first yes/no. Are PS genes overrepresented in this data?
#function. REQUIRES A DATA.TABLE
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


#data
positiveSet <- fread("~/Documents/Project/Programming/DataForAll/PositiveList_2.csv")
positiveSetIdCol <- "Ensembl_accession"
positiveSetNonNeefjes <- positiveSet %>% filter(Neefjes_accession == "")
allEnsemblGenes <- fread("~/Documents/Project/Programming/DataForAll/EnsemblBioMartIDs_GRCh38.p10.txt", sep = ",", header = T)
allEnsemblGenesIDCol <- "Gene stable ID"
ZscoreIDCol <- "ensembl_gene_id"
ZscoreValueCol <- "HighScore"
finalFinalDfZscores <- as.data.table(finalFinalDfZscores)
names(finalFinalDfZscores)
ZscoreYesNoEnrichment <- nrich_nalyse_yesno(positiveSet = as.data.table(positiveSetNonNeefjes),
                                            psIDCol = positiveSetIdCol, dataSet = finalFinalDfZscores,
                                            dsValueCol = ZscoreValueCol, dsIDCol = ZscoreIDCol,
                                            allEnsemblGenes = allEnsemblGenes, aegIDCol = allEnsemblGenesIDCol)
ZscoreYesNoEnrichment

#for raw Zscore Data
ZscoreRescreenOnlyYesNoEnrichment <- nrich_nalyse_yesno(positiveSet = as.data.table(positiveSetNonNeefjes),
                                                        psIDCol = positiveSetIdCol, dataSet = as.data.table(ZscoresRescreenOnly),
                                                        dsValueCol = ZscoreValueCol, dsIDCol = ZscoreIDCol,
                                                        allEnsemblGenes = allEnsemblGenes, aegIDCol = allEnsemblGenesIDCol)

ZscoreRescreenOnlyYesNoEnrichment
ZscoreYesNoEnrichment
range(ZscoresRescreenOnly$HighScore); sort(ZscoresRescreenOnly$HighScore)
#so yes, they are overrepresented. Let's see how the bins look.
range(finalFinalDfZscores$HighScore)
breaksZscore <- c (-3, 0, 10)
breaksZscore2 <- c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

#need to subset positive set without those I got from the Neefjes paper in the first place!
#that would be wrong. So use that.


PositiveSetZscores <- finalFinalDfZscores[finalFinalDfZscores$ensembl_gene_id %in% positiveSetNonNeefjes$Ensembl_accession,]$HighScore
PositiveSetZscores
NegativeSetZscores <- finalFinalDfZscores[finalFinalDfZscores$ensembl_gene_id %in% allEnsemblGenesNonPositive$`Gene stable ID`,]$HighScore
NegativeSetZscores

nrow(allEnsemblGenes)
nrow(finalFinalDfZscores)
allEnsemblGenesNonPositive <- allEnsemblGenes %>% filter(`Gene stable ID` %ni% positiveSetNonNeefjes$Ensembl_accession)
nrow(allEnsemblGenesNonPositive)




dataForPlottingZscoresWithoutWGD <- data.frame(finalFinalDfZscores$HighScore) %>%
  mutate(set = if_else(finalFinalDfZscores$ensembl_gene_id %in% positiveSetNonNeefjes$Ensembl_accession, "PS", "NS")) %>%
  rename(Zscore = finalFinalDfZscores.HighScore) %>%
  mutate(EnsemblID = finalFinalDfZscores$ensembl_gene_id) %>% arrange(EnsemblID)
head(dataForPlottingZscoresWithoutWGD)
tail(dataForPlottingZscoresWithoutWGD)
unique(dataForPlottingZscoresWithoutWGD$set)
dataForPlottingZscoresWithoutWGD[dataForPlottingZscoresWithoutWGD$set == "PS",]
dataForPlottingZscoresWithoutWGD
tail(dataForPlottingZscoresWithoutWGD)
dataForPlottingZscoresWithoutWGD$Zscore <- as.numeric(dataForPlottingZscoresWithoutWGD$Zscore)

#save this data for the final Bayesian classifier table
dataForBayesianZScores <- dataForPlottingZscoresWithoutWGD %>% dplyr::select(-set) %>% rename(ZscoreHighScore = Zscore)
#get also the additional info behind this
extraInfoZscores <- finalFinalDfZscores %>% arrange(ensembl_gene_id) %>% dplyr::select(ensembl_gene_id, Antibody, `L243-Cy3`, L243EffectDuplexAmount, `CLIP-Cy5`, CLIPEffectDuplexAmount) %>%
  rename(EnsemblID = ensembl_gene_id)
extraInfoZscores$L243EffectDuplexAmount <- as.numeric(extraInfoZscores$L243EffectDuplexAmount)
class(extraInfoZscores$L243EffectDuplexAmount); class(extraInfoZscores$CLIPEffectDuplexAmount)
class(extraInfoZscores$`L243-Cy3`)
extraInfoZscores <- extraInfoZscores %>% rename(`L243.Cy3` = `L243-Cy3`, `CLIP.Cy5` = `CLIP-Cy5`)
head(extraInfoZscores, 50); tail(extraInfoZscores, 50)
dataForBayesianZScores <- dataForBayesianZScores %>% inner_join (extraInfoZscores, by = "EnsemblID")
head(dataForBayesianZScores)
noZscoresData <- data.frame(EnsemblID = allEnsemblGenes$`Gene stable ID`[allEnsemblGenes$`Gene stable ID` %ni% dataForBayesianZScores$EnsemblID],
                            ZscoreHighScore = NA, Antibody = NA, `L243-Cy3` = NA, L243EffectDuplexAmount = NA, `CLIP-Cy5` = NA, CLIPEffectDuplexAmount = NA)
noZscoresData$ZscoreHighScore <- as.numeric(noZscoresData$ZscoreHighScore)
noZscoresData$Antibody = as.character(noZscoresData$Antibody); noZscoresData$`L243.Cy3` = as.character(noZscoresData$`L243.Cy3`)
noZscoresData$L243EffectDuplexAmount = as.numeric(noZscoresData$L243EffectDuplexAmount); noZscoresData$CLIPEffectDuplexAmount = as.numeric(noZscoresData$CLIPEffectDuplexAmount)
noZscoresData$CLIP.Cy5 = as.character(noZscoresData$CLIP.Cy5)
dataForBayesianZScores <- dataForBayesianZScores %>% full_join(noZscoresData)  %>% arrange(EnsemblID)
head(dataForBayesianZScores)
tail(dataForBayesianZScores)
dataForBayesianZScores[!is.na(dataForBayesianZScores$ZscoreHighScore),]
nrow(dataForBayesianZScores[!is.na(dataForBayesianZScores$ZscoreHighScore),])
###that's correct. Save it.
dataForBayesianZScores$ZscoreHighScore <- as.numeric(dataForBayesianZScores$ZscoreHighScore)

dataForBayesianZScores %>%
  write_csv(path = "~/Documents/Project/Programming/MHC2loadingandSurfaceExpressionNeefjes/Data/HighConfidenceZscore143GenesForBayesian.csv")



#only the data itself, without adding in all the unmeasured ensemblIDs
totalPlotZscoresBins <- dataForPlottingZscoresWithoutWGD %>%
  ggplot(aes(x = Zscore, fill = set, y = ..count..)) + geom_histogram(breaks = breaksZscore2, colour = "black") 
totalPlotZscoresBins



#get at the density data for calculating enrichments
dataTotalHistZscores <- ggplot_build(totalPlotZscoresBins)
dataTotalHistZscores$data

getScoresAdapted <- function(plotBuild) {
  
  logScores <- numeric()
  for(i in seq(2, nrow(plotBuild$data[[1]]), by = 2)) {
    
    logScores = c(logScores, log2((plotBuild$data[[1]][(i-1),]$count/nrow(positiveSetNonNeefjes))/(plotBuild$data[[1]][i,]$count/nrow(allEnsemblGenesNonPositive))))
    
  }
  newPlot <- plotBuild$plot + annotate("text", x = unique(plotBuild$data[[1]]$x), y = 100, label = round(logScores, digits = 3), angle = 90)
  
  list(plot = newPlot, scores = logScores)
  
}
test <- getScoresAdapted(dataTotalHistZscores)
test$plot
test$scores

finalZscoreLog2DataAllGenes <- dataForPlottingZscoresWithoutWGD %>% group_by(gr = cut(Zscore, breaksZscore)) %>%
  mutate(log2Scores = test$scores[gr])
finalZscoreLog2DataAllGenes
unique(finalZscoreLog2DataAllGenes$log2Scores)
unique(finalZscoreLog2DataAllGenes$gr)

write_csv(finalZscoreLog2DataAllGenes, path = "~/Documents/Project/Programming/MHC2loadingandSurfaceExpressionNeefjes/Data/Log2ZscoreDataBins.csv")



#####################################################
#####################################################
###                                               ###
###     SAME AS ABOVE BUT FOR RAW RESCREEN DATA   ###
###                                               ###
#####################################################
#####################################################
#This means that it is not checked for actual expression in immune cells and off-target effects


ZscoresRescreenOnly <- as.data.table(ZscoresRescreenOnly)
ZscoreRescreenOnlyYesNoEnrichment <- nrich_nalyse_yesno(positiveSet = as.data.table(positiveSetNonNeefjes),
                                                        psIDCol = positiveSetIdCol, dataSet = ZscoresRescreenOnly,
                                                        dsValueCol = ZscoreValueCol, dsIDCol = ZscoreIDCol,
                                                        allEnsemblGenes = allEnsemblGenes, aegIDCol = allEnsemblGenesIDCol)


head(ZscoresRescreenOnly, 30)


#write a .csv for the BayesianClassifier table

BayesianClassifierZscoresRescreenOnly <- ZscoresRescreenOnly %>% dplyr::select(ensembl_gene_id, HighScore, Antibody) %>%
  rename(RescreenOnlyZscores = HighScore, EnsemblID = ensembl_gene_id)
head(BayesianClassifierZscoresRescreenOnly)
genesNotInZScoresRescreenOnly <- data.frame(EnsemblID = allEnsemblGenes$`Gene stable ID`[allEnsemblGenes$`Gene stable ID` %ni% BayesianClassifierZscoresRescreenOnly$EnsemblID],
                                            RescreenOnlyZscores = as.numeric(NA), Antibody = as.character(NA))
head(genesNotInZScoresRescreenOnly); class(genesNotInZScoresRescreenOnly$RescreenOnlyZscores)
BayesianClassifierZscoresRescreenOnly <- BayesianClassifierZscoresRescreenOnly %>% full_join(genesNotInZScoresRescreenOnly) %>% arrange(EnsemblID)
BayesianClassifierZscoresRescreenOnly[BayesianClassifierZscoresRescreenOnly$EnsemblID == "ENSG00000204257",] <- c("ENSG00000204257", siDMStuffCLIP$MedianCLIP, "CLIP")
BayesianClassifierZscoresRescreenOnly %>%
  write_csv(path = "~/Documents/Project/Programming/MHC2loadingandSurfaceExpressionNeefjes/Data/ZscoresRescreenOnly314ProteinCodingGenesForBayesian.csv")

#is HLA-DMA in there?
BayesianClassifierZscoresRescreenOnly[BayesianClassifierZscoresRescreenOnly$EnsemblID == siDMStuffCLIP$ensembl_gene_id, ]
#no. So I added it in above. 


ZscoreRescreenOnlyYesNoEnrichment
range(ZscoresRescreenOnly$HighScore); sort(ZscoresRescreenOnly$HighScore)
#so yes, they are overrepresented. Let's see how the bins look.
range(ZscoresRescreenOnly$HighScore)
breaksZscore <- c (-4, 0, 10)
breaksZscore2 <- c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

#need to subset positive set without those I got from the Neefjes paper in the first place!
#that would be wrong. So use that.


PositiveSetZscores <- ZscoresRescreenOnly[ZscoresRescreenOnly$ensembl_gene_id %in% positiveSetNonNeefjes$Ensembl_accession,]$HighScore
PositiveSetZscores
OtherZscores       <- ZscoresRescreenOnly[ZscoresRescreenOnly$ensembl_gene_id %ni% positiveSetNonNeefjes$Ensembl_accession,]$HighScore
OtherZscores

nrow(allEnsemblGenes)
nrow(ZscoresRescreenOnly) + nrow(ZscoreNotInData)
#difference of 1...why?
sum(duplicated(ZscoresRescreenOnly$ensembl_gene_id))
#not anymore --> error in S2 of Neefjes paper


dataForPlottingZscoresRescreenOnlyWithoutWGD <- data.frame(ZscoresRescreenOnly$HighScore) %>%
  mutate(set = if_else(ZscoresRescreenOnly$ensembl_gene_id %in% positiveSetNonNeefjes$Ensembl_accession, "PS", "NS")) %>%
  rename(Zscore = ZscoresRescreenOnly.HighScore) %>%
  mutate(EnsemblID = ZscoresRescreenOnly$ensembl_gene_id) %>% arrange(EnsemblID)
head(dataForPlottingZscoresRescreenOnly)
tail(dataForPlottingZscoresRescreenOnly)
unique(dataForPlottingZscoresRescreenOnly$set)
nrow(dataForPlottingZscoresRescreenOnly[dataForPlottingZscoresRescreenOnly$set == "PS",])

#only the data itself, without adding in all the unmeasured ensemblIDs
subPlotZscoresBins <- dataForPlottingZscoresRescreenOnlyWithoutWGD %>%
  ggplot(aes(x = Zscore, fill = set, y = ..count..)) + geom_histogram(breaks = breaksZscore, colour = "black") 
subPlotZscoresBins



#get at the density data for calculating enrichments
dataTotalHistZscores <- ggplot_build(subPlotZscoresBins)

##
ZscoresRescreenOnlyPlot <- getScoresAdapted(dataTotalHistZscores)
ZscoresRescreenOnlyPlot$plot

finalZscoreRescreenOnlyLog2DataAllGenes <- dataForPlottingZscoresRescreenOnly %>% group_by(gr = cut(Zscore, breaksZscore)) %>%
  mutate(log2Scores = ZscoresRescreenOnlyPlot$scores[gr])
finalZscoreRescreenOnlyLog2DataAllGenes
unique(finalZscoreRescreenOnlyLog2DataAllGenes$log2Scores)
unique(finalZscoreRescreenOnlyLog2DataAllGenes$gr)

write_csv(finalZscoreLog2DataAllGenes, path = "~/Documents/Project/Programming/MHC2loadingandSurfaceExpressionNeefjes/Data/Log2ZscoreDataBins.csv")




##########
########## Orignal Z-scores (not the rescreen, not corrected)
##########
#RETRY: read in all the Z-scores there are, and work with that data


#things to do:
#1. there are inconsistencies in how siHLA-DM, siHLA-DO, etc. are labelled. Make sure that the same symbol is in both GeneID and GeneSymbol column
#2. remove all wt (not interested in that), collapse all DO and DM to one value (these are 2 positive set genes, need their value, but just once),
#         remove siControl (nonsense siRNA added)
#3. calculate the median (na.rm = true). Perhaps just use group_by GeneSymbol (since that will automatically take care of the above as well)
#4.translate Gene names (see below) to EnsemblIDS. If I go by GeneID, I can ignore all XM (models are not protein_coding, as I found out in
#working with the rescreen data)

#things to remember:
#1. for one of the screens, there are no GeneIDs for some of the genes. WIll have to go by Gene names.

originalScreenL243 <- fread("~/Documents/Project/Programming/MHC2loadingandSurfaceExpressionNeefjes/Data/OriginalZscores/Overview_zscores_GWS_L243_addedTotalDataSheetCsv.csv")
originalScreenL243 %>% head()
#are gene symbols present for everything(i.e. no NA)?
originalScreenL243[duplicated(originalScreenL243$GeneSymbol),] %>% `$`(GeneSymbol) %>% unique()

originalScreenCLIP <- fread("~/Documents/Project/Programming/MHC2loadingandSurfaceExpressionNeefjes/Data/OriginalZscores/Overview_zscores_GWS_CLIP_AddedTotalDataSheet.csv")
head(originalScreenCLIP)

#are gene symbols present for everything(i.e. no NA)?
originalScreenCLIP[duplicated(originalScreenCLIP$GeneSymbol),] %>% `$`(GeneSymbol) %>% unique()

originalScreenCLIP %<>% mutate(Screen = "CLIP")
originalScreenL243 %<>% mutate(Screen = "L243")

originalScreenL243ToBind <- originalScreenL243 %>% dplyr::select(GeneSymbol, GeneID, Screen, normalized_r1_ch1, normalized_r2_ch1, normalized_r3_ch1)

joinedTableCLIPL243OrgScreen <- originalScreenCLIP %>% dplyr::select(GeneSymbol, GeneID, Screen, normalized_r1_ch1, normalized_r2_ch1, normalized_r3_ch1) %>%
  rbind(originalScreenL243ToBind)

#check
head(joinedTableCLIPL243OrgScreen)
tail(joinedTableCLIPL243OrgScreen)

#geneSymbols are present for everything, so that is what I will group by. First remove rows that have no data at all, i.e. where 3 columns are NA


totalSignalOriginalZScreen <-joinedTableCLIPL243OrgScreen %>%
  filter(rowSums(is.na(.[c("normalized_r1_ch1","normalized_r2_ch1", "normalized_r3_ch1")])) != 3) %>%
  group_by(Screen, GeneSymbol) %>%
  summarise(Median = median(c(normalized_r1_ch1, normalized_r2_ch1, normalized_r3_ch1), na.rm = TRUE),
            StDev = sd(c(normalized_r1_ch1, normalized_r2_ch1, normalized_r3_ch1), na.rm = TRUE) )
head(as.data.frame(totalSignalOriginalZScreen), 50)
#so note that we have here the classic error that creeps into .xls files: gene names that have been turned into dates! How lovely.
#This was already wrong in the source files, so they will be filtered out automatically in the translation to EnsemblIDs
#(there is no gene by the name of 2-Sep or 3-Apr)

#now, split up, join together by GeneSymbol, and select what is highest. In theory, because I have a lot of data, I could make
#distributions again and fit them. In practice, since I have little Time, I will just do binning on this data as well.
#Since I am interested in any effect on MHC II, I will once again select the highest value of the two screens and go with that.
#Additionally, I take the absolute. I am interested solely in size of effect, not in direction
CLIPMedians <- totalSignalOriginalZScreen %>% ungroup() %>% filter(Screen == "CLIP") %>% mutate(MedianCLIP = Median) %>% dplyr::select(GeneSymbol, MedianCLIP)
L243Medians <- totalSignalOriginalZScreen %>% ungroup() %>% filter(Screen == "L243") %>% mutate(MedianL243 = Median) %>% dplyr::select(GeneSymbol, MedianL243)
head(CLIPMedians)  
joinedMedians <- CLIPMedians %>% left_join(L243Medians, by = "GeneSymbol")
head(joinedMedians) ; tail(joinedMedians)
#are there indeed no NAs
sum(is.na(joinedMedians$MedianCLIP)); sum(is.na(joinedMedians$MedianL243))
joinedMedians %<>% rowwise() %>% mutate(selectedMedian = max(abs(MedianCLIP), abs(MedianL243)))

#filter out wt, change siDM and remove DO (since that was overexpression, which is also lowering DM apparently) --> I don't know against which of the two
#(HLA-DMA or HLA-DMB). the siRNA was done. Since all other siRNAs are against
#one gene, I will simply assign them to HLA-DMA and I will remove DO.
head(joinedMedians)
joinedMedians[joinedMedians$GeneSymbol == "siDM",]$GeneSymbol = "HLA-DMA"
joinedMedians %<>% filter(GeneSymbol != "DO", GeneSymbol != "wt")
nrow(joinedMedians)

#well, filter out things that are either dates or just numbers, biomart doesn't accept those


#only numbers = str_match(joinedMedians$GeneSymbol, "^[^A-Za-z- .&!@#$%^*()\\s,]+$")
#dates        = str_match(joinedMedians$GeneSymbol, "^[0-9]{1,2}-\\w{3}$")
joinedMedians <- joinedMedians[-grep("^[0-9]{1,2}-\\w{3}$",joinedMedians$GeneSymbol ),]
joinedMedians <- joinedMedians[-grep("^\\d+$",joinedMedians$GeneSymbol ),]
nrow(joinedMedians)
#Done. Now translate to Ensembl Ids


originalDataAttrQuery   <- c("ensembl_gene_id", "hgnc_symbol", "description")

originalDataFilterQuery <- c("hgnc_symbol", "biotype")

originalDatavalues      <- list(joinedMedians$GeneSymbol[1:20], "protein_coding")

originalDataZscoresEnsemblGenes <- getBM(originalDataAttrQuery,originalDataFilterQuery,
                                         originalDatavalues, ensembl)

tools::showNonASCII(joinedMedians$GeneSymbol)
max(nchar(joinedMedians$GeneSymbol))
which(nchar(joinedMedians$GeneSymbol)==0)

#doesn't work because there are symbols there which arent recognised. It does not exclude those, apparently. I use the Biomart web interface.
#write file to data and see what happens then
joinedMedians %>% write_csv("~/Documents/Project/Programming/MHC2loadingandSurfaceExpressionNeefjes/Data/OriginalZscores/geneSymbolListFortesting.csv")

#Worked in ENsembl. Lost about half of the genes. That is consistent with the fact that many 'genes' were just gene models.

ensemblAvailableGenesOriginalZscreen <- fread("~/Documents/Project/Programming/MHC2loadingandSurfaceExpressionNeefjes/Data/OriginalZscores/BioMartLookupOfGeneSymbolsOriginalZscores.txt")
head(ensemblAvailableGenesOriginalZscreen)
nrow(ensemblAvailableGenesOriginalZscreen)
length(unique(ensemblAvailableGenesOriginalZscreen$`Gene stable ID`))
#they are already unique. Just match and we are done.

joinedMedians %<>% left_join(ensemblAvailableGenesOriginalZscreen, by = c("GeneSymbol" = "HGNC symbol"))
head(joinedMedians)


ensemblGenesOriginalZscores <- joinedMedians %>% filter(!is.na(`Gene stable ID`)) %>% rename("absoluteOfSelectedMedian" = "selectedMedian")
ensemblGenesOriginalZscores %<>% arrange(desc(absoluteOfSelectedMedian))
ensemblGenesOriginalZscores

#for some reason, HLA-DMA was mapped to all manner of EnsemblIDs. I will restrict that to the EnsemblID that I use: ENSG00000204257

ensemblGenesOriginalZscores <- ensemblGenesOriginalZscores[-c(1:4,6:8),] 
ensemblGenesOriginalZscores

#annotate by sets and do some binning

positiveSet <- fread("~/Documents/Project/Programming/DataForAll/PositiveList_2.csv")
positiveSetNonNeefjes <- positiveSet[positiveSet$Neefjes_accession == "",]

ensemblGenesOriginalZscores$Set <- ifelse(ensemblGenesOriginalZscores$`Gene stable ID` %in% positiveSetNonNeefjes$Ensembl_accession, "PS", "NS")
head(ensemblGenesOriginalZscores)
ensemblGenesOriginalZscores %>% filter(Set == "PS")


#add in NS genes and PS genes that aren't there at the moment
#ensemblGenesOriginalZscores %>% 
  
  
allEnsemblIds <- fread(allEnsemblIdsFile, sep = ",", header = T)


totalTableZscoresAllEnsemblIDs <- allEnsemblIds %>% left_join(ensemblGenesOriginalZscores, by = "Gene stable ID") %>%
  mutate(Set = ifelse(totalTableZscoresAllEnsemblIDs$`Gene stable ID` %in% positiveSetNonNeefjes$Ensembl_accession, "PS", "NS"))

countsperBinEnsemblGenesOriginalZScores <- totalTableZscoresAllEnsemblIDs %>% mutate(gr = cut(absoluteOfSelectedMedian, breaks = c(0,2,4,18))) %>%
  group_by(gr, Set) %>% summarize(Count = n())

#function below is a shameless grab from the ScriptForRunningTheBayesian.R
positiveSetGenesZScore <- totalTableZscoresAllEnsemblIDs %>% filter(Set == "PS") %>% pull(`Gene stable ID`)
negativeSetGenesZScore <- totalTableZscoresAllEnsemblIDs %>% filter(Set == "NS") %>% pull(`Gene stable ID`)
getLog2Scores <- function(summarizedData, positiveSetGenes, negativeSetGenes) {
  
  pSSize <- length(positiveSetGenes)
  nSSize <- length(negativeSetGenes)
  groupsToCycleThrough <- length(unique(summarizedData$gr))
  #print(groupsToCycleThrough)
  
  scoreVector <- numeric(length = groupsToCycleThrough)
  for (i in seq_len(groupsToCycleThrough)) {
    
    if(!is.na(summarizedData$gr[i]))relevantSubset <- summarizedData %>% filter(gr == unique(summarizedData$gr)[i]) 
    if(i == groupsToCycleThrough) relevantSubset <- summarizedData %>% filter(gr %ni% unique(summarizedData$gr)[seq_len(i-1)])
    #print(head(relevantSubset))
    fractionPositive <- relevantSubset %>% filter(Set == "PS") %>% `$`(Count)/pSSize
    print(fractionPositive)
    fractionNegative <- relevantSubset %>% filter(Set == "NS") %>% `$`(Count)/nSSize
    score <- log2(fractionPositive / fractionNegative)
    #print(score)
    scoreVector[i] <- score
  }
  
  scoreVector
  
}

log2ScoresAllEnsemblGenesZscore <- getLog2Scores(countsperBinEnsemblGenesOriginalZScores, positiveSetGenesZScore, negativeSetGenesZScore)
countsperBinEnsemblGenesOriginalZScores$log2Scores <- rep(log2ScoresAllEnsemblGenesZscore, each = 2)

countsperBinEnsemblGenesOriginalZScores

#check what the ROC does
roc(totalTableZscoresAllEnsemblIDs[,"Set"], totalTableZscoresAllEnsemblIDs$absoluteOfSelectedMedian,
    ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,
    plot = TRUE, grid=TRUE,
    print.auc=TRUE, show.thres=TRUE, main = "ZscoreTotalSet")


#shameless grab from the ScriptForRunningTheBayesian again, function below generates scores for bins
generateBinningPlots <- function(data, title = "Title", xlabel = "Bins", ylabel = "Density of scores in bin per set", log2ScoresColumn = 4, groupColumn = 1) {
  
  listReturn <- vector("list", length = 2)
  
  dataWithDensity <- data %>% group_by(Set) %>% mutate(densityPerSet = Count/sum(Count))
  
  
  #print(dataWithDensity[,groupColumn])
  #print(dataWithDensity[,log2ScoresColumn])
  
  plot <- ggplot(dataWithDensity, aes(x = as.factor(gr), y = densityPerSet, fill = Set)) +
    geom_bar(stat = "identity", colour = "black", width = 1, position = "dodge") +
    theme_bw() + theme(panel.grid.major.x = element_blank()) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_discrete(expand   = c(0,0),   name = xlabel) +
    scale_y_continuous(expand = c(0,0),   name = ylabel, limits = c(0,1), breaks = seq(0,1,0.1)) +
    annotate("text", x = as.character(dataWithDensity[,groupColumn][[1]][seq(1,nrow(dataWithDensity), by = 2)]), y = 0.8,
             label = round(dataWithDensity[,log2ScoresColumn][[1]][seq(1,nrow(dataWithDensity), by = 2)], digits = 2), angle = 90, size = 7) +
    ggtitle(label = title)
  
  listReturn[[1]] = dataWithDensity
  listReturn[[2]] = plot
  listReturn
  
}



countsperBinEnsemblGenesOriginalZScores %>% generateBinningPlots(log2ScoresColumn = 4, title = "Total Z Score screen")

ggplot(data = ensemblGenesOriginalZscores, aes(x= cut(absoluteOfSelectedMedian, breaks = c(0,2,4,18)), fill = Set)) + geom_bar(stat = "count", position = "dodge")


#write the new Zscores to a file for incorporation in the Bayesian Table
totalTableZscoresAllEnsemblIDs %>% write_csv(path = "~/Documents/Project/Programming/MHC2loadingandSurfaceExpressionNeefjes/Data/OriginalZscores/originalZscoresForInclusionBayesian.csv")


