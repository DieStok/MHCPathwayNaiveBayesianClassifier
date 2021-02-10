###create a .csv for the Bayesian classifier

######################################################################################################
######################################################################################################
######                                                                                          ######
######                                Loading Packages and Options                              ######
######                                                                                          ######
######################################################################################################
######################################################################################################

install.packages("pacman")
library("pacman")
pacman::p_load(ggplot2, dplyr, stringr, purrr, readr, LaF, ff, data.table, R.utils, base, utils, devtools, microbenchmark, reshape2, magrittr, tidyr)
options(stringsAsFactors = FALSE)





#positiveSet and all ENsemblIDs for checking
positiveSetFile       <- "~/Documents/Project/Programming/DataForAll/PositiveList_2.csv"
allEnsemblIdsFile     <- "~/Documents/Project/Programming/DataForAll/EnsemblBioMartIDs_GRCh38.p10.txt"
positiveSet <- fread(positiveSetFile)
allEnsemblIds <- fread("~/Documents/Project/Programming/DataForAll/EnsemblBioMartIDs_GRCh38.p10.txt", sep = ",", header = T)


#Changes for allowing different positive sets (19-9-2018): 
#make read-in of TFBSScore a function
#also make the final generation of the table a function that takes TFBSScores as an argument

#TFBSScores
TFBSScoreFileLocationMHCIAndII <- "/home/dieter/Documents/Project/Programming/Transcription Factor Binding/HumanWG_motifs/Data/singleLineFinalDf and Unified table/DataForBayesianCSVMHCIandIIPositiveSet2018-09-19_15:52:36.csv"
TFBSScoreFileLocationMHCI      <- "/home/dieter/Documents/Project/Programming/Transcription Factor Binding/HumanWG_motifs/Data/singleLineFinalDf and Unified table/DataForBayesianCSVMHCIPositiveSet2018-09-19_15:54:17.csv"
TFBSScoreFileLocationMHCII     <- "/home/dieter/Documents/Project/Programming/Transcription Factor Binding/HumanWG_motifs/Data/singleLineFinalDf and Unified table/DataForBayesianCSVMHCIIPositiveSet2018-09-19_15:56:19.csv"

readTFBSScores <- function(TFBSScoreFileLocation) {
  
  TFBSScores <- fread(TFBSScoreFileLocation) %>%
    arrange(EnsemblId) %>% rename(EnsemblID = EnsemblId, totalTFBSScore = totalScore)
  colnames(TFBSScores)[3:ncol(TFBSScores)] <- paste0("ScoreForTFBS_", colnames(TFBSScores)[3:ncol(TFBSScores)])
  print(head(TFBSScores,12))
  TFBSScores
  
}

MHCIAndIIPositiveSetTFBSScores <- readTFBSScores(TFBSScoreFileLocationMHCIAndII)
MHCIPositiveSetTFBSScores <- readTFBSScores(TFBSScoreFileLocationMHCI)
MHCIIPositiveSetTFBSScores <- readTFBSScores(TFBSScoreFileLocationMHCII)

sum(TFBSScores$EnsemblID %in% positiveSet$Ensembl_accession)
nrow(TFBSScores)

#Viral ppi
#will add both union and intersect columns
viralPPIData <- fread("~/Documents/Project/Programming/Ppi/Data/ViralPPIForBayesianCSV.csv") %>% rename(EnsemblID = ensembl_ID)
head(viralPPIData)
nrow(viralPPIData)
sum(viralPPIData$EnsemblID %in% positiveSet$Ensembl_accession)
sum(positiveSet$Ensembl_accession %in% viralPPIData$EnsemblID)


#Zscoresdata for disrupting peptide loading or correct surface expression of MHC II
ZscoreData143HighConfidence <- fread("~/Documents/Project/Programming/MHC2loadingandSurfaceExpressionNeefjes/Data/HighConfidenceZscore143GenesForBayesian.csv") %>%
  select(EnsemblID, everything())
head(ZscoreData143HighConfidence)
nrow(ZscoreData143HighConfidence)
sum(ZscoreData143HighConfidence$EnsemblID %in% positiveSet$Ensembl_accession)
#discrepancy with paper values is because they also have non-protein coding genes
ZscoreData314RescreenOnly   <- fread("~/Documents/Project/Programming/MHC2loadingandSurfaceExpressionNeefjes/Data/ZscoresRescreenOnly314ProteinCodingGenesForBayesian.csv") %>%
  select(EnsemblID, everything())
head(ZscoreData314RescreenOnly)
nrow(ZscoreData314RescreenOnly)
sum(ZscoreData314RescreenOnly$EnsemblID %in% positiveSet$Ensembl_accession)

#Get also the Zscores from the original screen
ZscoreDataOriginalScreenProteinCodingGenes <- fread("~/Documents/Project/Programming/MHC2loadingandSurfaceExpressionNeefjes/Data/OriginalZscores/originalZscoresForInclusionBayesian.csv")
ZscoreDataOriginalScreenProteinCodingGenes %<>%
  rename("EnsemblID" = "Gene stable ID", "OriginalZscoreHighestMedian" = "absoluteOfSelectedMedian") %>%
  select(EnsemblID, OriginalZscoreHighestMedian)

#Human tissue enrichment
HumanTissueData <- fread("~/Documents/Project/Programming/Human tissue enrichment/Data/TissueExpressionValuesForBayesianCSV.csv") %>%
  select(EnsemblID, EstimatedRankDifference, everything()) %>% rename(MannWhitneyUEstRankDiff = EstimatedRankDifference)
head(HumanTissueData, 20)
nrow(HumanTissueData)
range(HumanTissueData$MannWhitneyUEstRankDiff, na.rm = TRUE)
sum(HumanTissueData$EnsemblID %in% positiveSet$Ensembl_accession)

#Human macrophage activation profiles
HumanMacroPhageActivationData <- fread("~/Documents/Project/Programming/Human macrophage activation/Data/MacrophageActivationForBayesianCSV.csv")
head(HumanMacroPhageActivationData,20)
nrow(HumanMacroPhageActivationData)


#combine for final table
# FinalTableBayesian <- viralPPIData %>% left_join(TFBSScores, by = "EnsemblID") %>%
#   left_join(ZscoreData143HighConfidence, by = "EnsemblID") %>%
#   left_join(ZscoreData314RescreenOnly, by = "EnsemblID") %>%
#   left_join(HumanTissueData) %>%
#   left_join(HumanMacroPhageActivationData) %>%
#   left_join(ZscoreDataOriginalScreenProteinCodingGenes) %>%
#   rename(AntibodyHighConfidence = Antibody.x, AntibodyRescreenOnly = Antibody.y)
# head(FinalTableBayesian)
# nrow(FinalTableBayesian)

combineForFinalTable <- function(TFBSScores) {
  
  FinalTableBayesian <- viralPPIData %>% left_join(TFBSScores, by = "EnsemblID") %>%
    left_join(ZscoreData143HighConfidence, by = "EnsemblID") %>%
    left_join(ZscoreData314RescreenOnly, by = "EnsemblID") %>%
    left_join(HumanTissueData) %>%
    left_join(HumanMacroPhageActivationData) %>%
    left_join(ZscoreDataOriginalScreenProteinCodingGenes) %>%
    rename(AntibodyHighConfidence = Antibody.x, AntibodyRescreenOnly = Antibody.y)
  head(FinalTableBayesian) %>% print()
  nrow(FinalTableBayesian) %>% print()
  FinalTableBayesian
}

finalTableBayesianMHCIAndII <- combineForFinalTable(MHCIAndIIPositiveSetTFBSScores)
finalTableBayesianMHCI      <- combineForFinalTable(MHCIPositiveSetTFBSScores)
finalTableBayesianMHCII     <- combineForFinalTable(MHCIIPositiveSetTFBSScores)

write_csv(finalTableBayesianMHCIAndII, path = "~/Documents/Project/Programming/FinalBayesianClassifier/Data/BayesianClassifierCSVFileMHCIAndII.csv")
write_csv(finalTableBayesianMHCI     , path = "~/Documents/Project/Programming/FinalBayesianClassifier/Data/BayesianClassifierCSVFileMHCI.csv")
write_csv(finalTableBayesianMHCII    , path = "~/Documents/Project/Programming/FinalBayesianClassifier/Data/BayesianClassifierCSVFileMHCII.csv")
