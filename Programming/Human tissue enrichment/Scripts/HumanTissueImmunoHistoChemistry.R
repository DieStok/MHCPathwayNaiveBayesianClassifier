#Human tissue immunohistochemistry data
#select on lymphoid tissue, tonsil, thymus, bone marrow, lung macrophages
#see distribution of values. Author: Dieter Stoker

###########################################
###########################################
##                                       ##
##      Loading packages and data        ##
##                                       ##
###########################################
###########################################

library("pacman")
pacman::p_load(ggplot2, dplyr, stringr, purrr, readr, LaF, ff, data.table, R.utils, base, utils, devtools, microbenchmark, reshape2, magrittr, tidyr)
options(stringsAsFactors = FALSE)


tissueHistScores <- fread("/home/dieter/Documents/Project/Programming/Human tissue enrichment/Data/normal_tissue.tsv")
nrow(tissueHistScores)
head(tissueHistScores)


##############Data manipulation#############



tissueHistScores <- tissueHistScores %>% mutate(CellTypeTissue = paste(Tissue, `Cell type`)) %>%
  mutate(NumericScores = (as.numeric(factor(Level, levels = c("Not detected", "Low", "Medium", "High")))-1)) %>%
  filter(Reliability != "Uncertain")
print(unique(tissueHistScores$CellTypeTissue))
unique(tissueHistScores$`Cell type`)
length(unique(tissueHistScores$Tissue))
head(tissueHistScores)
tissueHistScores <-  tissueHistScores %>% transmute(Gene, NumericScores, CellTypeTissue) %>% tidyr::spread(CellTypeTissue, NumericScores)
head(tissueHistScores)
immuneScores <- tissueHistScores %>% select( `appendix lymphoid tissue` , `bone marrow hematopoietic cells` ,
                                            `lung macrophages` , `lymph node germinal center cells` ,
                                            `lymph node non-germinal center cells` , `skin 1 Langerhans` ,
                                            `spleen cells in white pulp` , `tonsil germinal center cells` ,
                                            `tonsil non-germinal center cells` , `thymus medullary cells` ,
                                            `thymus cortical cells`, `spleen cells in red pulp`, `tonsil squamous epithelial cells`)
head(immuneScores)
nonImmuneScores <- tissueHistScores %>% select(-`appendix lymphoid tissue` , -`bone marrow hematopoietic cells` ,
                                               -`lung macrophages` , -`lymph node germinal center cells` ,
                                               -`lymph node non-germinal center cells` , -`skin 1 Langerhans` ,
                                               -`spleen cells in white pulp` , -`tonsil germinal center cells` ,
                                               -`tonsil non-germinal center cells` , -`thymus medullary cells` ,
                                               -`thymus cortical cells`, -Gene, -`spleen cells in red pulp`, -`tonsil squamous epithelial cells`)
head(nonImmuneScores)



#use non-parametric test to test for sig. differences in protein stainings


mannWhitneyUTest <- wilcox.test(as.numeric(immuneScores[1,]), as.numeric(nonImmuneScores[1,]), conf.int = TRUE)
mannWhitneyUTest
mannWhitneyUTest$null.value
mannWhitneyUTest$p.value
mannWhitneyUTest$estimate
mannWhitneyUTest$statistic[[1]]

mWUTDataFrame <- as_tibble(data.frame("Ensembl-ID" = tissueHistScores$Gene, "Wstatistic" = 0, "p-value" = 0,
                                         "lower95%CI" = 0, "upper95%CI" = 0, "EstimatedRankDifference" = 0))
head(mWUTDataFrame); names(mWUTDataFrame)

for(row in seq_len(nrow(tissueHistScores))) {
  
  ensemblID             <- tissueHistScores[row,]$Gene
  
  #if test cannot be performed because no differences in rank
  if(sum(as.numeric(immuneScores[row,]), na.rm = TRUE) == sum(as.numeric(nonImmuneScores[row,]), na.rm = TRUE)) {
    
    mWUTDataFrame[row,]  <- c(ensemblID, NA, NA, NA, NA, NA)
    
  } else {
  
  mWUT                  <- wilcox.test(as.numeric(immuneScores[row,]), as.numeric(nonImmuneScores[row,]), conf.int = TRUE)
  confIntLow            <- as.numeric(mWUT$conf.int[[1]])
  confIntHigh           <- as.numeric(mWUT$conf.int[[2]])
  mWUTStat              <- as.numeric(mWUT$statistic[[1]])
  mWUTsig               <- as.numeric(mWUT$p.value)
  mWUTEstimate          <- as.numeric(mWUT$estimate)
  mWUTDataFrame[row,]  <- c(ensemblID, mWUTStat, mWUTsig, confIntLow, confIntHigh, mWUTEstimate) 
  }
  
}
head(mWUTDataFrame)
mWUTDataFrame$Wstatistic <- as.numeric(mWUTDataFrame$Wstatistic)
mWUTDataFrame$p.value <- as.numeric(mWUTDataFrame$p.value)
mWUTDataFrame$lower95.CI <- as.numeric(mWUTDataFrame$lower95.CI); mWUTDataFrame$upper95.CI <- as.numeric(mWUTDataFrame$upper95.CI)
mWUTDataFrame$EstimatedRankDifference <- as.numeric(mWUTDataFrame$EstimatedRankDifference)
head(mWUTDataFrame)
mWUTDataFrame$p.value.BH = p.adjust(mWUTDataFrame$p.value, method = "BH")
head(mWUTDataFrame)
nrow(mWUTDataFrame)


#####
#####OUTPUT FILE FOR BAYESIAN
colnames(immuneScores) <- paste0("ImmuneTissues_", gsub(" ", "_", colnames(immuneScores), fixed = TRUE))
head(immuneScores)
colnames(nonImmuneScores) <- paste0("NonImmuneTissues_", gsub(" ", "_", colnames(nonImmuneScores), fixed = TRUE))
head(nonImmuneScores)
nonImmuneScoresForBayesian <- nonImmuneScores %>% mutate(EnsemblID = tissueHistScores$Gene)
immuneScoresForBayesian <- immuneScores %>% mutate(EnsemblID = tissueHistScores$Gene)

totalTissueHistScoresForBayesian <- immuneScoresForBayesian %>% left_join(nonImmuneScoresForBayesian) %>% select(EnsemblID, everything()) 
ensemblGenesNotInTissueHistData <- allEnsemblGenes$`Gene stable ID`[allEnsemblGenes$`Gene stable ID` %ni% totalTissueHistScoresForBayesian$EnsemblID] %>%
  data.frame()
colnames(ensemblGenesNotInTissueHistData) <- c("EnsemblID")
head(ensemblGenesNotInTissueHistData)
totalTissueHistScoresForBayesian <- totalTissueHistScoresForBayesian %>% full_join(ensemblGenesNotInTissueHistData) %>% arrange(EnsemblID)

#add in also a column for mWUTscores
mWUTEstimatedRankDifferenceOnly <- mWUTDataFrame %>% select(Ensembl.ID, EstimatedRankDifference) %>%
  rename(EnsemblID = Ensembl.ID)
totalTissueHistScoresForBayesian <- totalTissueHistScoresForBayesian %>% full_join(mWUTEstimatedRankDifferenceOnly) %>% arrange(EnsemblID)
head(totalTissueHistScoresForBayesian)

write_csv(totalTissueHistScoresForBayesian, path = "~/Documents/Project/Programming/Human tissue enrichment/Data/TissueExpressionValuesForBayesianCSV.csv")


######





















sigMoreInImmune <- mWUTDataFrame %>% filter(p.value.BH <= 0.05, sign(EstimatedRankDifference) != -1)
head(sigMoreInImmune)
nrow(sigMoreInImmune)

####
####    Enrichment Stuff | Note that this is handled (with cross-validations) in ScriptForRunningTheBayesian.R, this is manual checking only
####    Therefore also not updated for additional (disparate) positive sets
####

positiveSet <- fread("~/Documents/Project/Programming/DataForAll/PositiveList_2.csv")
positiveSetIdCol <- "Ensembl_accession"
allEnsemblGenes <- fread("~/Documents/Project/Programming/DataForAll/EnsemblBioMartIDs_GRCh38.p10.txt", sep = ",", header = T)
allEnsemblGenesIDCol <- "Gene stable ID"
allEnsemblGenesNonPositiveTissueData <- allEnsemblGenes %>% filter(`Gene stable ID` %ni% positiveSet$Ensembl_accession) 
tissueIDCol <- "Ensembl.ID"
tissueValueCol <- "sampleEstDifference"
higherExpressedImmuneTissues <- as.data.table(sigMoreInImmune)


PositiveSetTissue <- higherExpressedImmuneTissues[higherExpressedImmuneTissues$Ensembl.ID %in% positiveSet$Ensembl_accession,]$EstimatedRankDifference
PositiveSetTissue
length(PositiveSetTissue)
OtherTissue       <- higherExpressedImmuneTissues[higherExpressedImmuneTissues$Ensembl.ID %ni% positiveSet$Ensembl_accession,]$EstimatedRankDifference
OtherTissue


dataForPlottingTissue <- data.frame(higherExpressedImmuneTissues$EstimatedRankDifference) %>%
  mutate(set = if_else(higherExpressedImmuneTissues$Ensembl.ID %in% positiveSet$Ensembl_accession, "PS", "NS")) %>%
  rename(EstimatedRankDifference = higherExpressedImmuneTissues.EstimatedRankDifference) %>%
  mutate(EnsemblID = higherExpressedImmuneTissues$Ensembl.ID)  %>%
  arrange(EnsemblID)

head(dataForPlottingTissue)
range(dataForPlottingTissue$EstimatedRankDifference)
startBreaks <- seq(0, 3, 0.5)
furtherBreaks <- c(-9, 0.0, 0.5, 1.0, 1.5, 3.0)


totalPlotTissueBins <- dataForPlottingTissue %>%
  ggplot(aes(x = EstimatedRankDifference, fill = set, y = ..count..)) + geom_histogram(breaks = furtherBreaks, colour = "black") 
totalPlotTissueBins



dataTotalHistTissue <- ggplot_build(totalPlotTissueBins)

getScoresAdapted <- function(plotBuild) {
  
  logScores <- numeric()
  for(i in seq(2, nrow(plotBuild$data[[1]]), by = 2)) {
    
    logScores = c(logScores, log2((plotBuild$data[[1]][(i-1),]$count/nrow(positiveSet))/(plotBuild$data[[1]][i,]$count/nrow(allEnsemblGenesNonPositiveTissueData))))
    
  }
  newPlot <- plotBuild$plot + annotate("text", x = unique(plotBuild$data[[1]]$x), y = 100, label = round(logScores, digits = 3), angle = 90)
  
  list(plot = newPlot, scores = logScores)
  
}


scoresTissue <- getScoresAdapted(dataTotalHistTissue)


scoresTissue$plot

