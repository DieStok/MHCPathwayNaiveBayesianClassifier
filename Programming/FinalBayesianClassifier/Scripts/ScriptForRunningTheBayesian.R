###############################
##Final script Bayesian Table##
###############################

# Wat ik wil doen:
# 
#   Definieer 10 subsets van de data, dus steeds 90% van de PS en 90% van de NS, met de overige 10% waar je niet op traint.
#   Loop over elke kolom van de Bayesian tabel heen. Bereken de scores voor alle genen. Houd per kolom de ranges
#   Van scores bij, of de log2 likelihoods voor de bins, etc. Je wil ook een boxplot met per keer dat de classifier runt
#   de range van scores van positieve en negatieve set.
#   Uiteindelijk krijg je een totale score en kan je de genen sorteren van hoog naar laag.
#   Vervolgens Geef je ze een rank (van 1 naar 22357). Dan selecteer je uit de zo verkregen tabel, waarbij voor
#   elke kolom de scores zijn berekend (bij TFBS op 2 manieren ook), de rijen van de testset (waarop niet getraind is).
#   Dat doe je dus bij de CrossVal heel de tijd, en uiteindelijk heb je dan je totale tabel met totalscores en ranks.
  
#Noot: check eerst nog of alle EnsemblIDS in de positive set kloppen, dat deed raar. En check ook of de ENsemblIDs
#van de human protein tissue atlas wel kloppen






######################################################################################################
######################################################################################################
######                                                                                          ######
######                                Loading Packages and Options                              ######
######                                                                                          ######
######################################################################################################
######################################################################################################


options(stringsAsFactors = FALSE)
install.packages("pacman", repos = "https://cloud.r-project.org/")
library("pacman")
pacman::p_load(tidyverse, stringr, purrr, readr, LaF, ff, data.table, R.utils, base, utils, reshape2, optparse, cowplot, ggplot2,
               dplyr, devtools, microbenchmark, magrittr, pROC, plotROC, lvplot, ggbeeswarm, RColorBrewer, superheat, rlang)
options(datatable.verbose = FALSE)


`%ni%` <- Negate(`%in%`)
set.seed(1234567890)

######################################################################################################
######################################################################################################
######                                                                                          ######
######                                  USER-DEFINED VARIABLES                                  ######
######                                                                                          ######
######################################################################################################
######################################################################################################

#input files
positiveSetFile            <- "~/Documents/Project/Programming/DataForAll/PositiveList_2.csv"
positiveSetFileMHCIOnly    <- "/home/dieter/Documents/Project/Programming/DataForAll/PositiveList_MHCIOnly.csv"
positiveSetFileMHCIIOnly   <- "/home/dieter/Documents/Project/Programming/DataForAll/PositiveList_MHCIIOnly.csv"
bayesianTableCSVTotal      <- "~/Documents/Project/Programming/FinalBayesianClassifier/Data/BayesianClassifierCSVFileMHCIAndII.csv"
bayesianTableCSVMHCIOnly   <- "~/Documents/Project/Programming/FinalBayesianClassifier/Data/BayesianClassifierCSVFileMHCI.csv"
bayesianTableCSVMHCIIOnly  <- "~/Documents/Project/Programming/FinalBayesianClassifier/Data/BayesianClassifierCSVFileMHCII.csv"
allEnsemblIdsFile          <- "~/Documents/Project/Programming/DataForAll/EnsemblBioMartIDs_GRCh38.p10.txt"

######################################################################################################
######################################################################################################
######                                                                                          ######
######                            READING IN FILES + SANITY CHECKS                              ######
######                                                                                          ######
######################################################################################################
######################################################################################################

bayesianTableTotal <- fread(bayesianTableCSVTotal)
bayesianTableCore  <- bayesianTableTotal %>% select(EnsemblID, UnionMeasuredPPI, IntersectMeasuredPPI,
                                                    totalTFBSScore, MannWhitneyUEstRankDiff, Profile, 
                                                    ZscoreHighScore, RescreenOnlyZscores, OriginalZscoreHighestMedian)

#check that HLA-DMA has a score in the ZscoreHighscore and RescreenOnly data
bayesianTableCore[bayesianTableCore$EnsemblID == "ENSG00000204257",]

positiveSet        <- fread(positiveSetFile)
allEnsemblIds      <- fread(allEnsemblIdsFile, sep = ",", header = T)


nrow(bayesianTableTotal); nrow(bayesianTableCore)
head(bayesianTableTotal, 20); tail(bayesianTableTotal, 20)
head(positiveSet); tail(positiveSet, 20)
head(allEnsemblIds)
sum(positiveSet$Ensembl_accession %in% bayesianTableTotal$EnsemblID)
sum(positiveSet$Ensembl_accession %in% allEnsemblIds$`Gene stable ID`)
#everything is in there

#note: the TFBS scores are only present for 18990 genes. This is because the ENCODE TFBS motifs were mapped
#on hg19, which is equal to GRCh37. Those motifs were assigned to EnsemblIDs based on that, and the resulting
#IDs were translated to Grch38.p10 (I HOPE, CHECK THIS!!!!)


#Add in information about positive set and negative set
bayesianTableCore <- bayesianTableCore %>% mutate(Set = ifelse(bayesianTableCore$EnsemblID %in% positiveSet$Ensembl_accession, "PS", "NS")) %>%
  select(EnsemblID, Set, everything())

#change Zscores to absolute
bayesianTableCore %<>% mutate(ZscoreHighScore = abs(ZscoreHighScore), RescreenOnlyZscores = abs(RescreenOnlyZscores))

bayesianTableTotal <- bayesianTableTotal %>% mutate(Set = ifelse(bayesianTableTotal$EnsemblID %in% positiveSet$Ensembl_accession, "PS", "NS")) %>%
  select(EnsemblID, Set, everything())


####rewrite data entry: make functions to ease read-in for all positive sets

readInBayesianCSVFiles <- function(positiveSetFile = positiveSetFile,
                                   bayesianCSVFile = bayesianTableCSVTotal,
                                   allEnsemblIDFile = allEnsemblIdsFile) {
  
  totalBayesianTable <- fread(bayesianCSVFile)
  coreBayesianTable  <- totalBayesianTable %>% dplyr::select(EnsemblID, UnionMeasuredPPI, IntersectMeasuredPPI,
                                                      totalTFBSScore, MannWhitneyUEstRankDiff, Profile, 
                                                      ZscoreHighScore, RescreenOnlyZscores, OriginalZscoreHighestMedian)
  
  #sanity checks
  #check that HLA-DMA has a score in the ZscoreHighscore and RescreenOnly data
  print(coreBayesianTable[coreBayesianTable$EnsemblID == "ENSG00000204257",])
  
  positiveSet        <- fread(positiveSetFile)
  allEnsemblIds      <- fread(allEnsemblIDFile, sep = ",", header = T)
  
  
  print(nrow(totalBayesianTable)); print(nrow(coreBayesianTable))
  print(head(totalBayesianTable, 20)); print(tail(coreBayesianTable, 20))
  print(head(positiveSet)); print(tail(positiveSet, 20))
  print(head(allEnsemblIds))
  print(sum(positiveSet$Ensembl_accession %in% totalBayesianTable$EnsemblID))
  print(sum(positiveSet$Ensembl_accession %in% allEnsemblIds$`Gene stable ID`))
  
  
  #Add in information about positive set and negative set
  coreBayesianTable %<>% mutate(Set = ifelse(coreBayesianTable$EnsemblID %in% positiveSet$Ensembl_accession, "PS", "NS")) %>%
    select(EnsemblID, Set, everything())
  
  #change Zscores to absolutes (only interested in maximum perturbation, not sign (more or less expression))
  coreBayesianTable %<>% mutate(ZscoreHighScore = abs(ZscoreHighScore), RescreenOnlyZscores = abs(RescreenOnlyZscores))
  
  totalBayesianTable %<>% mutate(Set = ifelse(totalBayesianTable$EnsemblID %in% positiveSet$Ensembl_accession, "PS", "NS")) %>%
    select(EnsemblID, Set, everything())
 
  returnList <- list(totalBayesianTable = totalBayesianTable,
                     coreBayesianTable  = coreBayesianTable)
   
}



MHCIAndIIBayesianTable <- readInBayesianCSVFiles(positiveSetFile = positiveSetFile,
                                                 bayesianCSVFile = bayesianTableCSVTotal,
                                                 allEnsemblIDFile = allEnsemblIdsFile)

MHCIBayesianTable      <- readInBayesianCSVFiles(positiveSetFile = positiveSetFileMHCIOnly,
                                                 bayesianCSVFile = bayesianTableCSVMHCIOnly,
                                                 allEnsemblIDFile = allEnsemblIdsFile)

MHCIIBayesianTable     <- readInBayesianCSVFiles(positiveSetFile = positiveSetFileMHCIIOnly,
                                                 bayesianCSVFile = bayesianTableCSVMHCIIOnly,
                                                 allEnsemblIDFile = allEnsemblIdsFile)

######################################################################################################
######################################################################################################
######                                                                                          ######
######                                Data Manipulation Functions                               ######
######                                                                                          ######
######################################################################################################
######################################################################################################




#function to generate the subsets of the table for cross-validation

createCrossVal <- function(fold, positiveSetIDs, negativeSetIDs, bayesianTable) {
  
  
  if(! is.numeric(fold) ) stop("Fold should be a number")
  if(fold <= 0 | fold > length(positiveSetIDs))  stop("Fold should be a number larger than 0 and less than the total size of the positive set")
  if(! is.character(positiveSetIDs) | ! is.character(negativeSetIDs)) stop("Both positiveSetIDs and negativeSetIDs need to be vectors of ensembl IDs")
  
  positiveSetNotTrainedOn <- vector("list", length = fold)
  positiveSetTrainedOn    <- vector("list", length = fold)
  #negative sets
  negativeSetNotTrainedOn <- vector("list", length = fold)
  negativeSetTrainedOn    <- vector("list", length = fold)
  #total Set
  totalSetTrainedOn       <- vector("list", length = fold)
  totalSetNotTrainedOn    <- vector("list", length = fold)
  fullTotalSet            <- vector("list", length = fold)
  #finalScoreTable
  totalScoreTable         <- vector("list", length = fold)
  extraInfoPerCrossval    <- vector("list", length = fold)
  
  
  #function that generates groups according to the fold crossval
  shuffleIdsIntoGroups    <- function(IdVector, groups = fold) {
    
    IdVector   <- sample(IdVector)
    safe_split <- quietly(split)
    returnList <- safe_split(IdVector, rep(1:groups, ceiling(length(IdVector)/groups)))
    returnList
    
  }
  
  #list of shuffled Ids, length fold, where each entry contains set size/fold IDs
  positiveSetShuffleList <- shuffleIdsIntoGroups(positiveSetIDs, groups = fold)
  negativeSetShuffleList <- shuffleIdsIntoGroups(negativeSetIDs, groups = fold)
  
  #generate framework for the final score table per crossval (i.e, the total table is filled in each time,
  #but only examples that were not trained on in that crossval are kept in the final table)
  finalBayesianScoreTableFramework <- bayesianTable %>% select(EnsemblID, Set) %>% mutate(PPIScoreUnion     = 0,
                                                                                          PPIScoreIntersect = 0,
                                                                                          TFBSScoreBins     = 0,
                                                                                          TFBSScoreFit      = 0,
                                                                                          ZscoreHighConf    = 0,
                                                                                          ZscoreLowerConf   = 0,
                                                                                          ImmuneTissueScore = 0,
                                                                                          MacrophageExpProf = 0
                                                                                          )
  
  
  #assign this to the lists defined previously
  if(fold == 1) {
    
    for(i in seq_len(fold)) {
      
      positiveSetTrainedOn[[i]]        <- bayesianTable %>% filter(EnsemblID %in% positiveSetShuffleList$result[[i]])   
      positiveSetNotTrainedOn[[i]]           <- bayesianTable %>% filter(EnsemblID %ni% positiveSetShuffleList$result[[i]], EnsemblID %ni% negativeSetIDs) 
      negativeSetTrainedOn[[i]]        <- bayesianTable %>% filter(EnsemblID %in% negativeSetShuffleList$result[[i]])
      negativeSetNotTrainedOn[[i]]           <- bayesianTable %>% filter(EnsemblID %ni% negativeSetShuffleList$result[[i]], EnsemblID %ni% positiveSetIDs)
      totalSetTrainedOn[[i]]              <- rbind(positiveSetTrainedOn[[i]], negativeSetTrainedOn[[i]]) %>% arrange(EnsemblID)
      totalSetNotTrainedOn[[i]]           <- rbind(positiveSetNotTrainedOn[[i]], negativeSetNotTrainedOn[[i]]) %>% arrange(EnsemblID)
      fullTotalSet[[i]]                   <- rbind(totalSetTrainedOn[[i]], totalSetNotTrainedOn[[i]]) %>% arrange(EnsemblID)
      totalScoreTable[[i]]                <- finalBayesianScoreTableFramework
    }
    
  }
  
  if(fold != 1) {
  for(i in seq_len(fold)) {
    
    positiveSetNotTrainedOn[[i]]        <- bayesianTable %>% filter(EnsemblID %in% positiveSetShuffleList$result[[i]])   
    positiveSetTrainedOn[[i]]           <- bayesianTable %>% filter(EnsemblID %ni% positiveSetShuffleList$result[[i]], EnsemblID %ni% negativeSetIDs) 
    negativeSetNotTrainedOn[[i]]        <- bayesianTable %>% filter(EnsemblID %in% negativeSetShuffleList$result[[i]])
    negativeSetTrainedOn[[i]]           <- bayesianTable %>% filter(EnsemblID %ni% negativeSetShuffleList$result[[i]], EnsemblID %ni% positiveSetIDs)
    totalSetTrainedOn[[i]]              <- rbind(positiveSetTrainedOn[[i]], negativeSetTrainedOn[[i]]) %>% arrange(EnsemblID)
    totalSetNotTrainedOn[[i]]           <- rbind(positiveSetNotTrainedOn[[i]], negativeSetNotTrainedOn[[i]]) %>% arrange(EnsemblID)
    fullTotalSet[[i]]                   <- rbind(totalSetTrainedOn[[i]], totalSetNotTrainedOn[[i]]) %>% arrange(EnsemblID)
    totalScoreTable[[i]]                <- finalBayesianScoreTableFramework
  }
  }
  
  
  list(
    positiveSetNotTrainedOn = positiveSetNotTrainedOn,
    positiveSetTrainedOn    = positiveSetTrainedOn,
    negativeSetNotTrainedOn = negativeSetNotTrainedOn,
    negativeSetTrainedOn    = negativeSetTrainedOn,
    totalSetTrainedOn       = totalSetTrainedOn,
    totalSetNotTrainedOn    = totalSetNotTrainedOn,
    fullTotalSet            = fullTotalSet,
    totalScoreTable         = totalScoreTable,
    extraInfoPerCrossval    = extraInfoPerCrossval
    
  )
  
  
}


#function that can take binned scores that are calculated in all functions that use binning
#and draw a nice bar chart with log 2 enrichment scores and density of positive and negative set per bin


generateBinningPlots <- function(data, title = "Title", xlabel = "Bins", ylabel = "Density of scores in bin per set", log2ScoresColumn = 4, groupColumn = 1, plotLabels = TRUE) {
  
  listReturn <- vector("list", length = 2)
  
  dataWithDensity <- data %>% group_by(Set) %>% mutate(densityPerSet = Count/sum(Count))
  dataWithDensity <- rbind(dataWithDensity[(nrow(dataWithDensity)-1):nrow(dataWithDensity),], dataWithDensity[1:(nrow(dataWithDensity)-2),])
  dataWithDensity$gr <- as.character(dataWithDensity$gr)
  dataWithDensity$gr[is.na(dataWithDensity$gr)] <- "None"
  print(dataWithDensity$gr)
  
  #print(dataWithDensity[,groupColumn])
  #print(dataWithDensity[,log2ScoresColumn])
  
  #print(unique(as.character(dataWithDensity$gr))[!is.na(unique(as.character(dataWithDensity$gr)))])
  #realGroups = unique(as.character(dataWithDensity$gr))[!is.na(unique(as.character(dataWithDensity$gr)))]
  #xFactor = factor(dataWithDensity$gr, exclude = NULL, levels = c(NA, realGroups), ordered = TRUE)
  #print(xFactor)
  #ordered so that NA is in front
  #print(xFactor)
  
  
  #breaks = unique(as.character(xFactor))
  
  plot <- ggplot(dataWithDensity, aes(x = factor(gr, levels = c("None", unique(gr)[unique(gr) != "None"])), y = densityPerSet, fill = Set)) +
    geom_bar(stat = "identity", colour = "black", width = 1, position = "dodge") +
    theme_bw() + theme(panel.grid.major.x = element_blank()) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_discrete(expand   = c(0,0),   name = xlabel) +
    scale_y_continuous(expand = c(0,0),   name = ylabel, limits = c(0,1), breaks = seq(0,1,0.1)) +
    ggtitle(label = title) +
    scale_fill_manual(breaks = c("PS","NS"),
                      values = c("#686caa", "#eabc4f")) #+
    #scale_x_manual
  
  if(plotLabels == TRUE) {
    
    plot <- plot +
      annotate("text", x = as.character(dataWithDensity[,groupColumn][[1]][seq(1,nrow(dataWithDensity), by = 2)]), y = 0.8,
                            label = round(dataWithDensity[,log2ScoresColumn][[1]][seq(1,nrow(dataWithDensity), by = 2)], digits = 2), angle = 90, size = 7) 
  }

  
  listReturn[[1]] = dataWithDensity
  listReturn[[2]] = plot
  listReturn
  
}

computeTFBSScoresBins <- function(bins, crossValOutput) {
  
  
  #do everything for the amount of times crossval is performed
  for (crossValSet in seq_len(length(crossValOutput$positiveSetNotTrainedOn))) {
  
    
    
  #define positive set and negative set genes that are used for training in this crossVal
  positiveSetGenes <- crossValOutput$positiveSetTrainedOn[[crossValSet]]$EnsemblID
  negativeSetGenes <- crossValOutput$negativeSetTrainedOn[[crossValSet]]$EnsemblID
  
  #compute binned scores. Always also include binning intervals of 1 to compare what manual binning changes
  
  #custom bins
  customBinScoresNotSummarizedTrain <- crossValOutput$totalSetTrainedOn[[crossValSet]] %>% select(EnsemblID, Set, totalTFBSScore) %>%
    group_by(gr = cut(totalTFBSScore, breaks = bins), Set)
  customBinScoresNotSummarizedTotal <- crossValOutput$fullTotalSet[[crossValSet]]      %>% select(EnsemblID, Set, totalTFBSScore) %>%
    group_by(gr = cut(totalTFBSScore, breaks = bins), Set)
  customBinScoresSummarized         <- customBinScoresNotSummarizedTrain %>% summarize(Count = n())
  
  
  #bins of 1 width
  oneWidthBinScoresNotSummarizedTrain <- crossValOutput$totalSetTrainedOn[[crossValSet]] %>% select(EnsemblID, Set, totalTFBSScore) %>%
    group_by(gr = cut(totalTFBSScore, breaks = seq(floor(min(totalTFBSScore, na.rm = TRUE)), ceiling(max(totalTFBSScore, na.rm = TRUE)))), Set)
  
  oneWidthBinScoresNotSummarizedTotal <- crossValOutput$fullTotalSet[[crossValSet]] %>% select(EnsemblID, Set, totalTFBSScore) %>%
    group_by(gr = cut(totalTFBSScore, breaks = seq(floor(min(totalTFBSScore, na.rm = TRUE)), ceiling(max(totalTFBSScore, na.rm = TRUE)))), Set)
  oneWidthBinScoresSummarized <- oneWidthBinScoresNotSummarizedTrain  %>% summarize(Count = n())
  
  #problem: not every bin contains both PS and NS scores. For score calculations, that is necessary. Need to add these in.
  #hence, function below:
  #for a grouped data frame of bins and the counts for positive and negative set genes, adds in 
  #an entry with the same bin, the missing set, and a count of 0
  addMissingSetCounts <- function(summarizedData) {
    
    setVector <- c("PS", "NS")
    for(set in seq_len(length(setVector))) {
      for(i in seq_len(length(unique(summarizedData$gr)))) {
        if(!is.na(unique(summarizedData$gr)[i])) {
          if(summarizedData %>% filter(gr == unique(summarizedData$gr)[i], Set == setVector[set]) %>% `$`(Count) %>% length() == 0 ) {
            summarizedData <- summarizedData %>% ungroup %>% add_row(gr = unique(summarizedData$gr)[i], Set = setVector[set], Count = 0) 
          }
        }
      }
      #for NA
      if(summarizedData %>% filter(is.na(gr) == TRUE, Set == setVector[set]) %>% `$`(Count) %>% length() == 0 ) {
        summarizedData <- summarizedData %>% ungroup %>% add_row(gr = NA, Set = setVector[set], Count = 0) }
    }
    summarizedData <- summarizedData %>% group_by(gr) %>% arrange(gr)
    summarizedData
  }
  
  customBinScoresSummarized   <- addMissingSetCounts(customBinScoresSummarized)
  oneWidthBinScoresSummarized <- addMissingSetCounts(oneWidthBinScoresSummarized)
  
  #supply a vector of positive Set and of negative Set genes. Insert the summarized bin data from above.
  #per bin, divide PS counts by PSsize, and divide NS counst by NSsize. Then divide those two by each other and take log2
  #special approach was taken to include NAs. (seq_len(i-1) selects the group which is NA)
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
  
  
  log2scoresCustomBins         <- getLog2Scores(customBinScoresSummarized, positiveSetGenes, negativeSetGenes)
  customBinScoresSummarized$log2ScoresTFBSBins   <- rep(log2scoresCustomBins, each = 2)
  log2scoresOneWidthBins       <- getLog2Scores(oneWidthBinScoresSummarized, positiveSetGenes, negativeSetGenes)
  oneWidthBinScoresSummarized$log2ScoresTFBSBins <- rep(log2scoresOneWidthBins, each = 2)
  
  #head(customBinScoresSummarized); tail(customBinScoresSummarized)
  #head(oneWidthBinScoresSummarized); tail(oneWidthBinScoresSummarized)
  
  
  #now give everything a score (traind on + not trained on)
  #only not trained on is kept from this crossVal run, but scores are assigned to everything to make a rank across crossVals.
  #Note that at present, NA are given a score. These are genes that were not in the TFBS dataset. SInce only 3 of the positive
  #set were not in the dataset, and many from the negative set, these are given a negative score. I need to discuss with John whether
  #this is really pertinent or not.
  
  
  
  assignScoresBinsToAllGenes          <- function(summarizedBinScores, notSummarizedBinScores) {
    notSummarizedBinScores$log2ScoresTFBSBins = 0
    print(head(notSummarizedBinScores))
    for (group in seq_len(length(unique(notSummarizedBinScores$gr)))) {
      #print(group)
      #print(is.na(unique(notSummarizedBinScores$gr)[group]))
      
      scoreToPlace <- summarizedBinScores %>%
        filter(gr %in% unique(notSummarizedBinScores$gr)[group], Set == "PS") %>% select(log2ScoresTFBSBins)  %>% as.data.frame() %>% `[`(1,2, drop = TRUE)
      #print(scoreToPlace)
      #print(head(notSummarizedBinScores[notSummarizedBinScores$gr %in% unique(notSummarizedBinScores$gr)[group],]))
      notSummarizedBinScores[notSummarizedBinScores$gr %in% unique(notSummarizedBinScores$gr)[group],]$log2ScoresTFBSBins = scoreToPlace
      #print(notSummarizedBinScores[notSummarizedBinScores$gr %in% unique(notSummarizedBinScores$gr)[group],])
      
    }
    
    notSummarizedBinScores %>% arrange(EnsemblID)
    
  }
  customBinScoresNotSummarizedTotal   <- assignScoresBinsToAllGenes(customBinScoresSummarized  , customBinScoresNotSummarizedTotal  )
  oneWidthBinScoresNotSummarizedTotal <- assignScoresBinsToAllGenes(oneWidthBinScoresSummarized, oneWidthBinScoresNotSummarizedTotal)
  
  # head(oneWidthBinScoresNotSummarizedTotal);                  tail(oneWidthBinScoresNotSummarizedTotal)
  # head(as.data.frame(customBinScoresNotSummarizedTotal), 50); tail(as.data.frame(customBinScoresNotSummarizedTotal), 50)
  # 
  
  
  print(customBinScoresSummarized)
  pieter <<- customBinScoresSummarized
  
  #make plots and add density values
  customBinScoresTFBSResult   <- (generateBinningPlots(customBinScoresSummarized,
                                                             title = "Custom bins TFBS score",
                                                             xlabel = "Binned TFBS scores"))
  
  oneWidthBinScoresTFBSResult <- (generateBinningPlots(oneWidthBinScoresSummarized,
                                                             title = "One width bins TFBS score",
                                                             xlabel = "Binned TFBS scores"))
  
  #Use the customBinScores, keep the oneWidthBinScores for reference. Keep also the summarizedTables.
  
  crossValOutput$totalScoreTable[[crossValSet]]$TFBSScoreBins <- customBinScoresNotSummarizedTotal$log2ScoresTFBSBins
  crossValOutput$extraInfoPerCrossval[[crossValSet]]          <- list(customBinScoresTFBSCountsPerBin   = customBinScoresTFBSResult[[1]],
                                                                      customBinScoresTFBSPlot           = customBinScoresTFBSResult[[2]],
                                                                      oneWidthBinScoresTFBSCountsPerBin = oneWidthBinScoresTFBSResult[[1]],
                                                                      oneWidthBinScoresTFBSPlot         = oneWidthBinScoresTFBSResult[[2]],
                                                                      oneWidthBinScoresTFBS             = oneWidthBinScoresNotSummarizedTotal)
  
  }
  crossValOutput
  
}

computeTFBSScoresPolynomial <- function(crossValOutput, smoothBandWidth = 2.5) {
  
  
  #calculate density of negative and positive training set,
  #make a model (20th order polynomial),
  #predict values for genes in each run, but ultimately keep only those from the test set per crossVal run
  
  for(crossValSet in seq_len(length(crossValOutput$positiveSetNotTrainedOn))) {
    
    print(paste0("Working on crossVal ", crossValSet, " out of ", length(crossValOutput$positiveSetNotTrainedOn),  "."))
    
    positiveSetTrainTFBSScores <- crossValOutput$positiveSetTrainedOn[[crossValSet]]$totalTFBSScore
    negativeSetTrainTFBSScores <- crossValOutput$negativeSetTrainedOn[[crossValSet]]$totalTFBSScore
    
    #use kernel smoothing to get a smoothened density
    densityPositiveTrain       <- density(positiveSetTrainTFBSScores, bw = smoothBandWidth, na.rm = TRUE)
    densityNegativeTrain       <- density(negativeSetTrainTFBSScores, bw = smoothBandWidth, na.rm = TRUE)
    
    #make a 20th order polynomial fit of the predictions (i.e. a very smooth line through the points)
    dfPositive                 <- data.frame(x = densityPositiveTrain$x, y = densityPositiveTrain$y)
    polynomialPositive         <- lm(lm(y ~ poly(x, 20), data = dfPositive))
    dfNegative                 <- data.frame(x=densityNegativeTrain$x, y=densityNegativeTrain$y)
    polynomialNegative         <- lm(y ~ poly(x, 20), data = dfNegative)
    
    
    #Predict stuff
    positivePolyScoresAllGenes <- predict(polynomialPositive, data.frame(x = crossValOutput$fullTotalSet[[crossValSet]]$totalTFBSScore))
    negativePolyScoresAllGenes <- predict(polynomialNegative, data.frame(x = crossValOutput$fullTotalSet[[crossValSet]]$totalTFBSScore))
    
    #get log2 overrepresentation scores
    log2ScoresPoly            <- log2(positivePolyScoresAllGenes/negativePolyScoresAllGenes)
    #print("Head and tail of polynomial log2 scores:")
    #print(head(log2ScoresPoly, 20)); print(tail(log2ScoresPoly, 20))
    
    
    #plot(densityNegativeTrain$y ~ densityNegativeTrain$x, col = "red", cex = 0.2)
    #lines(predict(polynomialNegative) ~ densityNegativeTrain$x, col = "darkred")
    #points(densityPositiveTrain$y ~ densityPositiveTrain$x, col = "green", cex = 0.2)
    #lines(predict(polynomialPositive) ~ densityPositiveTrain$x, col = "darkgreen")
   
    #save a ggplot object of the scores
    ggPlotDataFrame <- data.frame(y = c(densityNegativeTrain$y, densityPositiveTrain$y), x = c(densityNegativeTrain$x, densityPositiveTrain$x),
                                  posNeg = c(rep("negative", length(densityNegativeTrain$y)), rep("positive", length(densityPositiveTrain$y))))
    ggPlotDataFramePosOnly <- ggPlotDataFrame %>% filter(posNeg == "positive") %>% mutate(prediction = predict(polynomialPositive))
    ggPlotDataFrameNegOnly <- ggPlotDataFrame %>% filter(posNeg == "negative") %>% mutate(prediction = predict(polynomialNegative))
    head(ggPlotDataFrame); tail(ggPlotDataFrame, 30)
    polynomialDensityPlot <- ggplot(data = ggPlotDataFrame, aes(y = y, x = x, fill = factor(posNeg, levels = c("negative", "positive")))) +
      geom_area(data = ggPlotDataFramePosOnly, aes(y = prediction, x = x),  alpha = 0.8) +
      geom_area(data = ggPlotDataFrameNegOnly, aes(y = prediction, x = x),  alpha = 0.65) +
      geom_point(data = ggPlotDataFrameNegOnly, aes(y = y), colour = "#686caa", size = 0.5) +
      geom_point(data = ggPlotDataFramePosOnly, aes(y = y), colour = "#eabc4f", size = 0.5) +
      geom_line(data = ggPlotDataFrameNegOnly, aes(y = prediction), colour = "#686caa", size = 0.5) +
      geom_line(data = ggPlotDataFramePosOnly, aes(y = prediction), colour = "#eabc4f", size = 0.5) +
      scale_fill_manual(breaks = c("negative", "positive"), values = c("positive" = "#eabc4f", "negative" = "#686caa")) +
      guides(fill = guide_legend(title = NULL)) +
      scale_x_continuous(name = "total TFBS score", expand = c(0,0)) +
      scale_y_continuous(name = "density of values", expand = c(0,0)) +
      coord_cartesian(ylim = c(0,0.12)) +
      #ggtitle("Smoothened distributions of scores for TFBS")  +
      theme_bw() +  #theme(plot.title = element_text(hjust = 0.5), legend.justification = c(1.5,1.5), legend.position = c(0.92,0.5), legend.text = element_text(size = 14)) +
      theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) 
      #annotate("text", scoreSeq, 0.09, label = round(log2Seq, 3), angle = 90, size = 8) + annotate("text", 0, 0.11, label = "Log2 Enrichment Scores", size = 9) + panel_border(remove = TRUE) 
      polynomialDensityPlot
    
      #scale_fill_manual(breaks = c("PS","NS"),
      #                  values = c("#686caa", "#eabc4f"))
      
    
    #I now take all the data together, and bound the data on a minimum and maximum. SInce the positive set and negative set
    #don't completely overlap, the fit performs badly on the fringes. e.g., positive set genes with very high scores get
    #extremely high log likelihoods. To constrain that, I limit this to the region where the behaviour of the division of
    #the distributions by each other is still normal. I find a minimum and maximum within these fringes, and bound on that.
    #This still constitutes a sort of binning, but there is nothing for it: it must be done.
    dFScores <- data.frame(ID                      = crossValOutput$fullTotalSet[[crossValSet]]$EnsemblID,
                           TFBSScore               = crossValOutput$fullTotalSet[[crossValSet]]$totalTFBSScore,
                           TFBSPredictScore        = log2ScoresPoly,
                           TFBSPredictScoreBounded = log2ScoresPoly,
                           PositivePredicted       = positivePolyScoresAllGenes,
                           NegativePredicted       = negativePolyScoresAllGenes)
    
    
    #24-10-2018: De bound is bij de positive set bestaande uit MHC I en MHC II handmatig vastgesteld, door te kijken naar de unbounded scores en waar die
    #helemaal uit de hand beginnen te lopen (i.e. enorm snel rijzen aan de randen van de distributies). Dat is onhoudbaar als er verschillende
    #positive sets zijn dus ik moet iets anders bedenken. Als de density minder wordt dan 2.5 * 10^-3 dan bound ik de waardes. Dus check op
    #PositivePredicted en NegativePredicted, en zodra een lager wordt dan dat getal bound je alle waardes van dan af.
    
    
    print(head(dFScores))
    print(tail(dFScores))
    print(min(dFScores$TFBSPredictScore, na.rm = TRUE))
    print(max(dFScores$TFBSPredictScore, na.rm = TRUE))
    #Sys.sleep(5)
    minimumBound          <- min(dFScores[dFScores$TFBSScore >= -15 & dFScores$TFBSScore <= -10,]$TFBSPredictScore, na.rm = TRUE)
    minimumBoundTFBSScore <- dFScores[dFScores$TFBSPredictScore %in% minimumBound,]$TFBSScore
    maximumBound          <- max(dFScores[dFScores$TFBSScore >=  10 & dFScores$TFBSScore <=  15,]$TFBSPredictScore, na.rm = TRUE)
    maximumBoundTFBSScore <- dFScores[dFScores$TFBSPredictScore %in% maximumBound,]$TFBSScore
    dFScores[dFScores$TFBSScore > maximumBoundTFBSScore & !is.na(dFScores$TFBSScore), ]$TFBSPredictScoreBounded <- maximumBound
    dFScores[dFScores$TFBSScore < minimumBoundTFBSScore & !is.na(dFScores$TFBSScore), ]$TFBSPredictScoreBounded <- minimumBound
    
    #Bounding functions correctly. Output the minimum and maximum bound for checking. Also say how many values are affected by that bound
    #In addition, output plots.
    numericVectorBounding <- c(minimumBound = minimumBound, maximumBound = maximumBound,
                               minBoundedValueCount = nrow(dFScores[dFScores$TFBSScore < minimumBoundTFBSScore & !is.na(dFScores$TFBSScore), ]),
                               maxBoundedValueCount = nrow(dFScores[dFScores$TFBSScore > maximumBoundTFBSScore & !is.na(dFScores$TFBSScore), ]))
    unalteredPlot <- ggplot(data = dFScores, aes(y = TFBSPredictScore, x = TFBSScore )) + geom_point(alpha = 0.5 ) +
      scale_x_continuous(breaks = c(-15, -10, -5, 0, 5, 10, 15, 20)) +
      scale_y_continuous(breaks = c(-5, 0, 5, 10), limits = c(-5,10))
    alteredPlot <-   ggplot(data = dFScores, aes(y = TFBSPredictScoreBounded, x = TFBSScore )) + geom_point(alpha = 0.5 ) +
      scale_x_continuous(breaks = c(-15, -10, -5, 0, 5, 10, 15, 20)) +
      scale_y_continuous(breaks = c(-5, 0, 5, 10), limits = c(-5,10))
    #doesn't work perfectly, still some values are perhaps lower than necessary.
    #nevertheless, only ~12 values are affected
    
    crossValOutput$totalScoreTable[[crossValSet]]$TFBSScoreFit <- dFScores$TFBSPredictScoreBounded
    #this is growing a list, which is not nice. Perhaps in the final incarnation I could create the list beforehand and fill the items as needed.
    crossValOutput$extraInfoPerCrossval[[crossValSet]][[length(crossValOutput$extraInfoPerCrossval[[crossValSet]]) + 1]]    <- polynomialDensityPlot ; names(crossValOutput$extraInfoPerCrossval[[crossValSet]])[[length(crossValOutput$extraInfoPerCrossval[[crossValSet]])]] <- "smoothenedDistributionPlot"
    crossValOutput$extraInfoPerCrossval[[crossValSet]][[length(crossValOutput$extraInfoPerCrossval[[crossValSet]]) + 1]]    <- unalteredPlot         ; names(crossValOutput$extraInfoPerCrossval[[crossValSet]])[[length(crossValOutput$extraInfoPerCrossval[[crossValSet]])]] <- "unalteredTFBSSmoothenedLog2ScoresPlot"
    crossValOutput$extraInfoPerCrossval[[crossValSet]][[length(crossValOutput$extraInfoPerCrossval[[crossValSet]]) + 1]]    <- alteredPlot           ; names(crossValOutput$extraInfoPerCrossval[[crossValSet]])[[length(crossValOutput$extraInfoPerCrossval[[crossValSet]])]] <- "boundedTFBSSmoothenedLog2ScoresPlot"
    crossValOutput$extraInfoPerCrossval[[crossValSet]][[length(crossValOutput$extraInfoPerCrossval[[crossValSet]]) + 1]]    <- numericVectorBounding ; names(crossValOutput$extraInfoPerCrossval[[crossValSet]])[[length(crossValOutput$extraInfoPerCrossval[[crossValSet]])]] <- "statisticsBoundedScores"
    
  }
  
  print("Finished")
  print("Head of scoretable of crossVal 1:"); print (head(crossValOutput$totalScoreTable[[1]]))
  print(paste("Head of scoretable of crossVal", length(crossValOutput$positiveSetNotTrainedOn), ":")); print (head(crossValOutput$totalScoreTable[[length(crossValOutput$positiveSetNotTrainedOn)]]))
  
  crossValOutput
 
}

#Now: a function that computes viral PPI Scores. Should calculate simultaneously for Union and Intersect.
#Takes the yes and no for found in sets, calculates enrichments, and gives genes with yes the enrichment score, and
#genes with no the non-enrichment score.

computeViralPPIScores <- function(crossValOutput) {
  
  #do everything for the amount of times crossval is performed
  for (crossValSet in seq_len(length(crossValOutput$positiveSetNotTrainedOn))) {
    
    
    
    #define positive set and negative set genes that are used for training in this crossVal
    positiveSetGenes <- crossValOutput$positiveSetTrainedOn[[crossValSet]]$EnsemblID
    negativeSetGenes <- crossValOutput$negativeSetTrainedOn[[crossValSet]]$EnsemblID
    
    
    #load viral PPI scores, both Union and Intersect
    oneWidthBinScoresNotSummarizedTrain <- crossValOutput$totalSetTrainedOn[[crossValSet]] %>% select(EnsemblID, Set, UnionMeasuredPPI) %>%
      group_by(gr = UnionMeasuredPPI, Set)
    
    oneWidthBinScoresNotSummarizedTotal <- crossValOutput$fullTotalSet[[crossValSet]] %>% select(EnsemblID, Set, UnionMeasuredPPI) %>%
      group_by(gr = UnionMeasuredPPI, Set)
    oneWidthBinScoresSummarized <- oneWidthBinScoresNotSummarizedTrain  %>% summarize(Count = n())
    
    
    
    
    oneWidthBinScoresNotSummarizedTrainIntersect <- crossValOutput$totalSetTrainedOn[[crossValSet]] %>% select(EnsemblID, Set, IntersectMeasuredPPI) %>%
      group_by(gr = IntersectMeasuredPPI, Set)
    
    oneWidthBinScoresNotSummarizedTotalIntersect <- crossValOutput$fullTotalSet[[crossValSet]] %>% select(EnsemblID, Set, IntersectMeasuredPPI) %>%
      group_by(gr = IntersectMeasuredPPI, Set)
    oneWidthBinScoresSummarizedIntersect <- oneWidthBinScoresNotSummarizedTrainIntersect  %>% summarize(Count = n())
    
    
    
    #for a grouped data frame of profiles and the counts for positive and negative set genes, adds in 
    #an entry with the same profile, the missing set, and a count of 0
    #not necessary here.
    addMissingSetCounts <- function(summarizedData) {
      
      setVector <- c("PS", "NS")
      for(set in seq_len(length(setVector))) {
        for(i in seq_len(length(unique(summarizedData$gr)))) {
          if(!is.na(unique(summarizedData$gr)[i])) {
            if(summarizedData %>% filter(gr == unique(summarizedData$gr)[i], Set == setVector[set]) %>% `$`(Count) %>% length() == 0 ) {
              summarizedData <- summarizedData %>% ungroup %>% add_row(gr = unique(summarizedData$gr)[i], Set = setVector[set], Count = 0) 
            }
          }
        }
        #for NA
        if(summarizedData %>% filter(is.na(gr) == TRUE, Set == setVector[set]) %>% `$`(Count) %>% length() == 0 ) {
          summarizedData <- summarizedData %>% ungroup %>% add_row(gr = NA, Set = setVector[set], Count = 0) }
      }
      summarizedData <- summarizedData %>% group_by(gr) %>% arrange(gr)
      summarizedData
    }
    
    oneWidthBinScoresSummarized <- addMissingSetCounts(oneWidthBinScoresSummarized)
    oneWidthBinScoresSummarizedIntersect <- addMissingSetCounts(oneWidthBinScoresSummarizedIntersect)
    
    #supply a vector of positive Set and of negative Set genes. Insert the summarized bin data from above.
    #per bin, divide PS counts by PSsize, and divide NS counst by NSsize. Then divide those two by each other and take log2
    #special approach was taken to include NAs. (seq_len(i-1) selects the group which is NA)
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
    
    log2scoresOneWidthBins       <- getLog2Scores(oneWidthBinScoresSummarized, positiveSetGenes, negativeSetGenes)
    oneWidthBinScoresSummarized$log2ScoresViralPPI  <- rep(log2scoresOneWidthBins, each = 2)
    
    log2scoresOneWidthBinsIntersect       <- getLog2Scores(oneWidthBinScoresSummarizedIntersect, positiveSetGenes, negativeSetGenes)
    oneWidthBinScoresSummarizedIntersect$log2ScoresViralPPI  <- rep(log2scoresOneWidthBinsIntersect, each = 2)
    
    head(oneWidthBinScoresSummarized); tail(oneWidthBinScoresSummarized)
    
    
    #now give everything a score (traind on + not trained on)
    #only not trained on is kept from this crossVal run, but scores are assigned to everything to make a rank across crossVals.
    #Note that at present, NA are given a score. These are genes that were not in the Human Tissue Score dataset. 
    #I need to discuss with John whether
    #this is really pertinent or not.
    
    
    
    assignScoresBinsToAllGenes          <- function(summarizedBinScores, notSummarizedBinScores) {
      notSummarizedBinScores$log2ScoresViralPPI = 0
      print(head(notSummarizedBinScores))
      for (group in seq_len(length(unique(notSummarizedBinScores$gr)))) {
        #print(group)
        #print(is.na(unique(notSummarizedBinScores$gr)[group]))
        
        scoreToPlace <- summarizedBinScores %>%
          filter(gr %in% unique(notSummarizedBinScores$gr)[group], Set == "PS") %>% select(log2ScoresViralPPI)  %>% as.data.frame() %>% `[`(1,2, drop = TRUE)
        #print(scoreToPlace)
        #print(head(notSummarizedBinScores[notSummarizedBinScores$gr %in% unique(notSummarizedBinScores$gr)[group],]))
        notSummarizedBinScores[notSummarizedBinScores$gr %in% unique(notSummarizedBinScores$gr)[group],]$log2ScoresViralPPI = scoreToPlace
        #print(notSummarizedBinScores[notSummarizedBinScores$gr %in% unique(notSummarizedBinScores$gr)[group],])
        
      }
      
      notSummarizedBinScores %>% arrange(EnsemblID)
      
    }
    
    oneWidthBinScoresNotSummarizedTotal <- assignScoresBinsToAllGenes(oneWidthBinScoresSummarized, oneWidthBinScoresNotSummarizedTotal)
    oneWidthBinScoresNotSummarizedTotalIntersect <- assignScoresBinsToAllGenes(oneWidthBinScoresSummarizedIntersect, oneWidthBinScoresNotSummarizedTotalIntersect)
    
    
    print(head(oneWidthBinScoresNotSummarizedTotal));                  print(tail(oneWidthBinScoresNotSummarizedTotal))
    
    
    #make a plot for visual inspection
    #function generateBinningPlots creates plots and adds in the densityPerSet column to the summarized data
    
    viralPPIUnionBinScoresResult   <- (generateBinningPlots(oneWidthBinScoresSummarized[!is.na(oneWidthBinScoresSummarized$gr),],
                                                      title = "Human-viral protein-protein interactions measured (Union)",
                                                      xlabel = "Bin",
                                                      ylabel = "Density of genes in each profile per set"))
    
    viralPPIIntersectBinScoresResult <- (generateBinningPlots(oneWidthBinScoresSummarizedIntersect[!is.na(oneWidthBinScoresSummarizedIntersect$gr),],
                                                              title = "Human-viral protein-protein interactions measured (Intersect)",
                                                              xlabel = "Bin",
                                                              ylabel = "Density of genes in each profile per set"))
    
    #Use the customBinScores, keep the oneWidthBinScores for reference. Keep also the summarizedTables.
    
    crossValOutput$totalScoreTable[[crossValSet]]$PPIScoreUnion     <- oneWidthBinScoresNotSummarizedTotal$log2ScoresViralPPI
    crossValOutput$totalScoreTable[[crossValSet]]$PPIScoreIntersect <- oneWidthBinScoresNotSummarizedTotalIntersect$log2ScoresViralPPI
    
    
    crossValOutput$extraInfoPerCrossval[[crossValSet]][[length(crossValOutput$extraInfoPerCrossval[[crossValSet]]) + 1]]    <- viralPPIUnionBinScoresResult[[1]]
    names(crossValOutput$extraInfoPerCrossval[[crossValSet]])[[length(crossValOutput$extraInfoPerCrossval[[crossValSet]])]] <- "scoreCountsViralPPIUnion"
    crossValOutput$extraInfoPerCrossval[[crossValSet]][[length(crossValOutput$extraInfoPerCrossval[[crossValSet]]) + 1]]    <- viralPPIUnionBinScoresResult[[2]]
    names(crossValOutput$extraInfoPerCrossval[[crossValSet]])[[length(crossValOutput$extraInfoPerCrossval[[crossValSet]])]] <- "viralPPIUnionPlot"
    
    crossValOutput$extraInfoPerCrossval[[crossValSet]][[length(crossValOutput$extraInfoPerCrossval[[crossValSet]]) + 1]]    <- viralPPIIntersectBinScoresResult[[1]]
    names(crossValOutput$extraInfoPerCrossval[[crossValSet]])[[length(crossValOutput$extraInfoPerCrossval[[crossValSet]])]] <- "scoreCountsViralPPIIntersect"
    crossValOutput$extraInfoPerCrossval[[crossValSet]][[length(crossValOutput$extraInfoPerCrossval[[crossValSet]]) + 1]]    <- viralPPIIntersectBinScoresResult[[2]]
    names(crossValOutput$extraInfoPerCrossval[[crossValSet]])[[length(crossValOutput$extraInfoPerCrossval[[crossValSet]])]] <- "viralPPIIntersectPlot"
    
    
  }
  crossValOutput
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  #general function that outputs enrichments for 1/0 scenarios such as here (either you have viral ppi as a gene, or you haven't)
  # nrich_nalyse_yesno <- function(positiveSet, psIDCol = "Ensembl_accession", dataSet, dsValueCol, dsIDCol, allEnsemblGenes, aegIDCol) {
  #   
  #   
  #   pctCoverageDataSet          <- sum(allEnsemblGenes[, get(aegIDCol)] %in% dataSet[, get(dsIDCol)])/length(allEnsemblGenes[, get(aegIDCol)])*100
  #   print(paste0("In total, ", pctCoverageDataSet, "% of the genes are covered in this dataset."))
  #   
  #   
  #   presentPSGenes              <- positiveSet[, get(psIDCol)][positiveSet[, get(psIDCol)] %in% dataSet[, get(dsIDCol)]]
  #   presentDfPSGenes            <- positiveSet[positiveSet[, get(psIDCol)] %in% dataSet[, get(dsIDCol)]]
  #   
  #   absentPSGenes               <- positiveSet[, get(psIDCol)][positiveSet[, get(psIDCol)] %ni% dataSet[, get(dsIDCol)]]
  #   absentDfPSGenes             <- positiveSet[positiveSet[, get(psIDCol)] %ni% dataSet[, get(dsIDCol)]]
  #   amountPresentPSGenes        <- nrow(presentDfPSGenes)
  #   totalPSGenes                <- length(positiveSet[, get(psIDCol)])
  #   totalWGGenes                <- length(allEnsemblGenes[, get(aegIDCol)])
  #   pctPresentPSGenes           <- amountPresentPSGenes/totalPSGenes*100
  #   cat(paste0(pctPresentPSGenes, "% of the positive set genes are present in the dataset, i.e. ",
  #              amountPresentPSGenes, " out of ", totalPSGenes, " positive set genes."))
  #   
  #   
  #   presentDataSet             <- allEnsemblGenes[, get(aegIDCol)][allEnsemblGenes[, get(aegIDCol)] %in% dataSet[, get(dsIDCol)]]
  #   presentDfDataSet           <- allEnsemblGenes[allEnsemblGenes[,get(aegIDCol)] %in% presentDataSet, ]
  #   absentDataSet              <- allEnsemblGenes[, get(aegIDCol)][allEnsemblGenes[, get(aegIDCol)] %ni% dataSet[, get(dsIDCol)]]
  #   absentDfDataSet            <- allEnsemblGenes[allEnsemblGenes[,get(aegIDCol)] %ni% presentDataSet, ]
  #   
  #   
  #   allEnsemblGenesDf <- data.table(ENSEMBL_ID = allEnsemblGenes[, get(aegIDCol)],
  #                                   Value = "no")
  #   head(allEnsemblGenesDf)
  #   tail(allEnsemblGenesDf)
  #   allEnsemblGenesDf[ENSEMBL_ID %in% presentDataSet, ][, "Value"] <- "yes"
  #   
  #   #analysis
  #   
  #   bins <- c("no", "yes")
  #   
  #   
  #   
  #   fisherDataFrame <- data.frame(bins = character(),
  #                                 pValue=numeric(),
  #                                 sampleEstimate=numeric(),
  #                                 lowerConfidenceInterval=numeric(),
  #                                 upperConfidenceInterval=numeric())
  #   
  #   
  #   
  #   
  #   PSyes         <- filter(allEnsemblGenesDf, Value == "yes")
  #   PSyes         <- PSyes[ PSyes$ENSEMBL_ID %in% positiveSet[, get(psIDCol)],]
  #   
  #   PSno          <- filter(allEnsemblGenesDf, Value == "no")
  #   PSno          <- PSno[ PSno$ENSEMBL_ID %in% positiveSet[, get(psIDCol)],]
  #   
  #   WGyes         <- filter(allEnsemblGenesDf, Value == "yes")
  #   WGyes         <- WGyes[ WGyes$ENSEMBL_ID %ni% positiveSet[, get(psIDCol)],]
  #   
  #   WGno          <- filter(allEnsemblGenesDf, Value == "no")
  #   WGno          <- WGno[ WGno$ENSEMBL_ID %ni% positiveSet[, get(psIDCol)],]
  #   
  #   
  #   PSGenes       <- dataSet[dataSet[, get(dsIDCol)] %in% presentPSGenes, ]
  #   WGGenes       <- dataSet[dataSet[, get(dsIDCol)] %ni% presentPSGenes, ]
  #   
  #   
  #   listMatrices <- list()
  #   fisherMatrix <- matrix(c(nrow(PSyes),
  #                            totalPSGenes - nrow(PSyes),
  #                            nrow(WGyes),
  #                            totalWGGenes - nrow(WGyes)),
  #                          nrow = 2,
  #                          dimnames = list(Bin = c("yes", "no"),
  #                                          Set = c("Positive set", "Whole genome")))
  #   fisherTestResult <- fisher.test(fisherMatrix)
  #   
  #   pValue = fisherTestResult["p.value"]
  #   lowerConfInt <- unlist(fisherTestResult["conf.int"])[1]
  #   upperConfInt <- unlist(fisherTestResult["conf.int"])[2]
  #   sampleEst <- fisherTestResult["estimate"]
  #   
  #   listMatrices <- append(listMatrices,fisherMatrix)
  #   addVector <- c("viral PPI",
  #                  pValue,
  #                  sampleEst,
  #                  lowerConfInt,
  #                  upperConfInt)
  #   names(addVector) <- c("bins", "pValue", "sampleEstimate",
  #                         "lowerConfidenceInterval", "upperConfidenceInterval")
  #   #add to the Df the results of the fisherTest
  #   fisherDataFrame <- rbind(fisherDataFrame, addVector, make.row.names = FALSE)
  #   fisherDataFrame
  # }  
  # 
  # for (crossValSet in seq_len(length(crossValOutput$positiveSetNotTrainedOn))) {
  #   
  #   print(paste("Working on crossVal", crossValSet, "out of", length(crossValOutput$positiveSetNotTrainedOn), "."))
  # 
  #   
  #   #what follows below are definitions of arguments for the function above. It outputs enrichment of viral ppi interactions.
  # positiveSet          <- crossValOutput$positiveSetTrainedOn[[crossValSet]] %>% as.data.table()
  # positiveSetIDCol     <- "EnsemblID"
  # dataSetUnion         <- crossValOutput$totalSetTrainedOn[[crossValSet]] %>% filter(UnionMeasuredPPI == "yes") %>% as.data.table()
  # dsValueColUnion      <- "UnionMeasuredPPI"
  # dataSetIntersect     <- crossValOutput$totalSetTrainedOn[[crossValSet]] %>% filter(IntersectMeasuredPPI == "yes") %>% as.data.table()
  # dsValueColIntersect  <- "IntersectMeasuredPPI"
  # dsIDCol              <- "EnsemblID"
  # allEnsemblGenes      <- crossValOutput$totalSetTrainedOn[[crossValSet]] %>% as.data.table()
  # allEnsemblGenesIDCol <- "EnsemblID"
  # 
  # unionScore           <- nrich_nalyse_yesno(positiveSet, positiveSetIDCol, dataSetUnion, dsValueColUnion, dsIDCol, allEnsemblGenes, allEnsemblGenesIDCol)
  # log2UnionScore       <- log2(unionScore$sampleEstimate)
  # 
  # intersectScore       <- nrich_nalyse_yesno(positiveSet, positiveSetIDCol, dataSetIntersect, dsValueColIntersect, dsIDCol, allEnsemblGenes, allEnsemblGenesIDCol)
  # log2IntersectScore   <- log2(intersectScore$sampleEstimate)
  # 
  # #assign the scores
  # log2UnionScoresForTable     <- ifelse(crossValOutput$fullTotalSet[[crossValSet]]$UnionMeasuredPPI     == "yes", log2UnionScore    , 0)
  # log2IntersectScoresForTable <- ifelse(crossValOutput$fullTotalSet[[crossValSet]]$IntersectMeasuredPPI == "yes", log2IntersectScore, 0)
  # 
  # crossValOutput$totalScoreTable[[crossValSet]]$PPIScoreUnion     <- log2UnionScoresForTable
  # crossValOutput$totalScoreTable[[crossValSet]]$PPIScoreIntersect <- log2IntersectScoresForTable
  # 
  # #additional details. Include the enrichment details in the final output.
  # unionScore$bins     <- "Union yes"
  # intersectScore$bins <- "Intersect yes"
  # dFData              <- rbind(unionScore, intersectScore)
  # dFData              <- dFData %>% as_tibble()
  # 
  # #the below is growing a list, which is wrong. Once the code is finished, I should initialise the list with enough slots, and just fill them.
  # crossValOutput$extraInfoPerCrossval[[crossValSet]][[length(crossValOutput$extraInfoPerCrossval[[crossValSet]]) + 1]]    <- dFData ; names(crossValOutput$extraInfoPerCrossval[[crossValSet]])[[length(crossValOutput$extraInfoPerCrossval[[crossValSet]])]] <- "viralPPIEnrichmentDetails"
  # 
  # }
  # 
  # crossValOutput
}

#function to calculate the enrichment per Z-score bin of positive set genes.
#orgPositiveSetFile is the Dataframe directly obtained by loading in Positivelist_2.csv, which contains features
#of known MHC genes (the positive set). It is used to exclude positive set genes that were added based on the Neefjes
#paper to be used in determining the score for positive set genes in the Neefjes paper (would be circular). 

#Note that no Zscores are between [-0.9, 0.9], because they are already only the sig. different Zscores (in the case of the rescreened). 
binsZscore         <- c(0, 2, 10)
binsOriginalZscore <- c(0, 2, 4, 18)

#function uses positiveSet where Neefjes not in there (i.e. it removes those positive set genes which were added to the
#positive set based on the Neefjes MHC II surface expr. perturbation screen)

computeZscoreScores <- function(binsRescreen, binsOriginal, crossValOutput, orgPositiveSetFile) {
  
  
  
  #do everything for the amount of times crossval is performed
  for (crossValSet in seq_len(length(crossValOutput$positiveSetNotTrainedOn))) {
    
    
    
    #define positive set and negative set genes that are used for training in this crossVal
    positiveSetGenes <- crossValOutput$positiveSetTrainedOn[[crossValSet]]$EnsemblID[crossValOutput$positiveSetTrainedOn[[crossValSet]]$EnsemblID %in% orgPositiveSetFile[orgPositiveSetFile$Neefjes_accession == "",]$Ensembl_accession]
    #correct positive set for the fact that these genes should not have been added based on the Neefjes paper!
    negativeSetGenes <- crossValOutput$negativeSetTrainedOn[[crossValSet]]$EnsemblID
    
    
    #compute binned scores. Always also include binning intervals of 1 to compare
    
    
    #custom bins ZscoreHighscore
    customBinScoresNotSummarizedZscoreHighscoreTrain <- crossValOutput$totalSetTrainedOn[[crossValSet]] %>% select(EnsemblID, Set, ZscoreHighScore) %>%
      filter(EnsemblID %ni% orgPositiveSetFile[orgPositiveSetFile$Neefjes_accession != "",]) %>%
      group_by(gr = cut(ZscoreHighScore, breaks = binsRescreen), Set)
    customBinScoresNotSummarizedZscoreHighscoreTotal <- crossValOutput$fullTotalSet[[crossValSet]]      %>% select(EnsemblID, Set, ZscoreHighScore) %>%
      filter(EnsemblID %ni% orgPositiveSetFile[orgPositiveSetFile$Neefjes_accession != "",]) %>%
      group_by(gr = cut(ZscoreHighScore, breaks = binsRescreen), Set)
    customBinScoresSummarizedZscoreHighscore         <- customBinScoresNotSummarizedZscoreHighscoreTrain %>% summarize(Count = n())
    
    
    #bins of 1 width ZscoreHighscore
    oneWidthBinScoresNotSummarizedZscoreHighscoreTrain <- crossValOutput$totalSetTrainedOn[[crossValSet]] %>% select(EnsemblID, Set, ZscoreHighScore) %>%
      filter(EnsemblID %ni% orgPositiveSetFile[orgPositiveSetFile$Neefjes_accession != "",]) %>%
      group_by(gr = cut(ZscoreHighScore, breaks = seq(floor(min(ZscoreHighScore, na.rm = TRUE)), ceiling(max(ZscoreHighScore, na.rm = TRUE)))), Set)
      
    oneWidthBinScoresNotSummarizedZscoreHighscoreTotal <- crossValOutput$fullTotalSet[[crossValSet]] %>% select(EnsemblID, Set, ZscoreHighScore) %>%
      filter(EnsemblID %ni% orgPositiveSetFile[orgPositiveSetFile$Neefjes_accession != "",]) %>%
      group_by(gr = cut(ZscoreHighScore, breaks = seq(floor(min(ZscoreHighScore, na.rm = TRUE)), ceiling(max(ZscoreHighScore, na.rm = TRUE)))), Set)
    oneWidthBinScoresSummarizedZscoreHighscore <- oneWidthBinScoresNotSummarizedZscoreHighscoreTrain  %>% summarize(Count = n())
    
    #custom bins ZscoreRescreenOnly
    customBinScoresNotSummarizedZscoreRescreenOnlyTrain <- crossValOutput$totalSetTrainedOn[[crossValSet]] %>% select(EnsemblID, Set, RescreenOnlyZscores) %>%
      filter(EnsemblID %ni% orgPositiveSetFile[orgPositiveSetFile$Neefjes_accession != "",]) %>%
      group_by(gr = cut(RescreenOnlyZscores, breaks = binsRescreen), Set)
    customBinScoresNotSummarizedZscoreRescreenOnlyTotal <- crossValOutput$fullTotalSet[[crossValSet]]      %>% select(EnsemblID, Set, RescreenOnlyZscores) %>%
      filter(EnsemblID %ni% orgPositiveSetFile[orgPositiveSetFile$Neefjes_accession != "",]) %>%
      group_by(gr = cut(RescreenOnlyZscores, breaks = binsRescreen), Set)
    customBinScoresSummarizedZscoreRescreenOnly         <- customBinScoresNotSummarizedZscoreRescreenOnlyTrain %>% summarize(Count = n())
    
    #bins of 1 width ZscoreRescreenOnly
    oneWidthBinScoresNotSummarizedZscoreRescreenOnlyTrain <- crossValOutput$totalSetTrainedOn[[crossValSet]] %>% select(EnsemblID, Set, RescreenOnlyZscores) %>%
      filter(EnsemblID %ni% orgPositiveSetFile[orgPositiveSetFile$Neefjes_accession != "",]) %>%
      group_by(gr = cut(RescreenOnlyZscores, breaks = seq(floor(min(RescreenOnlyZscores, na.rm = TRUE)), ceiling(max(RescreenOnlyZscores, na.rm = TRUE)))), Set)
    
    oneWidthBinScoresNotSummarizedZscoreRescreenOnlyTotal <- crossValOutput$fullTotalSet[[crossValSet]] %>% select(EnsemblID, Set, RescreenOnlyZscores) %>%
      filter(EnsemblID %ni% orgPositiveSetFile[orgPositiveSetFile$Neefjes_accession != "",]) %>%
      group_by(gr = cut(RescreenOnlyZscores, breaks = seq(floor(min(RescreenOnlyZscores, na.rm = TRUE)), ceiling(max(RescreenOnlyZscores, na.rm = TRUE)))), Set)
    oneWidthBinScoresSummarizedZscoreRescreenOnly <- oneWidthBinScoresNotSummarizedZscoreRescreenOnlyTrain  %>% summarize(Count = n())
    
    
    #custom bins Original Zscore
    customBinScoresNotSummarizedZscoreOriginalTrain <- crossValOutput$totalSetTrainedOn[[crossValSet]] %>% select(EnsemblID, Set, OriginalZscoreHighestMedian) %>%
      filter(EnsemblID %ni% orgPositiveSetFile[orgPositiveSetFile$Neefjes_accession != "",]) %>%
      group_by(gr = cut(OriginalZscoreHighestMedian, breaks = binsOriginal), Set)
    customBinScoresNotSummarizedZscoreOriginalTotal <- crossValOutput$fullTotalSet[[crossValSet]]      %>% select(EnsemblID, Set, OriginalZscoreHighestMedian) %>%
      filter(EnsemblID %ni% orgPositiveSetFile[orgPositiveSetFile$Neefjes_accession != "",]) %>%
      group_by(gr = cut(OriginalZscoreHighestMedian, breaks = binsOriginal), Set)
    customBinScoresSummarizedZscoreOriginal         <- customBinScoresNotSummarizedZscoreOriginalTrain %>% summarize(Count = n())
    
    #bins of 1 width OriginalZScore
    oneWidthBinScoresNotSummarizedZscoreOriginalTrain <- crossValOutput$totalSetTrainedOn[[crossValSet]] %>% select(EnsemblID, Set, OriginalZscoreHighestMedian) %>%
      filter(EnsemblID %ni% orgPositiveSetFile[orgPositiveSetFile$Neefjes_accession != "",]) %>%
      group_by(gr = cut(OriginalZscoreHighestMedian, breaks = seq(floor(min(OriginalZscoreHighestMedian, na.rm = TRUE)), ceiling(max(OriginalZscoreHighestMedian, na.rm = TRUE)))), Set)
    
    oneWidthBinScoresNotSummarizedZscoreOriginalTotal <- crossValOutput$fullTotalSet[[crossValSet]] %>% select(EnsemblID, Set, OriginalZscoreHighestMedian) %>%
      filter(EnsemblID %ni% orgPositiveSetFile[orgPositiveSetFile$Neefjes_accession != "",]) %>%
      group_by(gr = cut(OriginalZscoreHighestMedian, breaks = seq(floor(min(OriginalZscoreHighestMedian, na.rm = TRUE)), ceiling(max(OriginalZscoreHighestMedian, na.rm = TRUE)))), Set)
    oneWidthBinScoresSummarizedZscoreOriginal <- oneWidthBinScoresNotSummarizedZscoreOriginalTrain  %>% summarize(Count = n())
    
    
    
    
    
    #problem: not every bin contains both PS and NS scores. For score calculations, that is necessary. Need to add these in.
    #hence, function below:
    #for a grouped data frame of bins and the counts for positive and negative set genes, adds in 
    #an entry with the same bin, the missing set, and a count of 0
    addMissingSetCounts <- function(summarizedData) {
      
      setVector <- c("PS", "NS")
      for(set in seq_len(length(setVector))) {
        for(i in seq_len(length(unique(summarizedData$gr)))) {
          if(!is.na(unique(summarizedData$gr)[i])) {
            if(summarizedData %>% filter(gr == unique(summarizedData$gr)[i], Set == setVector[set]) %>% `$`(Count) %>% length() == 0 ) {
              summarizedData <- summarizedData %>% ungroup %>% add_row(gr = unique(summarizedData$gr)[i], Set = setVector[set], Count = 0) 
            }
          }
        }
        #for NA
        if(summarizedData %>% filter(is.na(gr) == TRUE, Set == setVector[set]) %>% `$`(Count) %>% length() == 0 ) {
          summarizedData <- summarizedData %>% ungroup %>% add_row(gr = NA, Set = setVector[set], Count = 0) }
      }
      summarizedData <- summarizedData %>% group_by(gr) %>% arrange(gr)
      summarizedData
    }
    
    customBinScoresSummarizedZscoreHighscore      <- addMissingSetCounts(customBinScoresSummarizedZscoreHighscore)
    oneWidthBinScoresSummarizedZscoreHighscore    <- addMissingSetCounts(oneWidthBinScoresSummarizedZscoreHighscore)
    
    customBinScoresSummarizedZscoreRescreenOnly   <- addMissingSetCounts(customBinScoresSummarizedZscoreRescreenOnly)
    oneWidthBinScoresSummarizedZscoreRescreenOnly <- addMissingSetCounts(oneWidthBinScoresSummarizedZscoreRescreenOnly)
    
    customBinScoresSummarizedZscoreOriginal       <- addMissingSetCounts(customBinScoresSummarizedZscoreOriginal)
    oneWidthBinScoresSummarizedZscoreOriginal     <- addMissingSetCounts(oneWidthBinScoresSummarizedZscoreOriginal)
    
    #supply a vector of positive Set and of negative Set genes. Insert the summarized bin data from above.
    #per bin, divide PS counts by PSsize, and divide NS counst by NSsize. Then divide those two by each other and take log2
    #special approach was taken to include NAs. (seq_len(i-1) selects the group which is NA)
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
    
    
    #Zscores from the Rescreen, selected for actual expression in immune cells
    log2scoresCustomBinsZscoreHighscore                             <- getLog2Scores(customBinScoresSummarizedZscoreHighscore,
                                                                                     positiveSetGenes, negativeSetGenes)
    customBinScoresSummarizedZscoreHighscore$log2ScoresZscoreBins   <- rep(log2scoresCustomBinsZscoreHighscore, each = 2)
    log2scoresOneWidthBinsZscoreHighscore                           <- getLog2Scores(oneWidthBinScoresSummarizedZscoreHighscore,
                                                                                     positiveSetGenes, negativeSetGenes)
    oneWidthBinScoresSummarizedZscoreHighscore$log2ScoresZscoreBins <- rep(log2scoresOneWidthBinsZscoreHighscore , each = 2)
    
    #rescreen Zscores, not selected for actual expression in immune cells via microarray verification
    log2scoresCustomBinsZscoreRescreenOnly         <- getLog2Scores(customBinScoresSummarizedZscoreRescreenOnly , positiveSetGenes, negativeSetGenes)
    customBinScoresSummarizedZscoreRescreenOnly$log2ScoresZscoreBins   <- rep(log2scoresCustomBinsZscoreRescreenOnly, each = 2)
    log2scoresOneWidthBinsZscoreRescreenOnly       <- getLog2Scores(oneWidthBinScoresSummarizedZscoreRescreenOnly , positiveSetGenes, negativeSetGenes)
    oneWidthBinScoresSummarizedZscoreRescreenOnly$log2ScoresZscoreBins <- rep(log2scoresOneWidthBinsZscoreRescreenOnly, each = 2)
    
    #original Zscores
    log2scoresCustomBinsZscoreOriginal         <- getLog2Scores(customBinScoresSummarizedZscoreOriginal , positiveSetGenes, negativeSetGenes)
    customBinScoresSummarizedZscoreOriginal$log2ScoresZscoreBins   <- rep(log2scoresCustomBinsZscoreOriginal, each = 2)
    log2scoresOneWidthBinsZscoreOriginal       <- getLog2Scores(oneWidthBinScoresSummarizedZscoreOriginal , positiveSetGenes, negativeSetGenes)
    oneWidthBinScoresSummarizedZscoreOriginal$log2ScoresZscoreBins <- rep(log2scoresOneWidthBinsZscoreOriginal, each = 2)
    
    #head(customBinScoresSummarized); tail(customBinScoresSummarized)
    #head(oneWidthBinScoresSummarized); tail(oneWidthBinScoresSummarized)
    
    
    #now give everything a score (traind on + not trained on)
    #only not trained on is kept from this crossVal run, but scores are assigned to everything to make a rank across crossVals.
    #Note that at present, NA are given a score. These are genes that were not in the Z score dataset. SInce only 3 of the positive
    #set were not in the dataset, and many from the negative set, these are given a negative score. I need to discuss with John whether
    #this is really pertinent or not.
    
    
    
    assignScoresBinsToAllGenes          <- function(summarizedBinScores, notSummarizedBinScores) {
      notSummarizedBinScores$log2ScoresZscoreBins = 0
      print(head(notSummarizedBinScores))
      for (group in seq_len(length(unique(notSummarizedBinScores$gr)))) {
        #print(group)
        #print(is.na(unique(notSummarizedBinScores$gr)[group]))
        
        scoreToPlace <- summarizedBinScores %>%
          filter(gr %in% unique(notSummarizedBinScores$gr)[group], Set == "PS") %>% select(log2ScoresZscoreBins)  %>% as.data.frame() %>% `[`(1,2, drop = TRUE)
        #print(scoreToPlace)
        #print(head(notSummarizedBinScores[notSummarizedBinScores$gr %in% unique(notSummarizedBinScores$gr)[group],]))
        notSummarizedBinScores[notSummarizedBinScores$gr %in% unique(notSummarizedBinScores$gr)[group],]$log2ScoresZscoreBins = scoreToPlace
        #print(notSummarizedBinScores[notSummarizedBinScores$gr %in% unique(notSummarizedBinScores$gr)[group],])
        
      }
      
      notSummarizedBinScores %>% arrange(EnsemblID)
      
    }
    
    customBinScoresNotSummarizedTotalZscoreHighscore      <- assignScoresBinsToAllGenes(customBinScoresSummarizedZscoreHighscore  , customBinScoresNotSummarizedZscoreHighscoreTotal  )
    oneWidthBinScoresNotSummarizedTotalZscoreHighscore    <- assignScoresBinsToAllGenes(oneWidthBinScoresSummarizedZscoreHighscore, oneWidthBinScoresNotSummarizedZscoreHighscoreTotal)
    
    customBinScoresNotSummarizedTotalZscoreRescreenOnly   <- assignScoresBinsToAllGenes(customBinScoresSummarizedZscoreRescreenOnly  , customBinScoresNotSummarizedZscoreRescreenOnlyTotal  )
    oneWidthBinScoresNotSummarizedTotalZscoreRescreenOnly <- assignScoresBinsToAllGenes(oneWidthBinScoresSummarizedZscoreRescreenOnly, oneWidthBinScoresNotSummarizedZscoreRescreenOnlyTotal)
    
    customBinScoresNotSummarizedZscoreOriginalTotal       <- assignScoresBinsToAllGenes(customBinScoresSummarizedZscoreOriginal  , customBinScoresNotSummarizedZscoreOriginalTotal  )
    oneWidthBinScoresNotSummarizedZscoreOriginalTotal     <- assignScoresBinsToAllGenes(oneWidthBinScoresSummarizedZscoreOriginal, oneWidthBinScoresNotSummarizedZscoreOriginalTotal)
    
    
    
    #head(oneWidthBinScoresNotSummarizedTotal);                  tail(oneWidthBinScoresNotSummarizedTotal)
   # head(as.data.frame(customBinScoresNotSummarizedTotal), 50); tail(as.data.frame(customBinScoresNotSummarizedTotal), 50)
    
    
    #make a plot for visual inspection
    #function generateBinningPlots creates plots and adds in the densityPerSet column to the summarized data
    
    customBinScoresZscoreHighResult       <- (generateBinningPlots(customBinScoresSummarizedZscoreHighscore,
                                                               title = "Custom bins Z score high confidence set",
                                                               xlabel = "Binned Z scores Neefjes screen"))
    
    oneWidthBinScoresZscoreHighResult     <- (generateBinningPlots(oneWidthBinScoresSummarizedZscoreHighscore,
                                                     title = "One width bins Z score high confidence set",
                                                     xlabel = "Binned Z scores Neefjes screen"))
    
    customBinScoresZscoreRescreenResult   <- (generateBinningPlots(customBinScoresSummarizedZscoreRescreenOnly,
                                                               title = "Custom bins Z score lower confidence set",
                                                               xlabel = "Binned Z scores Neefjes screen"))
    
    oneWidthBinScoresZscoreRescreenResult <- (generateBinningPlots(oneWidthBinScoresSummarizedZscoreRescreenOnly,
                                                               title = "One width bins Z score lower confidence set",
                                                               xlabel = "Binned Z scores Neefjes screen"))
    
    customBinScoresZscoreOriginalResult   <- (generateBinningPlots(customBinScoresSummarizedZscoreOriginal,
                                                                   title = "Custom bins Z score original screen",
                                                                   xlabel = "Binned Z scores Neefjes screen"))
    
    oneWidthBinScoresZscoreOriginalResult <- (generateBinningPlots(oneWidthBinScoresSummarizedZscoreOriginal,
                                                                   title = "One width bins Z score original screen",
                                                                   xlabel = "Binned Z scores Neefjes screen"))
    
    

    
    #Use the customBinScores, keep the oneWidthBinScores for reference. Keep also the summarizedTables.
    
    crossValOutput$totalScoreTable[[crossValSet]]$ZscoreHighConf       <- customBinScoresNotSummarizedTotalZscoreHighscore$log2ScoresZscoreBins
    crossValOutput$totalScoreTable[[crossValSet]]$ZscoreLowerConf      <- customBinScoresNotSummarizedTotalZscoreRescreenOnly$log2ScoresZscoreBins
    crossValOutput$totalScoreTable[[crossValSet]]$ZscoreOriginalScreen <- customBinScoresNotSummarizedZscoreOriginalTotal$log2ScoresZscoreBins
    
    #high confidence, rescreen + expressed in immune cells
    crossValOutput$extraInfoPerCrossval[[crossValSet]][[length(crossValOutput$extraInfoPerCrossval[[crossValSet]]) + 1]]    <- customBinScoresZscoreHighResult[[1]] ; names(crossValOutput$extraInfoPerCrossval[[crossValSet]])[[length(crossValOutput$extraInfoPerCrossval[[crossValSet]])]] <- "customBinScoresCountsPerBinZscoreHighscore"
    crossValOutput$extraInfoPerCrossval[[crossValSet]][[length(crossValOutput$extraInfoPerCrossval[[crossValSet]]) + 1]]    <- customBinScoresZscoreHighResult[[2]] ; names(crossValOutput$extraInfoPerCrossval[[crossValSet]])[[length(crossValOutput$extraInfoPerCrossval[[crossValSet]])]] <- "customBinSZscoreHighConfPlot"
    
    crossValOutput$extraInfoPerCrossval[[crossValSet]][[length(crossValOutput$extraInfoPerCrossval[[crossValSet]]) + 1]]    <- oneWidthBinScoresZscoreHighResult[[1]] ; names(crossValOutput$extraInfoPerCrossval[[crossValSet]])[[length(crossValOutput$extraInfoPerCrossval[[crossValSet]])]] <- "oneWidthBinScoresCountsPerBinZscoreHighscore"
    crossValOutput$extraInfoPerCrossval[[crossValSet]][[length(crossValOutput$extraInfoPerCrossval[[crossValSet]]) + 1]]    <- oneWidthBinScoresZscoreHighResult[[2]] ; names(crossValOutput$extraInfoPerCrossval[[crossValSet]])[[length(crossValOutput$extraInfoPerCrossval[[crossValSet]])]] <- "oneWidthBinSZscoreHighConfPlot"
    
    crossValOutput$extraInfoPerCrossval[[crossValSet]][[length(crossValOutput$extraInfoPerCrossval[[crossValSet]]) + 1]]    <- oneWidthBinScoresNotSummarizedTotalZscoreHighscore ; names(crossValOutput$extraInfoPerCrossval[[crossValSet]])[[length(crossValOutput$extraInfoPerCrossval[[crossValSet]])]] <- "oneWidthBinScoresNotSummarizedZscoreHighscore"
    
    #Rescreen only
    crossValOutput$extraInfoPerCrossval[[crossValSet]][[length(crossValOutput$extraInfoPerCrossval[[crossValSet]]) + 1]]    <- customBinScoresZscoreRescreenResult[[1]] ; names(crossValOutput$extraInfoPerCrossval[[crossValSet]])[[length(crossValOutput$extraInfoPerCrossval[[crossValSet]])]] <- "customBinScoresCountsPerBinZscoreRescreenOnly"
    crossValOutput$extraInfoPerCrossval[[crossValSet]][[length(crossValOutput$extraInfoPerCrossval[[crossValSet]]) + 1]]    <- customBinScoresZscoreRescreenResult[[2]] ; names(crossValOutput$extraInfoPerCrossval[[crossValSet]])[[length(crossValOutput$extraInfoPerCrossval[[crossValSet]])]] <- "customBinSZscoreLowerConfPlot"
    
    crossValOutput$extraInfoPerCrossval[[crossValSet]][[length(crossValOutput$extraInfoPerCrossval[[crossValSet]]) + 1]]    <- oneWidthBinScoresZscoreRescreenResult[[1]] ; names(crossValOutput$extraInfoPerCrossval[[crossValSet]])[[length(crossValOutput$extraInfoPerCrossval[[crossValSet]])]] <- "oneWidthBinScoresCountsPerBinZscoreRescreenOnly"
    crossValOutput$extraInfoPerCrossval[[crossValSet]][[length(crossValOutput$extraInfoPerCrossval[[crossValSet]]) + 1]]    <- oneWidthBinScoresZscoreRescreenResult[[2]] ; names(crossValOutput$extraInfoPerCrossval[[crossValSet]])[[length(crossValOutput$extraInfoPerCrossval[[crossValSet]])]] <- "oneWidthBinSZscoreLowerConfPlot"
    
    crossValOutput$extraInfoPerCrossval[[crossValSet]][[length(crossValOutput$extraInfoPerCrossval[[crossValSet]]) + 1]]    <- oneWidthBinScoresNotSummarizedTotalZscoreRescreenOnly ; names(crossValOutput$extraInfoPerCrossval[[crossValSet]])[[length(crossValOutput$extraInfoPerCrossval[[crossValSet]])]] <- "oneWidthBinScoresNotSummarizedZscoreRescreenOnly"
    
    #original Z scores
    crossValOutput$extraInfoPerCrossval[[crossValSet]][[length(crossValOutput$extraInfoPerCrossval[[crossValSet]]) + 1]]    <- customBinScoresZscoreOriginalResult[[1]] ; names(crossValOutput$extraInfoPerCrossval[[crossValSet]])[[length(crossValOutput$extraInfoPerCrossval[[crossValSet]])]] <- "customBinScoresCountsPerBinOriginalZscores"
    crossValOutput$extraInfoPerCrossval[[crossValSet]][[length(crossValOutput$extraInfoPerCrossval[[crossValSet]]) + 1]]    <- customBinScoresZscoreOriginalResult[[2]] ; names(crossValOutput$extraInfoPerCrossval[[crossValSet]])[[length(crossValOutput$extraInfoPerCrossval[[crossValSet]])]] <- "customBinSOriginalZscoresPlot"
    
    crossValOutput$extraInfoPerCrossval[[crossValSet]][[length(crossValOutput$extraInfoPerCrossval[[crossValSet]]) + 1]]    <- oneWidthBinScoresZscoreOriginalResult[[1]] ; names(crossValOutput$extraInfoPerCrossval[[crossValSet]])[[length(crossValOutput$extraInfoPerCrossval[[crossValSet]])]] <- "oneWidthBinScoresCountsPerBinOriginalZscores"
    crossValOutput$extraInfoPerCrossval[[crossValSet]][[length(crossValOutput$extraInfoPerCrossval[[crossValSet]]) + 1]]    <- oneWidthBinScoresZscoreOriginalResult[[2]] ; names(crossValOutput$extraInfoPerCrossval[[crossValSet]])[[length(crossValOutput$extraInfoPerCrossval[[crossValSet]])]] <- "oneWidthBinsOriginalZscoresPlot"
    
    crossValOutput$extraInfoPerCrossval[[crossValSet]][[length(crossValOutput$extraInfoPerCrossval[[crossValSet]]) + 1]]    <- oneWidthBinScoresNotSummarizedZscoreOriginalTotal ; names(crossValOutput$extraInfoPerCrossval[[crossValSet]])[[length(crossValOutput$extraInfoPerCrossval[[crossValSet]])]] <- "oneWidthBinScoresNotSummarizedOriginalZscores"
    
    
  }
  crossValOutput
  
  
  
  
  
  
  
}




#Immune tissue score
#note that in the original analysis, I disregarded estRankDifferences below 0, since I am only interested in genes sig. more expressed in
#immune tissues. Also, these scores have not been multiple-testing corrected (for example with BH methodology), since
#the enrichment scores will still yield a signal even if the differences are non-significant. 



computeImmuneTissueHistScoreBins <- function(bins, crossValOutput) {
  
  
  #do everything for the amount of times crossval is performed
  for (crossValSet in seq_len(length(crossValOutput$positiveSetNotTrainedOn))) {
    
    
    
    #define positive set and negative set genes that are used for training in this crossVal
    positiveSetGenes <- crossValOutput$positiveSetTrainedOn[[crossValSet]]$EnsemblID
    negativeSetGenes <- crossValOutput$negativeSetTrainedOn[[crossValSet]]$EnsemblID
    
    #compute binned scores. Always also include binning intervals of 1 to compare
    
    #custom bins
    customBinScoresNotSummarizedTrain <- crossValOutput$totalSetTrainedOn[[crossValSet]] %>% select(EnsemblID, Set, MannWhitneyUEstRankDiff) %>%
      group_by(gr = cut(MannWhitneyUEstRankDiff, breaks = bins), Set)
    customBinScoresNotSummarizedTotal <- crossValOutput$fullTotalSet[[crossValSet]]      %>% select(EnsemblID, Set, MannWhitneyUEstRankDiff) %>%
      group_by(gr = cut(MannWhitneyUEstRankDiff, breaks = bins), Set)
    customBinScoresSummarized         <- customBinScoresNotSummarizedTrain %>% summarize(Count = n())
    
    
    #bins of 1 width
    oneWidthBinScoresNotSummarizedTrain <- crossValOutput$totalSetTrainedOn[[crossValSet]] %>% select(EnsemblID, Set, MannWhitneyUEstRankDiff) %>%
      group_by(gr = cut(MannWhitneyUEstRankDiff, breaks = seq(floor(min(MannWhitneyUEstRankDiff, na.rm = TRUE)), ceiling(max(MannWhitneyUEstRankDiff, na.rm = TRUE)))), Set)
    
    oneWidthBinScoresNotSummarizedTotal <- crossValOutput$fullTotalSet[[crossValSet]] %>% select(EnsemblID, Set, MannWhitneyUEstRankDiff) %>%
      group_by(gr = cut(MannWhitneyUEstRankDiff, breaks = seq(floor(min(MannWhitneyUEstRankDiff, na.rm = TRUE)), ceiling(max(MannWhitneyUEstRankDiff, na.rm = TRUE)))), Set)
    oneWidthBinScoresSummarized <- oneWidthBinScoresNotSummarizedTrain  %>% summarize(Count = n())
    
    #prepare data for plots
    plottingMannWhitneyUData <- crossValOutput$totalSetTrainedOn[[crossValSet]] %>% select(EnsemblID, Set, MannWhitneyUEstRankDiff)
    
    #problem: not every bin contains both PS and NS scores. For score calculations, that is necessary. Need to add these in.
    #hence, function below:
    #for a grouped data frame of bins and the counts for positive and negative set genes, adds in 
    #an entry with the same bin, the missing set, and a count of 0
    addMissingSetCounts <- function(summarizedData) {
      
      setVector <- c("PS", "NS")
      for(set in seq_len(length(setVector))) {
        for(i in seq_len(length(unique(summarizedData$gr)))) {
          if(!is.na(unique(summarizedData$gr)[i])) {
            if(summarizedData %>% filter(gr == unique(summarizedData$gr)[i], Set == setVector[set]) %>% `$`(Count) %>% length() == 0 ) {
              summarizedData <- summarizedData %>% ungroup %>% add_row(gr = unique(summarizedData$gr)[i], Set = setVector[set], Count = 0) 
            }
          }
        }
        #for NA
        if(summarizedData %>% filter(is.na(gr) == TRUE, Set == setVector[set]) %>% `$`(Count) %>% length() == 0 ) {
          summarizedData <- summarizedData %>% ungroup %>% add_row(gr = NA, Set = setVector[set], Count = 0) }
      }
      summarizedData <- summarizedData %>% group_by(gr) %>% arrange(gr)
      summarizedData
    }
    
    customBinScoresSummarized   <- addMissingSetCounts(customBinScoresSummarized)
    oneWidthBinScoresSummarized <- addMissingSetCounts(oneWidthBinScoresSummarized)
    
    #supply a vector of positive Set and of negative Set genes. Insert the summarized bin data from above.
    #per bin, divide PS counts by PSsize, and divide NS counst by NSsize. Then divide those two by each other and take log2
    #special approach was taken to include NAs. (seq_len(i-1) selects the group which is NA)
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
    
    
    log2scoresCustomBins         <- getLog2Scores(customBinScoresSummarized, positiveSetGenes, negativeSetGenes)
    customBinScoresSummarized$log2ScoresHumanTissueHistoScoresBins   <- rep(log2scoresCustomBins, each = 2)
    log2scoresOneWidthBins       <- getLog2Scores(oneWidthBinScoresSummarized, positiveSetGenes, negativeSetGenes)
    oneWidthBinScoresSummarized$log2ScoresHumanTissueHistoScoresBins  <- rep(log2scoresOneWidthBins, each = 2)
    
    print(head(customBinScoresSummarized)); print(tail(customBinScoresSummarized))
    head(oneWidthBinScoresSummarized); tail(oneWidthBinScoresSummarized)
   
    
    #now give everything a score (traind on + not trained on)
    #only not trained on is kept from this crossVal run, but scores are assigned to everything to make a rank across crossVals.
    #Note that at present, NA are given a score. These are genes that were not in the Human Tissue Score dataset. 
    #I need to discuss with John whether
    #this is really pertinent or not.
    
    
    
    assignScoresBinsToAllGenes          <- function(summarizedBinScores, notSummarizedBinScores) {
      notSummarizedBinScores$log2ScoresHumanTissueHistoScoresBins = 0
      print(head(notSummarizedBinScores))
      for (group in seq_len(length(unique(notSummarizedBinScores$gr)))) {
        #print(group)
        #print(is.na(unique(notSummarizedBinScores$gr)[group]))
        
        scoreToPlace <- summarizedBinScores %>%
          filter(gr %in% unique(notSummarizedBinScores$gr)[group], Set == "PS") %>% select(log2ScoresHumanTissueHistoScoresBins)  %>% as.data.frame() %>% `[`(1,2, drop = TRUE)
        #print(scoreToPlace)
        #print(head(notSummarizedBinScores[notSummarizedBinScores$gr %in% unique(notSummarizedBinScores$gr)[group],]))
        notSummarizedBinScores[notSummarizedBinScores$gr %in% unique(notSummarizedBinScores$gr)[group],]$log2ScoresHumanTissueHistoScoresBins = scoreToPlace
        #print(notSummarizedBinScores[notSummarizedBinScores$gr %in% unique(notSummarizedBinScores$gr)[group],])
        
      }
      
      notSummarizedBinScores %>% arrange(EnsemblID)
      
    }
    customBinScoresNotSummarizedTotal   <- assignScoresBinsToAllGenes(customBinScoresSummarized  , customBinScoresNotSummarizedTotal  )
    oneWidthBinScoresNotSummarizedTotal <- assignScoresBinsToAllGenes(oneWidthBinScoresSummarized, oneWidthBinScoresNotSummarizedTotal)
    
    print(head(oneWidthBinScoresNotSummarizedTotal));                  print(tail(oneWidthBinScoresNotSummarizedTotal))
    print(head(as.data.frame(customBinScoresNotSummarizedTotal), 50)); print(tail(as.data.frame(customBinScoresNotSummarizedTotal), 50))
    
    #Use the customBinScores, keep the oneWidthBinScores for reference. Keep also the summarizedTables.
    
    crossValOutput$totalScoreTable[[crossValSet]]$ImmuneTissueScore <- customBinScoresNotSummarizedTotal$log2ScoresHumanTissueHistoScoresBins
    
    
    
    #make a plot for visual inspection
    #function generateBinningPlots creates plots and adds in the densityPerSet column to the summarized data
    
    customBinScoresResult   <- (generateBinningPlots(customBinScoresSummarized,
                                                     title = "Custom bins MannWhitneyU Rank Differences",
                                                     xlabel = "Binned rank difference immune tissues versus non-immune tissues"))
    oneWidthBinScoresResult <- (generateBinningPlots(oneWidthBinScoresSummarized,
                                                     title = "One width bins MannWhitneyU Rank Differences",
                                                     xlabel = "Binned rank difference immune tissues versus non-immune tissues"))

    
      
   
    
    
    
    crossValOutput$extraInfoPerCrossval[[crossValSet]][[length(crossValOutput$extraInfoPerCrossval[[crossValSet]]) + 1]]    <- customBinScoresResult[[1]] ; names(crossValOutput$extraInfoPerCrossval[[crossValSet]])[[length(crossValOutput$extraInfoPerCrossval[[crossValSet]])]]           <- "customBinScoresCountsPerBinImmuneTissueHistMannWhitneyUEstRankDiff"
    crossValOutput$extraInfoPerCrossval[[crossValSet]][[length(crossValOutput$extraInfoPerCrossval[[crossValSet]]) + 1]]    <- customBinScoresResult[[2]] ; names(crossValOutput$extraInfoPerCrossval[[crossValSet]])[[length(crossValOutput$extraInfoPerCrossval[[crossValSet]])]]           <- "customBinsPlotMannWhitneyU"
    
    crossValOutput$extraInfoPerCrossval[[crossValSet]][[length(crossValOutput$extraInfoPerCrossval[[crossValSet]]) + 1]]    <- oneWidthBinScoresResult[[1]] ; names(crossValOutput$extraInfoPerCrossval[[crossValSet]])[[length(crossValOutput$extraInfoPerCrossval[[crossValSet]])]]         <- "oneWidthBinScoresCountsPerBinImmuneTissueHistMannWhitneyUEstRankDiff"
    crossValOutput$extraInfoPerCrossval[[crossValSet]][[length(crossValOutput$extraInfoPerCrossval[[crossValSet]]) + 1]]    <- oneWidthBinScoresResult[[2]] ; names(crossValOutput$extraInfoPerCrossval[[crossValSet]])[[length(crossValOutput$extraInfoPerCrossval[[crossValSet]])]]         <- "oneWidthBinsPlotMannWhitneyU"
    crossValOutput$extraInfoPerCrossval[[crossValSet]][[length(crossValOutput$extraInfoPerCrossval[[crossValSet]]) + 1]]    <- oneWidthBinScoresNotSummarizedTotal ; names(crossValOutput$extraInfoPerCrossval[[crossValSet]])[[length(crossValOutput$extraInfoPerCrossval[[crossValSet]])]] <- "oneWidthBinScoresNotSummarizedImmuneTissueHistMannWhitneyUEstRankDiff"
    
    
  }
  crossValOutput
  
}




#think about what to do with -Inf or +Inf scores, as well as NA scores (right now, scores are calculated for the NA category.)


#function that calculates enrichment of positive set genes in certain gene expression trajectory profiles
#and scores genes accordingly. Similar to binning approaches, but without the need to bin, since all
#genes are pre-assigned to the profiles. NA here has a meaning: these genes do not change.
computeGeneExpressionProfileMacrophageActivationScores <- function(crossValOutput) {
  
  
  #do everything for the amount of times crossval is performed
  for (crossValSet in seq_len(length(crossValOutput$positiveSetNotTrainedOn))) {
    
    
    
    #define positive set and negative set genes that are used for training in this crossVal
    positiveSetGenes <- crossValOutput$positiveSetTrainedOn[[crossValSet]]$EnsemblID
    negativeSetGenes <- crossValOutput$negativeSetTrainedOn[[crossValSet]]$EnsemblID
    
    
    #load profile scores
    oneWidthBinScoresNotSummarizedTrain <- crossValOutput$totalSetTrainedOn[[crossValSet]] %>% select(EnsemblID, Set, Profile) %>%
      group_by(gr = Profile, Set)
    
    oneWidthBinScoresNotSummarizedTotal <- crossValOutput$fullTotalSet[[crossValSet]] %>% select(EnsemblID, Set, Profile) %>%
      group_by(gr = Profile, Set)
    oneWidthBinScoresSummarized <- oneWidthBinScoresNotSummarizedTrain  %>% summarize(Count = n())
    
    #problem: not every Profile contains both PS and NS scores. For score calculations, that is necessary. Need to add these in.
    #hence, function below:
    #for a grouped data frame of profiles and the counts for positive and negative set genes, adds in 
    #an entry with the same profile, the missing set, and a count of 0
    addMissingSetCounts <- function(summarizedData) {
      
      setVector <- c("PS", "NS")
      for(set in seq_len(length(setVector))) {
        for(i in seq_len(length(unique(summarizedData$gr)))) {
          if(!is.na(unique(summarizedData$gr)[i])) {
            if(summarizedData %>% filter(gr == unique(summarizedData$gr)[i], Set == setVector[set]) %>% `$`(Count) %>% length() == 0 ) {
              summarizedData <- summarizedData %>% ungroup %>% add_row(gr = unique(summarizedData$gr)[i], Set = setVector[set], Count = 0) 
            }
          }
        }
        #for NA
        if(summarizedData %>% filter(is.na(gr) == TRUE, Set == setVector[set]) %>% `$`(Count) %>% length() == 0 ) {
          summarizedData <- summarizedData %>% ungroup %>% add_row(gr = NA, Set = setVector[set], Count = 0) }
      }
      summarizedData <- summarizedData %>% group_by(gr) %>% arrange(gr)
      summarizedData
    }
    
    oneWidthBinScoresSummarized <- addMissingSetCounts(oneWidthBinScoresSummarized)
    
    #supply a vector of positive Set and of negative Set genes. Insert the summarized bin data from above.
    #per bin, divide PS counts by PSsize, and divide NS counst by NSsize. Then divide those two by each other and take log2
    #special approach was taken to include NAs. (seq_len(i-1) selects the group which is NA)
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
    
    log2scoresOneWidthBins       <- getLog2Scores(oneWidthBinScoresSummarized, positiveSetGenes, negativeSetGenes)
    oneWidthBinScoresSummarized$log2ScoresGeneExprProfileMacrophage  <- rep(log2scoresOneWidthBins, each = 2)
    
    head(oneWidthBinScoresSummarized); tail(oneWidthBinScoresSummarized)
    
    
    #now give everything a score (traind on + not trained on)
    #only not trained on is kept from this crossVal run, but scores are assigned to everything to make a rank across crossVals.
    #Note that at present, NA are given a score. These are genes that were not in the Human Tissue Score dataset. 
    #I need to discuss with John whether
    #this is really pertinent or not.
    
    
    
    assignScoresBinsToAllGenes          <- function(summarizedBinScores, notSummarizedBinScores) {
      notSummarizedBinScores$log2ScoresGeneExprProfileMacrophage = 0
      print(head(notSummarizedBinScores))
      for (group in seq_len(length(unique(notSummarizedBinScores$gr)))) {
        #print(group)
        #print(is.na(unique(notSummarizedBinScores$gr)[group]))
        
        scoreToPlace <- summarizedBinScores %>%
          filter(gr %in% unique(notSummarizedBinScores$gr)[group], Set == "PS") %>% select(log2ScoresGeneExprProfileMacrophage)  %>% as.data.frame() %>% `[`(1,2, drop = TRUE)
        #print(scoreToPlace)
        #print(head(notSummarizedBinScores[notSummarizedBinScores$gr %in% unique(notSummarizedBinScores$gr)[group],]))
        notSummarizedBinScores[notSummarizedBinScores$gr %in% unique(notSummarizedBinScores$gr)[group],]$log2ScoresGeneExprProfileMacrophage = scoreToPlace
        #print(notSummarizedBinScores[notSummarizedBinScores$gr %in% unique(notSummarizedBinScores$gr)[group],])
        
      }
      
      notSummarizedBinScores %>% arrange(EnsemblID)
      
    }
    
    oneWidthBinScoresNotSummarizedTotal <- assignScoresBinsToAllGenes(oneWidthBinScoresSummarized, oneWidthBinScoresNotSummarizedTotal)
    
    print(head(oneWidthBinScoresNotSummarizedTotal));                  print(tail(oneWidthBinScoresNotSummarizedTotal))
    
    
    #make a plot for visual inspection
    #function generateBinningPlots creates plots and adds in the densityPerSet column to the summarized data
    
    ProfileBinScoresResult   <- (generateBinningPlots(oneWidthBinScoresSummarized,
                                                     title = "Short Time-series Expression Miner Gene Expression Profiles In Macrophages",
                                                     xlabel = "Gene Expression Profile",
                                                     ylabel = "Density of genes in each profile per set"))
    
    #Use the customBinScores, keep the oneWidthBinScores for reference. Keep also the summarizedTables.
    
    crossValOutput$totalScoreTable[[crossValSet]]$MacrophageExpProf <- oneWidthBinScoresNotSummarizedTotal$log2ScoresGeneExprProfileMacrophage
    
    
    crossValOutput$extraInfoPerCrossval[[crossValSet]][[length(crossValOutput$extraInfoPerCrossval[[crossValSet]]) + 1]]    <- ProfileBinScoresResult[[1]]
    names(crossValOutput$extraInfoPerCrossval[[crossValSet]])[[length(crossValOutput$extraInfoPerCrossval[[crossValSet]])]] <- "scoreCountsPerProfileMacrophageGeneExpressionProfile"
    
    crossValOutput$extraInfoPerCrossval[[crossValSet]][[length(crossValOutput$extraInfoPerCrossval[[crossValSet]]) + 1]]    <- ProfileBinScoresResult[[2]]
    names(crossValOutput$extraInfoPerCrossval[[crossValSet]])[[length(crossValOutput$extraInfoPerCrossval[[crossValSet]])]] <- "MacrophageGeneExpressionProfilePlot"
    
    
  }
  crossValOutput
  

}





###########################################################
###########################################################
####
####
####  Generation of Cross-validation sets and calculation of enrichment scores
####
####
###########################################################
###########################################################


#generates cross-validation sets, calculates enrichment scores for all five datasets (depending on binning parameters). debug = TRUE gives
#information on structure, heads, tails, and scores per dataset that scores are calculated for. Function runs other functions which are
#defined at the top of this script.
generateBayesianScores <- function(bayesianTableList, nrCrossVals = 10,
                                   customBinsTFBS = c(-14, -4, -2, 2, 4, 6, 8, 14, 22),
                                   polynomialSmoothBandWidthTFBS = 2.5,
                                   binsZscore = c(0, 2, 10), binsOriginalZscore = c(0, 2, 4, 18),
                                   binsImmuneTissueHistScore = c(-2, 0, 1.0, 1.5, 3.0),
                                   debug = FALSE) {
  
  
  bayesianTableCore  <- bayesianTableList$coreBayesianTable
  bayesianTableTotal <- bayesianTableList$totalBayesianTable
  
  positiveSetIDs <- bayesianTableCore$EnsemblID[bayesianTableCore$Set == "PS"]
  negativeSetIDs <- bayesianTableCore$EnsemblID[bayesianTableCore$EnsemblID %ni% positiveSetIDs]
  
  crossValOutput <- createCrossVal(nrCrossVals, positiveSetIDs, negativeSetIDs, bayesianTableCore)
  

  print(str(crossValOutput))  
  
  
  #Compute the TFBS Score in 2 ways. 1: with predefined bins 2: with polynomial
  #want to put in:
  #vector defining the bins
  #list of sets created by createCrossVal
  #want as output:
  #list of diagnostics per crossval: what is the range of the scores, a boxplot of scores perhaps. How many values have the min/max scores. How many
  #     genes of positive and negative set per bin/score, etc.
  #per crossval set 2 columns of output: binned and polynomial scores
  
  
  #calculate binned TFBS scores. Add to total table.
  
  crossValOutput <- computeTFBSScoresBins(customBinsTFBS, crossValOutput)
  

  #for debugging and testing. Commented out for now.
  if(debug == TRUE) {
    print("Output of first two crossVals after TFBS Scoring (check whether differences are there):")
  crossValOutput$extraInfoPerCrossval[[1]]$oneWidthBinScoresCountsPerBin %>% print()
  crossValOutput$extraInfoPerCrossval[[2]]$oneWidthBinScoresCountsPerBin %>% print()
  
  Sys.sleep(2)
  crossValOutput$totalScoreTable[[1]][crossValOutput$totalScoreTable[[1]]$EnsemblID %in% crossValOutput$positiveSetNotTrainedOn[[1]]$EnsemblID,] %>% print()
  crossValOutput$totalScoreTable[[2]][crossValOutput$totalScoreTable[[2]]$EnsemblID %in% crossValOutput$positiveSetNotTrainedOn[[2]]$EnsemblID,] %>% print()
  Sys.sleep(2)

  }
  
  #calculate smoothened, polynomial-fitted TFBS scores

  crossValOutput <- computeTFBSScoresPolynomial(crossValOutput, smoothBandWidth = polynomialSmoothBandWidthTFBS)
 
  
  if(debug == TRUE) {
    print("Output after smoothing and polynomial-fitting of scores:")
  crossValOutput$extraInfoPerCrossval[[1]]$statisticsBoundedScores %>% print()
  crossValOutput$extraInfoPerCrossval[[2]]$statisticsBoundedScores %>% print()
  crossValOutput$extraInfoPerCrossval[[7]]$statisticsBoundedScores %>% print()
  Sys.sleep(2)
  
  crossValOutput$totalScoreTable[[1]] %>% head(20) %>% print(); crossValOutput$totalScoreTable[[5]] %>% head(20) %>% print()
  Sys.sleep(2)
  
  print("Showing plots with bounded scores for first three crossVals ")
  crossValOutput$extraInfoPerCrossval[[1]]$boundedTFBSSmoothenedLog2ScoresPlot 
  crossValOutput$extraInfoPerCrossval[[2]]$boundedTFBSSmoothenedLog2ScoresPlot
  crossValOutput$extraInfoPerCrossval[[3]]$boundedTFBSSmoothenedLog2ScoresPlot
  }
  
  
  
  #viral PPI enrichment
  crossValOutput <- computeViralPPIScores(crossValOutput)
  
  #viral PPI check
  if(debug == TRUE) {
  crossValOutput$extraInfoPerCrossval[[1]]$viralPPIEnrichmentDetails %>% print()
  crossValOutput$extraInfoPerCrossval[[2]]$viralPPIEnrichmentDetails %>% print()
  
  head(crossValOutput$totalScoreTable[[1]])
  head(crossValOutput$totalScoreTable[[2]])
  #Yes, differences in different crossvals.
  }
  
  
  #Zscores


  crossValOutput <- computeZscoreScores(binsZscore, binsOriginalZscore, crossValOutput, positiveSet)
  #Check ZscoreCalculations
  if(debug == TRUE) {
    print("Output of MHC II surface expression Zscore enrichment calculations:")
    
    crossValOutput$totalScoreTable[[1]] %>% head() %>% print()
    crossValOutput$totalScoreTable[[1]][crossValOutput$totalScoreTable[[1]]$ZscoreHighConf != crossValOutput$totalScoreTable[[1]][1,7],] %>% print()
    Sys.sleep(2)
    crossValOutput$extraInfoPerCrossval[[1]]$oneWidthBinScoresCountsPerBinZscoreHighscore %>% print()
    crossValOutput$extraInfoPerCrossval[[1]]$customBinScoresCountsPerBinZscoreHighscore %>% print()
    Sys.sleep(2)
    
    print("Showing head of 3 different crossValSet score tables")
    head(crossValOutput$totalScoreTable[[2]]) %>% print()
    head(crossValOutput$totalScoreTable[[1]]) %>% print()
    head(crossValOutput$totalScoreTable[[7]]) %>% print()
    Sys.sleep(2)
    
  }
  
  
  
  #immuneTissueHistScore enrichment --> is the gene enriched in antibody stainings of immune tissues (as measured with
  #Mann-Whitney-U sig. rank difference)?

  


  
  
  crossValOutput <- computeImmuneTissueHistScoreBins(binsImmuneTissueHistScore, crossValOutput)
  
  
  #Checking whether that worked
  if(debug == TRUE) {
  #range of values of MWU-test in 1st crossVal:
  print("range of values of MWU-test in 1st crossVal for immune tissue histology scores:")
  range(crossValOutput$fullTotalSet[[1]]$MannWhitneyUEstRankDiff, na.rm = TRUE) %>% print()
  print("heads and tails of score tables:")
  crossValOutput$totalScoreTable[[1]] %>% head(20) %>% print()
  crossValOutput$totalScoreTable[[1]] %>% tail(20) %>% print()
  crossValOutput$totalScoreTable[[2]] %>% head(20) %>% print()
  Sys.sleep(2)
  print("Information on bins and scores (custom Bins first, then one-width bins):")
  crossValOutput$extraInfoPerCrossval[[1]]$customBinScoresCountsPerBinImmuneTissueHistMannWhitneyUEstRankDiff %>% print()
  crossValOutput$extraInfoPerCrossval[[1]]$oneWidthBinScoresCountsPerBinImmuneTissueHistMannWhitneyUEstRankDiff %>% print()
  Sys.sleep(1)
  crossValOutput$extraInfoPerCrossval[[2]]$customBinScoresCountsPerBinImmuneTissueHistMannWhitneyUEstRankDiff
  crossValOutput$extraInfoPerCrossval[[2]]$oneWidthBinScoresCountsPerBinImmuneTissueHistMannWhitneyUEstRankDiff
  Sys.sleep(1)
}
  
  #geneExpressionProfilesMacroPhages. Profiles made with STEM (http://www.cs.cmu.edu/~jernst/stem/), grouping gene expression in human
  #activated macrophages into profiles. Enrichment to see whether certain profiles are informative to being an MHC gene.
  crossValOutput <- computeGeneExpressionProfileMacrophageActivationScores(crossValOutput)
  #check whether this worked
  if(debug == TRUE) {
  crossValOutput$totalScoreTable[[1]]
  

  }
  
  crossValOutput
}
  
  
  

#Resultant scores are sometimes +Inf or -Inf. Many are also NA, though NA scores carry a proper meaning in, I think, all but one datasets.
#The following function rectifies this, bounding infinite scores to the nearest measured upper or lower bound in the data.

#The following discusses, per dataset, the meaning of a NA value:

#TFBS --> pipeline could have identified TFBS in promoters of genes, but didn't. So NAs means: sampled, but nothing found. That is a result, not real missing data.
#though you could argue that you expect at least a TFBS, and it is thus a failure of the methods. This goes for binned.
#TFBS fitting --> If there is no data here that should simply be the case. So these NAs remain.
#virus --> missing data is encoded as 'no'. So NA is absence. Hence no NAs in the data.
#Zscore --> means something, all genes were sampled with siRNAs, those without data were either 1. unconfirmed in follow-up
#or 2. not sig. different. 
#THE ABOVE HAS CHANGED now that the original Zscores are in. There, not being seen just means you are not a protein_coding gene.
#Specifically for the originalZscores, the meaning of NA is thus different.
#Human Tissue Score --> doesn't mean anything. They are created by grouping the expression for every gene 
#in the original data by tissue and cell type, and some genes were simply not measured in some tissues.
#Hence, NAs here should be disregarded.
#Macrophage expression (Profiles) --> does have a meaning, these are genes that have no notable gene trajectories, i.e. are
#mostly unchanging. A group in its own right, so NA-scores are kept.

#For now, NA values are given a score in all datasets, although in the Human Tissue Score, we should disregard them. 

correctBayesianScores <- function(generateBayesianScoresOutput,
                                  debug = FALSE) {
  
  minMaxToCorrectTo <- data.frame()
  minMaxBeforeCorrection <- data.frame()
  for(entries in seq_len(length(generateBayesianScoresOutput$totalScoreTable))) {
    
    #dirty rbind is dirty. Could make an empty dataframe with an nrow and then delete empty rows. Would be better.
    minMaxToCorrectTo      <- rbind(minMaxToCorrectTo, generateBayesianScoresOutput$totalScoreTable[[entries]] %>% 
                                      select(-EnsemblID, -Set) %>%
                                      lapply(FUN = range, finite = TRUE, na.rm = TRUE))
    minMaxBeforeCorrection <- rbind(minMaxBeforeCorrection, generateBayesianScoresOutput$totalScoreTable[[entries]] %>% 
                                      select(-EnsemblID, -Set) %>%
                                      lapply(FUN = range, finite = FALSE, na.rm = FALSE))
    
  }
  #I got both the minimum and maximum. Now disentangle, and write a for-loop to replace them in the relevant datastructures
  maximaToCorrectTo <- minMaxToCorrectTo %>% slice(seq(2,20, by = 2))
  minimaToCorrectTo <- minMaxToCorrectTo %>% slice(seq(1,19, by = 2))
  
  maximaUncorrected <- minMaxBeforeCorrection %>% slice(seq(2,20, by = 2))
  minimaUncorrected <- minMaxBeforeCorrection %>% slice(seq(1,19, by = 2))
  
  if(debug == TRUE) {
    
    print("Showing heads and tails of 1) data values before correction, 2) data to which they should be corrected:")
    minMaxBeforeCorrection %>% head(20) %>% print()
    minMaxToCorrectTo      %>% head(20) %>% print()
    Sys.sleep(2)
    print("Showing heads and tails of 1) uncorrected maxima and minima, 2) maxima and minima to correct to:")
    maximaUncorrected %>% print()
    minimaUncorrected %>% print()
    maximaToCorrectTo %>% print()
    minimaToCorrectTo %>% print()
    Sys.sleep(2)
  }
  
  
  #select only those columns that I want to change (i.e. exclude those that need not be bounded or corrected)
  maximaUncorrectedToChange <- maximaUncorrected %>% select(-PPIScoreUnion, -PPIScoreIntersect, -TFBSScoreFit)
  minimaUncorrectedToChange <- minimaUncorrected %>% select(-PPIScoreUnion, -PPIScoreIntersect, -TFBSScoreFit)
  minimaToCorrectToToChange <- minimaToCorrectTo %>% select(-PPIScoreUnion, -PPIScoreIntersect, -TFBSScoreFit)
  maximaToCorrectToToChange <- maximaToCorrectTo %>% select(-PPIScoreUnion, -PPIScoreIntersect, -TFBSScoreFit)
  
  
  #minimaToCorrectToToChange
  #range(crossValOutput$totalScoreTable[[1]]$TFBSScoreBins, finite = TRUE)
  
  
  #per crossVal, get the scores, find those uncorrected, and bound the minima and maxima
  #NOTE: in the current crossVal, this is only really necessary for the minima, the maxima are already correct.
  #Nevertheless, programmed such that it would work, if ever the maxima do change.
  
  for(crossVal in seq_len(length(generateBayesianScoresOutput$positiveSetNotTrainedOn))) {
    
    columnsToBeChanged <- generateBayesianScoresOutput$totalScoreTable[[crossVal]] %>% select(-EnsemblID, -Set, -PPIScoreUnion, -PPIScoreIntersect, -TFBSScoreFit)
    columnsToKeep      <- generateBayesianScoresOutput$totalScoreTable[[crossVal]] %>% select(EnsemblID, Set, PPIScoreUnion, PPIScoreIntersect, TFBSScoreFit)
    generateBayesianScoresOutput$totalScoreTableCorrected[[crossVal]] <- columnsToBeChanged
    
    #print(crossVal)
    
    for (column in seq_len(ncol(generateBayesianScoresOutput$totalScoreTableCorrected[[crossVal]]))) {
      
      #print(column)
      
      generateBayesianScoresOutput$totalScoreTableCorrected[[crossVal]][generateBayesianScoresOutput$totalScoreTableCorrected[[crossVal]][,column, drop = FALSE] == minimaUncorrectedToChange[crossVal,column][[1]],][,column] <- minimaToCorrectToToChange[crossVal,column][[1]] 
      generateBayesianScoresOutput$totalScoreTableCorrected[[crossVal]][generateBayesianScoresOutput$totalScoreTableCorrected[[crossVal]][,column, drop = FALSE] == maximaUncorrectedToChange[crossVal,column][[1]],][,column] <- maximaToCorrectToToChange[crossVal,column][[1]] 
      
    }
    
    print("Reached the join part")
    generateBayesianScoresOutput$totalScoreTableCorrected[[crossVal]] <- cbind(columnsToKeep, generateBayesianScoresOutput$totalScoreTableCorrected[[crossVal]]) %>%
      select(EnsemblID, Set, PPIScoreUnion, PPIScoreIntersect, TFBSScoreBins, TFBSScoreFit, everything())
    
    
  }
  
  
  generateBayesianScoresOutput
  
  
}




#THe datasets have different variants. This function takes all those variants, and calculates the totalscore for each gene, 
#and its rank (a comparable measure of likeliness to be MHC-related) for all combinations of dataset variants. It then outputs those
#to another table, for inspection and, later, plotting.

#all these different combinations were made to see what was most predictive, and whether continuous and binned enrichment calculations
#differed. Additionally, scores are also calculated for 4 datasets only, leaving 1 out every time. This to show how every dataset
#influences the final score. In the event that, in the future, other dataset versions should be chosen (bins versus fitted, or ppiIntersect
#vs. union), these should be changed.

calculateRanksAndScores <- function(correctBayesianScoresOutput,
                                    debug = FALSE) {

toSum_Union_Fit_LowerConf     <- c("PPIScoreUnion", "TFBSScoreFit",  "ZscoreLowerConf", "ImmuneTissueScore", "MacrophageExpProf")
toSum_Union_Fit_HighConf      <- c("PPIScoreUnion", "TFBSScoreFit",  "ZscoreHighConf",  "ImmuneTissueScore", "MacrophageExpProf")
toSum_Union_Bin_LowerConf     <- c("PPIScoreUnion", "TFBSScoreBins", "ZscoreLowerConf", "ImmuneTissueScore", "MacrophageExpProf")
toSum_Union_Bin_HighConf      <- c("PPIScoreUnion", "TFBSScoreBins", "ZscoreHighConf",  "ImmuneTissueScore", "MacrophageExpProf")
toSum_Intersect_Fit_LowerConf <- c("PPIScoreIntersect", "TFBSScoreFit",  "ZscoreLowerConf", "ImmuneTissueScore", "MacrophageExpProf")
toSum_Intersect_Fit_HighConf  <- c("PPIScoreIntersect", "TFBSScoreFit",  "ZscoreHighConf", "ImmuneTissueScore", "MacrophageExpProf")
toSum_Intersect_Bin_LowerConf <- c("PPIScoreIntersect", "TFBSScoreBins",  "ZscoreLowerConf", "ImmuneTissueScore", "MacrophageExpProf")
toSum_Intersect_Bin_HighConf  <- c("PPIScoreIntersect", "TFBSScoreBins",  "ZscoreHighConf", "ImmuneTissueScore", "MacrophageExpProf")
toSum_Union_Fit_OriginalSet   <- c("PPIScoreUnion", "TFBSScoreFit",  "ZscoreOriginalScreen",  "ImmuneTissueScore", "MacrophageExpProf")

#All variations of 4 datasets with one left out, used to judge the impact of every dataset on final classifier score.
ToSum_NotViralPPI             <- c( "TFBSScoreFit",  "ZscoreHighConf",  "ImmuneTissueScore", "MacrophageExpProf")
ToSum_NotTFBS                 <- c("PPIScoreUnion",  "ZscoreHighConf",  "ImmuneTissueScore", "MacrophageExpProf")
ToSum_NotZscore               <- c("PPIScoreUnion", "TFBSScoreFit",    "ImmuneTissueScore", "MacrophageExpProf")
ToSum_NotImmuneTissueScore    <- c("PPIScoreUnion", "TFBSScoreFit",  "ZscoreHighConf", "MacrophageExpProf")
ToSum_notMacrophageProf      <- c("PPIScoreUnion", "TFBSScoreFit",  "ZscoreHighConf",  "ImmuneTissueScore")


toSumList <- list(Union_Fit_LowerConf = toSum_Union_Fit_LowerConf, Union_Fit_HighConf = toSum_Union_Fit_HighConf,
                  Union_Bin_LowerConf = toSum_Union_Bin_LowerConf, Union_Bin_HighConf = toSum_Union_Bin_HighConf,
                  Intersect_Fit_LowerConf = toSum_Intersect_Fit_LowerConf, Intersect_Fit_HighConf = toSum_Intersect_Fit_HighConf,
                  Intersect_Bin_LowerConf = toSum_Intersect_Bin_LowerConf, Intersect_Bin_HighConf = toSum_Intersect_Bin_HighConf,
                  Union_Fit_OriginalZ     = toSum_Union_Fit_OriginalSet,
                  notViralPPI = ToSum_NotViralPPI, notTFBS = ToSum_NotTFBS, notZScore = ToSum_NotZscore,
                  notImmuneTissueScore = ToSum_NotImmuneTissueScore, notMacrophageProfiles = ToSum_notMacrophageProf)

for(crossVal in seq_len(length(correctBayesianScoresOutput$positiveSetNotTrainedOn))) {
  
  for(scoreMethod in seq_len(length(toSumList))) {
    
    #arrange_at and one_of because of dplyr syntax which works with quosures and nse which I barely understand. This works.
    nameForScoreColumn <- paste0("totalScore_",names(toSumList[scoreMethod]))
    nameForRankColumn  <- paste0("rank_",names(toSumList[scoreMethod]))
    
    totalScore <- correctBayesianScoresOutput$totalScoreTableCorrected[[crossVal]] %>% select( one_of(toSumList[[scoreMethod]]) ) %>% rowSums(na.rm = TRUE)
    correctBayesianScoresOutput$totalScoreTableCorrected[[crossVal]][nameForScoreColumn] <- totalScore
    correctBayesianScoresOutput$totalScoreTableCorrected[[crossVal]] <- correctBayesianScoresOutput$totalScoreTableCorrected[[crossVal]] %>% arrange_at(nameForScoreColumn, desc)
    correctBayesianScoresOutput$totalScoreTableCorrected[[crossVal]][nameForRankColumn] <- seq(1, nrow(correctBayesianScoresOutput$totalScoreTableCorrected[[crossVal]]), by = 1)
    correctBayesianScoresOutput$totalScoreTableCorrected[[crossVal]] <- correctBayesianScoresOutput$totalScoreTableCorrected[[crossVal]] %>% arrange(EnsemblID)
    
    
    
    
    
  }
  
  correctBayesianScoresOutput$totalScoreTableCorrected[[crossVal]] <- correctBayesianScoresOutput$totalScoreTableCorrected[[crossVal]] %>% select(EnsemblID, Set, starts_with("rank"),
                                                                                                                        PPIScoreUnion, PPIScoreIntersect,
                                                                                                                        TFBSScoreBins, TFBSScoreFit, ZscoreHighConf,
                                                                                                                        ZscoreLowerConf, ZscoreOriginalScreen,
                                                                                                                        ImmuneTissueScore, MacrophageExpProf,
                                                                                                                        starts_with("totalScore_"))
  
}

correctBayesianScoresOutput

}


#function below generates a final score table (taking only the scors per gene
#obtained in the cross-validation in which it was not included for calculating the enrichments)
#and adds the prior to the total score (this does not change the ranks, as the prior is the same for every gene)

generateFinalBayesianScoreTable <- function(calculateRanksAndScoresOutput,
                                            expectedNPositiveSetGenes = 200,
                                            debug = FALSE) {
  
  
  
  testSetList <- vector("list", length = length(calculateRanksAndScoresOutput$positiveSetNotTrainedOn))
  
  
  for(crossVal in seq_len(length(calculateRanksAndScoresOutput$positiveSetNotTrainedOn))) {
    
    testSetList[[crossVal]] <- calculateRanksAndScoresOutput$totalScoreTableCorrected[[crossVal]] %>%
      filter(EnsemblID %in% calculateRanksAndScoresOutput$totalSetNotTrainedOn[[crossVal]]$EnsemblID )
    
  }
  
  totalScoreTableBayesianScores <- do.call(rbind, args = testSetList)
  
  totalScoreTableBayesianScores <- totalScoreTableBayesianScores %>% arrange(EnsemblID)
  
  
  
  #add in the prior
  totalEnsemblGenes <- nrow(totalScoreTableBayesianScores)
  ##NOOT: in het bestand in het begin zijn het er 22387. Wat gaat er in de tussentijd mis? Waar verlies ik 470 genen?
  prior = log2(expectedNPositiveSetGenes/totalEnsemblGenes)
  
  #add the prior to all the total scores (note: this leaves the rankings intact)
  totalScoreTableBayesianScores[, startsWith(colnames(totalScoreTableBayesianScores), "totalScore")] <- totalScoreTableBayesianScores[, startsWith(colnames(totalScoreTableBayesianScores), "totalScore")] + prior
  
  calculateRanksAndScoresOutput$totalScoreTableBayesianScores <- totalScoreTableBayesianScores
  calculateRanksAndScoresOutput
  
}
  



#         RUN THESE FUNCTIONS ON ALL POSITIVE SETS

finalScoresBayesianTableMHCIAndII <- MHCIAndIIBayesianTable %>%
  generateBayesianScores() %>%  correctBayesianScores() %>% calculateRanksAndScores() %>% generateFinalBayesianScoreTable()

finalScoresBayesianTableMHCI <- MHCIBayesianTable %>%
  generateBayesianScores(debug = TRUE)# %>%  correctBayesianScores(debug = TRUE) %>% calculateRanksAndScores() %>% generateFinalBayesianScoreTable()

finalScoresBayesianTableMHCII <- MHCIIBayesianTable %>%
  generateBayesianScores() %>%  correctBayesianScores() %>% calculateRanksAndScores() %>% generateFinalBayesianScoreTable()




##
##
## NEW PLOTS
##
##


#function to plot plots for every positive set.

plotBayesianTablePlots <- function(finalScoresBayesianTable,
                                   vectorComparisonEnsemblIDs = c("~/Documents/Project/Programming/FinalBayesianClassifier/Interpretation/Data/Toll-like_receptor_pathway_genes.txt",
                                                                   "~/Documents/Project/Programming/FinalBayesianClassifier/Interpretation/Data/Chemokine signalling pathway _Biomart_genes.txt"),
                                   kaas = "braap") {
  
  plottingTable <- finalScoresBayesianTable %>% select(EnsemblID, Set, starts_with("totalScore_"))
  plottingTable <- plottingTable %>% gather(scoreMethod, Score, 3:ncol(plottingTable))
  plottingTable <- plottingTable %>% left_join(positiveSet, by = c("EnsemblID" = "Ensembl_accession"))
  MHCData <- plottingTable %>% filter(scoreMethod == "totalScore_Union_Fit_HighConf", Set == "PS")
  MHCDataFullSet <- MHCData %>% mutate(`Subdivision types` = "Total positive set")
  
  
  
  
  
  
}


##2 plots: 1 with negative set, toll-like pathway, chemokine signalling pathway, positive set
#1 that compares MHCI, MHC II, and proteasome (subset of MHC I)

plottingTable <- totalScoreTableBayesianScores %>% select(EnsemblID, Set, starts_with("totalScore_"))
plottingTable <- plottingTable %>% gather(scoreMethod, Score, 3:ncol(plottingTable))
plottingTable <- plottingTable %>% left_join(positiveSet, by = c("EnsemblID" = "Ensembl_accession"))
MHCData <- plottingTable %>% filter(scoreMethod == "totalScore_Union_Fit_HighConf", Set == "PS")
MHCDataFullSet <- MHCData %>% mutate(`Subdivision types` = "Total positive set")

violinPlotProteasomeMHCIMHCII <- MHCData %>% ggplot(aes(x = `Subdivision types`, y = Score, fill = `Subdivision types`)) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75), bw = 0.75, width = 0.9) +
  geom_beeswarm(priority = "density", alpha = 0.8, shape = 1, fill = "black") +
  geom_violin(data = MHCDataFullSet, draw_quantiles = c(0.25,0.5,0.75), bw = 0.75, width = 0.9) +
  geom_beeswarm(data = MHCDataFullSet, priority = "density", alpha = 0.8, shape = 1, fill = "black") +
  theme(title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
violinPlotProteasomeMHCIMHCII

boxPlotLocationOfAction <- MHCData %>% ggplot(aes(x = `Location of action`, y = Score, fill = `Location of action`)) +
  geom_boxplot(draw_quantiles = c(0.25,0.5,0.75), bw = 0.75, width = 0.9) +
  geom_beeswarm(priority = "density", alpha = 0.8, shape = 1, fill = "black") +
  geom_boxplot(data = MHCDataFullSet, draw_quantiles = c(0.25,0.5,0.75), bw = 0.75, width = 0.9) +
  geom_beeswarm(data = MHCDataFullSet, priority = "density", alpha = 0.8, shape = 1, fill = "black") +
  theme(title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
boxPlotLocationOfAction


#All sets together

#read in data for Chemokine and Toll-like
TollLike              <- fread("~/Documents/Project/Programming/FinalBayesianClassifier/Interpretation/Data/Toll-like_receptor_pathway_genes.txt")
chemokineSignallingPW <- fread("~/Documents/Project/Programming/FinalBayesianClassifier/Interpretation/Data/Chemokine signalling pathway _Biomart_genes.txt")


fullDataPlotting <- plottingTable %>% filter(scoreMethod == "totalScore_Union_Fit_HighConf")
fullDataPlotting <- fullDataPlotting %>% group_by(scoreMethod, Set) %>% mutate(outlierDivided = Score > median(Score) + IQR(Score) * 1.5) %>% ungroup()
tollLikeDataPlotting <- fullDataPlotting %>%
  filter(EnsemblID %in% TollLike$`Gene stable ID`) %>%
  mutate(Set = "Toll-like pathway")
chemokineDataPlotting <- fullDataPlotting %>% filter(EnsemblID %in% chemokineSignallingPW$`Gene stable ID`) %>% mutate(Set = "Chemokine signalling pathway")

#get negative set top scorers
negativeSetTopScorers <- fullDataPlotting %>% filter(Set == "NS") %>%
  arrange(desc(outlierDivided), desc(Score)) %>%
  slice(1:(n()/200))


fullDataPlotting %>% ggplot(aes(x = Set, y = Score, fill = Set)) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75), bw = 0.75, width = 0.9) +
  geom_beeswarm(data = negativeSetTopScorers, priority = "density", alpha = 0.8, shape = 21, fill = "yellow") +
  geom_beeswarm(data = filter(fullDataPlotting, Set == "PS"),  priority = "density", alpha = 0.8, shape = 1, fill = "black") +
  geom_violin(data = tollLikeDataPlotting, draw_quantiles = c(0.25,0.5,0.75), bw = 0.75, width = 0.9) +
  geom_beeswarm(data = tollLikeDataPlotting, priority = "density", alpha = 0.8, shape = 1, fill = "black") +
  geom_violin(data = chemokineDataPlotting, draw_quantiles = c(0.25,0.5,0.75), bw = 0.75, width = 0.9) +
  geom_beeswarm(data = chemokineDataPlotting, priority = "density", alpha = 0.8, shape = 1, fill = "black") +
  scale_x_discrete(limits = c("NS", "Chemokine signalling pathway", "Toll-like pathway", "PS")) +
  scale_y_continuous(name = "Log 2 likelihood score") +
  scale_fill_manual(limits = c("NS", "Chemokine signalling pathway", "Toll-like pathway", "PS"),
                    values = c("#f8766d", "green", "orange1", "#00bfc4")) +
  theme_bw() +
  theme(title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) 


#barchart with all the scores in order
plottingBarChartAllGenes <- totalScoreTableBayesianScores %>%
  arrange(rank_Union_Fit_HighConf) %>%
  select(EnsemblID, Set, rank_Union_Fit_HighConf,totalScore_Union_Fit_HighConf) %>%
  mutate(barheight = 1, ranksWithoutAllowingTies = seq(1, nrow(plottingBarChartAllGenes)))

subsetPlottingBarChart <- plottingBarChartAllGenes %>% slice(1:500)
head(plottingBarChartAllGenes)

ggplot(data = subsetPlottingBarChart, aes(x = ranksWithoutAllowingTies,
                                          fill = Set,
                                          y = barheight)) +
  geom_bar(stat = "identity", position = "identity", width = 1) +
  scale_y_continuous(limits = c(0,1), expand = c (0,0), name = "", breaks = NULL) +
  scale_x_continuous(expand = c(0,0), name = "ranks", breaks = c(1,100,200,300,400,500))



##TOp scoring candidates from the negative set!




#write the negative set top Scorers to a file
#these are the top 0.5% genes in the negative set.
negativeSetTopScorers %>% write_csv("~/Documents/Project/Programming/FinalBayesianClassifier/Interpretation/NStopScores.csv")







voorCan <-  totalScoreTableBayesianScores %>% arrange(rank_Intersect_Bin_HighConf)
voorCan <- voorCan %>% select(EnsemblID, Set, rank_Intersect_Bin_HighConf) %>% head(50)

voorCan %>% write_csv(path = "~/Documents/Project/Programming/FinalBayesianClassifier/canBespreking.csv")


#ROC curves

#example
data(aSAH)
head(aSAH)



# Build a ROC object and compute the AUC
roc(aSAH$outcome, aSAH$s100b)
roc(outcome ~ s100b, aSAH)
# Smooth ROC curve
roc(outcome ~ s100b, aSAH, smooth=TRUE)
# more options, CI and plotting
roc1 <- roc(aSAH$outcome,
            aSAH$s100b, percent=TRUE,
            # arguments for auc
            partial.auc=c(100, 90), partial.auc.correct=TRUE,
            partial.auc.focus="sens",
            # arguments for ci
            ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,
            # arguments for plot
            plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
            print.auc=TRUE, show.thres=TRUE)

pROC::

head(totalScoreTableBayesianScores)
roc(Set ~ rank_Union_Fit_HighConf, data = totalScoreTableBayesianScores,
    ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,
    plot = TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
    print.auc=TRUE, show.thres=TRUE)

columnsToGoOver <- names(totalScoreTableBayesianScores[,3:11])
setRepeat <- rep("Set",9)

rocList <- vector("list", length = 9)
#do this for all 9 possible metrics
for(i in seq(1,9,1)) {
  
  column = names(totalScoreTableBayesianScores[,3:11])[i]
  print(column)
  
 rocList[[i]] <- roc(totalScoreTableBayesianScores[,"Set"], totalScoreTableBayesianScores[,column],
      ci=TRUE, boot.n=100, ci.alpha=0.9,  stratified=FALSE,
      plot = TRUE, grid=TRUE,
      print.auc=TRUE, show.thres=TRUE, main = column)
 names(rocList)[i] <- column
  
  
}


#Aha, so using the original Z-scores lowers the quality of the final classifier.


roc(totalScoreTableBayesianScores[,"Set"], )


#compare with ROCs for individual data sets
a <- roc(bayesianTableCore[,"Set"], factor(bayesianTableCore$UnionMeasuredPPI, levels = c("no", "yes"), ordered = TRUE),
    ci=TRUE, boot.n=100, ci.alpha=0.9, ci.type = "bars", stratified=FALSE,
    plot = TRUE, grid=TRUE,
    print.auc=TRUE, show.thres=TRUE, main = "UnionPPI") ;a

b <- roc(bayesianTableCore[,"Set"], factor(bayesianTableCore$IntersectMeasuredPPI, levels = c("no", "yes"), ordered = TRUE),
    ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,
    plot = TRUE, grid=TRUE,
    print.auc=TRUE, show.thres=TRUE, main = "IntersectPPI") ; b

roc(bayesianTableCore[,"Set"],bayesianTableCore$totalTFBSScore,
    ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,
    plot = TRUE, grid=TRUE,
    print.auc=TRUE, show.thres=TRUE, main = "TFBSScore")
# could also perhaps compare with a yes/no approach to TFBS




#MannWhitneyU: more negative --> significantly less found in immune Tissues
#more positive --> significantly more in immune tissues 
roc(bayesianTableCore[,"Set"], bayesianTableCore$MannWhitneyUEstRankDiff ,
    ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,
    plot = TRUE, grid=TRUE,
    print.auc=TRUE, show.thres=TRUE, main = "MannWhitneyUEstRankDiff")

#for profiles, you need to put the most high scoring profile first, so that at first you get
#the best profile for being MHC, then the 2nd best, etc.
#in principle, there are 15 profiles, so the profile that is the best MHC-predictive
#should get the number 15, 2nd best 14, etc.
#So, check over all crossVals the predictive value of each profile, filter -Inf and make a ranking, and change the real scores to that ranking
totalProfileDataFrame <- data.frame( gr = character(), Set = character(), Count = numeric(), log2ScoresGeneExprProfileMacrophage = numeric(), densityPerSet = numeric())
for (i in seq_len(length(crossValOutput$positiveSetNotTrainedOn))) {
  
  totalProfileDataFrame <- full_join(totalProfileDataFrame, crossValOutput$extraInfoPerCrossval[[i]]$scoreCountsPerProfileMacrophageGeneExpressionProfile)
}

ranking <- totalProfileDataFrame %>% 
  group_by(gr) %>%
  filter(is.finite(log2ScoresGeneExprProfileMacrophage)) %>%
  summarize(meanScore = mean(log2ScoresGeneExprProfileMacrophage), sdMean = sd(log2ScoresGeneExprProfileMacrophage)  ) %>%
  arrange(desc(meanScore))
ranking

sort(ranking$gr)

#profiles 0, 3 and 12 never have any positive set genes in any of the crossVals. Hence, the scores are always -Inf. I will give those 
#a rank equal to the lowest possible score (that of profile 7)
finalRankingProfiles <- c(ranking$gr, 0, 3, 12)
scoresProfile <- seq(length(finalRankingProfiles), by = -1)
scoresProfile[15:17] <- 4

ranksProfile <- vector('numeric', length = nrow(bayesianTableCore))
for(values in seq_len(length(bayesianTableCore$Profile))) {
  #print(values)
  #print(which(finalRankingProfiles %in% bayesianTableCore$Profile[values]))
  ranksProfile[values] <- scoresProfile[which(as.numeric(finalRankingProfiles) %in% bayesianTableCore$Profile[values])]
  #print(ranksProfile)

}
head(ranksProfile)
ranksProfile

roc(bayesianTableCore[,"Set"], ranksProfile ,
    ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,
    plot = TRUE, grid=TRUE,
    print.auc=TRUE, show.thres=TRUE, main = "GeneExpressionProfile")


roc(bayesianTableCore[,"Set"], bayesianTableCore$Profile ,
    ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,
    plot = TRUE, grid=TRUE,
    print.auc=TRUE, show.thres=TRUE, main = "GeneExpressionProfile")

roc(bayesianTableCore[!is.na(bayesianTableCore$ZscoreHighScore),"Set"], abs(bayesianTableCore$ZscoreHighScore[!is.na(bayesianTableCore$ZscoreHighScore)]) ,
    ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,
    plot = TRUE, grid=TRUE,
    print.auc=TRUE, show.thres=TRUE, main = "ZscoreHighScore")

roc(bayesianTableCore[,"Set"], bayesianTableCore$RescreenOnlyZscores ,
    ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,
    plot = TRUE, grid=TRUE,
    print.auc=TRUE, show.thres=TRUE, main = "ZscoreRescreenOnly")

roc(bayesianTableCore[,"Set"], bayesianTableCore$OriginalZscoreHighestMedian ,
    ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,
    plot = TRUE, grid=TRUE,
    print.auc=TRUE, show.thres=TRUE, main = "ZscoresOriginal")


#why does the Z-score perform so badly?

ZscoreHighScoreBadPerformanceDetective <- bayesianTableCore %>%
                                          select(EnsemblID, Set, ZscoreHighScore, RescreenOnlyZscores) %>%
                                          mutate(ZscoreHighScore = abs(ZscoreHighScore), RescreenOnlyZscores = abs(RescreenOnlyZscores)) %>%
                                          left_join(positiveSet, by = c("EnsemblID" = "Ensembl_accession")) %>%
                                          select(EnsemblID, Ensembl_name, Set, ZscoreHighScore, RescreenOnlyZscores) %>%
                                          arrange(desc(RescreenOnlyZscores))
ZscoreHighScoreBadPerformanceDetective %>% head(30)

#the ROC curve ranks on genes that I disallow in the enrichment because they came from the Neefjes Set. Therefore, I should
#change those to NS in using it to predict things. But that would make it worse? Well, I will try anyways.

#So what was done below: If I were to use only this dataset as a classifier, I need to do that based only on
#the positive set I use in my classifier. I use, for this dataset, in my classifier, a positive set
#without the candidates that came from results in the Neefjes paper (otherwise: circularity)
#So the below code selects those genes from the Neefjes paper, sets their 'Set' to NS, 
#and generates a new table. It also makes the Z-scores absolute (to function in the ROC curve, since
#both a very low Z-score and a very high Z-score are equally predictive of involvement with MHC II surface expression,
#though in different directions (increase/decrease)). 

positiveSetNonNeefjes <- positiveSet %>% filter(Neefjes_accession == "")


changedPSNeefjesToNS <- bayesianTableCore %>%
  filter(Set == "PS", EnsemblID %ni% positiveSetNonNeefjes$Ensembl_accession) %>%
  mutate(Set = "NS")

bayesianTableCoreWithoutNeefjesInPositiveSet <- bayesianTableCore %>%
  filter(Set == "PS", EnsemblID %in% positiveSetNonNeefjes$Ensembl_accession) %>%
  rbind(bayesianTableCore[bayesianTableCore$Set == "NS",]) %>%
  rbind(changedPSNeefjesToNS) %>%
  mutate(ZscoreHighScore = abs(ZscoreHighScore), RescreenOnlyZscores = abs(RescreenOnlyZscores))

bayesianTableCoreWithoutNeefjesInPositiveSet %<>% arrange(desc(ZscoreHighScore))


#try again on this new set
roc(bayesianTableCoreWithoutNeefjesInPositiveSet[,"Set"], bayesianTableCoreWithoutNeefjesInPositiveSet$ZscoreHighScore,
    ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,
    plot = TRUE, grid=TRUE,
    print.auc=TRUE, show.thres=TRUE, main = "ZscoreHighScore")

roc(bayesianTableCoreWithoutNeefjesInPositiveSet[,"Set"], bayesianTableCoreWithoutNeefjesInPositiveSet$RescreenOnlyZscores,
    ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,
    plot = TRUE, grid=TRUE,
    print.auc=TRUE, show.thres=TRUE, main = "RescreenScores")  

roc(bayesianTableCoreWithoutNeefjesInPositiveSet[,"Set"], bayesianTableCoreWithoutNeefjesInPositiveSet$OriginalZscoreHighestMedian,
    ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,
    plot = TRUE, grid=TRUE,
    print.auc=TRUE, show.thres=TRUE, main = "OriginalZscore")  

#The ROC plots already look nicer than those above. Nevertheless, I think I threw away positive controls in cleaning the data,
#and these positive controls would up the predictive value of the set. Let me see whether I can get those back.


#above plots with geom_ROC
rescreenZscoreDF <- data.frame(Set = bayesianTableCoreWithoutNeefjesInPositiveSet[,"Set"],
                         ZscoreHighConfidence = bayesianTableCoreWithoutNeefjesInPositiveSet$ZscoreHighScore,
                         ZscoreRescreenOnly   = bayesianTableCoreWithoutNeefjesInPositiveSet$RescreenOnlyZscores,
                         OriginalZScores      = bayesianTableCoreWithoutNeefjesInPositiveSet$OriginalZscoreHighestMedian)

rescreenZscoreLongFormat <- rescreenZscoreDF %>% melt(id.vars = c("Set"), measure.vars = seq(2, length(rescreenZscoreDF)),
                                                          variable.name = "Data Type")

ggplot(rescreenZscoreLongFormat, aes(d = Set, m = value, linetype = factor(`Data Type`), colour = factor(`Data Type`))) +
  geom_roc(labels = FALSE, n.cuts = 0, size = 1.5) + theme(legend.position = 'none') +
  geom_abline(intercept = c(0,0), slope = 1, colour = "gray") +
  scale_linetype_discrete(name = "Data Type",
                          breaks = c("ZscoreHighConfidence", "ZscoreRescreenOnly", "OriginalZScores"
                                     ),
                          labels = c("High confidence candidates", "Rescreened candidates", "Scores original screen" 
                                     )) +
  scale_colour_discrete(name = "Data Type",
                        breaks = c("ZscoreHighConfidence", "ZscoreRescreenOnly", "OriginalZScores"
                        ),
                        labels = c("High confidence candidates", "Rescreened candidates", "Scores original screen"
                                   ))



#misschien heb ik in verwerking data weggegooid (HLA-DO, HLA-DM?)
#filteren positiveSet op wat niet in NEefjes paper mogelijk reden? --> check hoe filteren gebeurt bij het binnen


kaas <- totalScoreTableBayesianScores %>% arrange(rank_Union_Fit_HighConf)
freek <- cbind(bayesianTableCore, ranksProfile)


totalROCPlottingData <- bayesianTableCore %>% mutate(UnionMeasuredPPI = as.numeric(factor(UnionMeasuredPPI, levels = c("no", "yes"), ordered = TRUE )),
                                                     IntersectMeasuredPPI = as.numeric(factor(IntersectMeasuredPPI, levels = c("no", "yes"), ordered = TRUE )),
                                                     BayesianScore_rank_Union_Fit_HighConf = abs(22357-totalScoreTableBayesianScores$rank_Union_Fit_HighConf),
                                                     ProfileRank = ranksProfile)




#filter the non-ranked profiles
totalROCPlottingData <- totalROCPlottingData %>% select(-Profile)

head(totalROCPlottingData)

totalROCPlottingDataLongFormat <- totalROCPlottingData %>% melt(id.vars = c("Set"), measure.vars = seq(3, length(totalROCPlottingData)),
                                                                variable.name = "Data Type")
totalROCPlottingDataLongFormat

levels = c("BayesianScore_rank_Union_Fit_HighConf",
           "totalTFBSScore",
           "MannWhitneyUEstRankDiff",
           "UnionMeasuredPPI",
           "ZscoreHighScore",
           "ProfileRank",
           "IntersectMeasuredPPI",
           "RescreenOnlyZscores")


ggplot(totalROCPlottingDataLongFormat, aes(d = Set, m = value, linetype = factor(`Data Type`), colour = factor(`Data Type`))) +
  geom_roc(labels = FALSE, n.cuts = 0) +
  geom_abline(intercept = c(0,0), slope = 1, colour = "gray") +
  scale_linetype_discrete(breaks = c("BayesianScore_rank_Union_Fit_HighConf", "MannWhitneyUEstRankDiff", "totalTFBSScore",
                                     "UnionMeasuredPPI", "IntersectMeasuredPPI", "ProfileRank", "ZscoreHighScore",
                                     "RescreenOnlyZscores")) +
  scale_colour_discrete(breaks = c("BayesianScore_rank_Union_Fit_HighConf", "MannWhitneyUEstRankDiff", "totalTFBSScore",
                                   "UnionMeasuredPPI", "IntersectMeasuredPPI", "ProfileRank", "ZscoreHighScore",
                                   "RescreenOnlyZscores"))


totalROCPlottingDataLongFormat %>% head()
totalROCPlottingDataLongFormat %<>% filter(`Data Type` != "RescreenOnlyZscores", `Data Type` != "IntersectMeasuredPPI")

#presentation BioSB figures
pipo <- totalROCPlottingDataLongFormat %>% filter(`Data Type` != "OriginalZscoreHighestMedian")

ggplot(pipo, aes(d = Set, m = value, linetype = factor(`Data Type`), colour = factor(`Data Type`))) +
  geom_roc(labels = FALSE, n.cuts = 0, size = 1.5) + theme(legend.position = 'none') +
  geom_abline(intercept = c(0,0), slope = 1, colour = "gray") +
  scale_linetype_discrete(name = "Data Type",
    breaks = c("BayesianScore_rank_Union_Fit_HighConf", "MannWhitneyUEstRankDiff", "totalTFBSScore",
                                     "UnionMeasuredPPI", "ProfileRank", "ZscoreHighScore"),
                          labels = c("Final Bayesian Classifier", "Immune tissue scores only", "TFBS score only", 
                                     "PPI only", "Co-expression in macrophage only", "MHC II surface expression perturbance only")) +
  scale_colour_discrete(name = "Data Type",
                        breaks = c("BayesianScore_rank_Union_Fit_HighConf", "MannWhitneyUEstRankDiff", "totalTFBSScore",
                                   "UnionMeasuredPPI", "ProfileRank", "ZscoreHighScore"
                                   ),
                        labels = c("Final Bayesian Classifier", "Immune tissue scores only", "TFBS score only", 
                                   "PPI only", "Co-expression in macrophage only", "MHC II surface expression perturbance only"))


ggplot( data = bayesianTableCore, aes(d = Set)) + 
  #geom_roc(aes(m = as.numeric(factor(UnionMeasuredPPI, levels = c("no", "yes"), ordered = TRUE ))), labels = FALSE) +
  #geom_roc(aes(m = as.numeric(factor(IntersectMeasuredPPI, levels = c("no", "yes"), ordered = TRUE ))), labels = FALSE) +
  #geom_roc(aes(m = totalTFBSScore), labels = FALSE) +
  geom_roc(aes(m = MannWhitneyUEstRankDiff), labels = TRUE, cutoffs.at = c(3,2,1,0,-1,-2)) +
  #geom_roc( data = freek, aes(m = ranksProfile), labels = TRUE) +
  #geom_roc(aes(m = ZscoreHighScore), labels = FALSE) +
  #geom_roc(aes(m = RescreenOnlyZscores), labels = FALSE) +
  geom_abline(intercept = c(0,0), slope = 1, colour = "gray") +
  geom_roc(data = kaas,  aes(d = Set, m = abs(22357-rank_Union_Fit_HighConf)), colour = "green", labels = TRUE, n.cuts = 20) + theme_bw() +theme(legend.position = 'none') #+ style_roc()


##make ROCs of partial Bayesian classifier scores

#make alternative total scores, without one data set every time
plottingDataComparativeROCLeaveOneOutDataset <- bayesianTableCore %>% mutate(BayesianScore_rank_Union_Fit_HighConf = abs(22357-totalScoreTableBayesianScores$rank_Union_Fit_HighConf),
                                                                             BayesianWithoutViralPPI = abs(22357-totalScoreTableBayesianScores$rank_notViralPPI),
                                                                             BayesianWithoutTFBS = abs(22357-totalScoreTableBayesianScores$rank_notTFBS),
                                                                             BayesianWithoutZscore = abs(22357-totalScoreTableBayesianScores$rank_notZScore),
                                                                             BayesianWithoutImmuneScore = abs(22357-totalScoreTableBayesianScores$rank_notImmuneTissueScore),
                                                                             BayesianWithoutMacrophageProfs = abs(22357-totalScoreTableBayesianScores$rank_notMacrophageProfiles)
)


longFormatComparisonDataSetData <- plottingDataComparativeROCLeaveOneOutDataset %>% melt(id.vars = c("Set"), measure.vars = seq(3, length(plottingDataComparativeROCLeaveOneOutDataset)),
                                                                variable.name = "Data Type")


longFormatComparisonDataSetData %>% head() 
unique(longFormatComparisonDataSetData$`Data Type`)
longFormatComparisonROCDataToUse <- c("BayesianScore_rank_Union_Fit_HighConf", "BayesianWithoutViralPPI", "BayesianWithoutTFBS",
                                     "BayesianWithoutZscore", "BayesianWithoutImmuneScore", "BayesianWithoutMacrophageProfs")
longFormatComparisonDataSetData %<>% filter(`Data Type` %in% longFormatComparisonROCDataToUse)
unique(longFormatComparisonDataSetData$`Data Type`)
head(longFormatComparisonDataSetData)

ggplot(longFormatComparisonDataSetData, aes(d = Set, m = as.numeric(value), linetype = factor(`Data Type`, levels = longFormatComparisonROCDataToUse), colour = factor(`Data Type`, levels = longFormatComparisonROCDataToUse))) +
  geom_roc(labels = FALSE, n.cuts = 0, size = 1.5) + theme(legend.position = 'none') +
  geom_abline(intercept = c(0,0), slope = 1, colour = "gray") +
  scale_linetype_discrete(name = "Data Type",
                          breaks = c("BayesianScore_rank_Union_Fit_HighConf", "BayesianWithoutViralPPI", "BayesianWithoutTFBS",
                                     "BayesianWithoutZscore", "BayesianWithoutImmuneScore", "BayesianWithoutMacrophageProfs"),
                          labels = c("Final Bayesian Classifier", "no viral PPI data", "no TFBS data", 
                                     "no MHC II perturbation data", "no immune tissue overrepresentation data", "no macrophage expr. profile data")) +
  scale_colour_discrete(name = "Data Type",
                        breaks = c("BayesianScore_rank_Union_Fit_HighConf", "BayesianWithoutViralPPI", "BayesianWithoutTFBS",
                                   "BayesianWithoutZscore", "BayesianWithoutImmuneScore", "BayesianWithoutMacrophageProfs"),
                        labels = c("Final Bayesian Classifier", "no viral PPI data", "no TFBS data", 
                                   "no MHC II perturbation data", "no immune tissue overrepresentation data", "no macrophage expr. profile data"))






#Extract the data for all cross-validations for plotting

outputFigures <- function(crossValOutput, directory) {
  
  classVector = character()
  for(entries in crossValOutput$extraInfoPerCrossval[[1]])
  {classVector <- c(classVector, class(entries)[1])}
  
  plotsToExport <- names(crossValOutput$extraInfoPerCrossval[[1]])[classVector == 'gg']
  
  
  for(crossValNr in seq_len(length(crossValOutput$extraInfoPerCrossval))) {
    
    for(plot in plotsToExport) {
      
      toExport <- crossValOutput$extraInfoPerCrossval[[crossValNr]][[plot]]
      
      if(plot == "oneWidthBinScoresTFBSPlot" ) {
        toExport <- toExport + theme(axis.text.x = element_text(angle = 90, size = 14))
      } else if (plot == "unalteredTFBSSmoothenedLog2ScoresPlot") {
        toExport <- toExport + ggtitle("unaltered continuous TFBS score")
      } else if (plot == "boundedTFBSSmoothenedLog2ScoresPlot") {
        toExport <- toExport + ggtitle("bounded continuous TFBS score")
      }
      
      
      toExport <- toExport + ggtitle(paste0(toExport[['labels']][['title']], " -", crossValNr, "-"))
      
      
      ggsave(paste0(directory, "/", toExport[['labels']][['title']], ".svg"),
             plot = toExport, width = 264.583333, height = 228.6, units = "mm")
  }
  
  
  }
}

outputFigures(crossValOutput, directory = "~//Documents//Project//writing//FIgures//crossValFigures")


#ranges of the scores per dataset. Removes -Inf, doesn't account for NAs.
outputScoresCrossVal <- function(crossValOutput, directory) {
  
  classVector = character()
  for(entries in crossValOutput$extraInfoPerCrossval[[1]])
  {classVector <- c(classVector, class(entries)[1])}
  
 scoresToExport <- names(crossValOutput$extraInfoPerCrossval[[1]])[classVector == 'grouped_df']
 outputList <- list()
  
  
  for(scoreType in scoresToExport) {
    
    #skip the extra Info that has one-width bin scores for all genes based on TFBS
    if(scoreType == "oneWidthBinScoresTFBS") {next()}
    if(grepl("NotSummarized", scoreType, fixed = TRUE)) {next()}
    

    
    for(crossValNr in seq_len(length(crossValOutput$extraInfoPerCrossval))) {
      
      if(crossValNr == 1) {
      totalTable <- crossValOutput$extraInfoPerCrossval[[1]][[scoreType]]
      totalTable$ordering <- seq(1,nrow(totalTable))
      } else {
        toBind <- crossValOutput$extraInfoPerCrossval[[crossValNr]][[scoreType]]
        toBind$ordering <- seq(1, nrow(toBind))
        totalTable <- rbind(totalTable, toBind)
        
      }
    }
    piet <<- as.data.frame(totalTable)
    print(as.data.frame(totalTable))
    #Sys.sleep(5)
    scoreToGet <- sym(colnames(totalTable)[grepl( "log2", colnames(totalTable), fixed = TRUE)])
    tableInfiniteRemoved <- totalTable
    tableInfiniteRemoved <- do.call(data.frame, lapply(tableInfiniteRemoved, function(x) {
      replace(x, is.infinite(x), NA)
    })
    )
    michel <<- print(tableInfiniteRemoved)
    #Sys.sleep(10)
    
    
    summaryTable <- tableInfiniteRemoved %>%
                    group_by(gr, Set) %>%
                    summarise(medianLog2Score = median(!!scoreToGet, na.rm = TRUE),
                              sdLog2Score  = sd(!!scoreToGet, na.rm = TRUE),
                              maxLog2Score = range(!!scoreToGet, finite = TRUE, na.rm = TRUE)[2],
                              minLog2Score = range(!!scoreToGet, finite = TRUE, na.rm = TRUE)[1],
                              medianDensity = median(densityPerSet, na.rm = TRUE),
                              sdDensity = sd(densityPerSet, na.rm = TRUE),
                              maxDensity = range(densityPerSet, finite = TRUE, na.rm = TRUE)[2],
                              minDensity = range(densityPerSet, finite = TRUE, na.rm = TRUE)[1],
                              ordering = mean(ordering))
    
    #print(summaryTable)

    outputList[[scoreType]] <- summaryTable %>% arrange(ordering) %>% select(-ordering)
    
  }
 outputList
}


statisticsOnScoresInCrossVals <- outputScoresCrossVal(crossValOutput, "")

for(i in seq_len(length(statisticsOnScoresInCrossVals))) {
  
  write_csv(statisticsOnScoresInCrossVals[[i]],
            path = paste0("~//Documents//Project//writing//FIgures//crossValFigures","/",
                          "scoreSummary", names(statisticsOnScoresInCrossVals)[i], ".csv" ))
  
}


#get also the bounding statistics
boundingTable <- data.frame(minimumBound = numeric(),
                            maximumBound = numeric(),
                            minBoundedValueCount = numeric(),
                            maxBoundedValueCount = numeric())
for (crossValNr in seq_len(length(crossValOutput))) {
  
  boundingTable <- rbind(boundingTable, crossValOutput$extraInfoPerCrossval[[crossValNr]]$statisticsBoundedScores)
  
}
colnames(boundingTable) <- names(crossValOutput$extraInfoPerCrossval[[2]]$statisticsBoundedScores)
boundingTable$crossValidation <- seq(1, nrow(boundingTable))
boundingTable %<>% select(crossValidation, everything())
boundingTable
boundingTable %>% write_csv(path = paste0("~//Documents//Project//writing//FIgures//crossValFigures","/",
                                          "BoundingStats.csv" ))






#generate figure 2 for report from statistics on scores in Crossvals
generatePlotsWithCrossValConfInt <- function(data, plotLabels = TRUE) {
  
  listPlots <- vector("list", length = length(data))
  counter   <- 1
  for(entry in data) {



plot <- ggplot(entry, aes(x = factor(gr, levels = c("None", unique(gr)[unique(gr) != "None"])), y = medianDensity, fill = Set)) +
  geom_bar(stat = "identity", colour = "black", width = 0.75, position = "dodge") +
  geom_errorbar(aes(x = factor(gr, levels = c("None", unique(gr)[unique(gr) != "None"])),
                    ymin = medianDensity - 1.96 * sdDensity,
                    ymax = medianDensity + 1.96 * sdDensity), position = position_dodge(width = 0.75), width = 0.6) +
  theme_bw() + theme(panel.grid.major.x = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(expand   = c(0,0), name = NULL) +
  scale_y_continuous(expand = c(0,0),   name = "fraction of genes", limits = c(0,1), breaks = seq(0,1,0.1)) +
  ggtitle(label = names(data)[counter]) +
  scale_fill_manual(breaks = c("PS","NS"),
                    values = c("#686caa", "#eabc4f")) +
  theme(axis.text = element_text(size = 20), plot.title = element_text(size = 20),
        legend.text = element_text(size = 20), legend.title = element_text(size = 20),
        axis.title = element_text(size = 20))
#+
#scale_x_manual

if(plotLabels == TRUE) {
  
  if(names(data)[counter] != "oneWidthBinScoresTFBSCountsPerBin") {
  plot <- plot +
    annotate("text", x = as.character(entry[,1][[1]][seq(1,nrow(entry), by = 2)]), y = 0.8,
             label = round(entry[,3][[1]][seq(1,nrow(entry), by = 2)], digits = 2), angle = 90, size = 12) 
  } else {
    
    plot = plot + theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
      scale_y_continuous(expand = c(0,0),   name = "fraction of genes", limits = c(0,0.4), breaks = seq(0,1,0.1)) +
      #annotate("text", x = as.character(entry[,1][[1]][seq(1,nrow(entry), by = 4)]), y = 0.3,
      #         label = round(entry[,3][[1]][seq(1,nrow(entry), by = 4)], digits = 2), angle = 90, size = 12) +
      #annotate("text", x = as.character(entry[,1][[1]][seq(3,nrow(entry), by = 4)]), y = 0.24,
      #         label = round(entry[,3][[1]][seq(3,nrow(entry), by = 4)], digits = 2), angle = 90, size = 12) +
      scale_x_discrete(expand   = c(0.0,0.0), name = NULL) +
      annotate("text", x = as.character(entry[,1][[1]][seq(1,nrow(entry), by = 2)]), y = 0.22,
               label = round(entry[,3][[1]][seq(1,nrow(entry), by = 2)], digits = 2), angle = 90, size = 8)
    
  }
}


listPlots[[counter]] <- plot
names(listPlots)[counter] <- names(data)[counter]
counter <- counter + 1

  }
  
  
listPlots  
}




crossvalConfIntPlots <- generatePlotsWithCrossValConfInt(statisticsOnScoresInCrossVals)
for (plot in seq_len(length(crossvalConfIntPlots))) {
  
  ggsave(paste0("~//Documents//Project//writing//FIgures//crossValFigures//newerFigures","/", names(crossvalConfIntPlots)[plot], ".svg"),
         plot = crossvalConfIntPlots[[plot]],  width = 338.66666667, height = 254, units = "mm")
  
  
  
}






















#generate correlation plots with actual values rather than 0/1
toCorrelate <- bayesianTableCore %>% select(UnionMeasuredPPI, totalTFBSScore, MannWhitneyUEstRankDiff, Profile, ZscoreHighScore)



#generate correlation plots. See per data set which genes are predictive and check whether that coincides


crossValOutput$extraInfoPerCrossval[[1]]$customBinScoresTFBSPlot
#TFBS --> score of 2 onwards =1, otherwise 0
crossValOutput$extraInfoPerCrossval[[1]]$viralPPIUnionPlot
#viralPPiunion --> yes = 1, otherwise 0. Same for Intersect
crossValOutput$extraInfoPerCrossval[[1]]$customBinSZscoreHighConfPlot
#anything with a Z-score is 1, the rest is 0
crossValOutput$extraInfoPerCrossval[[1]]$customBinsPlotMannWhitneyU
#immune tissues overrep. 0-3 is 1, NA and <0 is 0
crossValOutput$extraInfoPerCrossval[[1]]$MacrophageGeneExpressionProfilePlot
crossValOutput$extraInfoPerCrossval[[2]]$MacrophageGeneExpressionProfilePlot
crossValOutput$extraInfoPerCrossval[[4]]$MacrophageGeneExpressionProfilePlot
crossValOutput$extraInfoPerCrossval[[5]]$MacrophageGeneExpressionProfilePlot
#more difficult. Looking at variable ranking from above, where I calculated the mean enrichment score over all 10 crossVals:
ranking
#I will go for everything that is definitely above 0, even without stDev
selectedForUseInCorrelationProfile <- ranking %>% mutate(meanMinusSD = meanScore-sdMean)
selectedForUseInCorrelationProfile$use = with(selectedForUseInCorrelationProfile, ifelse(meanMinusSD > 0, "yes", "no"))
selectedForUseInCorrelationProfile %<>% arrange(desc(use))
selectedForUseInCorrelationProfile
#so, use 11, 10, 15, 6, 9 and 4 to say 1. Other profiles and NA are 0. Note that perhaps I should only use the top three profiles or something..
profilesToUseForOne <- selectedForUseInCorrelationProfile %>% filter(use == "yes") %>% pull(gr)

correlationTable <- bayesianTableCore %>% mutate(ViralPPIOneZero            = ifelse(UnionMeasuredPPI == "yes", 1, 0),
                                                 TFBSScoreOneZero           = ifelse(totalTFBSScore >= 2 & !is.na(totalTFBSScore), 1, 0),
                                                 ImmuneTissueOverRepOneZero = ifelse(MannWhitneyUEstRankDiff >= 0 & !is.na(MannWhitneyUEstRankDiff), 1, 0),
                                                 ProfileMacrophageOneZero   = ifelse(Profile %in% profilesToUseForOne, 1, 0),
                                                 ZscoreHighScoreOneZero     = ifelse(!is.na(OriginalZscoreHighestMedian), 1, 0)) %>%
  select(contains("OneZero"))

correlationList <- vector("list", length = ncol(correlationTable))

for(columns in seq_len(ncol(correlationTable))) {
  
  correlationVectorColumn <- vector("numeric", ncol(correlationTable))
  for(dataSetsToCorrelateAgainst in seq_len(ncol(correlationTable))) {
    
    correlationVectorColumn[dataSetsToCorrelateAgainst] <- cor.test(correlationTable[,columns],
                                                                    correlationTable[,dataSetsToCorrelateAgainst],
                                                                    method = "spearman") %>% `$`(estimate)
    names(correlationVectorColumn)[dataSetsToCorrelateAgainst] <- paste0(colnames(correlationTable)[columns], " against ",
                                                                         colnames(correlationTable)[dataSetsToCorrelateAgainst])
    
    
    
    
  }
  correlationList[[columns]] <- correlationVectorColumn
  
  
}

orderOfDataSets = c("ViralPPI", "TFBS", "ImmuneTissueOverRepresentation", "MacrophageGeneExpressionProfiles", "MHCII surface expression")

correlationResultsTable <- do.call(rbind, correlationList)
rownames(correlationResultsTable) <- orderOfDataSets
colnames(correlationResultsTable) <- orderOfDataSets
correlationResultsTable
#correlationResultsTable <- log(correlationResultsTable)
correlationResultsTable %>% heatmap(keep.dendro = FALSE, symm = TRUE, scale = "row", Rowv = NA, Colv = NA, revC = TRUE)
correlationResultsTable %>% superheat(X.text = round(as.matrix(correlationResultsTable), 3), bottom.label.text.angle = 90, bottom.label.size = 1, bottom.label.text.size = 4,
                                      left.label.text.size = 4, left.label.size = 1)


#make a total figure of datasets

TFBSFigureTalk <- ggplot_build(crossValOutput$extraInfoPerCrossval[[1]]$customBinScoresTFBSPlot)
mannWhitneyUTestFigureTalk <- ggplot_build(crossValOutput$extraInfoPerCrossval[[1]]$customBinsPlotMannWhitneyU)
ViralPPIFigureTalk         <- ggplot_build(crossValOutput$extraInfoPerCrossval[[1]]$viralPPIUnionPlot)
ZscoreFigureTalk           <- ggplot_build(crossValOutput$extraInfoPerCrossval[[1]]$customBinSZscoreHighConfPlot)
macrophageFigureTalk       <- ggplot_build(crossValOutput$extraInfoPerCrossval[[1]]$MacrophageGeneExpressionProfilePlot)
macrophageFigureTalk







###########################
##
##    CODE GRAVEYARD
##
###########################

#OLD CODE FOR MANUALLY ADDING IN SCORES FOR GENES NOT IN CERTAIN BINS
#for both cases, add in any missing rows to allow for +Inf or -Inf scores in output, instead of missing scores
for(i in 1:length(unique(freek2$gr))) {
  if(!is.na(unique(freek2$gr)[i])) {
    if(freek2 %>% filter(gr == unique(freek2$gr)[i], Set == "PS") %>% `$`(Count) %>% length() == 0 ) {
      freek2 <- freek2 %>% ungroup %>% add_row(gr = unique(freek2$gr)[i], Set = "PS", Count = 0) 
    }
    if(freek2 %>% filter(gr == unique(freek2$gr)[i], Set == "NS") %>% `$`(Count) %>% length() == 0 ) {
      freek2 <- freek2 %>% ungroup %>% add_row(gr = unique(freek2$gr)[i], Set = "NS", Count = 0) 
    }
  }
}
#for NA
if(freek2 %>% filter(is.na(gr) == TRUE, Set == "PS") %>% `$`(Count) %>% length() == 0 ) {
  freek2 <- freek2 %>% ungroup %>% add_row(gr = NA, Set = "PS", Count = 0) }
if(freek2 %>% filter(is.na(gr) == TRUE, Set == "NS") %>% `$`(Count) %>% length() == 0 ) {
  freek2 <- freek2 %>% ungroup %>% add_row(gr = NA, Set = "NS", Count = 0) }

}
freek2 <- freek2 %>% group_by(gr) %>% arrange(gr)


for (group in seq_len(length(unique(customBinScoresNotSummarized$gr)))) {
  #print(group)
  #print(is.na(unique(customBinScoresNotSummarized$gr)[group]))
  if(!is.na(unique(customBinScoresNotSummarized$gr)[group])) {
    scoreToPlace <- customBinScoresSummarized %>%
      filter(gr == unique(customBinScoresNotSummarized$gr)[group], Set == "PS") %>% select(log2scores)  %>% as.data.frame() %>% `[`(1,2, drop = TRUE)
    #print(scoreToPlace)
    #print(head(customBinScoresNotSummarized[customBinScoresNotSummarized$gr %in% unique(customBinScoresNotSummarized$gr)[group],]))
    customBinScoresNotSummarized[customBinScoresNotSummarized$gr %in% unique(customBinScoresNotSummarized$gr)[group],]$log2scores = scoreToPlace
    #print(customBinScoresNotSummarized[customBinScoresNotSummarized$gr %in% unique(customBinScoresNotSummarized$gr)[group],])
  } else {
    scoreToPlace <- customBinScoresSummarized %>%
      filter(is.na(unique(customBinScoresNotSummarized$gr)[group]) == TRUE, Set == "PS") %>% select(log2scores) %>% as.data.frame() %>% `[`(1,2, drop = TRUE)
    #print(scoreToPlace)
    #print(group)
    #print(is.na(unique(customBinScoresNotSummarized$gr)[group]))
    customBinScoresNotSummarized[is.na(customBinScoresNotSummarized$gr) == TRUE,]$log2scores = scoreToPlace
  }
}



#old not yet a function binscorecalculator
#

computeTFBSScoresBins <- function(bins, crossValOutput, negativeSetGenes, positiveSetGenes) {
  
  #compute binned scores. Always also include binning intervals of 1 to compare
  
  #custom bins
  customBinScoresNotSummarized <- crossValOutput$totalSetTrainedOn[[1]] %>% select(EnsemblID, Set, totalTFBSScore) %>%
    group_by(gr = cut(totalTFBSScore, breaks = bins), Set)
  customBinScoresSummarized    <- customBinScoresNotSummarized %>% summarize(Count = n())
  
  
  #bins of 1 width
  oneWidthBinScoresNotSummarized <- crossValOutput$totalSetTrainedOn[[1]] %>% select(EnsemblID, Set, totalTFBSScore) %>%
    group_by(gr = cut(totalTFBSScore, breaks = seq(floor(min(totalTFBSScore, na.rm = TRUE)), ceiling(max(totalTFBSScore, na.rm = TRUE)))), Set)
  oneWidthBinScoresSummarized <- oneWidthBinScoresNotSummarized  %>% summarize(Count = n())
  
  #problem: not every bin contains both PS and NS scores. For score calculations, that is necessary. Need to add these in.
  #hence, function below:
  #for a grouped data frame of bins and the counts for positive and negative set genes, adds in 
  #an entry with the same bin, the missing set, and a count of 0
  addMissingSetCounts <- function(summarizedData) {
    
    setVector <- c("PS", "NS")
    for(set in seq_len(length(setVector))) {
      for(i in 1:length(unique(summarizedData$gr))) {
        if(!is.na(unique(summarizedData$gr)[i])) {
          if(summarizedData %>% filter(gr == unique(summarizedData$gr)[i], Set == setVector[set]) %>% `$`(Count) %>% length() == 0 ) {
            summarizedData <- summarizedData %>% ungroup %>% add_row(gr = unique(summarizedData$gr)[i], Set = setVector[set], Count = 0) 
          }
        }
      }
      #for NA
      if(summarizedData %>% filter(is.na(gr) == TRUE, Set == setVector[set]) %>% `$`(Count) %>% length() == 0 ) {
        summarizedData <- summarizedData %>% ungroup %>% add_row(gr = NA, Set = setVector[set], Count = 0) }
    }
    summarizedData <- summarizedData %>% group_by(gr) %>% arrange(gr)
    summarizedData
  }
  
  customBinScoresSummarized   <- addMissingSetCounts(customBinScoresSummarized)
  oneWidthBinScoresSummarized <- addMissingSetCounts(oneWidthBinScoresSummarized)
  
  #supply a vector of positive Set and of negative Set genes. Insert the summarized bin data from above.
  #per bin, divide PS counts by PSsize, and divide NS counst by NSsize. Then divide those two by each other and take log2
  #special approach was taken to include NAs. (seq_len(i-1) selects the group which is NA)
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
  
  
  log2scoresCustomBins         <- getLog2Scores(customBinScoresSummarized, positiveSetTestVector, negativeSetTestVector)
  customBinScoresSummarized$log2ScoresTFBSBins   <- rep(log2scoresCustomBins, each = 2)
  log2scoresOneWidthBins       <- getLog2Scores(oneWidthBinScoresSummarized, positiveSetTestVector, negativeSetTestVector)
  oneWidthBinScoresSummarized$log2ScoresTFBSBins <- rep(log2scoresOneWidthBins, each = 2)
  
  #head(customBinScoresSummarized); tail(customBinScoresSummarized)
  #head(oneWidthBinScoresSummarized); tail(oneWidthBinScoresSummarized)
  
  
  #now give everything a score (traind on + not trained on)
  #only not trained on is kept from this crossVal run
  #Note that at present, NA are given a score. These are genes that were not in the TFBS dataset. SInce only 3 of the positive
  #set were not in the dataset, and many from the negative set, these are given a negative score. I need to discuss with John whether
  #this is really pertinent or not.
  
  
  
  assignScoresBinsToAllGenes     <- function(summarizedBinScores, notSummarizedBinScores) {
    notSummarizedBinScores$log2ScoresTFBSBins = 0
    for (group in seq_len(length(unique(notSummarizedBinScores$gr)))) {
      #print(group)
      #print(is.na(unique(notSummarizedBinScores$gr)[group]))
      
      scoreToPlace <- summarizedBinScores %>%
        filter(gr %in% unique(notSummarizedBinScores$gr)[group], Set == "PS") %>% select(log2ScoresTFBSBins)  %>% as.data.frame() %>% `[`(1,2, drop = TRUE)
      #print(scoreToPlace)
      #print(head(notSummarizedBinScores[notSummarizedBinScores$gr %in% unique(notSummarizedBinScores$gr)[group],]))
      notSummarizedBinScores[notSummarizedBinScores$gr %in% unique(notSummarizedBinScores$gr)[group],]$log2ScoresTFBSBins = scoreToPlace
      #print(notSummarizedBinScores[notSummarizedBinScores$gr %in% unique(notSummarizedBinScores$gr)[group],])
      
    }
    
    notSummarizedBinScores %>% arrange(EnsemblID)
    
  }
  customBinScoresNotSummarized   <- assignScoresBinsToAllGenes(customBinScoresSummarized, customBinScoresNotSummarized)
  oneWidthBinScoresNotSummarized <- assignScoresBinsToAllGenes(oneWidthBinScoresSummarized, oneWidthBinScoresNotSummarized)
  
  head(oneWidthBinScoresNotSummarized); tail(oneWidthBinScoresNotSummarized)
  head(as.data.frame(customBinScoresNotSummarized), 50); tail(as.data.frame(customBinScoresNotSummarized), 50)
  
  
  #Use the customBinScores, keep the oneWidthBinScores for reference. Keep also the summarizedTables.
  
  crossValOutput$totalScoreTable[[1]]$TFBSScoreBins <- customBinScoresNotSummarized$log2ScoresTFBSBins
  crossValOutput$extraInfoPerCrossval[[1]] <- list(customBinScoresCountsPerBin   = customBinScoresSummarized,
                                          oneWidthBinScoresCountsPerBin = oneWidthBinScoresSummarized,
                                          oneWidthBinScores             = oneWidthBinScoresNotSummarized)
  
}


#old not yet a function polynomialscoreCalculator
computeTFBSScoresPolynomial <- function(crossValOutput, smoothBandWidth = 2.5) {
  
  
  #calculate density of negative and positive training set,
  #make a model (20th order polynomial),
  #predict values for genes in each run, but ultimately keep only those from the test set per crossVal run
  
  for(crossValSet in seq_len(length(crossValOutput$positiveSetNotTrainedOn))) {
    
    positiveSetTrainTFBSScores <- crossValOutput$positiveSetTrainedOn[[1]]$totalTFBSScore
    negativeSetTrainTFBSScores <- crossValOutput$negativeSetTrainedOn[[1]]$totalTFBSScore
    
    #use kernel smoothing to get a smoothened density
    densityPositiveTrain       <- density(positiveSetTrainTFBSScores, bw = smoothBandWidth, na.rm = TRUE)
    densityNegativeTrain       <- density(negativeSetTrainTFBSScores, bw = smoothBandWidth, na.rm = TRUE)
    
    #make a 20th order polynomial fit of the predictions (i.e. a very smooth line through the points)
    dfPositive                 <- data.frame(x = densityPositiveTrain$x, y = densityPositiveTrain$y)
    polynomialPositive         <- lm(lm(y ~ poly(x, 20), data = dfPositive))
    dfNegative                 <- data.frame(x=densityNegativeTrain$x, y=densityNegativeTrain$y)
    polynomialNegative         <- lm(y ~ poly(x, 20), data = dfNegative)
    
    
    #Predict stuff
    positivePolyScoresAllGenes <- predict(polynomialPositive, data.frame(x = crossValOutput$fullTotalSet[[1]]$totalTFBSScore))
    negativePolyScoresAllGenes <- predict(polynomialNegative, data.frame(x = crossValOutput$fullTotalSet[[1]]$totalTFBSScore))
    
    #get log2 overrepresentation scores
    log2ScoresPoly            <- log2(positivePolyScoresAllGenes/negativePolyScoresAllGenes)
    print("Head and tail of polynomial log2 scores:")
    print(head(log2ScoresPoly, 20)); print(tail(log2ScoresPoly, 20))
    
    
    #plot(densityNegativeTrain$y ~ densityNegativeTrain$x, col = "red", cex = 0.2)
    #lines(predict(polynomialNegative) ~ densityNegativeTrain$x, col = "darkred")
    #points(densityPositiveTrain$y ~ densityPositiveTrain$x, col = "green", cex = 0.2)
    #lines(predict(polynomialPositive) ~ densityPositiveTrain$x, col = "darkgreen")
    
    #save a ggplot object of the scores
    ggPlotDataFrame <- data.frame(y = c(densityNegativeTrain$y, densityPositiveTrain$y), x = c(densityNegativeTrain$x, densityPositiveTrain$x),
                                  posNeg = c(rep("negative", length(densityNegativeTrain$y)), rep("positive", length(densityPositiveTrain$y))))
    ggPlotDataFramePosOnly <- ggPlotDataFrame %>% filter(posNeg == "positive") %>% mutate(prediction = predict(polynomialPositive))
    ggPlotDataFrameNegOnly <- ggPlotDataFrame %>% filter(posNeg == "negative") %>% mutate(prediction = predict(polynomialNegative))
    head(ggPlotDataFrame); tail(ggPlotDataFrame, 30)
    polynomialDensityPlot <- ggplot(data = ggPlotDataFrame, aes(y = y, x = x, fill = factor(posNeg, levels = c("negative", "positive")))) +
      geom_area(data = ggPlotDataFramePosOnly, aes(y = prediction, x = x),  alpha = 0.8) +
      geom_area(data = ggPlotDataFrameNegOnly, aes(y = prediction, x = x),  alpha = 0.65) +
      geom_point(data = ggPlotDataFrameNegOnly, aes(y = y), colour = "darkorchid4", size = 0.5) +
      geom_point(data = ggPlotDataFramePosOnly, aes(y = y), colour = "goldenrod1", size = 0.5) +
      geom_line(data = ggPlotDataFrameNegOnly, aes(y = prediction), colour = "darkorchid4", size = 0.5) +
      geom_line(data = ggPlotDataFramePosOnly, aes(y = prediction), colour = "goldenrod1", size = 0.5) +
      scale_fill_manual(breaks = c("negative", "positive"), values = c("positive" = "goldenrod1", "negative" = "darkorchid4")) +
      guides(fill = guide_legend(title = NULL)) +
      scale_x_continuous(name = "total TFBS score", expand = c(0,0)) +
      scale_y_continuous(name = "density of values", expand = c(0,0)) +
      coord_cartesian(ylim = c(0,0.12)) +
      #ggtitle("Smoothened distributions of scores for TFBS")  +
      theme_bw() +  #theme(plot.title = element_text(hjust = 0.5), legend.justification = c(1.5,1.5), legend.position = c(0.92,0.5), legend.text = element_text(size = 14)) +
      theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) 
    #annotate("text", scoreSeq, 0.09, label = round(log2Seq, 3), angle = 90, size = 8) + annotate("text", 0, 0.11, label = "Log2 Enrichment Scores", size = 9) + panel_border(remove = TRUE) 
    polynomialDensityPlot
    
    
    #I now take all the data together, and bound the data on a minimum and maximum. SInce the positive set and negative set
    #don't completely overlap, the fit performs badly on the fringes. e.g., positive set genes with very high scores get
    #extremely high log likelihoods. To constrain that, I limit this to the region where the behaviour of the division of
    #the distributions by each other is still normal. I find a minimum and maximum within these fringes, and bound on that.
    #This still constitutes a sort of binning, but there is nothing for it: it must be done.
    dFScores <- data.frame(ID                      = crossValOutput$fullTotalSet[[1]]$EnsemblID,
                           TFBSScore               = crossValOutput$fullTotalSet[[1]]$totalTFBSScore,
                           TFBSPredictScore        = log2ScoresPoly,
                           TFBSPredictScoreBounded = log2ScoresPoly)
    minimumBound          <- min(dFScores[dFScores$TFBSScore >= -15 & dFScores$TFBSScore <= -10,]$TFBSPredictScore, na.rm = TRUE)
    minimumBoundTFBSScore <- dFScores[dFScores$TFBSPredictScore %in% minimumBound,]$TFBSScore
    maximumBound          <- max(dFScores[dFScores$TFBSScore >=  10 & dFScores$TFBSScore <=  15,]$TFBSPredictScore, na.rm = TRUE)
    maximumBoundTFBSScore <- dFScores[dFScores$TFBSPredictScore %in% maximumBound,]$TFBSScore
    dFScores[dFScores$TFBSScore > maximumBoundTFBSScore & !is.na(dFScores$TFBSScore), ]$TFBSPredictScoreBounded <- maximumBound
    dFScores[dFScores$TFBSScore < minimumBoundTFBSScore & !is.na(dFScores$TFBSScore), ]$TFBSPredictScoreBounded <- minimumBound
    
    #Bounding functions correctly. Output the minimum and maximum bound for checking. Also say how many values are affected by that bound
    #In addition, output plots.
    numericVectorBounding <- c(minimumBound = minimumBound, maximumBound = maximumBound,
                               minBoundedValueCount = nrow(dFScores[dFScores$TFBSScore < minimumBoundTFBSScore & !is.na(dFScores$TFBSScore), ]),
                               maxBoundedValueCount = nrow(dFScores[dFScores$TFBSScore > maximumBoundTFBSScore & !is.na(dFScores$TFBSScore), ]))
    unalteredPlot <- ggplot(data = dFScores, aes(y = TFBSPredictScore, x = TFBSScore )) + geom_point(alpha = 0.5 ) +
      scale_x_continuous(breaks = c(-15, -10, -5, 0, 5, 10, 15, 20)) +
      scale_y_continuous(breaks = c(-5, 0, 5, 10), limits = c(-5,10))
    alteredPlot <-   ggplot(data = dFScores, aes(y = TFBSPredictScoreBounded, x = TFBSScore )) + geom_point(alpha = 0.5 ) +
      scale_x_continuous(breaks = c(-15, -10, -5, 0, 5, 10, 15, 20)) +
      scale_y_continuous(breaks = c(-5, 0, 5, 10), limits = c(-5,10))
    #doesn't work perfectly, still some values are perhaps lower than necessary.
    #nevertheless, only ~12 values are affected
    
    crossValOutput$totalScoreTable[[1]]$TFBSScoreFit <- dFScores$TFBSPredictScoreBounded
    crossValOutput$extraInfoPerCrossval[[1]][[4]]    <- polynomialDensityPlot ; names(crossValOutput$extraInfoPerCrossval[[1]])[[4]] <- "smoothenedDistributionPlot"
    crossValOutput$extraInfoPerCrossval[[1]][[5]]    <- unalteredPlot         ; names(crossValOutput$extraInfoPerCrossval[[1]])[[5]] <- "unalteredTFBSSmoothenedLog2ScoresPlot"
    crossValOutput$extraInfoPerCrossval[[1]][[6]]    <- alteredPlot           ; names(crossValOutput$extraInfoPerCrossval[[1]])[[6]] <- "boundedTFBSSmoothenedLog2ScoresPlot"
    crossValOutput$extraInfoPerCrossval[[1]][[7]]    <- numericVectorBounding ; names(crossValOutput$extraInfoPerCrossval[[1]])[[7]] <- "statisticsBoundedScores"
    
    
    
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
}

#1. get density data
# customBinScoresSummarized   <- customBinScoresSummarized  %>% group_by(Set) %>% mutate(densityPerSet = Count/sum(Count))
# oneWidthBinScoresSummarized <- oneWidthBinScoresSummarized %>% group_by(Set) %>% mutate(densityPerSet = Count/sum(Count))
# #2. make plot
# customBinsPlot <- ggplot(customBinScoresSummarized, aes(x = gr, y = densityPerSet, fill = Set)) +
#   geom_bar(stat = "identity", colour = "black", width = 1) +
#   theme_bw() + theme(panel.grid.major.x = element_blank()) +
#   scale_x_discrete(expand   = c(0,0),   name = "Binned rank difference immune tissues versus non-immune tissues") +
#   scale_y_continuous(expand = c(0,0),   name = "Density of scores in bin per set", limits = c(0,1), breaks = seq(0,1,0.1)) +
#   annotate("text", x = unique(customBinScoresSummarized[,1][[1]]), y = 0.8,
#            label = round(unique(customBinScoresSummarized[,4][[1]]), digits = 2), angle = 90, size = 8)
# 
# oneWidthBinsPlot <- ggplot(oneWidthBinScoresSummarized, aes(x = gr, y = densityPerSet, fill = Set)) +
#   geom_bar(stat = "identity", colour = "black", width = 1) +
#   theme_bw() + theme(panel.grid.major.x = element_blank()) +
#   scale_x_discrete(expand   = c(0,0),   name = "Binned rank difference immune tissues versus non-immune tissues") +
#   scale_y_continuous(expand = c(0,0),   name = "Density of scores in bin per set", limits = c(0,1), breaks = seq(0,1,0.1)) +
#   annotate("text", x = unique(oneWidthBinScoresSummarized[,1][[1]]), y = 0.8,
#            label = round(unique(oneWidthBinScoresSummarized[,4][[1]]), digits = 2), angle = 90, size = 8)





ggplot(braam, aes(x = MannWhitneyUEstRankDiff, fill = Set)) +stat_bin(breaks = furtherBreaks, colour = "black",aes(y = ..count../sum(..count..))) + theme_bw()



  
bram2 <- braam %>% mutate(grCustom   = cut(MannWhitneyUEstRankDiff, breaks = binsImmuneTissueHistScore),
                            grStandard = cut(MannWhitneyUEstRankDiff, breaks = seq(floor(min(MannWhitneyUEstRankDiff, na.rm = TRUE)), ceiling(max(MannWhitneyUEstRankDiff, na.rm = TRUE)))))

braam3 <- bram2 %>% group_by(Set, grCustom) %>%
                mutate(density)



#code graveyard entry: Z-scores without adding in the function for the Full Z-score data (~6000 genes, initial screen values)
computeZscoreScores <- function(bins, crossValOutput, orgPositiveSetFile) {
  
  
  
  #do everything for the amount of times crossval is performed
  for (crossValSet in seq_len(length(crossValOutput$positiveSetNotTrainedOn))) {
    
    
    
    #define positive set and negative set genes that are used for training in this crossVal
    positiveSetGenes <- crossValOutput$positiveSetTrainedOn[[crossValSet]]$EnsemblID[crossValOutput$positiveSetTrainedOn[[crossValSet]]$EnsemblID %in% orgPositiveSetFile[orgPositiveSetFile$Neefjes_accession == "",]$Ensembl_accession]
    #correct positive set for the fact that these genes should not have been added based on the Neefjes paper!
    negativeSetGenes <- crossValOutput$negativeSetTrainedOn[[crossValSet]]$EnsemblID
    
    
    #compute binned scores. Always also include binning intervals of 1 to compare
    
    
    #custom bins ZscoreHighscore
    customBinScoresNotSummarizedZscoreHighscoreTrain <- crossValOutput$totalSetTrainedOn[[crossValSet]] %>% select(EnsemblID, Set, ZscoreHighScore) %>%
      filter(EnsemblID %ni% orgPositiveSetFile[orgPositiveSetFile$Neefjes_accession != "",]) %>%
      group_by(gr = cut(ZscoreHighScore, breaks = bins), Set)
    customBinScoresNotSummarizedZscoreHighscoreTotal <- crossValOutput$fullTotalSet[[crossValSet]]      %>% select(EnsemblID, Set, ZscoreHighScore) %>%
      filter(EnsemblID %ni% orgPositiveSetFile[orgPositiveSetFile$Neefjes_accession != "",]) %>%
      group_by(gr = cut(ZscoreHighScore, breaks = bins), Set)
    customBinScoresSummarizedZscoreHighscore         <- customBinScoresNotSummarizedZscoreHighscoreTrain %>% summarize(Count = n())
    
    
    #bins of 1 width ZscoreHighscore
    oneWidthBinScoresNotSummarizedZscoreHighscoreTrain <- crossValOutput$totalSetTrainedOn[[crossValSet]] %>% select(EnsemblID, Set, ZscoreHighScore) %>%
      filter(EnsemblID %ni% orgPositiveSetFile[orgPositiveSetFile$Neefjes_accession != "",]) %>%
      group_by(gr = cut(ZscoreHighScore, breaks = seq(floor(min(ZscoreHighScore, na.rm = TRUE)), ceiling(max(ZscoreHighScore, na.rm = TRUE)))), Set)
    
    oneWidthBinScoresNotSummarizedZscoreHighscoreTotal <- crossValOutput$fullTotalSet[[crossValSet]] %>% select(EnsemblID, Set, ZscoreHighScore) %>%
      filter(EnsemblID %ni% orgPositiveSetFile[orgPositiveSetFile$Neefjes_accession != "",]) %>%
      group_by(gr = cut(ZscoreHighScore, breaks = seq(floor(min(ZscoreHighScore, na.rm = TRUE)), ceiling(max(ZscoreHighScore, na.rm = TRUE)))), Set)
    oneWidthBinScoresSummarizedZscoreHighscore <- oneWidthBinScoresNotSummarizedZscoreHighscoreTrain  %>% summarize(Count = n())
    
    #custom bins ZscoreRescreenOnly
    customBinScoresNotSummarizedZscoreRescreenOnlyTrain <- crossValOutput$totalSetTrainedOn[[crossValSet]] %>% select(EnsemblID, Set, RescreenOnlyZscores) %>%
      filter(EnsemblID %ni% orgPositiveSetFile[orgPositiveSetFile$Neefjes_accession != "",]) %>%
      group_by(gr = cut(RescreenOnlyZscores, breaks = bins), Set)
    customBinScoresNotSummarizedZscoreRescreenOnlyTotal <- crossValOutput$fullTotalSet[[crossValSet]]      %>% select(EnsemblID, Set, RescreenOnlyZscores) %>%
      filter(EnsemblID %ni% orgPositiveSetFile[orgPositiveSetFile$Neefjes_accession != "",]) %>%
      group_by(gr = cut(RescreenOnlyZscores, breaks = bins), Set)
    customBinScoresSummarizedZscoreRescreenOnly         <- customBinScoresNotSummarizedZscoreRescreenOnlyTrain %>% summarize(Count = n())
    
    #bins of 1 width ZscoreRescreenOnly
    oneWidthBinScoresNotSummarizedZscoreRescreenOnlyTrain <- crossValOutput$totalSetTrainedOn[[crossValSet]] %>% select(EnsemblID, Set, RescreenOnlyZscores) %>%
      filter(EnsemblID %ni% orgPositiveSetFile[orgPositiveSetFile$Neefjes_accession != "",]) %>%
      group_by(gr = cut(RescreenOnlyZscores, breaks = seq(floor(min(RescreenOnlyZscores, na.rm = TRUE)), ceiling(max(RescreenOnlyZscores, na.rm = TRUE)))), Set)
    
    oneWidthBinScoresNotSummarizedZscoreRescreenOnlyTotal <- crossValOutput$fullTotalSet[[crossValSet]] %>% select(EnsemblID, Set, RescreenOnlyZscores) %>%
      filter(EnsemblID %ni% orgPositiveSetFile[orgPositiveSetFile$Neefjes_accession != "",]) %>%
      group_by(gr = cut(RescreenOnlyZscores, breaks = seq(floor(min(RescreenOnlyZscores, na.rm = TRUE)), ceiling(max(RescreenOnlyZscores, na.rm = TRUE)))), Set)
    oneWidthBinScoresSummarizedZscoreRescreenOnly <- oneWidthBinScoresNotSummarizedZscoreRescreenOnlyTrain  %>% summarize(Count = n())
    
    
    #problem: not every bin contains both PS and NS scores. For score calculations, that is necessary. Need to add these in.
    #hence, function below:
    #for a grouped data frame of bins and the counts for positive and negative set genes, adds in 
    #an entry with the same bin, the missing set, and a count of 0
    addMissingSetCounts <- function(summarizedData) {
      
      setVector <- c("PS", "NS")
      for(set in seq_len(length(setVector))) {
        for(i in seq_len(length(unique(summarizedData$gr)))) {
          if(!is.na(unique(summarizedData$gr)[i])) {
            if(summarizedData %>% filter(gr == unique(summarizedData$gr)[i], Set == setVector[set]) %>% `$`(Count) %>% length() == 0 ) {
              summarizedData <- summarizedData %>% ungroup %>% add_row(gr = unique(summarizedData$gr)[i], Set = setVector[set], Count = 0) 
            }
          }
        }
        #for NA
        if(summarizedData %>% filter(is.na(gr) == TRUE, Set == setVector[set]) %>% `$`(Count) %>% length() == 0 ) {
          summarizedData <- summarizedData %>% ungroup %>% add_row(gr = NA, Set = setVector[set], Count = 0) }
      }
      summarizedData <- summarizedData %>% group_by(gr) %>% arrange(gr)
      summarizedData
    }
    
    customBinScoresSummarizedZscoreHighscore      <- addMissingSetCounts(customBinScoresSummarizedZscoreHighscore)
    oneWidthBinScoresSummarizedZscoreHighscore    <- addMissingSetCounts(oneWidthBinScoresSummarizedZscoreHighscore)
    
    customBinScoresSummarizedZscoreRescreenOnly   <- addMissingSetCounts(customBinScoresSummarizedZscoreRescreenOnly)
    oneWidthBinScoresSummarizedZscoreRescreenOnly <- addMissingSetCounts(oneWidthBinScoresSummarizedZscoreRescreenOnly)
    
    #supply a vector of positive Set and of negative Set genes. Insert the summarized bin data from above.
    #per bin, divide PS counts by PSsize, and divide NS counst by NSsize. Then divide those two by each other and take log2
    #special approach was taken to include NAs. (seq_len(i-1) selects the group which is NA)
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
    
    
    log2scoresCustomBinsZscoreHighscore                             <- getLog2Scores(customBinScoresSummarizedZscoreHighscore,
                                                                                     positiveSetGenes, negativeSetGenes)
    customBinScoresSummarizedZscoreHighscore$log2ScoresZscoreBins   <- rep(log2scoresCustomBinsZscoreHighscore, each = 2)
    log2scoresOneWidthBinsZscoreHighscore                           <- getLog2Scores(oneWidthBinScoresSummarizedZscoreHighscore,
                                                                                     positiveSetGenes, negativeSetGenes)
    oneWidthBinScoresSummarizedZscoreHighscore$log2ScoresZscoreBins <- rep(log2scoresOneWidthBinsZscoreHighscore , each = 2)
    
    log2scoresCustomBinsZscoreRescreenOnly         <- getLog2Scores(customBinScoresSummarizedZscoreRescreenOnly , positiveSetGenes, negativeSetGenes)
    customBinScoresSummarizedZscoreRescreenOnly$log2ScoresZscoreBins   <- rep(log2scoresCustomBinsZscoreRescreenOnly, each = 2)
    log2scoresOneWidthBinsZscoreRescreenOnly       <- getLog2Scores(oneWidthBinScoresSummarizedZscoreRescreenOnly , positiveSetGenes, negativeSetGenes)
    oneWidthBinScoresSummarizedZscoreRescreenOnly$log2ScoresZscoreBins <- rep(log2scoresOneWidthBinsZscoreRescreenOnly, each = 2)
    
    #head(customBinScoresSummarized); tail(customBinScoresSummarized)
    #head(oneWidthBinScoresSummarized); tail(oneWidthBinScoresSummarized)
    
    
    #now give everything a score (traind on + not trained on)
    #only not trained on is kept from this crossVal run, but scores are assigned to everything to make a rank across crossVals.
    #Note that at present, NA are given a score. These are genes that were not in the Z score dataset. SInce only 3 of the positive
    #set were not in the dataset, and many from the negative set, these are given a negative score. I need to discuss with John whether
    #this is really pertinent or not.
    
    
    
    assignScoresBinsToAllGenes          <- function(summarizedBinScores, notSummarizedBinScores) {
      notSummarizedBinScores$log2ScoresZscoreBins = 0
      print(head(notSummarizedBinScores))
      for (group in seq_len(length(unique(notSummarizedBinScores$gr)))) {
        #print(group)
        #print(is.na(unique(notSummarizedBinScores$gr)[group]))
        
        scoreToPlace <- summarizedBinScores %>%
          filter(gr %in% unique(notSummarizedBinScores$gr)[group], Set == "PS") %>% select(log2ScoresZscoreBins)  %>% as.data.frame() %>% `[`(1,2, drop = TRUE)
        #print(scoreToPlace)
        #print(head(notSummarizedBinScores[notSummarizedBinScores$gr %in% unique(notSummarizedBinScores$gr)[group],]))
        notSummarizedBinScores[notSummarizedBinScores$gr %in% unique(notSummarizedBinScores$gr)[group],]$log2ScoresZscoreBins = scoreToPlace
        #print(notSummarizedBinScores[notSummarizedBinScores$gr %in% unique(notSummarizedBinScores$gr)[group],])
        
      }
      
      notSummarizedBinScores %>% arrange(EnsemblID)
      
    }
    
    customBinScoresNotSummarizedTotalZscoreHighscore      <- assignScoresBinsToAllGenes(customBinScoresSummarizedZscoreHighscore  , customBinScoresNotSummarizedZscoreHighscoreTotal  )
    oneWidthBinScoresNotSummarizedTotalZscoreHighscore    <- assignScoresBinsToAllGenes(oneWidthBinScoresSummarizedZscoreHighscore, oneWidthBinScoresNotSummarizedZscoreHighscoreTotal)
    
    customBinScoresNotSummarizedTotalZscoreRescreenOnly   <- assignScoresBinsToAllGenes(customBinScoresSummarizedZscoreRescreenOnly  , customBinScoresNotSummarizedZscoreRescreenOnlyTotal  )
    oneWidthBinScoresNotSummarizedTotalZscoreRescreenOnly <- assignScoresBinsToAllGenes(oneWidthBinScoresSummarizedZscoreRescreenOnly, oneWidthBinScoresNotSummarizedZscoreRescreenOnlyTotal)
    
    
    #head(oneWidthBinScoresNotSummarizedTotal);                  tail(oneWidthBinScoresNotSummarizedTotal)
    # head(as.data.frame(customBinScoresNotSummarizedTotal), 50); tail(as.data.frame(customBinScoresNotSummarizedTotal), 50)
    
    
    #make a plot for visual inspection
    #function generateBinningPlots creates plots and adds in the densityPerSet column to the summarized data
    
    customBinScoresZscoreHighResult   <- (generateBinningPlots(customBinScoresSummarizedZscoreHighscore,
                                                               title = "Custom bins Z score high confidence set",
                                                               xlabel = "Binned Z scores Neefjens screen"))
    
    oneWidthBinScoresZscoreHighResult <- (generateBinningPlots(oneWidthBinScoresSummarizedZscoreHighscore,
                                                               title = "One width bins Z score high confidence set",
                                                               xlabel = "Binned Z scores Neefjens screen"))
    
    customBinScoresZscoreRescreenResult   <- (generateBinningPlots(customBinScoresSummarizedZscoreRescreenOnly,
                                                                   title = "Custom bins Z score lower confidence set",
                                                                   xlabel = "Binned Z scores Neefjens screen"))
    
    oneWidthBinScoresZscoreRescreenResult <- (generateBinningPlots(oneWidthBinScoresSummarizedZscoreRescreenOnly,
                                                                   title = "One width bins Z score lower confidence set",
                                                                   xlabel = "Binned Z scores Neefjens screen"))
    
    
    
    #Use the customBinScores, keep the oneWidthBinScores for reference. Keep also the summarizedTables.
    
    crossValOutput$totalScoreTable[[crossValSet]]$ZscoreHighConf <- customBinScoresNotSummarizedTotalZscoreHighscore$log2ScoresZscoreBins
    crossValOutput$totalScoreTable[[crossValSet]]$ZscoreLowerConf  <- customBinScoresNotSummarizedTotalZscoreRescreenOnly$log2ScoresZscoreBins
    
    
    crossValOutput$extraInfoPerCrossval[[crossValSet]][[length(crossValOutput$extraInfoPerCrossval[[crossValSet]]) + 1]]    <- customBinScoresZscoreHighResult[[1]] ; names(crossValOutput$extraInfoPerCrossval[[crossValSet]])[[length(crossValOutput$extraInfoPerCrossval[[crossValSet]])]] <- "customBinScoresCountsPerBinZscoreHighscore"
    crossValOutput$extraInfoPerCrossval[[crossValSet]][[length(crossValOutput$extraInfoPerCrossval[[crossValSet]]) + 1]]    <- customBinScoresZscoreHighResult[[2]] ; names(crossValOutput$extraInfoPerCrossval[[crossValSet]])[[length(crossValOutput$extraInfoPerCrossval[[crossValSet]])]] <- "customBinSZscoreHighConfPlot"
    
    crossValOutput$extraInfoPerCrossval[[crossValSet]][[length(crossValOutput$extraInfoPerCrossval[[crossValSet]]) + 1]]    <- oneWidthBinScoresZscoreHighResult[[1]] ; names(crossValOutput$extraInfoPerCrossval[[crossValSet]])[[length(crossValOutput$extraInfoPerCrossval[[crossValSet]])]] <- "oneWidthBinScoresCountsPerBinZscoreHighscore"
    crossValOutput$extraInfoPerCrossval[[crossValSet]][[length(crossValOutput$extraInfoPerCrossval[[crossValSet]]) + 1]]    <- oneWidthBinScoresZscoreHighResult[[2]] ; names(crossValOutput$extraInfoPerCrossval[[crossValSet]])[[length(crossValOutput$extraInfoPerCrossval[[crossValSet]])]] <- "oneWidthBinSZscoreHighConfPlot"
    
    crossValOutput$extraInfoPerCrossval[[crossValSet]][[length(crossValOutput$extraInfoPerCrossval[[crossValSet]]) + 1]]    <- oneWidthBinScoresNotSummarizedTotalZscoreHighscore ; names(crossValOutput$extraInfoPerCrossval[[crossValSet]])[[length(crossValOutput$extraInfoPerCrossval[[crossValSet]])]] <- "oneWidthBinScoresNotSummarizedZscoreHighscore"
    
    
    crossValOutput$extraInfoPerCrossval[[crossValSet]][[length(crossValOutput$extraInfoPerCrossval[[crossValSet]]) + 1]]    <- customBinScoresZscoreRescreenResult[[1]] ; names(crossValOutput$extraInfoPerCrossval[[crossValSet]])[[length(crossValOutput$extraInfoPerCrossval[[crossValSet]])]] <- "customBinScoresCountsPerBinZscoreRescreenOnly"
    crossValOutput$extraInfoPerCrossval[[crossValSet]][[length(crossValOutput$extraInfoPerCrossval[[crossValSet]]) + 1]]    <- customBinScoresZscoreRescreenResult[[2]] ; names(crossValOutput$extraInfoPerCrossval[[crossValSet]])[[length(crossValOutput$extraInfoPerCrossval[[crossValSet]])]] <- "customBinSZscoreLowerConfPlot"
    
    crossValOutput$extraInfoPerCrossval[[crossValSet]][[length(crossValOutput$extraInfoPerCrossval[[crossValSet]]) + 1]]    <- oneWidthBinScoresZscoreRescreenResult[[1]] ; names(crossValOutput$extraInfoPerCrossval[[crossValSet]])[[length(crossValOutput$extraInfoPerCrossval[[crossValSet]])]] <- "oneWidthBinScoresCountsPerBinZscoreRescreenOnly"
    crossValOutput$extraInfoPerCrossval[[crossValSet]][[length(crossValOutput$extraInfoPerCrossval[[crossValSet]]) + 1]]    <- oneWidthBinScoresZscoreRescreenResult[[2]] ; names(crossValOutput$extraInfoPerCrossval[[crossValSet]])[[length(crossValOutput$extraInfoPerCrossval[[crossValSet]])]] <- "oneWidthBinSZscoreLowerConfPlot"
    
    crossValOutput$extraInfoPerCrossval[[crossValSet]][[length(crossValOutput$extraInfoPerCrossval[[crossValSet]]) + 1]]    <- oneWidthBinScoresNotSummarizedTotalZscoreRescreenOnly ; names(crossValOutput$extraInfoPerCrossval[[crossValSet]])[[length(crossValOutput$extraInfoPerCrossval[[crossValSet]])]] <- "oneWidthBinScoresNotSummarizedZscoreRescreenOnly"
    
    
  }
  crossValOutput
  
  
  
  
  
  
  
}


##
##
## OLD PLOTS
##
##
#boxplots of the distributions of highest scores, coloured by set

plottingTable <- totalScoreTableBayesianScores %>% select(EnsemblID, Set, starts_with("totalScore_"))

plottingTable <- plottingTable %>% gather(scoreMethod, Score, 3:ncol(plottingTable))
plottingTable <- plottingTable %>% group_by(scoreMethod, Set) %>% mutate(outlier = Score > median(Score) + IQR(Score) * 1.5) %>% ungroup()
head(plottingTable)
plottingTable <- plottingTable %>% left_join(positiveSet, by = c("EnsemblID" = "Ensembl_accession"))


plottingTable$MHC1OR2[plottingTable$MHC1OR2 == '1']     <- "MHC I only"
plottingTable$MHC1OR2[plottingTable$MHC1OR2 == '2']     <- "MHC II only"
plottingTable$MHC1OR2[plottingTable$MHC1OR2 == '2AND1'] <- "MHC I AND II"
plottingTable <- plottingTable %>% rename("Gene group" = MHC1OR2)
head(plottingTable)
unique(plottingTable$`Gene group`)

#load in the Toll-like receptor pathway
#Note: benefit of the doubt approach in mapping genes to pathway: TNF mapped to 3 ensembl IDs, used all (though could be TNF-a, TNF-b, and TNF-g together?)
TollLike <- fread("~/Documents/Project/Programming/FinalBayesianClassifier/Interpretation/Data/Toll-like_receptor_pathway_genes.txt")
TollLike$Group <- "Toll-like pathway"

plottingTable <- plottingTable %>% left_join(TollLike, by = c("EnsemblID" = "Gene stable ID"))
head(plottingTable); plottingTable
unique(plottingTable$`Gene group`)
plottingTable <- plottingTable %>% mutate(`Gene group` = replace(`Gene group`, which(Group == "Toll-like pathway"), "Toll-like pathway"))
plottingTable <- plottingTable %>% group_by(scoreMethod, `Gene group`) %>% mutate(outlierDivided = Score > median(Score) + IQR(Score) * 1.5) %>% ungroup()
#filter out the one MHC1and2 Gene
plottingTableWithoutMHC1AND2 <- plottingTable %>% filter(`Gene group` %ni% "MHC I AND II")

#Plotting all of the score methods for a comparison
plottingTableForPlottingSubsets <- plottingTable %>% select(EnsemblID, scoreMethod, `Gene group`)
plottingTableForPlottingPositiveNegativeSet <- plottingTable %>% select(EnsemblID, scoreMethod, Set)
plottingTableForPlottingSubsets

k <- ggplot(data = plottingTable, aes(y = Score, fill = Set)) + facet_wrap(~scoreMethod) +
  geom_violin(data = plottingTableWithoutMHC1AND2, aes(x = factor(`Gene group`, levels = c("Toll-like pathway", "MHC I only", "MHC II only"))), draw_quantiles = c(0.25,0.5,0.75), bw = 0.75, width = 0.9) +
  theme_bw() +
  geom_beeswarm(data = function(x) dplyr::filter_(x, ~ x$`Gene group` != "NS" & x$`Gene group` != "MHC I AND II"), aes(x = factor(`Gene group`, levels = c("Toll-like pathway", "MHC I only", "MHC II only"))), priority = 'density', alpha = 0.8, shape = 1) +
  geom_violin(data = plottingTable, aes(y = Score, x = Set, fill = Set), draw_quantiles = c(0.25,0.5,0.75), bw = 0.75, width = 0.9) +
  geom_beeswarm(data = function(x) dplyr::filter_(x, ~ x$`Set` != "NS"), aes(x = Set), priority = 'density', alpha = 0.8, shape = 1) +
  scale_x_discrete(limits = c("NS", "Toll-like pathway", "PS", "MHC I only", "MHC II only"  ), name = "Gene Set")
k


#plotting only the final picked one (Union_Fit_HighConfidenceScore)
finalViolinVisualisationTable <- plottingTable %>% ungroup() %>% filter(scoreMethod == "totalScore_Union_Fit_HighConf")
finalViolinVisualisationTableWithoutMHC1AND2 <- plottingTableWithoutMHC1AND2 %>% ungroup() %>% filter(scoreMethod == "totalScore_Union_Fit_HighConf")

#Load in the chemokine signalling pathway
chemokineSignallingPW <- fread("~/Documents/Project/Programming/FinalBayesianClassifier/Interpretation/Data/Chemokine signalling pathway _Biomart_genes.txt")
chemokineGenesToAlter <- chemokineSignallingPW %>% select(`Gene stable ID`)
chemokineSignallingPlotData <- finalViolinVisualisationTableWithoutMHC1AND2 %>% mutate(`Gene group` = replace(`Gene group`, which(EnsemblID %in% chemokineGenesToAlter$`Gene stable ID`), "Chemokine signalling pathway")) %>% mutate(Set = "NS")
#Some genes are in the positive set. For colouration purposes in the plot, I have labelled all as negative set.

#show the top 0.5% hits in the negative set.
negativeSetTopScorers <- finalViolinVisualisationTableWithoutMHC1AND2 %>% filter(Set == "NS") %>%
  arrange(desc(outlierDivided), desc(Score)) %>%
  slice(1:(n()/200))


finalPlot <- ggplot(data = finalViolinVisualisationTable, aes(y = Score, fill = Set)) +
  geom_violin(data = finalViolinVisualisationTableWithoutMHC1AND2, aes(x = factor(`Gene group`, levels = c("Toll-like pathway", "MHC I only", "MHC II only"))), draw_quantiles = c(0.25,0.5,0.75), bw = 0.75, width = 0.9) +
  theme_bw() +
  geom_beeswarm(data = function(x) dplyr::filter_(x, ~ x$`Gene group` != "NS" & x$`Gene group` != "MHC I AND II"), aes(x = factor(`Gene group`, levels = c("Toll-like pathway", "MHC I only", "MHC II only"))), priority = 'density', alpha = 0.8, shape = 1) +
  geom_violin(data   = finalViolinVisualisationTable, aes(y = Score, x = Set, fill = Set), draw_quantiles = c(0.25,0.5,0.75), bw = 0.75, width = 0.9) +
  geom_beeswarm(data = function(x) dplyr::filter_(x, ~ x$`Set` != "NS"), aes(x = Set), priority = 'density', alpha = 0.8, shape = 1) +
  geom_violin(data   = chemokineSignallingPlotData[chemokineSignallingPlotData$`Gene group` == "Chemokine signalling pathway",], aes(y = Score, x = factor(`Gene group`), fill = Set), draw_quantiles = c(0.25,0.5,0.75), bw = 0.75, width = 0.9) +
  geom_beeswarm(data = chemokineSignallingPlotData[chemokineSignallingPlotData$`Gene group` == "Chemokine signalling pathway",], aes(x = factor(`Gene group`)), priority = 'density', alpha = 0.8, shape = 1) +
  scale_x_discrete(limits = c("NS", "Toll-like pathway", "Chemokine signalling pathway", "PS", "MHC I only", "MHC II only"  ), name = "Gene Set") +
  theme(title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1, size = 15)) + ggtitle("Distributions of scores for positive, negative and known immune gene sets") +
  geom_abline(slope = 0, intercept = median(finalViolinVisualisationTable[finalViolinVisualisationTable$Set == "PS",]$Score), linetype = 5, alpha = 0.7, colour = "#00BFC4") +
  geom_abline(slope = 0, intercept = quantile(finalViolinVisualisationTable[finalViolinVisualisationTable$Set == "NS",]$Score, na.rm = TRUE, probs = 0.5), linetype = 5, alpha = 0.7, colour = "#F8766D") 
finalPlot

finalPlot <- finalPlot + geom_point(data = negativeSetTopScorers, aes(x = Set, y = Score), fill = "green", color = "black", shape = 21, alpha = 0.8, position = position_jitter(0.3, 0.01) )
finalPlot


####
####
#### Original code for generating final bayesian score table (not yet encapsulated in functions)
####

######################################################################################################
######################################################################################################
######                                                                                          ######
######                        Execution of the crossval score calculations                      ######
######                                                                                          ######
######################################################################################################
######################################################################################################



positiveSetIDs <- positiveSet$Ensembl_accession
negativeSetIDs <- bayesianTableCore$EnsemblID[bayesianTableCore$EnsemblID %ni% positiveSetIDs]
nrCrossVals    <- 10
#nrCrossVals    <- 10 #normal setting

crossValOutput <- createCrossVal(nrCrossVals, positiveSetIDs, negativeSetIDs, bayesianTableCore)

length(crossValOutput$positiveSetNotTrainedOn)
crossValOutput$positiveSetTrainedOn[[1]]
crossValOutput$positiveSetNotTrainedOn[[1]]
crossValOutput$negativeSetTrainedOn[[1]]
crossValOutput$negativeSetNotTrainedOn[[1]]
crossValOutput$totalSetTrainedOn[[1]]
crossValOutput$totalScoreTable[[1]]
crossValOutput$extraInfoPerCrossval[[1]]



#Compute the TFBS Score in 2 ways. 1: with predefined bins 2: with polynomial
#want to put in:
#vector defining the bins
#list of sets created by createCrossVal
#want as output:
#list of diagnostics per crossval: what is the range of the scores, a boxplot of scores perhaps. How many values have the min/max scores. How many
#     genes of positive and negative set per bin/score, etc.
#per crossval set 2 columns of output: binned and polynomial scores


#calculate binned TFBS scores




##Tests
customBins = c(-14, -4, -2, 2, 4, 6, 8, 14, 22)
#positiveSetTestVector <- crossValOutput$positiveSetTrainedOn[[1]]$EnsemblID
#negativeSetTestVector <- crossValOutput$negativeSetTrainedOn[[1]]$EnsemblID

crossValOutput <- computeTFBSScoresBins(customBins, crossValOutput)

#check whether this worked, i.e. is something really different between the crossVals?

crossValOutput$extraInfoPerCrossval[[1]]$oneWidthBinScoresCountsPerBin
crossValOutput$extraInfoPerCrossval[[2]]$oneWidthBinScoresCountsPerBin

crossValOutput$totalScoreTable[[1]][crossValOutput$totalScoreTable[[1]]$EnsemblID %in% crossValOutput$positiveSetNotTrainedOn[[1]]$EnsemblID,]
crossValOutput$totalScoreTable[[2]][crossValOutput$totalScoreTable[[2]]$EnsemblID %in% crossValOutput$positiveSetNotTrainedOn[[2]]$EnsemblID,]
#yes

#smoothened, polynomial-fitted TFBS scores
smoothBandWidth = 2.5
crossValOutput <- computeTFBSScoresPolynomial(crossValOutput, smoothBandWidth = smoothBandWidth)
#check whether this worked, i.e. for differences between crossVals
crossValOutput$extraInfoPerCrossval[[1]]$statisticsBoundedScores
crossValOutput$extraInfoPerCrossval[[2]]$statisticsBoundedScores
crossValOutput$extraInfoPerCrossval[[7]]$statisticsBoundedScores

crossValOutput$totalScoreTable[[1]]; crossValOutput$totalScoreTable[[5]]
#I notice that there is an NA in the score of gene ENSG00000460, while it does not have that when scores are binned.
#This is because there is an NA-bin, but no NAs are in the continuous scores.

#viral PPI enrichment
crossValOutput <- computeViralPPIScores(crossValOutput)
#viral PPI check
crossValOutput$extraInfoPerCrossval[[1]]$viralPPIEnrichmentDetails
crossValOutput$extraInfoPerCrossval[[2]]$viralPPIEnrichmentDetails

head(crossValOutput$totalScoreTable[[1]])
head(crossValOutput$totalScoreTable[[2]])
#Yes, differences in different crossvals.


#Zscores
#Note that no Zscores are between [-0.9, 0.9], because they are already only the sig. different Zscores. 
binsZscore         <- c(0, 2, 10)
binsOriginalZscore <- c(0, 2, 4, 18)
crossValOutput <- computeZscoreScores(binsZscore, binsOriginalZscore, crossValOutput, positiveSet)
#did this work?
crossValOutput$totalScoreTable[[1]]
crossValOutput$totalScoreTable[[1]][crossValOutput$totalScoreTable[[1]]$ZscoreHighConf != crossValOutput$totalScoreTable[[1]][1,7],]
crossValOutput$extraInfoPerCrossval[[1]]$oneWidthBinScoresCountsPerBinZscoreHighscore
crossValOutput$extraInfoPerCrossval[[1]]$customBinScoresCountsPerBinZscoreHighscore


head(crossValOutput$totalScoreTable[[2]])
head(crossValOutput$totalScoreTable[[1]])
head(crossValOutput$totalScoreTable[[7]])


#immuneTissueHistScore enrichment --> is it enriched in immune tissues?
binsImmuneTissueHistScore <- c(-2, 0, 1.0, 1.5, 3.0)

crossValOutput$positiveSetTrainedOn[[1]]
range(crossValOutput$fullTotalSet[[1]]$MannWhitneyUEstRankDiff, na.rm = TRUE)
#function is just another incarnation of the TFBSscore function with bins, and the principle is therefore the same.
crossValOutput <- computeImmuneTissueHistScoreBins(binsImmuneTissueHistScore, crossValOutput)


#Checking whether that worked
crossValOutput$totalScoreTable[[1]] %>% head(20)
crossValOutput$totalScoreTable[[1]] %>% tail(20)
crossValOutput$totalScoreTable[[2]] %>% head(20)
crossValOutput$extraInfoPerCrossval[[1]]$customBinScoresCountsPerBinImmuneTissueHistMannWhitneyUEstRankDiff
crossValOutput$extraInfoPerCrossval[[1]]$oneWidthBinScoresCountsPerBinImmuneTissueHistMannWhitneyUEstRankDiff

crossValOutput$extraInfoPerCrossval[[2]]$customBinScoresCountsPerBinImmuneTissueHistMannWhitneyUEstRankDiff
crossValOutput$extraInfoPerCrossval[[2]]$oneWidthBinScoresCountsPerBinImmuneTissueHistMannWhitneyUEstRankDiff


#geneExpressionProfilesMacroPhages
crossValOutput <- computeGeneExpressionProfileMacrophageActivationScores(crossValOutput)
#check whether this worked
crossValOutput$totalScoreTable[[1]]

crossValOutput$extraInfoPerCrossval[[1]]$boundedTFBSSmoothenedLog2ScoresPlot
crossValOutput$extraInfoPerCrossval[[2]]$boundedTFBSSmoothenedLog2ScoresPlot
crossValOutput$extraInfoPerCrossval[[3]]$boundedTFBSSmoothenedLog2ScoresPlot



######################################################################################################
######################################################################################################
######                                                                                          ######
######                                Quality Control and Bounding                              ######
######                                                                                          ######
######################################################################################################
######################################################################################################

#for every score, bound scores with -Inf and +Inf at the lowest and highest real score, respectively
#NAs mean something in almost every data set:
#TFBS --> pipeline could have identified TFBS in promoters of genes, but didn't. So NAs means: sampled, but nothing found. That is a result, not real missing data.
#though you could argue that you expect at least a TFBS, and it is thus a failure of the methods. THis goes for binned.
#TFBS fitting --> If there is no data here that should simply be the case. So these NAs remain.
#virus --> missing data is encoded as 'no'. So NA is absence. Hence no NAs in the data.
#Zscore --> means something, all genes were sampled with siRNAs, those without data were either 1. unconfirmed in follow-up
#or 2. not sig. different. 
#THE ABOVE HAS CHANGED now that the original Zscores are in. There, not being seen just means you are not a protein_coding gene.
#Human Tissue Score --> doesn't mean anything. THey are created by grouping the expression for every gene 
#in the original data by tissue and cell type, and some genes were simply not measured in some tissues.
#Hence, NAs here should be disregarded.
#Macrophage expression (Profiles) --> does have a meaning, these are genes that have no notable gene trajectories, i.e. are
#mostly unchanging. A group in its own right, so NA-scores are kept.








crossValOutput$totalScoreTable[[1]]
min(crossValOutput$totalScoreTable[[1]]$TFBSScoreBins[crossValOutput$totalScoreTable[[1]]$TFBSScoreBins != -Inf])


#Determine the maxima and minima of each of the 10 Cross Validations
#Correct the -Inf and other nonsense scores to the minima

minMaxToCorrectTo <- data.frame()
minMaxBeforeCorrection <- data.frame()
for(entries in seq_len(length(crossValOutput$totalScoreTable))) {
  
  minMaxToCorrectTo      <- rbind(minMaxToCorrectTo, crossValOutput$totalScoreTable[[entries]] %>% 
                                    select(-EnsemblID, -Set) %>%
                                    lapply(FUN = range, finite = TRUE, na.rm = TRUE))
  minMaxBeforeCorrection <- rbind(minMaxBeforeCorrection, crossValOutput$totalScoreTable[[entries]] %>% 
                                    select(-EnsemblID, -Set) %>%
                                    lapply(FUN = range, finite = FALSE, na.rm = FALSE))
  
  
}
#I got both the minimum and maximum. Now disentangle, and write a for-loop to replace them in the relevant datastructures
maximaToCorrectTo <- minMaxToCorrectTo %>% slice(seq(2,20, by = 2))
minimaToCorrectTo <- minMaxToCorrectTo %>% slice(seq(1,19, by = 2))
maximaToCorrectTo
minMaxBeforeCorrection
maximaUncorrected <- minMaxBeforeCorrection %>% slice(seq(2,20, by = 2))
minimaUncorrected <- minMaxBeforeCorrection %>% slice(seq(1,19, by = 2))

maximaUncorrected
minimaUncorrected
minimaToCorrectTo

#select only those columns that I want to change
maximaUncorrectedToChange <- maximaUncorrected %>% select(-PPIScoreUnion, -PPIScoreIntersect, -TFBSScoreFit)
minimaUncorrectedToChange <- minimaUncorrected %>% select(-PPIScoreUnion, -PPIScoreIntersect, -TFBSScoreFit)
minimaToCorrectToToChange <- minimaToCorrectTo %>% select(-PPIScoreUnion, -PPIScoreIntersect, -TFBSScoreFit)
maximaToCorrectToToChange <- maximaToCorrectTo %>% select(-PPIScoreUnion, -PPIScoreIntersect, -TFBSScoreFit)
minimaToCorrectToToChange
range(crossValOutput$totalScoreTable[[1]]$TFBSScoreBins, finite = TRUE)


#per crossVal, get the scores, find those uncorrected, and bound the minima and maxima
#NOTE: in the current crossVal, this is only really necessary for the minima, the maxima are already correct.
#Nevertheless, programmed such that it would work, if ever the maxima do change

for(crossVal in seq_len(length(crossValOutput$positiveSetNotTrainedOn))) {
  
  columnsToBeChanged <- crossValOutput$totalScoreTable[[crossVal]] %>% select(-EnsemblID, -Set, -PPIScoreUnion, -PPIScoreIntersect, -TFBSScoreFit)
  columnsToKeep      <- crossValOutput$totalScoreTable[[crossVal]] %>% select(EnsemblID, Set, PPIScoreUnion, PPIScoreIntersect, TFBSScoreFit)
  crossValOutput$totalScoreTableCorrected[[crossVal]] <- columnsToBeChanged
  
  print(crossVal)
  
  for (column in seq_len(ncol(crossValOutput$totalScoreTableCorrected[[crossVal]]))) {
    
    print(column)
    
    crossValOutput$totalScoreTableCorrected[[crossVal]][crossValOutput$totalScoreTableCorrected[[crossVal]][,column, drop = FALSE] == minimaUncorrectedToChange[crossVal,column][[1]],][,column] <- minimaToCorrectToToChange[crossVal,column][[1]] 
    crossValOutput$totalScoreTableCorrected[[crossVal]][crossValOutput$totalScoreTableCorrected[[crossVal]][,column, drop = FALSE] == maximaUncorrectedToChange[crossVal,column][[1]],][,column] <- maximaToCorrectToToChange[crossVal,column][[1]] 
    
  }
  
  print("Reached the join part")
  crossValOutput$totalScoreTableCorrected[[crossVal]] <- cbind(columnsToKeep, crossValOutput$totalScoreTableCorrected[[crossVal]]) %>%
    select(EnsemblID, Set, PPIScoreUnion, PPIScoreIntersect, TFBSScoreBins, TFBSScoreFit, everything())
  
  
}

#DONE. Now select pertinent columns and do rowSums to get ranks


toSum_Union_Fit_LowerConf     <- c("PPIScoreUnion", "TFBSScoreFit",  "ZscoreLowerConf", "ImmuneTissueScore", "MacrophageExpProf")
toSum_Union_Fit_HighConf      <- c("PPIScoreUnion", "TFBSScoreFit",  "ZscoreHighConf",  "ImmuneTissueScore", "MacrophageExpProf")
toSum_Union_Bin_LowerConf     <- c("PPIScoreUnion", "TFBSScoreBins", "ZscoreLowerConf", "ImmuneTissueScore", "MacrophageExpProf")
toSum_Union_Bin_HighConf      <- c("PPIScoreUnion", "TFBSScoreBins", "ZscoreHighConf",  "ImmuneTissueScore", "MacrophageExpProf")
toSum_Intersect_Fit_LowerConf <- c("PPIScoreIntersect", "TFBSScoreFit",  "ZscoreLowerConf", "ImmuneTissueScore", "MacrophageExpProf")
toSum_Intersect_Fit_HighConf  <- c("PPIScoreIntersect", "TFBSScoreFit",  "ZscoreHighConf", "ImmuneTissueScore", "MacrophageExpProf")
toSum_Intersect_Bin_LowerConf <- c("PPIScoreIntersect", "TFBSScoreBins",  "ZscoreLowerConf", "ImmuneTissueScore", "MacrophageExpProf")
toSum_Intersect_Bin_HighConf  <- c("PPIScoreIntersect", "TFBSScoreBins",  "ZscoreHighConf", "ImmuneTissueScore", "MacrophageExpProf")
toSum_Union_Fit_OriginalSet   <- c("PPIScoreUnion", "TFBSScoreFit",  "ZscoreOriginalScreen",  "ImmuneTissueScore", "MacrophageExpProf")

#add in scores that are like the final score but with 1 data set missing
ToSum_NotViralPPI             <- c( "TFBSScoreFit",  "ZscoreHighConf",  "ImmuneTissueScore", "MacrophageExpProf")
ToSum_NotTFBS                 <- c("PPIScoreUnion",  "ZscoreHighConf",  "ImmuneTissueScore", "MacrophageExpProf")
ToSum_NotZscore               <- c("PPIScoreUnion", "TFBSScoreFit",    "ImmuneTissueScore", "MacrophageExpProf")
ToSum_NotImmuneTissueScore    <- c("PPIScoreUnion", "TFBSScoreFit",  "ZscoreHighConf", "MacrophageExpProf")
ToSum_notMacrophageProf      <- c("PPIScoreUnion", "TFBSScoreFit",  "ZscoreHighConf",  "ImmuneTissueScore")


toSumList <- list(Union_Fit_LowerConf = toSum_Union_Fit_LowerConf, Union_Fit_HighConf = toSum_Union_Fit_HighConf,
                  Union_Bin_LowerConf = toSum_Union_Bin_LowerConf, Union_Bin_HighConf = toSum_Union_Bin_HighConf,
                  Intersect_Fit_LowerConf = toSum_Intersect_Fit_LowerConf, Intersect_Fit_HighConf = toSum_Intersect_Fit_HighConf,
                  Intersect_Bin_LowerConf = toSum_Intersect_Bin_LowerConf, Intersect_Bin_HighConf = toSum_Intersect_Bin_HighConf,
                  Union_Fit_OriginalZ     = toSum_Union_Fit_OriginalSet,
                  notViralPPI = ToSum_NotViralPPI, notTFBS = ToSum_NotTFBS, notZScore = ToSum_NotZscore,
                  notImmuneTissueScore = ToSum_NotImmuneTissueScore, notMacrophageProfiles = ToSum_notMacrophageProf)

for(crossVal in seq_len(length(crossValOutput$positiveSetNotTrainedOn))) {
  
  for(scoreMethod in seq_len(length(toSumList))) {
    
    #arrange_at and one_of because of dplyr syntax which works with quosures and nse which I barely understand. This works.
    nameForScoreColumn <- paste0("totalScore_",names(toSumList[scoreMethod]))
    nameForRankColumn  <- paste0("rank_",names(toSumList[scoreMethod]))
    
    totalScore <- crossValOutput$totalScoreTableCorrected[[crossVal]] %>% select( one_of(toSumList[[scoreMethod]]) ) %>% rowSums(na.rm = TRUE)
    crossValOutput$totalScoreTableCorrected[[crossVal]][nameForScoreColumn] <- totalScore
    crossValOutput$totalScoreTableCorrected[[crossVal]] <- crossValOutput$totalScoreTableCorrected[[crossVal]] %>% arrange_at(nameForScoreColumn, desc)
    crossValOutput$totalScoreTableCorrected[[crossVal]][nameForRankColumn] <- seq(1, nrow(crossValOutput$totalScoreTableCorrected[[crossVal]]), by = 1)
    crossValOutput$totalScoreTableCorrected[[crossVal]] <- crossValOutput$totalScoreTableCorrected[[crossVal]] %>% arrange(EnsemblID)
    
    
    
    
    
  }
  
  crossValOutput$totalScoreTableCorrected[[crossVal]] <- crossValOutput$totalScoreTableCorrected[[crossVal]] %>% select(EnsemblID, Set, starts_with("rank"),
                                                                                                                        PPIScoreUnion, PPIScoreIntersect,
                                                                                                                        TFBSScoreBins, TFBSScoreFit, ZscoreHighConf,
                                                                                                                        ZscoreLowerConf, ZscoreOriginalScreen,
                                                                                                                        ImmuneTissueScore, MacrophageExpProf,
                                                                                                                        starts_with("totalScore_"))
  
}






###Now give me one total table with all the scores and rankings for all the different combinations of scoring methods tested

testSetList <- vector("list", length = length(crossValOutput$positiveSetNotTrainedOn))


for(crossVal in seq_len(length(crossValOutput$positiveSetNotTrainedOn))) {
  
  testSetList[[crossVal]] <- crossValOutput$totalScoreTableCorrected[[crossVal]] %>%
    filter(EnsemblID %in% crossValOutput$totalSetNotTrainedOn[[crossVal]]$EnsemblID )
  
}

totalScoreTableBayesianScores <- do.call(rbind, args = testSetList)

totalScoreTableBayesianScores <- totalScoreTableBayesianScores %>% arrange(EnsemblID)



##I have my final score table
head(totalScoreTableBayesianScores)

#add in the prior
expectedNrMHCGenes <- 200
totalEnsemblGenes <- nrow(totalScoreTableBayesianScores)
##NOOT: in het bestand in het begin zijn het er 22387. Wat gaat er in de tussentijd mis? Waar verlies ik 470 genen?
prior = log2(expectedNrMHCGenes/totalEnsemblGenes)

#add the prior to all the total scores (note: this leaves the rankings intact)
totalScoreTableBayesianScores[, startsWith(colnames(totalScoreTableBayesianScores), "totalScore")] <- totalScoreTableBayesianScores[, startsWith(colnames(totalScoreTableBayesianScores), "totalScore")] + prior

totalScoreTableBayesianScores


#look up the score of specific candidates
totalScoreTableBayesianScores %>% filter(EnsemblID == "ENSG00000198821")


