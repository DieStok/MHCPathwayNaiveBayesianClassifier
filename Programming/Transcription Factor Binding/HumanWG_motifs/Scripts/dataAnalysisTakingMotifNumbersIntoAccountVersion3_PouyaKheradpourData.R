##Dieter Stoker


##########################################################
#Initial options
##########################################################

options(stringsAsFactors = FALSE)
install.packages("pacman", repos = "https://cloud.r-project.org/")
library("pacman")
pacman::p_load(ggplot2, dplyr, stringr, purrr, readr, LaF, ff, data.table, R.utils, base, utils, reshape2, tcltk)
options(datatable.verbose = FALSE)

oldwd <- getwd()
setwd("~/Documents/Project/Programming/Transcription Factor Binding/HumanWG_motifs/Data/singleLineFinalDf and Unified table/")
`%ni%` <- Negate(`%in%`)















##########################################################
#Function to load in the data and required libraries
#reloadFromSource = whether source files like positive list etc. should be selected manually
##########################################################


#####
##
##Testing stuff
##
####
# testTotalGenesMultipleLines <- data.frame(EnsemblId = c("ENSG0000024", "ENSG0000024", "ENSG0000024"  , "ENSG0000024", "ENSG0000048", "ENSG0000048", "ENSG0000048", "ENSG0000060", "ENSG0000060", "ENSG0000060", "ENSG0000060"),
#                                           TFBS =      c("Mef2_5"     , "Drakl3"     , "Frakki_known5", "Drakl2"     , "Drakl3"     , "Mokki"      , "TCF::Sjaak" , "Drakl3"     , "Misantroop" , "Vondelpark" , "Mef2_5")     ,
#                                           Counts =    c(2            ,  5           , 33              , 23            , 4            , 20           , 15            , 100            , 2            , 13            ,  100))
# 
# testTotalGenesOneLine       <- data.frame(EnsemblId = c("ENSG0000024"                                  , "ENSG0000048"                    , "ENSG0000060"),
#                                           TFBS      = c("Mef2_5-16|Drakl3-5|Frakki_known5-33|Drakl2-23", "Drakl3-4|Mokki-20|TCF::Sjaak-15", "Drakl3-100|Misantroop-2|Vondelpark-1|Mef2_5-100"))
# 
# testPositiveGeneset <- data.frame("ENSG0000060")
# colnames(testPositiveGeneset) <- c("EnsemblId")
# 
# 
# ####
####
###
####

#singleLineTableHumanMotifs <- fread("/home/dieter/Documents/Project/Programming/Transcription Factor Binding/HumanWG_motifs/analysisdata/combinedMotifFile.csv")



loadLibraryAndData <- function(reloadFromSource = FALSE, TenXCrossVal = TRUE, foldCrossVal = 10) {
  


  currentDir <- getwd()

  options(stringsAsFactors = FALSE)

  if(reloadFromSource == TRUE) {

    colNamesTFBSTables <- c("EnsemblId","TFBS")


    print("Choose the file that contains the ENCODE TFBS on one line")
    tFTableOneLine <- fread(tk_choose.files(default = "~/Documents/Project/Programming/Transcription Factor Binding/HumanWG_motifs/Data/singleLineFinalDf and Unified table/underscoresRemovedSingleLinePouya.csv",
                                                 caption = "Choose the .csv file that contains the ENCODE TFBS on one line",
                                                 multi = FALSE))
    colnames(tFTableOneLine) <- colNamesTFBSTables
    print(head(tFTableOneLine))
    print(paste0("rows original : ", nrow(tFTableOneLine)))
    
    #remove genes that have no motifs
    tFTableOneLine <- tFTableOneLine[tFTableOneLine$TFBS != "",]
    print(paste0("rows non-empty : ", nrow(tFTableOneLine)))

    tFTableMultipleLines <- fread(tk_choose.files(default = "~/Documents/Project/Programming/Transcription Factor Binding/HumanWG_motifs/Data/singleLineFinalDf and Unified table/underscoresRemovedMultipleLinesPouya.csv",
                                                       caption = "Choose the .csv file that contains the ENCODE TFBS on multiple lines",
                                                       multi = FALSE))
    colnames(tFTableMultipleLines) <- c(colNamesTFBSTables, "Count")
    print(head(tFTableMultipleLines))
    #remove genes without motifs/NA motifs
    print(paste0("rows original : ", nrow(tFTableMultipleLines)))
    tFTableMultipleLines <- tFTableMultipleLines[tFTableMultipleLines$Count > 0,]
    print(paste0("rows non-empty : ", nrow(tFTableMultipleLines)))
    
    
    MHCGeneTable <- read.csv(tk_choose.files(default = "~/Documents/Project/Programming/DataForAll/PositiveList_2.csv",
                                             caption = "Choose the .csv file that contains the positive MHC gene set",
                                             multi = FALSE),header = TRUE)

    names(MHCGeneTable)[7] <- "EnsemblId"
    print(head(MHCGeneTable))

    #save as .R data files

    saveRDS(tFTableOneLine, file = paste(currentDir, "/tFTableOneLinePouya.rds", sep=""))
    saveRDS(tFTableMultipleLines, file = paste(currentDir, "/tFTableMultipleLinesPouya.rds", sep=""))
    saveRDS(MHCGeneTable, file = paste(currentDir, "/MHCGeneTable.rds", sep=""))

  } else if (reloadFromSource == FALSE) {

    tFTableOneLine <- readRDS(file = paste(currentDir, "/tFTableOneLinePouya.rds", sep=""))
    tFTableMultipleLines <- readRDS(file = paste(currentDir, "/tFTableMultipleLinesPouya.rds", sep=""))
    MHCGeneTable <- readRDS(file = paste(currentDir, "/MHCGeneTable.rds", sep=""))

  }
  

  
  
  
  
  ########################
  ##       TEST         ##
  ########################
  
  #amountPositiveSetInList <- sum(testPositiveGeneset$EnsemblId %in% testTotalGenesMultipleLines$EnsemblId)
  #positiveSetTFBSDf <- testTotalGenesMultipleLines[testTotalGenesMultipleLines$EnsemblId %in% testPositiveGeneset$EnsemblId,]
  
  
  
  
  
  
  
  amountPositiveSetInList <- sum(MHCGeneTable$EnsemblId %in% tFTableMultipleLines$EnsemblId)
  positiveSetTFBSDf <- tFTableMultipleLines[tFTableMultipleLines$EnsemblId %in% MHCGeneTable$EnsemblId,]
  negativeSetTFBSDf <- tFTableMultipleLines[tFTableMultipleLines$EnsemblId %ni% MHCGeneTable$EnsemblId,]
  totalSetTFBSDf    <- tFTableMultipleLines
  
  print(paste("Amount of positive set genes in data set:", amountPositiveSetInList))
  print(head(positiveSetTFBSDf))
  #note: positiveSetTFBSDf thus contains the TFBS present in MHC pathway related genes
  
  
  if(TenXCrossVal == TRUE) {
    
    print("Generating sets for tenfold cross-validation, hold on...")
    #now generate a list of 10 sets to leave out per time
    ensemblIdsToLeaveOutListPositive <- vector("list", length = foldCrossVal)
    ensemblIdsToLeaveOutListNegative <- vector("list", length = foldCrossVal)
    #list to hold the 1/10th of the positive set that is left out each time
    positiveSetNotTrainedOn <- vector("list", length = foldCrossVal)
    positiveSetTrainedOn    <- vector("list", length = foldCrossVal)
    #negative sets
    negativeSetNotTrainedOn <- vector("list", length = foldCrossVal)
    negativeSetTrainedOn    <- vector("list", length = foldCrossVal)
    #total Set
    totalSetTrainedOn       <- vector("list", length = foldCrossVal)
    totalSetNotTrainedOn    <- vector("list", length = foldCrossVal)
    #total set subsets of whole genome data. Used in function scoresOverRows (crossval implementation of it)
    totalSetNotTrainedOnWGDOneLine <- vector("list", length = foldCrossVal)
    
    #random volgorde van de EnsemblIds in de positive set
    shuffledEnsemblIdsPositive <- sample(unique(positiveSetTFBSDf$EnsemblId))
    shuffledEnsemblIdsNegative <- sample(unique(negativeSetTFBSDf$EnsemblId))

    counter <- 1
    
     for(i in 1:length(shuffledEnsemblIdsPositive)) {
       
       ensemblIdsToLeaveOutListPositive[[counter]] <- c(ensemblIdsToLeaveOutListPositive[[counter]], shuffledEnsemblIdsPositive[[i]])
       counter <- counter + 1
       if(counter > foldCrossVal) {
         counter = 1
       }
      
     }
    
    counter <- 1
    
    for(i in 1:length(shuffledEnsemblIdsNegative)) {
      
      ensemblIdsToLeaveOutListNegative[[counter]] <- c(ensemblIdsToLeaveOutListNegative[[counter]], shuffledEnsemblIdsNegative[[i]])
      counter <- counter + 1
      if(counter > foldCrossVal) {
        counter = 1
      }
      
    }
    
    
    #get the ten sets to keep for testing post-training and the ten sets to train on pre-testing, both negative and positive, and combined
    for(i in 1:length(positiveSetNotTrainedOn)) {
      
    positiveSetNotTrainedOn[[i]]        <- positiveSetTFBSDf[positiveSetTFBSDf$EnsemblId %in% ensemblIdsToLeaveOutListPositive[[i]], ]
    positiveSetTrainedOn[[i]]           <- positiveSetTFBSDf[positiveSetTFBSDf$EnsemblId %ni% ensemblIdsToLeaveOutListPositive[[i]], ]
    negativeSetNotTrainedOn[[i]]        <- negativeSetTFBSDf[negativeSetTFBSDf$EnsemblId %in% ensemblIdsToLeaveOutListNegative[[i]], ]
    negativeSetTrainedOn[[i]]           <- negativeSetTFBSDf[negativeSetTFBSDf$EnsemblId %ni% ensemblIdsToLeaveOutListNegative[[i]], ]
    totalSetTrainedOn[[i]]              <- rbind(positiveSetTrainedOn[[i]], negativeSetTrainedOn[[i]])
    totalSetNotTrainedOn[[i]]           <- rbind(positiveSetNotTrainedOn[[i]], negativeSetNotTrainedOn[[i]])
    totalSetNotTrainedOnWGDOneLine[[i]] <- tFTableOneLine[tFTableOneLine$EnsemblId %in% totalSetNotTrainedOn[[i]]$EnsemblId, ]
    }
  }
  
  
  #count the totals for the positive set
  totalCountsMHCSet <- sum(positiveSetTFBSDf$Counts)
  
  #aggregate the amounts by TFBS for the positive set
  perMotifCountsMHCSet <- with(positiveSetTFBSDf,
                               aggregate(positiveSetTFBSDf,
                                         by = list(TFBS),
                                         FUN = function(x) { 
                                           if (is.double(x) | is.numeric(x)) {
                                             return(sum(x))
                                           }
                                           else {
                                             return(0)
                                           }}))
  perMotifCountsMHCSet <- perMotifCountsMHCSet[,c(1,4)]
  colnames(perMotifCountsMHCSet)[1] <- "TFBS"
  head(perMotifCountsMHCSet)
  
  #put that in a function for use in the cross-validation data
  generateMHCMotifCounts <- function(positiveSetTFBSDf) {
    
    perMotifCountsMHCSet <- with(positiveSetTFBSDf,
                                 aggregate(positiveSetTFBSDf,
                                           by = list(TFBS),
                                           FUN = function(x) { 
                                             if (is.double(x) | is.numeric(x)) {
                                               return(sum(x))
                                             }
                                             else {
                                               return(0)
                                             }}))
    perMotifCountsMHCSet <- perMotifCountsMHCSet[,c(1,4)]
    colnames(perMotifCountsMHCSet)[1] <- "TFBS"
    perMotifCountsMHCSet
    
  }
  
  #aggregate the motif amounts by TFBS for all the genes
  totalSetTFBSCounts <- with(tFTableMultipleLines,
                             aggregate(tFTableMultipleLines,
                                       by = list(TFBS),
                                       FUN = function(x) { 
                                         if (is.double(x) | is.numeric(x)) {
                                           return(sum(x))
                                           
                                         }
                                         else {
                                           return(0)
                                         }}))
  totalSetTFBSCounts <- totalSetTFBSCounts[,c(1,4)]
  colnames(totalSetTFBSCounts)[1] <- "TFBS"
  head(totalSetTFBSCounts)
  
  generateTotalSetMotifCounts <- function(totalSetPositiveAndNegativeTFBSDf) {
    
    totalSetTFBSCounts <- with(totalSetPositiveAndNegativeTFBSDf,
                               aggregate(totalSetPositiveAndNegativeTFBSDf,
                                         by = list(TFBS),
                                         FUN = function(x) { 
                                           if (is.double(x) | is.numeric(x)) {
                                             return(sum(x))
                                             
                                           }
                                           else {
                                             return(0)
                                           }}))
    totalSetTFBSCounts <- totalSetTFBSCounts[,c(1,4)]
    colnames(totalSetTFBSCounts)[1] <- "TFBS"
    head(totalSetTFBSCounts)
    totalSetTFBSCounts
    
  }
 
  
  
  
  totalCountsFullSet = sum(totalSetTFBSCounts$Counts)
  
  if(TenXCrossVal == TRUE) {
    
    print("Generating count data for tenfold cross-validation, hold on...")
    PerMotifCountsListPositive <- vector("list", length = foldCrossVal)
    PerMotifCountsListTotal    <- vector("list", length = foldCrossVal)
    
    for(i in 1:length(PerMotifCountsListPositive)) {
      
      PerMotifCountsListPositive[[i]] <- generateMHCMotifCounts(positiveSetTFBSDf = positiveSetTrainedOn[[i]])
      PerMotifCountsListTotal[[i]]    <- generateTotalSetMotifCounts(totalSetPositiveAndNegativeTFBSDf = totalSetTrainedOn[[i]])
      
    }
  }
  
  
  
  #calculate the odds that a random gene in the full set has a certain TFBS
  
  genomeOdds <- apply(totalSetTFBSCounts, FUN = function(x) {
    as.numeric(x[2])/totalCountsFullSet
  }, MARGIN = 1)
  names(genomeOdds) <- totalSetTFBSCounts$TFBS
  
  #calculate the odds that a gene in the MHC set has a certain TFBS
  
  MHCOdds <- apply(perMotifCountsMHCSet, FUN = function(x) {
    as.numeric(x[2])/totalCountsMHCSet
  }, MARGIN = 1)
  names(MHCOdds) <- perMotifCountsMHCSet$TFBS
  
  
  if(TenXCrossVal == TRUE) {
  list(                   tFTableOneLine                = tFTableOneLine,
                          tFTableMultipleLines          = tFTableMultipleLines,
                          MHCGeneTable                  = MHCGeneTable,
                          MHCGenesInDataset             = amountPositiveSetInList,
                          genomeOdds                    = genomeOdds,
                          MHCOdds                       = MHCOdds,
                          positiveSetTFBSCounts         = perMotifCountsMHCSet,
                          totalSetTFBSCounts            = totalSetTFBSCounts,
                          TFBSinMHCGenes                = positiveSetTFBSDf,
                          cvPositiveSetNotTrainedOnList = positiveSetNotTrainedOn,
                          cvPositiveSetTrainedOnList    = positiveSetTrainedOn,
                          cvNegativeSetTrainedOnList    = negativeSetTrainedOn,
                          cvNegativeSetNotTrainedOnList = negativeSetNotTrainedOn,
                          cvTotalSetTrainedOnList       = totalSetTrainedOn,
                          cvTotalSetNotTrainedOnList    = totalSetNotTrainedOn,
                          cvTotalSetNotTrainedOnWGDOL   = totalSetNotTrainedOnWGDOneLine,
                          cvPositiveSetTFBSCountsList   = PerMotifCountsListPositive,
                          cvTotalSetTFBSCountsList      = PerMotifCountsListTotal
  )
  } else {
    
    list(tFTableOneLine = tFTableOneLine,
         tFTableMultipleLines = tFTableMultipleLines,
         MHCGeneTable = MHCGeneTable,
         MHCGenesInDataset = amountPositiveSetInList,
         genomeOdds = genomeOdds,
         MHCOdds = MHCOdds,
         positiveSetTFBSCounts = perMotifCountsMHCSet,
         totalSetTFBSCounts = totalSetTFBSCounts,
         TFBSinMHCGenes = positiveSetTFBSDf)
    
  }
  
  }
  
  
  
  ##########################################################
  #computerFishersTest
  #function takes whole genome count data and positive set count data,
  #executes fisher's exact p on all relevant TFSB to find if they are enriched
  #in positive set
  ##########################################################
  
  computeFishersTest = function(wholeGenomeTFBSCountData, positiveSetTFBSCountData, run = FALSE, outputMatrices = TRUE) {
    
    #don't compute again unless explicitly told to. Otherwise, just open the .rds
    if(run == FALSE) {
      
      fisherDataFrame <- readRDS(file = paste(getwd(), "/humanMotifsFisherDataFramePouya.rds", sep = ""))
      return(fisherDataFrame)
      
    }
    else {
      
      print(head(wholeGenomeTFBSCountData))
      print(head(positiveSetTFBSCountData))
      #Get only those genes from the whole genome that have TFBS that are also found in the positive MHC set
      relevantWholeGenome    <- wholeGenomeTFBSCountData[wholeGenomeTFBSCountData$TFBS %in% positiveSetTFBSCountData$TFBS,] 
      totalCountsWholeGenome <- sum(wholeGenomeTFBSCountData$Count)
      totalCountsMHCGenes    <- sum(positiveSetTFBSCountData$Count)
      #print(paste("TotalCountsWholeGenome:", totalCountsWholeGenome))
      #print(paste("TotalCountsMHCGenes:", totalCountsMHCGenes))
      print(head(relevantWholeGenome))
      
      #listMatrices is for debug purposes, to check whether correct contingency tables are created
      listMatrices <- list()
      fisherDataFrame <- data.frame(TFBS = character(),
                                    pValue=numeric(),
                                    sampleEstimate=numeric(),
                                    lowerConfidenceInterval=numeric(),
                                    upperConfidenceInterval=numeric(),
                                    WholeGenomeCount = numeric(),
                                    MHCCount = numeric())
      
      #cycle over all the rows of 
      for(bindingSites in 1:nrow(positiveSetTFBSCountData)) {
        
        ##debugging stuff
        #print(totalCountsMHCGenes)
        #print(positiveSetTFBSCountData$Counts[bindingSites])
        
        ##Note: since I often forget exactly how this works I am writing in here what I am doing:
        ##You have the amount of this specific TFBS in the positive set.
        ##You have the amount of all other TFBS in the positive set (so, total - this specific one)
        ##You have the amount of this specific TFBS in the whole genome
        ##You have the amount of all other TFBS in the whole genome (so, whole genome total - this specific one in the whole genome data)
        ##What I'm worried about: is this in the right order (i.e. is relevantWholeGenome automatically in the right order if you ask %in%
        ##positiveSetCountData$TFBS. apparently, it is)
        fisherMatrix <- matrix(c(positiveSetTFBSCountData$Count[bindingSites],
                                 totalCountsMHCGenes - positiveSetTFBSCountData$Count[bindingSites],
                                 relevantWholeGenome$Count[bindingSites],
                                 totalCountsWholeGenome - relevantWholeGenome$Count[bindingSites]),
                               nrow = 2,
                               dimnames = list(TFBS = c(relevantWholeGenome$TFBS[bindingSites], "TotalTFBS"),
                                               Set = c("Positive set", "Whole genome")))
        #print(fisherMatrix)
        fisherTestResult <- fisher.test(fisherMatrix)
        pValue           <- fisherTestResult["p.value"]
        lowerConfInt     <- unlist(fisherTestResult["conf.int"])[1]
        upperConfInt     <- unlist(fisherTestResult["conf.int"])[2]
        sampleEst        <- fisherTestResult["estimate"]
        WholeGenomePresence <- fisherMatrix[1,2]
        MHCPresence <- fisherMatrix[1,1]
        
        listMatrices     <- list(listMatrices,list(fisherMatrix))
        addVector        <- c(as.character(relevantWholeGenome$TFBS[bindingSites]),
                              pValue,
                              sampleEst,
                              lowerConfInt,
                              upperConfInt,
                              WholeGenomePresence,
                              MHCPresence)
        names(addVector) <- c("TFBS", "pValue", "sampleEstimate",
                              "lowerConfidenceInterval", "upperConfidenceInterval",
                              "WholeGenomeCount", "MHCCount")
        #add to the Df the results of the fisherTest
        fisherDataFrame  <- rbind(fisherDataFrame, addVector, make.row.names = FALSE)
        
      }
      #adjust for multiple testing with Benjamini-Hochberg methodology (utilises FDR)
      fisherDataFrame$qValue <- p.adjust(fisherDataFrame$pValue, method="BH")
      saveRDS(fisherDataFrame, file = paste(getwd(), "/humanMotifsFisherDataFramePouya.rds", sep = ""))
      if(outputMatrices == TRUE) {
        return(list(matrices = listMatrices, fisherDataFrame = fisherDataFrame))
      } else if (outputMatrices == FALSE) {
        return(fisherDataFrame)
      }
    }
    
    
  }

  crossValFisherTests <- function(crossValTotalSetTFBSCounts, crossValPositiveTFBSCountList) {
  
  result <- vector("list", length(crossValPositiveTFBSCountList))
  for(i in 1:length(result)) {
    
    result[[i]] <- computeFishersTest(crossValTotalSetTFBSCounts[[i]],
                                      crossValPositiveTFBSCountList[[i]],
                                      run = TRUE,
                                      outputMatrices = FALSE)
  }
  
  
  result 
}

  
  ##########################################################
  #function scoresOverRows
  #takes the dataframe with Fisher tests and sample estimates as well as whole genome data (one line).
  #returns dataframe with score per TFBS and total score of all genes in the Pouya Kheradpout data.
  ##########################################################
  
  scoresOverRows <- function(wholeGenomeDataOneLine, fisherTests, PositiveSetGenesInData,
                             scoreCorrectionMethod = 3) {
  
  #scoreCorrectionMethods:
    #1: uncorrected for the amount of times the motif is present. That is to say:
    # a gene gets the score for presence of a motif whether it is present just once
    #or present 5 times.
    
    #2: corrected for the #times present in a target gene and the average presence in MHC genes
    #by assigning the score: log2(sampleEstimate) * countInTargetGene/AverageCountPerMHCGene
    #In this way, for example for NFkB, an 'average' MHC gene has ~3.2 NFkB. The full score
    #for a gene for this TFBS is only obtained when 3.2 NFkB's are present. If there are more
    #the score is increased somewhat. THis is biologically motivated in that more Motifs
    #in theory increased chances that something actually binds.
    #the problem is that there are also TFBS underrepresented in MHC. When dividing there
    #any gene that has even one of these TFBS will get very higly penalised. To combat this, 
    #see correction 3.
    
    #3: Same as correction 2, but given that we are more interested in TFBS that are overrepresented
    #in MHC genes rather than those that are underrepresented, the correction for the counts
    #is only applied to MHC genes that are overrepresented in MHC. i.e. if a TFBS is only present 4
    #times in all MHC genes, but much more often in whole genome, it might be that each MHC gene
    #on average has only 0.2 of that TFBS. Let's say the sampleEstimate of that TFBS is 0.4 for MHC.
    #If another gene in the whole genome now has 20 IRF motifs (present on average 10 times per MHC gene)
    #but two times that depleted TFBS, its score will be very low --> log2(0.4) * (2/0.2) = -26.4, while the
    #IRF only yields log2(1.42) * 20/10 = 1. That is an excessive disparity between the weights on depleted and enriched TFBS,
    #so in this option, no correction for average MHC amount is made for MHC-depleted TFBS.
    
    
    
  #DEBUGGINF    DEBUGGING         DEBUGGING
  
  #accountForTFBSCount             <- TRUE
  #wholeGenomeDataOneLine          <- testTotalGenesOneLine
  #amountTFBSMotifsForNominalScore <- 1
  
  #find TFBS that are significantly overrepresented in positive set genes
  #assign them a score. In this case, the point estimate from the Fisher tests.
  
  
  informativeTFBS            <- subset(fisherTests, qValue <= 0.10 & MHCCount > 3)
  informativeTFBS            <- informativeTFBS[order(informativeTFBS$qValue),]
  rownames(informativeTFBS)  <- seq(1:nrow(informativeTFBS))
  informativeTFBS$score      <- 0
  informativeTFBS$MHCAverage <- informativeTFBS$MHCCount/PositiveSetGenesInData
  
  informativeTFBS[informativeTFBS$sampleEstimate <1, ]$score = log2(informativeTFBS[informativeTFBS$sampleEstimate <1, ]$sampleEstimate)
  informativeTFBS[informativeTFBS$sampleEstimate >=1, ]$score = log2(informativeTFBS[informativeTFBS$sampleEstimate >=1, ]$sampleEstimate) 
  print(informativeTFBS[informativeTFBS$sampleEstimate <1, ])
  print( informativeTFBS[informativeTFBS$sampleEstimate >=1, ])
  
  #print(head(informativeTFBS))
  
  
  scoreCompendium  <- vector()
  dataFrameResults <- data.frame()
  #cycle over rows. Could be refined with an apply function
  for(rows in 1:nrow(wholeGenomeDataOneLine)) {
    
    #print(rows)
    
    #check whether any of the informative TFBS are present in a gene's promoter region
    currentRow                    <- wholeGenomeDataOneLine[rows,]
    #print(currentRow)
    currentRowTFs                 <- currentRow[1,2]
    #print(currentRowTFs)
    TFBSCountsAndPresenceThisGene <- str_match_all(currentRowTFs,"([^|]+)-(\\d+)")
    #print(TFBSCountsAndPresenceThisGene)
    #Sys.sleep(5)
    characterVectorTFs            <- TFBSCountsAndPresenceThisGene[[1]][,2]
    countVectorTFBS               <- as.numeric(TFBSCountsAndPresenceThisGene[[1]][,3])
    names(countVectorTFBS)        <- characterVectorTFs
    logicalPresence               <- informativeTFBS$TFBS %in% characterVectorTFs
    currentRowEnsemblId           <- currentRow[1,1]
    #print(currentRow)
    #print(currentRowEnsemblId)
    
    for(bindingSites in 1:length(logicalPresence)) {
      #for all informative TFBS, if they are present, add their score, otherwise, add 0
      
      
      if(logicalPresence[bindingSites] == TRUE) {
        
        #print("In the true part of logicalpresence")
        
        
        if(scoreCorrectionMethod == 1) {
          
          valueToAdd <- informativeTFBS$score[bindingSites] 

        } else if(scoreCorrectionMethod == 2) {
          
          TFBSCount  <- countVectorTFBS[informativeTFBS$TFBS[bindingSites]]
          valueToAdd <- informativeTFBS$score[bindingSites] * (TFBSCount / informativeTFBS$MHCAverage[bindingSites])
          
        } else if(scoreCorrectionMethod == 3) {
          
          if(informativeTFBS$sampleEstimate[bindingSites] < 1) {
            
          valueToAdd <- informativeTFBS$score[bindingSites]
          } else {
            
            TFBSCount  <- countVectorTFBS[informativeTFBS$TFBS[bindingSites]]
            valueToAdd <- informativeTFBS$score[bindingSites] * (TFBSCount / informativeTFBS$MHCAverage[bindingSites])
          }
  
        }
        
        
      } else if (logicalPresence[bindingSites] == FALSE) {
        
        valueToAdd <- 0
        
      }
      
      scoreCompendium <- append(scoreCompendium,valueToAdd)
      names(scoreCompendium)[bindingSites] <-  informativeTFBS$TFBS[bindingSites]
      
    }
    #set rownames as the EnsemblId of the gene from the whole genome
    dataFrameResults <- rbind(dataFrameResults, scoreCompendium)
    rownames(dataFrameResults) <- c(head(rownames(dataFrameResults),-1),currentRow[1,1])
    scoreCompendium <- vector()
    #print(nrow(dataFrameResults))
    
  }
  colnames(dataFrameResults) <- paste("Score for TFBS", informativeTFBS$TFBS)
  dataFrameResults$totalScore <- rowSums(dataFrameResults)
  dataFrameResults <- dataFrameResults[order(dataFrameResults$totalScore, decreasing = TRUE),]
  print("These are the informative TFBS:")
  print(informativeTFBS$TFBS)
  
  #DEBUGGINF    DEBUGGING         DEBUGGING #DEBUGGINF    DEBUGGING         DEBUGGING #DEBUGGINF    DEBUGGING         DEBUGGING

  
     return(list(
       scoresWholeGenome = dataFrameResults,
       informativeTFBS = informativeTFBS))
   }
  
  crossValScoresOverRows <- function(crossValWholeGenomeDataOneLineList, crossValFisherTestList, crossValPositiveSets) {
    
    
    
    result <- vector("list", length = length(crossValFisherTestList))
    
    for(i in 1: length(result)) {
      
      positiveSetGenesInData <- length(unique(crossValPositiveSets[[i]]$EnsemblId))
      result[[i]] <- scoresOverRows(crossValWholeGenomeDataOneLineList[[i]],
                                    crossValFisherTestList[[i]],
                                    positiveSetGenesInData,
                                    scoreCorrectionMethod = 3)

      
    }
    
    result
   
    
  }
  
  ##########################################################
  #actual execution of the programme
  #loads positive set, checks which TFBS are present there and how often
  #does the same for whole genome data (WGD)
  #determines which TFBS overrepresented in MHC via fisher tests
  #assigns each gene a score log(sample odds ratio estimate from Fisher's exact p/corrected p-value Fisher's exact p) for every TFBS present
  #concatenates these individual scores into a total gene score and saves that to a file.
  ##########################################################
  dataStructureTFBS <- loadLibraryAndData(reloadFromSource = TRUE)
  str(dataStructureTFBS)
  dataStructureTFBS$tFTableOneLine
  fisherTests <- computeFishersTest(dataStructureTFBS$perMotifCountsFullSet,dataStructureTFBS$perMotifCountsPositiveSet, run = TRUE, outputMatrices = FALSE)
  crossValFisherTestData <- crossValFisherTests(dataStructureTFBS$cvTotalSetTFBSCountsList, dataStructureTFBS$cvPositiveSetTFBSCountsList)
  
  str(fisherTests)
  informativeTFBS <- subset(fisherTests, qValue <= 0.10 & MHCCount > 3)
  informativeTFBS
  
  
  #switch for execution of the selection of all genes in the WGD that have the TFBS sig. overrepresented in positive set genes
  #set to TRUE when the positive set has changed and/or the analysis should be run again.
  #toggle which file to open (specific motifs versus more general motifs) below
  
  runselection = TRUE
  
  if (runselection == TRUE) {
    
    finalScoresWholeGenomePouyaCorrectionThree <- scoresOverRows(dataStructureTFBS$tFTableOneLine,fisherTests, dataStructureTFBS$MHCGenesInDataset, scoreCorrectionMethod = 3)
    saveRDS(finalScoresWholeGenomePouyaCorrectionThree, file = "/home/dieter/Documents/Project/Programming/Transcription Factor Binding/HumanWG_motifs/Data/singleLineFinalDf and Unified table/GeneralisedPouyaWholeGenomeScores.rds")
     }
  if (runselection == FALSE) {
    finalScoresWholeGenomePouyaCorrectionThree <- readRDS(file = "/home/dieter/Documents/Project/Programming/Transcription Factor Binding/HumanWG_motifs/Data/singleLineFinalDf and Unified table/GeneralisedPouyaWholeGenomeScores.rds")
      }
  
  str(finalScoresWholeGenomePouyaCorrectionThree)
  head(finalScoresWholeGenomePouyaCorrectionThree$scoresWholeGenome, 15)
  tail(finalScoresWholeGenomePouyaCorrectionThree$scoresWholeGenome, 15)
  
  
  #TenXCrossVal
  crossValScoresWholeGenome <- crossValScoresOverRows(crossValWholeGenomeDataOneLineList = dataStructureTFBS$cvTotalSetNotTrainedOnWGDOL,
                                                              crossValFisherTestList = crossValFisherTestData,
                                                              crossValPositiveSets = dataStructureTFBS$cvPositiveSetTrainedOnList)
  
  
  #get the crossval final results
  finalResults <- data.frame(totalScore = numeric())
  for (i in 1:length(crossValScoresWholeGenome)) {
    p <- as.data.frame(crossValScoresWholeGenome[[i]]$scoresWholeGenome$totalScore)
    rownames(p) <- rownames(crossValScoresWholeGenome[[i]]$scoresWholeGenome)
    finalResults <- rbind(finalResults, p)
  }
  
  
  
  #make some boxplots of the amount of TFBS each MHC gene has
  nrow(dataStructureTFBS$MHCGeneTable)
  dataStructureTFBS$MHCGenesInDataset
  plotDataBoxplot <- data.frame(TFBS = character(), Count = numeric())
  for(i in 1:nrow(informativeTFBS)) {
  TFBSToAdd <- dataStructureTFBS$TFBSinMHCGenes[dataStructureTFBS$TFBSinMHCGenes$TFBS == informativeTFBS$TFBS[i]][,c(2,3)]
  plotDataBoxplot <- rbind(plotDataBoxplot, TFBSToAdd, use.names = TRUE)
  }
  ggplot(data = plotDataBoxplot, aes(x = TFBS, y = Count)) + geom_boxplot() + theme_bw() +
    theme(axis.text.x=element_text(angle=60,hjust=1)) 
  
  ##changed so that zeroes count:
  plotDataBoxplotZeroesIncluded <- data.frame(TFBS = character(), Count = numeric())
  for(i in 1:nrow(informativeTFBS)) {
    TFBSToAdd <- dataStructureTFBS$TFBSinMHCGenes[dataStructureTFBS$TFBSinMHCGenes$TFBS == informativeTFBS$TFBS[i]][,c(2,3)]
    amountOfZeroes <- dataStructureTFBS$MHCGenesInDataset - nrow(TFBSToAdd)
    TFBScol <- rep(as.character(TFBSToAdd[1,1]), amountOfZeroes)
    Countcol <- rep(0, amountOfZeroes)
    zeroDataFrame <- data.frame(TFBS = TFBScol, Count = Countcol)
    TFBSToAdd <- rbind(TFBSToAdd, zeroDataFrame, use.names = TRUE)
    print(head(TFBSToAdd))
    print(tail(TFBSToAdd))
    plotDataBoxplotZeroesIncluded <- rbind(plotDataBoxplotZeroesIncluded, TFBSToAdd, use.names = TRUE)
  }
  head(plotDataBoxplotZeroesIncluded, 60)
  ggplot(data = plotDataBoxplotZeroesIncluded, aes(x = TFBS, y = Count)) + geom_boxplot() + theme_bw() +
    theme(axis.text.x=element_text(angle=60,hjust=1)) 
  
  #how much are there in MHC, and what are the means?
  aggregate(Count ~ TFBS, plotDataBoxplotZeroesIncluded, FUN = function(x) {k <- cbind(sum(x), mean(x))
  names(k) <- c("Amount", "MeanperGene")
  k})
  
  
  #Save the final scores of each gene with respect to the amount of MHC-associated TFBS it has
  
  ###################################################################################################
  #Two functions follow below. One takes the frequency data in the whole genome and positive set and 
  #creates a ggplot-compatible data frame. The other makes the actual plot.
  ###################################################################################################
  
  #step 1: create a ggplot compatible data frame:
  createGraphData <- function(TFBSdata = dataStructureTFBS, run = TRUE ) {
    
    infTFBS <- finalScoresWholeGenomePouyaCorrectionThree$informativeTFBS[finalScoresWholeGenomePouyaCorrectionThree$informativeTFBS$qValue <= 0.1 & finalScoresWholeGenomePouyaCorrectionThree$informativeTFBS$MHCCount > 3, ]$TFBS
    resultsDf <- data.frame(TF = character(), yesno = character(), genomeMHC = character(), proportion = numeric(), qValue = numeric(), tag = character())
    for (significantTFBS in 1: length(infTFBS))
    {
      
      
      totalTFBSInMHCGenes    <- sum(TFBSdata$TFBSinMHCGenes$Count)
      specificTFBSInMHCGenes <- TFBSdata$TFBSinMHCGenes$TFBS %in% infTFBS[significantTFBS]
      totalSpecificTFBSInMHC <- sum(TFBSdata$TFBSinMHCGenes[specificTFBSInMHCGenes,]$Count)
      
      totalTFBSInWGDGenes    <- sum(TFBSdata$totalSetTFBSCounts$Count)
      specificTFBSInWGDGenes <- TFBSdata$totalSetTFBSCounts$TFBS %in% infTFBS[significantTFBS]
      totalSpecificTFBSInWGD <- sum(TFBSdata$totalSetTFBSCounts[specificTFBSInWGDGenes,]$Count)
      
      
      proportionInMHC        <- totalSpecificTFBSInMHC/totalTFBSInMHCGenes
      proportionNotInMHC     <- 1 - proportionInMHC
      proportionInWGD        <- totalSpecificTFBSInWGD/totalTFBSInWGDGenes
      proportionNotInWGD     <- 1- totalSpecificTFBSInWGD/totalTFBSInWGDGenes 
      tempDf                 <- data.frame(TF = rep(infTFBS[significantTFBS],4),
                                           yesno = c("yes", "no", "yes", "no"),
                                           genomeMHC = c("MHC", "MHC", "genome", "genome"),
                                           proportion = c(proportionInMHC, proportionNotInMHC, proportionInWGD, proportionNotInWGD),
                                           qValue = rep(finalScoresWholeGenomePouyaCorrectionThree$informativeTFBS[finalScoresWholeGenomePouyaCorrectionThree$informativeTFBS == infTFBS[significantTFBS],]$qValue,4))
      #adds in three marks for significance for drawing
      if (unique(tempDf$qValue) < 0.01) {
        tempDf = data.frame(tempDf, data.frame(tag = rep("***",4)))
      } else if (unique(tempDf$qValue) < 0.05){
        tempDf = data.frame(tempDf, data.frame(tag = rep("**",4)))
      } else {
        tempDf = data.frame(tempDf, data.frame(tag = rep("*",4)))
      }
      resultsDf = rbind(resultsDf, tempDf)
      
    }
    return(resultsDf)
  }
  
  
  
  yesNoPlotData <- createGraphData(dataStructureTFBS, run = TRUE)
  #step 2:
  #Now plot this data.
  plotGraphData <- function(plottingData, run = TRUE) {
    
    if(run == TRUE) {
      green <- rgb(0, 140, 71, maxColorValue = 255)
      red   <- rgb(237, 45, 46, maxColorValue = 255)
      initialPlot <-ggplot(data = plottingData,
                           aes(x= factor(yesno, levels = c("yes","no")),
                               y = proportion, fill = factor(genomeMHC, levels = c ("MHC", "genome")))) + 
        facet_wrap(~TF, scales = "free") +
        geom_bar(stat = "identity", position = position_dodge(), colour = "black") +
        theme_minimal() +
        xlab ("TFBS present") +
        scale_fill_manual(breaks = c("MHC","genome"),
                          values = c(green, red)) +
        guides(fill = guide_legend(title = NULL)) +
        scale_x_discrete(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0)) +
        coord_cartesian(ylim = c(0,1))
      #Add significance values, lines, and asterisks indicating sig. qvalues (*** <0.01, ** <0.05, * < 0.1)
      annotatedPlot <- initialPlot +
        geom_text(aes(x= 1, y = 0.9, label = paste("qValue =\n",round(qValue,5))), size = 2.5)  +
        geom_text(aes(x= 1, y = 0.77, label = tag))                                             +
        geom_segment(aes(x = 0.75, y = 0.75, xend = 1.25, yend = 0.75)) 
      annotatedPlot
      ggsave(plot = annotatedPlot, filename = "/home/dieter/Documents/Project/Programming/Transcription Factor Binding/HumanWG_motifs/Data/singleLineFinalDf and Unified table/TFBSOverrepresentationPlotPouya_MoreGeneral.pdf", device = "pdf", dpi = 600)

      return(annotatedPlot)
    }
    
  }
  
  plotGraphData(yesNoPlotData, run = TRUE)
  
  #Create a .csv that outputs, per TFBS, the genes from the positive set that have this TFBS and their characteristics
  #wat ik wil: weten welke positive genes een bepaalde TFBS hebben
  #methode: pak de EnsemblIds van deze positive genes
  #pak de relevante rijen uit de transcriptionfactorTablemultiplelines
  #split ze --> TFs per gen
  #kijk in welk van deze positive genes nu een bepaalde TFBS zit.
  #doe dit alleen als ik iets veranderd heb, laat voor de rest dat bestand ongemoeid. Dus manueel op run = TRUE zetten.
  outputMHCGeneTFBS <- function(run = FALSE) {
    
    if (run == TRUE) {
      Takegenes <- dataStructureTFBS$tFTableMultipleLines[dataStructureTFBS$tFTableMultipleLines$EnsemblId %in% dataStructureTFBS$MHCGeneTable$EnsemblId,]
      print(head(Takegenes))
      
      #now, for every relevant transcription factor, return the genes that have it
      genelist = list()
      for(relevantTFBS in 1:length(finalScoresWholeGenome$informativeTFBS$TFBS)) {
        
        positiveGenesWithTFBS = Takegenes[finalScoresWholeGenome$informativeTFBS$TFBS[relevantTFBS] == as.vector(Takegenes$TFBS),]
        print(positiveGenesWithTFBS)
        print(positiveGenesWithTFBS)
        print(class(Takegenes$TranscriptionFactorBindingSites))
        EnsemblIds = positiveGenesWithTFBS$EnsemblId
        genelist[[relevantTFBS]] <- MHCGeneTable[MHCGeneTable$EnsemblId  %in%  EnsemblIds,]
        names(genelist)[relevantTFBS] <- finalScoresWholeGenome$informativeTFBS$TFBS[relevantTFBS]
        genelist[[relevantTFBS]] <- cbind(TFBS = rep(finalScoresWholeGenome$informativeTFBS$TFBS[relevantTFBS],length(genelist[[relevantTFBS]]$MHC1OR2)),
                                          genelist[[relevantTFBS]] )
        if (relevantTFBS == 1) {
          
          write.csv(genelist[[relevantTFBS]], "/home/dieter/Documents/Project/Programming/Transcription Factor Binding/HumanWG_motifs/Data/singleLineFinalDf and Unified table/PositiveSetGenesAndAssociatedGeneralisedTFBSPouyaSet.csv")
         
        }
        
        else if (relevantTFBS > 1) {
          
          write.table(genelist[[relevantTFBS]], "/home/dieter/Documents/Project/Programming/Transcription Factor Binding/HumanWG_motifs/Data/singleLineFinalDf and Unified table/PositiveSetGenesAndAssociatedGeneralisedTFBSPouyaSet.csv",
                      append = TRUE, sep = ",", col.names = FALSE)
          
          
        }
      }
    }
    
  }
  #manual switch for running this, set to true if the positive set changes in some form.
  outputMHCGeneTFBS(run = TRUE)
  #outputMHCGeneTFBS(run = FALSE)
  ###################################################################
  #Output every gene in the WGD that has a sig. TFBS to a .csv file
  #method:
  #Take all the EnsemblIds from the whole genome data that have an informative transcription factor
  #Output them as a csv with the TFBS as header and the EnsemblIDs below
  #Function below, run is a switch determining whether to run it
  ###################################################################
  outputEnsIdsWithTFBS <- function(run = FALSE) {
    
    if(run == TRUE) {
      listIds <- list()
      for(relevantTFBS in 1:length(finalScoresWholeGenome$informativeTFBS$TFBS)) {
        currentTFBS <- finalScoresWholeGenome$informativeTFBS$TFBS[relevantTFBS]
        WGD_relevantgenes <-  dataStructureTFBS$tFTableMultipleLines[dataStructureTFBS$tFTableMultipleLines$TranscriptionFactorBindingSites %in%  currentTFBS,]
        
        #get the EnsemblIds for the current transcription factor from the data
        #count the amount that have it per TFBS
        Ids <- WGD_relevantgenes$EnsemblId
        listIds[[currentTFBS]] <- Ids
        lengthIds <- length(Ids)
        names(lengthIds) <- currentTFBS
        listIds[["amount"]] <- append(listIds[["amount"]], lengthIds)
        
      }
      
      #now construct a dataframe with this knowledge
      
      neededLength <- max(listIds$amount)
      
      #loop through the transcription factors, get their EnsIDs, and add "" as necessary for dataframe construction
      for(relevantTFBS in 1:length(finalScoresWholeGenome$informativeTFBS$TFBS)) {
        currentTFBS <- finalScoresWholeGenome$informativeTFBS$TFBS[relevantTFBS]
        
        #repeat "" required amount of times (i.e. the needed length minus the length it has)
        #append that to the EnsemblID list
        reqReps <-  neededLength - listIds$amount[currentTFBS]
        listIds[[currentTFBS]] <- c(listIds[[currentTFBS]], rep("", reqReps)) 
        
        if (exists("dfWGDEnsmblIdsPerTFBS")) {
          
          dfWGDEnsmblIdsPerTFBS <- data.frame(dfWGDEnsmblIdsPerTFBS,listIds[[currentTFBS]])
          
        } else {
          
          dfWGDEnsmblIdsPerTFBS <- data.frame(listIds[[currentTFBS]])
          
        }
      }
      
      colnames(dfWGDEnsmblIdsPerTFBS) <- names(listIds$amount)
      #Once that is done, put all the columns into a dataframe, and save that as .csv
      write.csv(dfWGDEnsmblIdsPerTFBS, "/home/dieter/Documents/Project/Programming/Transcription Factor Binding/HumanWG_motifs/Data/singleLineFinalDf and Unified table/EnsemblIdsWithPositiveSetEnrichedGeneralisedTFBSPouya.csv")
      #switch for different filename
      #write.csv(dfWGDEnsmblIdsPerTFBS, "/home/dieter/Documents/Project/Programming/Transcription Factor Binding/HumanWG_motifs/FinalFilesForEnrichmentAnalysis/EnsemblIdsWithThePositiveSetEnrichedSpecificTFBS.csv")
      
      
      #return the df for checking in R
      return (dfWGDEnsmblIdsPerTFBS)
    }
    
  }
  
  outputEnsIdsDataFrame <- outputEnsIdsWithTFBS(run = TRUE)
  
  
  
  
  ##################################################################################################
  # This piece of code creates histograms of the positive and negative set for the scores of all TFBS
  #
  ##################################################################################################
  
  
  makePlotsDistributions <- function (WholeGenomeScoreTable = 3, tenXCrossVal = FALSE) {
  
  if (tenXCrossVal == FALSE) {
  positiveSetIDs <- dataStructureTFBS$MHCGeneTable$EnsemblId
  inPositiveSet <- rownames(listFinalScoresPouya[[WholeGenomeScoreTable]]$scoresWholeGenome) %in% positiveSetIDs
  notInPositiveSet <- rownames(listFinalScoresPouya[[WholeGenomeScoreTable]]$scoresWholeGenome) %ni% positiveSetIDs
  finalPlotDfHist <- data.frame(scores = numeric(), set = character(), ScoreFor = character())
  for (i in 1:ncol(listFinalScoresPouya[[WholeGenomeScoreTable]]$scoresWholeGenome)) {
    finalScoresPositiveSetTotal <- listFinalScoresPouya[[WholeGenomeScoreTable]]$scoresWholeGenome[inPositiveSet, i]
    finalScoresNegativeSetTotal <- listFinalScoresPouya[[WholeGenomeScoreTable]]$scoresWholeGenome[notInPositiveSet, i]
    plotDfHistPositive <- data.frame(scores = finalScoresPositiveSetTotal, set = "positive")
    plotDfHistNegative <- data.frame(scores = finalScoresNegativeSetTotal, set = "negative")
    totalPlotDfHist <- rbind(plotDfHistNegative, plotDfHistPositive)
    totalPlotDfHist$ScoreFor <- colnames(listFinalScoresPouya[[WholeGenomeScoreTable]]$scoresWholeGenome)[i]
    finalPlotDfHist <- rbind(finalPlotDfHist, totalPlotDfHist)
  }
  } else if (tenXCrossVal == TRUE) {
    
    
    positiveSetIDs <- dataStructureTFBS$MHCGeneTable$EnsemblId
    inPositiveSet <- rownames(finalResults) %in% positiveSetIDs
    notInPositiveSet <- rownames(finalResults) %ni% positiveSetIDs
    finalPlotDfHist <- data.frame(scores = numeric(), set = character(), ScoreFor = character())
    for (i in 1:ncol(finalResults)) {
      finalScoresPositiveSetTotal <- finalResults[inPositiveSet, i]
      finalScoresNegativeSetTotal <- finalResults[notInPositiveSet, i]
      plotDfHistPositive <- data.frame(scores = finalScoresPositiveSetTotal, set = "positive")
      plotDfHistNegative <- data.frame(scores = finalScoresNegativeSetTotal, set = "negative")
      totalPlotDfHist <- rbind(plotDfHistNegative, plotDfHistPositive)
      totalPlotDfHist$ScoreFor <- "totalScore"
      finalPlotDfHist <- rbind(finalPlotDfHist, totalPlotDfHist)
    
    
    
    
    }
  }
  

  
    
kaas <- subset(finalPlotDfHist, ScoreFor == "totalScore")
chicken <- ggplot(data = kaas, aes(x = scores)) +
  geom_histogram(aes(y = 4 * ..density.., fill = set), color = "black", position = "dodge", stat = "bin", binwidth = 4) +
  scale_y_continuous(breaks=waiver(), expand=c(0,0), trans = "log1p") +
  theme_bw()   
chicken

pieter <- ggplot(data = kaas, aes(x = scores)) +
  geom_density(aes(y = ..density.., fill = set), color = "black",  alpha = 0.3) +
  scale_y_continuous(breaks=waiver()) +
  #scale_x_continuous(limits = c(-100, 100)) +
  theme_bw() #+ 
  #facet_wrap(~ScoreFor, scales = "free")
pieter

  list(chicken, pieter, finalPlotDfHist)
    

  }
  
  
  makePlotsDistributions(3)
  j <- makePlotsDistributions(3, tenXCrossVal = TRUE)
  j
  #########################Further analysis of TFBS data##########################
  
  
  setwd(oldwd)
  ###################################################
  #                       END                       #
  ###################################################
  
  