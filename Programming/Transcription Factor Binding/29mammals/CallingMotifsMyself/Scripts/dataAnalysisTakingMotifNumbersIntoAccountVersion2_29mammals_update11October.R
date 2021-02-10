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
setwd("~/Documents/Project/Programming/Transcription Factor Binding/29mammals/CallingMotifsMyself/Data/FinalDfAndUnifiedMotifFile/")
















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



loadLibraryAndData <- function(reloadFromSource = FALSE) {
  


  currentDir <- getwd()

  options(stringsAsFactors = FALSE)

  if(reloadFromSource == TRUE) {

    colNamesTFBSTables <- c("EnsemblId","TFBS")


    print("Choose the file that contains the ENCODE TFBS on one line")
    tFTableOneLine <- fread(tk_choose.files(default = "~/Documents/Project/Programming/Transcription Factor Binding/29mammals/CallingMotifsMyself/Data/FinalDfAndUnifiedMotifFile/GeneralisedMotifs/underscoresRemovedSingleLine29mammals.csv",
                                                 caption = "Choose the .csv file that contains the ENCODE TFBS on one line",
                                                 multi = FALSE))
    colnames(tFTableOneLine) <- colNamesTFBSTables
    print(head(tFTableOneLine))
    print(paste0("rows original : ", nrow(tFTableOneLine)))
    
    #remove genes that have no motifs
    tFTableOneLine <- tFTableOneLine[tFTableOneLine$TFBS != "",]
    print(paste0("rows non-empty : ", nrow(tFTableOneLine)))

    tFTableMultipleLines <- fread(tk_choose.files(default = "~/Documents/Project/Programming/Transcription Factor Binding/29mammals/CallingMotifsMyself/Data/FinalDfAndUnifiedMotifFile/GeneralisedMotifs/underscoresRemovedMultipleLines29mammals.csv",
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

    saveRDS(tFTableOneLine, file = paste(currentDir, "/tFTableOneLine29mammals.rds", sep=""))
    saveRDS(tFTableMultipleLines, file = paste(currentDir, "/tFTableMultipleLines29mammals.rds", sep=""))
    saveRDS(MHCGeneTable, file = paste(currentDir, "/MHCGeneTable.rds", sep=""))

  } else if (reloadFromSource == FALSE) {

    tFTableOneLine <- readRDS(file = paste(currentDir, "/tFTableOneLine29mammals.rds", sep=""))
    tFTableMultipleLines <- readRDS(file = paste(currentDir, "/tFTableMultipleLines29mammals.rds", sep=""))
    MHCGeneTable <- readRDS(file = paste(currentDir, "/MHCGeneTable.rds", sep=""))

  }
  

  
  
  
  
  ########################
  ##       TEST         ##
  ########################
  
  #amountPositiveSetInList <- sum(testPositiveGeneset$EnsemblId %in% testTotalGenesMultipleLines$EnsemblId)
  #positiveSetTFBSDf <- testTotalGenesMultipleLines[testTotalGenesMultipleLines$EnsemblId %in% testPositiveGeneset$EnsemblId,]
  
  amountPositiveSetInList <- sum(MHCGeneTable$EnsemblId %in% tFTableMultipleLines$EnsemblId)
  positiveSetTFBSDf <- tFTableMultipleLines[tFTableMultipleLines$EnsemblId %in% MHCGeneTable$EnsemblId,]
  
  print(paste("Amount of positive set genes in data set:", amountPositiveSetInList))
  print(head(positiveSetTFBSDf))
  #note: positiveSetTFBSDf thus contains the TFBS present in MHC pathway related genes
  
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
  
  
  totalCountsFullSet = sum(totalSetTFBSCounts$Counts)
  
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
  
  list(tFTableOneLine = tFTableOneLine,
                          tFTableMultipleLines = tFTableMultipleLines,
                          MHCGeneTable = MHCGeneTable,
                          MHCGenesInDataset = amountPositiveSetInList,
                          genomeOdds = genomeOdds,
                          MHCOdds = MHCOdds,
                          positiveSetTFBSCounts = perMotifCountsMHCSet,
                          totalSetTFBSCounts = totalSetTFBSCounts,
                          TFBSinMHCGenes = positiveSetTFBSDf
  )
  
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
      
      fisherDataFrame <- readRDS(file = paste(getwd(), "/humanMotifsFisherDataFrame29mammals.rds", sep = ""))
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
        print(totalCountsMHCGenes)
        print(positiveSetTFBSCountData$Counts[bindingSites])
        
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
        print(fisherMatrix)
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
      saveRDS(fisherDataFrame, file = paste(getwd(), "/humanMotifsFisherDataFrame29mammals.rds", sep = ""))
      if(outputMatrices == TRUE) {
        return(list(matrices = listMatrices, fisherDataFrame = fisherDataFrame))
      } else if (outputMatrices == FALSE) {
        return(fisherDataFrame)
      }
    }
    
    
  }
  
  ##########################################################
  #function scoresOverRows
  #takes the dataframe with Fisher tests and sample estimates as well as whole genome data (one line).
  #returns dataframe with score per TFBS and total score of all genes in the Pouya Kheradpout data.
  ##########################################################
  
scoresOverRows <- function(wholeGenomeDataOneLine, fisherTests, PositiveSetGenesInData,
                           scoreCorrectionMethod = 1) {
  
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
  
  
  informativeTFBS           <- subset(fisherTests, qValue <= 0.10 & MHCCount > 3)
  informativeTFBS[order(informativeTFBS$qValue),]
  rownames(informativeTFBS) <- seq(1:nrow(informativeTFBS))
  informativeTFBS$score <- 0
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

  
  ##########################################################
  #actual execution of the programme
  #loads positive set, checks which TFBS are present there and how often
  #does the same for whole genome data (WGD)
  #determines which TFBS overrepresented in MHC via fisher tests
  #assigns each gene a score log(sample odds ratio estimate from Fisher's exact p/corrected p-value Fisher's exact p) for every TFBS present
  #concatenates these individual scores into a total gene score and saves that to a file.
  ##########################################################
  dataStructureTFBS <- loadLibraryAndData(TRUE)
  str(dataStructureTFBS)
  dataStructureTFBS$tFTableOneLine
  fisherTests <- computeFishersTest(dataStructureTFBS$totalSetTFBSCounts,dataStructureTFBS$positiveSetTFBSCounts, run = TRUE, outputMatrices = FALSE)
  str(fisherTests)
  informativeTFBS <- subset(fisherTests, qValue <= 0.10 & MHCCount > 3)
  informativeTFBS
  
  
  #analyse the outcome of the fishertests in greater detail
  
  #what about Irf?
  fisherTests[fisherTests$TFBS %in% grep("IRF", fisherTests$TFBS, value = TRUE, fixed = TRUE),]
  #what about NFKB?
  fisherTests[fisherTests$TFBS %in% grep("NFKB", fisherTests$TFBS, value = TRUE, fixed = TRUE),]
  
  #switch for execution of the selection of all genes in the WGD that have the TFBS sig. overrepresented in positive set genes
  #set to TRUE when the positive set has changed and/or the analysis should be run again.
  #toggle which file to open (specific motifs versus more general motifs) below
  
  runselection = TRUE
  
  
  if (runselection == TRUE) {
    finalScoresWholeGenomeTwentyNineMammalsCorrectionOne <- scoresOverRows(dataStructureTFBS$tFTableOneLine,fisherTests, dataStructureTFBS$MHCGenesInDataset, scoreCorrectionMethod = 1)
    finalScoresWholeGenomeTwentyNineMammalsCorrectionTwo <- scoresOverRows(dataStructureTFBS$tFTableOneLine,fisherTests, dataStructureTFBS$MHCGenesInDataset, scoreCorrectionMethod = 2)
    finalScoresWholeGenomeTwentyNineMammalsCorrectionThree <- scoresOverRows(dataStructureTFBS$tFTableOneLine,fisherTests, dataStructureTFBS$MHCGenesInDataset, scoreCorrectionMethod = 3)
    listFinalScoresTwentyNineMammals <- list(finalScoresWholeGenomeTwentyNineMammalsCorrectionOne, finalScoresWholeGenomeTwentyNineMammalsCorrectionTwo, finalScoresWholeGenomeTwentyNineMammalsCorrectionThree)
    saveRDS(listFinalScoresTwentyNineMammals, file = "/home/dieter/Documents/Project/Programming/Transcription Factor Binding/29mammals/CallingMotifsMyself/Data/FinalDfAndUnifiedMotifFile/Generalised29mammalsMotifsTotalScoreWholeGenome.rds")
  }
  if (runselection == FALSE) {
    listFinalScoresTwentyNineMammals <- readRDS(file = "/home/dieter/Documents/Project/Programming/Transcription Factor Binding/29mammals/CallingMotifsMyself/Data/FinalDfAndUnifiedMotifFile/Generalised29mammalsMotifsTotalScoreWholeGenome.rds")
    finalScoresWholeGenomeTwentyNineMammalsCorrectionOne <- listFinalScoresTwentyNineMammals[[1]]
    finalScoresWholeGenomeTwentyNineMammalsCorrectionTwo <- listFinalScoresTwentyNineMammals[[2]]
    finalScoresWholeGenomeTwentyNineMammalsCorrectionThree <- listFinalScoresTwentyNineMammals[[3]]
  }
  str(finalScoresWholeGenomeTwentyNineMammalsCorrectionOne)
  head(finalScoresWholeGenomeTwentyNineMammalsCorrectionOne$scoresWholeGenome, 15)
  tail(finalScoresWholeGenomeTwentyNineMammalsCorrectionOne$scoresWholeGenome, 15)
  
  str(finalScoresWholeGenomeTwentyNineMammalsCorrectionTwo)
  head(finalScoresWholeGenomeTwentyNineMammalsCorrectionTwo$scoresWholeGenome, 15)
  tail(finalScoresWholeGenomeTwentyNineMammalsCorrectionTwo$scoresWholeGenome, 15)
  
  str(finalScoresWholeGenomeTwentyNineMammalsCorrectionThree)
  head(finalScoresWholeGenomeTwentyNineMammalsCorrectionThree$scoresWholeGenome, 15)
  tail(finalScoresWholeGenomeTwentyNineMammalsCorrectionThree$scoresWholeGenome, 15)
  

  
  #draw the distributions of the different scores
  
  listHistogramData <- list(length = length(listFinalScoresPouya))
  for(i in 1:length(listFinalScoresPouya)) {
    listHistogramData[[i]] <- melt(listFinalScoresTwentyNineMammals[[i]]$scoresWholeGenome)
    
    
  }
  
  #plot the scores for all TFBS
  listPlots <- list()
  for (i in 1:length(listHistogramData)) {
    
    b <- ggplot(data = listHistogramData[[i]], aes(x = value)) + geom_histogram() + facet_wrap( ~ variable) +
      theme_bw()
    listPlots <- append(listPlots, list(b))
    
  }
  listPlots
  
  
  
  
  
  
  #plot how much each MHC gene has of certain TFBS (either when including or excluding 0)
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
    
    infTFBS <- finalScoresWholeGenome$informativeTFBS$TFBS
    plotList = list()
    resultsDf <- data.frame(TF = character(), yesno = character(), genomeMHC = character(), proportion = numeric(), qValue = numeric(), tag = character())
    for (significantTFBS in 1: length(infTFBS))
    {
      
      #rewrite to take the actual counts instead of the sum of TRUEs which is presence/absence
      #
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
                                           qValue = rep(finalScoresWholeGenome$informativeTFBS[finalScoresWholeGenome$informativeTFBS == infTFBS[significantTFBS],]$qValue,4))
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
  #note: switch in the function with commenting (manual)
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
      ggsave(plot = annotatedPlot, filename = "/home/dieter/Documents/Project/Programming/Transcription Factor Binding/29mammals/CallingMotifsMyself/Data/FinalDfAndUnifiedMotifFile/TFBSOverrepresentationPlot_MoreGeneral.pdf", device = "pdf", dpi = 600)
      #switch for other filenames
      #ggsave(plot = annotatedPlot, filename = "/home/dieter/Documents/Project/Programming/Transcription Factor Binding/HumanWG_motifs/FinalFilesForEnrichmentAnalysis/TFBSOverrepresentationPlot_MoreSpecific.pdf", device = "pdf", dpi = 600)
      
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
          
          write.csv(genelist[[relevantTFBS]], "/home/dieter/Documents/Project/Programming/Transcription Factor Binding/29mammals/CallingMotifsMyself/Data/FinalDfAndUnifiedMotifFile/PositiveSetGenesAndAssociatedGeneralisedTFBS.csv")
          #swith different filename
          #write.csv(genelist[[relevantTFBS]], "/home/dieter/Documents/Project/Programming/Transcription Factor Binding/HumanWG_motifs/FinalFilesForEnrichmentAnalysis/PositiveSetGenesAndAssociatedSpecificTFBS.csv")
          
        }
        
        else if (relevantTFBS > 1) {
          
          write.table(genelist[[relevantTFBS]], "/home/dieter/Documents/Project/Programming/Transcription Factor Binding/29mammals/CallingMotifsMyself/Data/FinalDfAndUnifiedMotifFile/PositiveSetGenesAndAssociatedGeneralisedTFBS.csv",
                      append = TRUE, sep = ",", col.names = FALSE)
          #switch different filename
          #write.table(genelist[[relevantTFBS]], "/home/dieter/Documents/Project/Programming/Transcription Factor Binding/HumanWG_motifs/FinalFilesForEnrichmentAnalysis/PositiveSetGenesAndAssociatedSpecificTFBS.csv",
                      #append = TRUE, sep = ",", col.names = FALSE)
          
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
          
        } else if(!exists("dfWGDEnsmblIdsPerTFBS")) {
          
          dfWGDEnsmblIdsPerTFBS <- data.frame(listIds[[currentTFBS]])
          
        }
      }
      
      colnames(dfWGDEnsmblIdsPerTFBS) <- names(listIds$amount)
      #Once that is done, put all the columns into a dataframe, and save that as .csv
      write.csv(dfWGDEnsmblIdsPerTFBS, "/home/dieter/Documents/Project/Programming/Transcription Factor Binding/29mammals/CallingMotifsMyself/Data/FinalDfAndUnifiedMotifFile/EnsemblIdsWithPositiveSetEnrichedGeneralisedTFBS29mammals.csv")
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
  
  makePlotsDistributions <- function (WholeGenomeScoreTable = 1) {
    
    positiveSetIDs <- dataStructureTFBS$MHCGeneTable$EnsemblId
    inPositiveSet <- rownames(listFinalScoresTwentyNineMammals[[WholeGenomeScoreTable]]$scoresWholeGenome) %in% positiveSetIDs
    notInPositiveSet <- rownames(listFinalScoresTwentyNineMammals[[WholeGenomeScoreTable]]$scoresWholeGenome) %ni% positiveSetIDs
    finalPlotDfHist <- data.frame(scores = numeric(), set = character(), ScoreFor = character())
    for (i in 1:ncol(listFinalScoresTwentyNineMammals[[WholeGenomeScoreTable]]$scoresWholeGenome)) {
      finalScoresPositiveSetTotal <- listFinalScoresTwentyNineMammals[[WholeGenomeScoreTable]]$scoresWholeGenome[inPositiveSet, i]
      finalScoresNegativeSetTotal <- listFinalScoresTwentyNineMammals[[WholeGenomeScoreTable]]$scoresWholeGenome[notInPositiveSet, i]
      plotDfHistPositive <- data.frame(scores = finalScoresPositiveSetTotal, set = "positive")
      plotDfHistNegative <- data.frame(scores = finalScoresNegativeSetTotal, set = "negative")
      totalPlotDfHist <- rbind(plotDfHistNegative, plotDfHistPositive)
      totalPlotDfHist$ScoreFor <- colnames(listFinalScoresTwentyNineMammals[[WholeGenomeScoreTable]]$scoresWholeGenome)[i]
      finalPlotDfHist <- rbind(finalPlotDfHist, totalPlotDfHist)
    }
    
    
    #   ggplot(data = finalPlotDfHist, aes(x = scores)) +
    #     geom_density(aes(y = ..density.., fill = set), color = "black",  alpha = 0.3) +
    #     scale_y_continuous(breaks=waiver()) +
    #     #scale_x_continuous(limits = c(-100, 100)) +
    #     theme_bw() + 
    #     facet_wrap(~ScoreFor, scales = "free")
    #   
    #   bwidth <- 4
    #   chicken <- ggplot(data = finalPlotDfHist, aes(x = scores)) +
    #     geom_histogram(aes(y = 4 * ..density.., fill = set), color = "black", position = "dodge", stat = "bin", binwidth = 4) +
    #     scale_y_continuous(breaks=waiver(), expand=c(0,0), trans = "log1p") +
    #     theme_bw()  + 
    #     facet_wrap(~ScoreFor, scales = "free") 
    #   
    # chicken  
    
    kaas <- subset(finalPlotDfHist, ScoreFor == "totalScore")
    chicken <- ggplot(data = kaas, aes(x = scores)) +
      geom_histogram(aes(y = 4 * ..density.., fill = set), color = "black", position = "dodge", stat = "bin", binwidth = 4) +
      scale_y_continuous(breaks=waiver(), expand=c(0,0), trans = "log1p") +
      theme_bw()   
    
    pieter <- ggplot(data = kaas, aes(x = scores)) +
      geom_density(aes(y = ..density.., fill = set), color = "black",  alpha = 0.3) +
      scale_y_continuous(breaks=waiver()) +
      #scale_x_continuous(limits = c(-100, 100)) +
      theme_bw() + 
      facet_wrap(~ScoreFor, scales = "free")
    list(chicken, pieter)
  }
  makePlotsDistributions(1)
  makePlotsDistributions(2)
  makePlotsDistributions(3)
  
  
  
  
  #########################Further analysis of TFBS data##########################
  #How often is GSC present in all genes?
  GSC_subset <- subset(dataStructureTFBS$tFTableMultipleLines, TFBS == "GSC")
  head(GSC_subset)
  print(sum(GSC_subset$Count))
  
  #How often HoxD13
  HOXD13_subset <- subset(dataStructureTFBS$tFTableMultipleLines, TFBS == "HOXD13")
  head(HOXD13_subset)
  print(sum(HOXD13_subset$Count))
  
  subsetCountAnalyse <- function(x) {
    print(x)
    subset <- subset(dataStructureTFBS$tFTableMultipleLines, TFBS == x)
    
    print(quote = FALSE, paste("In total data:",sum(subset$Count)))
    print(quote = FALSE, paste("In MHC:", dataStructureTFBS$positiveSetTFBSCounts[dataStructureTFBS$positiveSetTFBSCounts$TFBS == x, 2]))
    }
  
  for( i in informativeTFBS$TFBS) {
    subsetCountAnalyse(i)
  }
  
  setwd(oldwd)
  ###################################################
  #                       END                       #
  ###################################################
  
  