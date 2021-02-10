

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
pacman::p_load(tidyverse, stringr, purrr, readr, LaF, ff, data.table, R.utils, base, utils, reshape2, optparse)
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
positiveSetFile       <- "~/Documents/Project/Programming/DataForAll/PositiveList_2.csv"
multipleLinesTFBSFile <- "~/Documents/Project/Programming/Transcription Factor Binding/HumanWG_motifs/Data/singleLineFinalDf and Unified table/underscoresRemovedMultipleLinesPouya.csv"
singleLineTFBSFile    <- "~/Documents/Project/Programming/Transcription Factor Binding/HumanWG_motifs/Data/singleLineFinalDf and Unified table/underscoresRemovedSingleLinePouya.csv"

#if you do not wish to load data from source files, input the .rds file locations here
TFBSDataObjectRDS   <- NA
fisherDataObjectRDS <- NA


#working directory
oldwd <- getwd()
setwd("~/Documents/Project/Programming/Transcription Factor Binding/HumanWG_motifs/Data/singleLineFinalDf and Unified table/")

#directory to save output of the run
customDir   <- NA
if(is.na(customDir)) {
  dirForRun   <- paste0(getwd(), "/program_run_at_", gsub(" ", "_", Sys.time(), fixed = TRUE))
  dir.create(dirForRun)
} else {
  dirForRun <- customDir
  if (!dir.exists(customDir)) {
    
    dir.create(customDir)
  }
}
print(paste0("Directory for this run: ", dirForRun))

#initiate log file for this run
file.create(paste0(dirForRun, "/log.txt"))
fileConn <- file(paste0(dirForRun, "/log.txt"))
writeLines(paste0("Log for TFBSenrichr, run at ", Sys.time()), fileConn)



######################################################################################################
######################################################################################################
######                                                                                          ######
######                                Data Manipulation Functions                               ######
######                                                                                          ######
######################################################################################################
######################################################################################################

logWrite <- function(toWrite) {
  
  fileConn <- file(paste0(dirForRun, "/log.txt"))
  readr::write_lines(path = fileConn, x = toWrite, append = TRUE)
  
}

#loadingLibrary
loadLibraryAndData     <- function(foldCrossVal = 10, entityName = "TFBS", UI = FALSE, reloadFromSource = TRUE) {
  
  
  
  #write arguments to file
  logWrite("Running function loadLibraryAndData\n arguments:\n")
  logWrite(as.character(unlist(sys.call())))
  
  
  if(reloadFromSource == TRUE) {
    
    colNamesFeatureTables <- c("EnsemblId", entityName)
    
    
    
    if (UI == TRUE) {
      print(paste0("reading in file ", singleLineTFBSFile))
      tFTableOneLine <- fread(tk_choose.files(default = "~/Documents/Project/Programming/Transcription Factor Binding/HumanWG_motifs/Data/singleLineFinalDf and Unified table/underscoresRemovedSingleLinePouya.csv",
                                              caption = "Choose the .csv file that contains the ENCODE TFBS on one line",
                                              multi = FALSE))
      
      colnames(tFTableOneLine) <- colNamesFeatureTables
      print(head(tFTableOneLine))
      print(paste0("rows original : ", nrow(tFTableOneLine)))
      
      #remove genes that have no motifs
      tFTableOneLine <- tFTableOneLine[tFTableOneLine$TFBS != "",]
      print(paste0("rows non-empty : ", nrow(tFTableOneLine)))
      
      
      print(paste0("reading in file ", multipleLinesTFBSFile))
      tFTableMultipleLines <- fread(tk_choose.files(default = "~/Documents/Project/Programming/Transcription Factor Binding/HumanWG_motifs/Data/singleLineFinalDf and Unified table/underscoresRemovedMultipleLinesPouya.csv",
                                                    caption = "Choose the .csv file that contains the ENCODE TFBS on multiple lines",
                                                    multi = FALSE))
      
      colnames(tFTableMultipleLines) <- c(colNamesFeatureTables, "Count")
      print(head(tFTableMultipleLines))
      #remove genes without motifs/NA motifs
      print(paste0("rows original : ", nrow(tFTableMultipleLines)))
      tFTableMultipleLines <- tFTableMultipleLines[tFTableMultipleLines$Count > 0,]
      print(paste0("rows non-empty : ", nrow(tFTableMultipleLines)))
      
      print(paste0("reading in file ", positiveSetFile))
      positiveSetTable <- fread(tk_choose.files(default = "~/Documents/Project/Programming/DataForAll/PositiveList_2.csv",
                                                caption = "Choose the .csv file that contains the positive MHC gene set",
                                                multi = FALSE))
      names(positiveSetTable)[7] <- "EnsemblId"
      print(head(positiveSetTable))
      
    } else {
      
      print(paste0("reading in file ", singleLineTFBSFile))
      tFTableOneLine <- fread(singleLineTFBSFile)
      
      colnames(tFTableOneLine) <- colNamesFeatureTables
      print(head(tFTableOneLine))
      print(paste0("rows original : ", nrow(tFTableOneLine)))
      
      #remove genes that have no motifs
      tFTableOneLine <- tFTableOneLine[tFTableOneLine$TFBS != "",]
      print(paste0("rows non-empty : ", nrow(tFTableOneLine)))
      
      
      print(paste0("reading in file ", multipleLinesTFBSFile))
      tFTableMultipleLines <- fread(multipleLinesTFBSFile)
      
      colnames(tFTableMultipleLines) <- c(colNamesFeatureTables, "Count")
      print(head(tFTableMultipleLines))
      #remove genes without motifs/NA motifs
      print(paste0("rows original : ", nrow(tFTableMultipleLines)))
      tFTableMultipleLines <- tFTableMultipleLines[tFTableMultipleLines$Count > 0,]
      print(paste0("rows non-empty : ", nrow(tFTableMultipleLines)))
      
      print(paste0("reading in file ", positiveSetFile))
      positiveSetTable <- fread(positiveSetFile)
      names(positiveSetTable)[7] <- "EnsemblId"
      print(head(positiveSetTable))
      
    }
    
    
    amountPositiveSetInList <- sum(positiveSetTable$EnsemblId %in% tFTableMultipleLines$EnsemblId)
    positiveSetTFBSDf <- tFTableMultipleLines[tFTableMultipleLines$EnsemblId %in% positiveSetTable$EnsemblId,]
    negativeSetTFBSDf <- tFTableMultipleLines[tFTableMultipleLines$EnsemblId %ni% positiveSetTable$EnsemblId,]
    totalSetTFBSDf    <- tFTableMultipleLines
    
    print(paste("Amount of positive set genes in data set:", amountPositiveSetInList))
    print(head(positiveSetTFBSDf))
    #note: positiveSetTFBSDf thus contains the TFBS present in MHC pathway related genes
    
    
    
    
    print("Generating sets for tenfold cross-validation, hold on...")
    #now generate a list of 10 sets to leave out per time
    #ensemblIdsToLeaveOutListPositive <- vector("list", length = foldCrossVal)
    #ensemblIdsToLeaveOutListNegative <- vector("list", length = foldCrossVal)
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
    totalSetNotTrainedOnWGDOneLine   <- vector("list", length = foldCrossVal)
    
    #shuffle Ids of both positive and negative set for exclusion
    shuffleIdsIntoGroups    <- function(IdVector, groups = foldCrossVal) {
      
      IdVector <- sample(IdVector)
      safe_split <- quietly(split)
      returnList <- safe_split(IdVector, rep(1:10, ceiling(length(IdVector)/10)))
      returnList
      
    }
    shuffledIdsListPositive <- shuffleIdsIntoGroups(unique(positiveSetTFBSDf$EnsemblId))
    shuffledIdsListNegative <- shuffleIdsIntoGroups(unique(negativeSetTFBSDf$EnsemblId))
    
    
    #get the ten sets to keep for testing post-training and the ten sets to train on pre-testing, both negative and positive, and combined
    for(i in 1:length(positiveSetNotTrainedOn)) {
      
      positiveSetNotTrainedOn[[i]]        <- positiveSetTFBSDf[positiveSetTFBSDf$EnsemblId %in% shuffledIdsListPositive$result[[i]], ]
      positiveSetTrainedOn[[i]]           <- positiveSetTFBSDf[positiveSetTFBSDf$EnsemblId %ni% shuffledIdsListPositive$result[[i]], ]
      negativeSetNotTrainedOn[[i]]        <- negativeSetTFBSDf[negativeSetTFBSDf$EnsemblId %in% shuffledIdsListNegative$result[[i]], ]
      negativeSetTrainedOn[[i]]           <- negativeSetTFBSDf[negativeSetTFBSDf$EnsemblId %ni% shuffledIdsListNegative$result[[i]], ]
      totalSetTrainedOn[[i]]              <- rbind(positiveSetTrainedOn[[i]], negativeSetTrainedOn[[i]])
      totalSetNotTrainedOn[[i]]           <- rbind(positiveSetNotTrainedOn[[i]], negativeSetNotTrainedOn[[i]])
      totalSetNotTrainedOnWGDOneLine[[i]] <- tFTableOneLine[tFTableOneLine$EnsemblId %in% totalSetNotTrainedOn[[i]]$EnsemblId, ]
    }
    
    
    #aggregate the amounts by TFBS
    generateTotalMotifCounts <- function(setTFBSDf) {
      
      perMotifCountsPositiveSet <- setTFBSDf %>% group_by(TFBS) %>% summarise(Count = sum(Count))
      perMotifCountsPositiveSet
    }
    generateTotalMotifCountsList <- function(listOfSets) {
      
      resultList <- vector("list", length = length(listOfSets))
      
      for(i in 1:length(listOfSets)) {
        
        resultList[[i]] <- generateTotalMotifCounts(setTFBSDf = listOfSets[[i]])
        
      }
      resultList
    }
    
    totalCountsPositiveSet    <- sum(positiveSetTFBSDf$Count)
    perMotifCountsPositiveSet <- generateTotalMotifCounts(positiveSetTFBSDf)
    perMotifCountsFullSet     <- generateTotalMotifCounts(totalSetTFBSDf)
    totalCountsFullSet        <- sum(perMotifCountsFullSet$Count)
    
    #for Crossval
    print("Generating count data for tenfold cross-validation, hold on...")
    perMotifCountsPositiveSetListCrossVal <- lapply(positiveSetTrainedOn, FUN = function(x) { perMotifCountsPositiveSet <- x %>% group_by(TFBS) %>% summarise(Count = sum(Count)); perMotifCountsPositiveSet}  )
    PerMotifCountsListTotal               <- lapply(totalSetTrainedOn, FUN = function(x) { perMotifCountsPositiveSet <- x %>% group_by(TFBS) %>% summarise(Count = sum(Count)); perMotifCountsPositiveSet}  )
    
    #return
    
    returnList <- list(                   tFTableOneLine                = tFTableOneLine,
                                          tFTableMultipleLines          = tFTableMultipleLines,
                                          positiveSetTable              = positiveSetTable,
                                          MHCGenesInDataset             = amountPositiveSetInList,
                                          perMotifCountsPositiveSet     = perMotifCountsPositiveSet,
                                          perMotifCountsFullSet         = perMotifCountsFullSet,
                                          TFBSinMHCGenes                = positiveSetTFBSDf,
                                          cvPositiveSetNotTrainedOnList = positiveSetNotTrainedOn,
                                          cvPositiveSetTrainedOnList    = positiveSetTrainedOn,
                                          cvNegativeSetTrainedOnList    = negativeSetTrainedOn,
                                          cvNegativeSetNotTrainedOnList = negativeSetNotTrainedOn,
                                          cvTotalSetTrainedOnList       = totalSetTrainedOn,
                                          cvTotalSetNotTrainedOnList    = totalSetNotTrainedOn,
                                          cvTotalSetNotTrainedOnWGDOL   = totalSetNotTrainedOnWGDOneLine,
                                          cvPositiveSetTFBSCountsList   = perMotifCountsPositiveSetListCrossVal,
                                          cvTotalSetTFBSCountsList      = PerMotifCountsListTotal
    )
    print("Saving TFBS data object, hold on...")
    saveRDS(returnList, file = paste0(dirForRun, "/TFBSDataObject.rds"))
    
    
    
    
  } else if (reloadFromSource == FALSE) {
    
    if(is.na(TFBSDataObjectRDS)) {
      print("You have not supplied the data object used when not reloading all data from source. Specify it above in variable TFBSDataObjectRDS")
      print("Choose the file now via the file browser")
      choice <- tk_choose.files(default = "~/Documents/Project/Programming/Transcription Factor Binding/HumanWG_motifs/Data/singleLineFinalDf and Unified table",
                                caption = "Choose the .rds file that contains the processed TFBS data",
                                multi = FALSE)
      if(endsWith(choice, ".rds")) {
        returnList <- readRDS(choice)
      }
    } else {
      returnList <- readRDS(TFBSDataObjectRDS)
    }
    
  }
  
  return(returnList)
}

#fisherTestFunctions
computeFishersTest     <- function(wholeGenomeTFBSCountData, positiveSetTFBSCountData, runNewTest = TRUE, debug = FALSE) {
  
  #don't compute again unless explicitly told to. Otherwise, just open the .rds
  if(runNewTest == FALSE) {
    
    if(! is.na(fisherDataObjectRDS)) {
      
      fisherDataFrame <- readRDS(file = fisherDataObjectRDS) 
      
    } else {
      
      print("You did not specify a file from which to load the previously computed fisherScores")
      print("please do so now:")
      choice <- tcltk::tk_choose.files(default = "", caption = "select the .rds containing the previously computed fisher test scores", multi = FALSE)
      if(endsWith(choice, ".rds")) {
        fisherDataFrame <- readRDS(file = choice)
        
      } else {
        stop("That is not a correct .rds file, mister!")
      }
      
    }
    
    return(fisherDataFrame)
    
  }
  else {
    
    print(head(wholeGenomeTFBSCountData))
    print(head(positiveSetTFBSCountData))
    
    #Get only those genes from the whole genome that have TFBS that are also found in the positive MHC set
    relevantWholeGenome    <- wholeGenomeTFBSCountData[wholeGenomeTFBSCountData$TFBS %in% positiveSetTFBSCountData$TFBS,] 
    totalCountsWholeGenome <- sum(wholeGenomeTFBSCountData$Count)
    totalCountsMHCGenes    <- sum(positiveSetTFBSCountData$Count)
    if (debug == TRUE) {
      print(head(relevantWholeGenome))
    }
    
    #listMatrices is for debug purposes, to check whether correct contingency tables are created
    listMatrices <- vector("list", length = length(unique(relevantWholeGenome$TFBS)))
    fisherDataFrame <- as.data.table(matrix(ncol = 7, nrow = length(unique(relevantWholeGenome$TFBS))))
    colnames(fisherDataFrame) <- c("TFBS", "pValue", "sampleEstimate", "lowerConfidenceInterval", "upperConfidenceInterval", "wholeGenomeCount", "positiveSetCount")
    fisherDataFrame$TFBS <- rep("NA", nrow(fisherDataFrame))
    fisherDataFrame[, colnames(fisherDataFrame)[-1] := lapply(.SD, as.numeric), .SDcols = colnames(fisherDataFrame)[-1]]
    fisherDataFrame
    if(debug == TRUE) {
      print(class(fisherDataFrame$TFBS))
      print(class(fisherDataFrame$pValue))
    }
    
    #cycle over all the rows of 
    for(bindingSites in 1:nrow(positiveSetTFBSCountData)) {
      
      ##Note: since I often forget exactly how this works I am writing in here what I am doing:
      ##You have the amount of this specific TFBS in the positive set.
      ##You have the amount of all other TFBS in the positive set (so, total - this specific one)
      ##You have the amount of this specific TFBS in the whole genome
      ##You have the amount of all other TFBS in the whole genome (so, whole genome total - this specific one in the whole genome data)
      
      fisherMatrix <- matrix(c(positiveSetTFBSCountData$Count[bindingSites],
                               totalCountsMHCGenes - positiveSetTFBSCountData$Count[bindingSites],
                               relevantWholeGenome$Count[bindingSites],
                               totalCountsWholeGenome - relevantWholeGenome$Count[bindingSites]),
                             nrow = 2,
                             dimnames = list(TFBS = c(relevantWholeGenome$TFBS[bindingSites], "TotalTFBS"),
                                             Set = c("Positive set", "Whole genome")))
      #print(fisherMatrix)
      fisherTestResult    <- fisher.test(fisherMatrix)
      pValue              <- fisherTestResult["p.value"]
      lowerConfInt        <- unlist(fisherTestResult["conf.int"])[1]
      upperConfInt        <- unlist(fisherTestResult["conf.int"])[2]
      sampleEst           <- fisherTestResult["estimate"]
      wholeGenomePresence <- fisherMatrix[1,2]
      positiveSetPresence <- fisherMatrix[1,1]
      
      listMatrices[[bindingSites]] <- fisherMatrix
      addVector        <- c(as.character(relevantWholeGenome$TFBS[bindingSites]),
                            pValue,
                            sampleEst,
                            lowerConfInt,
                            upperConfInt,
                            wholeGenomePresence,
                            positiveSetPresence)
      names(addVector) <- c("TFBS", "pValue", "sampleEstimate",
                            "lowerConfidenceInterval", "upperConfidenceInterval",
                            "wholeGenomeCount", "positiveSetCount")
      #add to the Df the results of the fisherTest
      fisherDataFrame[bindingSites,  ]  <- addVector
      if(debug == TRUE) {
        print(addVector)
        print(head(fisherDataFrame))
      }
      
    }
    
    #adjust for multiple testing with Benjamini-Hochberg methodology (utilises FDR)
    fisherDataFrame$qValue <- p.adjust(fisherDataFrame$pValue, method="BH")
    
    
    if(debug == TRUE) {
      saveRDS(list(matrices = listMatrices, fisherDataFrame = fisherDataFrame), file = paste0(dirForRun, "/humanMotifsFisherDataFramePouya.rds"))
      return(list(matrices = listMatrices, fisherDataFrame = fisherDataFrame))
    } else {
      saveRDS(fisherDataFrame, file = paste0(dirForRun, "/humanMotifsFisherDataFramePouya.rds"))
      return(fisherDataFrame)
    }
  }
  
  
}

crossValFisherTests    <- function(crossValTotalSetTFBSCounts, crossValPositiveTFBSCountList) {
  
  result <- purrr::pmap(list(crossValTotalSetTFBSCounts,crossValPositiveTFBSCountList), .f = computeFishersTest)
}

#scoring over rows functions
scoresOverRows         <- function(wholeGenomeDataOneLine, fisherTests, positiveSetGenesInData, dataStrucTFBS = dataStructureTFBS,
                                   debug = FALSE, runBayesFactor = FALSE, includeCounts = FALSE, scoreStruc = newTFBSData) {
  
  #find TFBS that are significantly overrepresented in positive set genes
  #assign them a score. In this case, the point estimate from the Fisher tests.
  
  if (runBayesFactor == FALSE) {
  informativeTFBS                    <- subset(fisherTests, qValue <= 0.10 & positiveSetCount > 3)
  informativeTFBS                    <- informativeTFBS[order(informativeTFBS$qValue),]
  #rownames(informativeTFBS)         <- seq(1:nrow(informativeTFBS))
  informativeTFBS$score              <- log2(informativeTFBS$sampleEstimate)
  informativeTFBS$positiveSetAverage <- informativeTFBS$positiveSetCount/positiveSetGenesInData
  
  
  if(debug == TRUE) {
    print(head(informativeTFBS))
  }
  
  
  scoreCompendium  <- vector()
  dataFrameResults <- data.frame()
  #cycle over rows. Could be refined with an apply function
  
  
  
  
  
  wholeGenomeScores <- as.data.table(do.call(rbind, purrr::map2(wholeGenomeDataOneLine$EnsemblId,
                                                                wholeGenomeDataOneLine$TFBS,
                                                                .f = function(rowEnsemblId, rowTFBS) {
                                                                  
                                                                  TFBSCounts <- stringr::str_match_all(rowTFBS,"([^|]+)-(\\d+)")
                                                                  namesTFBS  <- TFBSCounts[[1]][,2]
                                                                  countsTFBS <- as.numeric(TFBSCounts[[1]][,3])
                                                                  names(countsTFBS) <- namesTFBS
                                                                  
                                                                  presenceOfInformativeTFBS <- informativeTFBS$TFBS %in% namesTFBS
                                                                  #initialse the score. Note the +1 for the totalScore column
                                                                  scoreCompendium           <- vector("numeric", length = (length(presenceOfInformativeTFBS) + 1))
                                                                  
                                                                  for(particularTFBS in 1:length(presenceOfInformativeTFBS)) {
                                                                    
                                                                    if(presenceOfInformativeTFBS[particularTFBS] == TRUE) {
                                                                      
                                                                      #if it is there but the score is negative, keep as is. If the score is positive, correct to the average positive set amount
                                                                      
                                                                      if(informativeTFBS$sampleEstimate[particularTFBS] < 1) {
                                                                        
                                                                        valueToAdd <- informativeTFBS$score[particularTFBS]
                                                                      } else {
                                                                        
                                                                        #als gemiddeld aantal keer deze TFBS >1 in MHC, corrigeer dan het te scoren gen voor het aantal keer dat het deze TFBS heeft
                                                                        if(informativeTFBS$positiveSetAverage[particularTFBS] > 1) {
                                                                          particularTFBSCount  <- countsTFBS[informativeTFBS$TFBS[particularTFBS]]
                                                                          valueToAdd <- informativeTFBS$score[particularTFBS] * (particularTFBSCount / informativeTFBS$positiveSetAverage[particularTFBS])
                                                                        } else {
                                                                          valueToAdd <- informativeTFBS$score[particularTFBS]
                                                                        }
                                                                        
                                                                      }
                                                                      
                                                                    } else {
                                                                      
                                                                      valueToAdd <- 0
                                                                      
                                                                    }
                                                                    
                                                                    scoreCompendium[particularTFBS] <- valueToAdd 
                                                                    
                                                                    
                                                                  }
                                                                  #add the totalScore
                                                                  
                                                                  names(scoreCompendium) <- c(informativeTFBS$TFBS, "totalScore")
                                                                  scoreCompendium
                                                                }
  )))
  wholeGenomeScores$totalScore <- rowSums(wholeGenomeScores)
  wholeGenomeScores$EnsemblId  <- wholeGenomeDataOneLine$EnsemblId
  
  print("These are the informative TFBS:")
  print(informativeTFBS$TFBS)
  
  list(
    scoresWholeGenome = wholeGenomeScores,
    informativeTFBS   = informativeTFBS)
  } else {
    
    if (debug == TRUE) {
    print("stepping into purrr now")
    }
    wholeGenomeScores <- as.data.table(do.call(rbind, purrr::map2(wholeGenomeDataOneLine$EnsemblId,
                                                                  wholeGenomeDataOneLine$TFBS,
                                                                  .f = function(rowEnsemblId, rowTFBS) {
                                                                    
                                                                    TFBSCounts <- stringr::str_match_all(rowTFBS,"([^|]+)-(\\d+)")
                                                                    namesTFBS  <- TFBSCounts[[1]][,2]
                                                                    countsTFBS <- as.numeric(TFBSCounts[[1]][,3])
                                                                    names(countsTFBS) <- namesTFBS
                                                                    
                                                                    ########################
                                                                    #changed to BayesFactor#
                                                                    ########################
                                                                  
                                                                    
                                                                    
                                                                    
                                                                    presenceOfInformativeTFBS <- unique(dataStrucTFBS$tFTableMultipleLines$TFBS) %in% namesTFBS
                                                                    print(presenceOfInformativeTFBS)
                                                                    #initialse the score. Note the +1 for the totalScore column
                                                                    scoreCompendium           <- vector("numeric", length = (length(presenceOfInformativeTFBS) + 1))
                                                                    
                                                                    for(particularTFBS in 1:length(presenceOfInformativeTFBS)) {
                                                                      
                                                                      nameOfCurrentTFBS <- unique(dataStrucTFBS$tFTableMultipleLines$TFBS)[particularTFBS]
                                                                      scoreRow <- scoreStruc %>% filter(TFBS == nameOfCurrentTFBS)
                                                                      
                                                                      if(debug == TRUE) {
                                                                      print(nameOfCurrentTFBS)
                                                                      }
                                                                      
                                                                      ##this is a ratio of (counts for this TFBS in all A genes / total counts of all TFBS in A)  / counts for this TFBS in all B genes / total counts of all TFBS in B
                                                                      ##so this would correct on the counts (versus yesno)
                                                                      #particularCountsA <- dataStrucTFBS$tFTableMultipleLines %>% filter(EnsemblId %in% dataStrucTFBS$MHCGeneTable$EnsemblId & TFBS == particularTFBS) %>% select(Count) %>% sum()
                                                                      #totalCountsA <- dataStrucTFBS$tFTableMultipleLines %>% filter(EnsemblId %in% dataStrucTFBS$MHCGeneTable$EnsemblId) %>% select(Count) %>% sum()
                                                                      #particularCountsB <- dataStrucTFBS$tFTableMultipleLines %>% filter(EnsemblId %ni% dataStrucTFBS$MHCGeneTable$EnsemblId & TFBS == particularTFBS) %>% select(Count) %>% sum()
                                                                      #totalCountsB <- dataStrucTFBS$tFTableMultipleLines %>% filter(EnsemblId %ni% dataStrucTFBS$MHCGeneTable$EnsemblId) %>% select(Count) %>% sum()
                                                                      
                                                                      #ratioA <- particularCountsA/totalCountsA
                                                                      #ratioB <- particularCountsB/totalCountsB
                                                                      
                                                                      #ratioAOverRatioB <- ratioA/ratioB
                                                                      
                                                                      ##this is a ratio that gives how much of the average MHC count a gene has
                                                                      
                                                                      if(presenceOfInformativeTFBS[particularTFBS] == TRUE) {
                                                                        
                                                                        #if it is there but the score is negative, keep as is. If the score is positive, correct to the average positive set amount
                                                                        
                                                                        
                                                                      
                                                                        #probleem --> zelfs voor iets als NfKB of IRF is dit niet of nauwelijks positief. Waarom: veel genen hebben 1 keer dat motief...
                                                                        
                                                                        if (includeCounts == FALSE) {
                                                                          valueToAdd <- log2(scoreRow$ratioPaiPbi)
                                                                        } else {
                                                                          
                                                                          
                                                                          
                                                                          
                                                                          
                                                                          if(scoreRow$avgCountsMHC < 1) {
                                                                          
                                                                          valueToAdd <- log2(scoreRow$ratioPaiPbi)
                                                                          } else {
                                                                            particularCountsThisGene <-  dataStrucTFBS$tFTableMultipleLines %>% filter(EnsemblId == rowEnsemblId & TFBS == nameOfCurrentTFBS) %>% select(Count) %>% `[[`(1)
                                                                            
                                                                            
                                                                          valueToAdd <- log2((scoreRow$ratioPaiPbi) * (particularCountsThisGene/scoreRow$avgCountsMHC))
                                                                          }}
                                                                          
                                                                        
                                                                        
                                                                          
                                                                          } else {
                                                                        
                                                                        valueToAdd <- log2(scoreRow$absentratioPaiPbi)
                                                                        
                                                                      }
                                                                      
                                                                      scoreCompendium[particularTFBS] <- valueToAdd 
                                                                      
                                                                      
                                                                    }
                                                                    #add the totalScore
                                                                    
                                                                    names(scoreCompendium) <- c(informativeTFBS$TFBS, "totalScore")
                                                                    scoreCompendium
                                                                  }
    )))
    wholeGenomeScores$totalScore <- rowSums(wholeGenomeScores)
    wholeGenomeScores$EnsemblId  <- wholeGenomeDataOneLine$EnsemblId
    
    print("These are the informative TFBS:")
    print(informativeTFBS$TFBS)
    
    wholeGenomeScores
    
    
    
    
  }
  
  
}

crossValScoresOverRows <- function(crossValWholeGenomeDataOneLineList, crossValFisherTestList, crossValPositiveSets) {
  
  positiveSetGenesInData <- lapply(crossValPositiveSets, FUN = function(x) {length(unique(x$EnsemblId))})
  result <- purrr::pmap(list(crossValWholeGenomeDataOneLineList, crossValFisherTestList, positiveSetGenesInData), .f = scoresOverRows)
  result
  
}

######################################################################################################
######################################################################################################
######                                                                                          ######
######                              Plotting/Visualisation Functions                            ######
######                                                                                          ######
######################################################################################################
######################################################################################################

plotTFBSPresencePositiveSetGenes <- function(dataStructureTFBS, informativeTFBS) {
  
  #see (for those genes that have the TFBS at least once) the distribution of the counts of TFBS in positive set Genes  
  plotDataBoxplot <- data.table(TFBS = character(), Count = numeric())
  for(i in 1:nrow(informativeTFBS)) {
    TFBSToAdd <- dataStructureTFBS$TFBSinMHCGenes[dataStructureTFBS$TFBSinMHCGenes$TFBS == informativeTFBS$TFBS[i]][,c(2,3)]
    plotDataBoxplot <- rbind(plotDataBoxplot, TFBSToAdd) 
  }
  
  plotWithoutZeroValues <- ggplot(data = plotDataBoxplot, aes(x = TFBS, y = Count)) + geom_boxplot(outlier.alpha = 0) +
    theme_bw() + geom_point(alpha = 0.3, position = position_jitter(width = 0.4, height = 0.2)) +
    theme(axis.text.x=element_text(angle=60,hjust=1)) + scale_y_continuous(breaks = seq(0,(max(plotDataBoxplot$Count) + 10), 2)) +
    theme(panel.grid.major.x = element_blank())
  
  
  
}

createBarPlotData                <- function(TFBSdata = dataStructureTFBS, informativeTFBS) {
  
  infTFBS <- informativeTFBS$TFBS
  resultsDf <- data.frame(TF = character(), yesno = character(), genomeMHC = character(), proportion = numeric(), qValue = numeric(), tag = character())
  for (significantTFBS in seq_along(infTFBS))
  {
    totalTFBSInPositiveSetGenes    <- sum(TFBSdata$TFBSinMHCGenes$Count)
    specificTFBSInPositiveSetGenes <- TFBSdata$TFBSinMHCGenes$TFBS %in% infTFBS[significantTFBS]
    totalSpecificTFBSInPositiveSet <- sum(TFBSdata$TFBSinMHCGenes[specificTFBSInPositiveSetGenes,]$Count)
    
    totalTFBSInWGDGenes            <- sum(TFBSdata$perMotifCountsFullSet$Count)
    specificTFBSInWGDGenes         <- TFBSdata$perMotifCountsFullSet$TFBS %in% infTFBS[significantTFBS]
    totalSpecificTFBSInWGD         <- sum(TFBSdata$perMotifCountsFullSet[specificTFBSInWGDGenes,]$Count)
    
    
    proportionInPositiveSet        <- totalSpecificTFBSInPositiveSet/totalTFBSInPositiveSetGenes
    proportionNotInPositiveSet     <- 1 - proportionInPositiveSet
    proportionInWGD                <- totalSpecificTFBSInWGD/totalTFBSInWGDGenes
    proportionNotInWGD             <- 1- totalSpecificTFBSInWGD/totalTFBSInWGDGenes 
    tempDf                         <- data.frame(TF = rep(infTFBS[significantTFBS],4),
                                                 yesno = c("yes", "no", "yes", "no"),
                                                 genomePositiveSet = c("positive set", "positive set", "genome", "genome"),
                                                 proportion = c(proportionInPositiveSet, proportionNotInPositiveSet, proportionInWGD, proportionNotInWGD),
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
  resultsDf
}

plotBarPlotData                  <- function(plottingData) {
  
  #takes as argument the output generated by createGraphData
  #generates and saves a bar plot
  
  
  green <- rgb(0, 140, 71, maxColorValue = 255)
  red   <- rgb(237, 45, 46, maxColorValue = 255)
  initialPlot <-ggplot(data = plottingData,
                       aes(x= factor(yesno, levels = c("yes","no")),
                           y = proportion, fill = factor(genomePositiveSet, levels = c ("positive set", "genome")))) + 
    facet_wrap(~TF, scales = "free") +
    geom_bar(stat = "identity", position = position_dodge(), colour = "black") +
    theme_minimal() +
    xlab ("TFBS present") +
    scale_fill_manual(breaks = c("positive set","genome"),
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
  ggsave(plot = annotatedPlot,
         filename = paste0(dirForRun, "/TFBSOverrepresentationPlotPouyaGeneralised.pdf"),
         device = "pdf", dpi = 600)
  
  annotatedPlot
}

plotScoreDistributions           <- function (dataStructureTFBS      = dataStructureTFBS,
                                              finalScoresWholeGenome = finalScoresWholeGenome$scoresWholeGenome,
                                              debug = FALSE) {
  
  #plots histograms of the scores obtained by scoresOverRows
  positiveSetIDs   <- dataStructureTFBS$positiveSetTable$EnsemblId
  inPositiveSet    <- finalScoresWholeGenome$EnsemblId %in% positiveSetIDs
  notInPositiveSet <- finalScoresWholeGenome$EnsemblId %ni% positiveSetIDs
  finalPlotDfHist  <- data.table(scores = numeric(), set = character(), ScoreFor = character())
  scoresOnlyDT     <- finalScoresWholeGenome[, !"EnsemblId", with=FALSE]
  
  for (TFBSScores in 1:ncol(scoresOnlyDT)) {
    if(debug == TRUE) {
      print(TFBSScores)
    }
    finalScoresPositiveSetTotal <- scoresOnlyDT[inPositiveSet, TFBSScores, with = FALSE]
    if(debug == TRUE) {
      print(finalScoresPositiveSetTotal)
      Sys.sleep(3) }
    finalScoresNegativeSetTotal <- scoresOnlyDT[notInPositiveSet, TFBSScores, with = FALSE]
    if(debug == TRUE) {
      print(finalScoresNegativeSetTotal)
      Sys.sleep(3) }
    plotDfHistPositive          <- data.table(scores = finalScoresPositiveSetTotal[[1]], set = "positive")
    if(debug == TRUE) {
      print(plotDfHistPositive)
      Sys.sleep(3) }
    plotDfHistNegative          <- data.table(scores = finalScoresNegativeSetTotal[[1]], set = "negative")
    if(debug == TRUE) {
      print(plotDfHistNegative)
      Sys.sleep(3) }
    totalPlotDfHist             <- rbind(plotDfHistNegative, plotDfHistPositive)
    totalPlotDfHist$ScoreFor    <- colnames(scoresOnlyDT)[TFBSScores]
    finalPlotDfHist             <- rbind(finalPlotDfHist, totalPlotDfHist, fill = TRUE)
    if(debug == TRUE) {
      print(finalPlotDfHist)
      Sys.sleep(3) }
  }
  
  totalScoreData <- finalPlotDfHist %>% dplyr::filter(ScoreFor == "totalScore")
  
  histogramPlotTotalScore  <- ggplot(data = totalScoreData, aes(x = scores)) +
    geom_histogram(aes(y = ..density.., fill = set), color = "black", position = "dodge", stat = "bin", binwidth = 1) +
    scale_y_continuous(breaks=waiver(), expand=c(0,0)) +
    theme_bw() + ggtitle("Distribution of total TFBS score in positive set (MHC genes)\n and negative set (all genes in the human genome)") +
    theme(plot.title = element_text(hjust = 0.5))  
  
  densityPlotTotalScore    <- ggplot(data = totalScoreData, aes(x = scores)) +
    geom_density(aes(y = ..density.., fill = set), color = "black",  alpha = 0.3) +
    scale_y_continuous(breaks=waiver()) +
    theme_bw() + ggtitle("Distribution of total TFBS score in positive set (MHC genes)\n and negative set (all genes in the human genome)") +
    theme(plot.title = element_text(hjust = 0.5))
  
  #now do this also for all the separate TFBS
  TFBSForWhichToScore <- unique(finalPlotDfHist$ScoreFor)
  listPlotsPerTFBS    <- vector("list", length(TFBSForWhichToScore))
  for(TFBS in seq_along(TFBSForWhichToScore)) {
    
    dataToUse     <- finalPlotDfHist %>% dplyr::filter(ScoreFor == TFBSForWhichToScore[TFBS])
    histogramPlot <- ggplot(data = dataToUse, aes(x = scores)) +
      geom_histogram(aes(y = ..density.., fill = set), color = "black", position = "dodge", stat = "bin", binwidth = 1) +
      scale_y_continuous(breaks=waiver(), expand=c(0,0)) +
      theme_bw() + ggtitle(paste0("Score distribution for TFBS ", TFBSForWhichToScore[TFBS])) +
      theme(plot.title = element_text(hjust = 0.5))
    
    densityPlot   <- ggplot(data = dataToUse, aes(x = scores)) +
      geom_density(aes(y = ..density.., fill = set), color = "black",  alpha = 0.3) +
      scale_y_continuous(breaks=waiver()) +
      theme_bw() + ggtitle(paste0("Score distribution for TFBS ", TFBSForWhichToScore[TFBS])) +
      theme(plot.title = element_text(hjust = 0.5))
    
    toAdd                         <- list(dataUsedForPlot = dataToUse, histogram = histogramPlot, densityPlot = densityPlot) 
    listPlotsPerTFBS[[TFBS]]      <- toAdd
    names(listPlotsPerTFBS)[TFBS] <- TFBSForWhichToScore[TFBS]
  }
  
  
  
  
  
  returnList     <- list(plottingData            = finalPlotDfHist,
                         histogramPlotTotalScore = histogramPlotTotalScore,
                         densityPlotTotalScore   = densityPlotTotalScore,
                         plotsAndDataPerTFBS     = listPlotsPerTFBS)
}

######################################################################################################
######################################################################################################
######                                                                                          ######
######                      Generate and save output data (gene lists etc.)                     ######
######                                                                                          ######
######################################################################################################
######################################################################################################


#Note: this function has not been streamlined so much as of yet
outputEnsIdsWithTFBS                               <- function(dataStructureTFBS = dataStructureTFBS,
                                                               finalScoresWholeGenome = finalScoresWholeGenome) {
  
  
  listIds <- vector("list", length = nrow(finalScoresWholeGenome$informativeTFBS))
  for(relevantTFBS in seq_along(finalScoresWholeGenome$informativeTFBS$TFBS)) {
    currentTFBS       <- finalScoresWholeGenome$informativeTFBS$TFBS[relevantTFBS]
    wGDRelevantGenes <-  dataStructureTFBS$tFTableMultipleLines[dataStructureTFBS$tFTableMultipleLines$TranscriptionFactorBindingSites %in%  currentTFBS,]
    
    #get the EnsemblIds for the current transcription factor from the data
    #count the amount that have it per TFBS
    Ids <-  wGDRelevantGenes$EnsemblId
    listIds[[currentTFBS]] <- Ids
    lengthIds <- length(Ids)
    names(lengthIds) <- currentTFBS
    listIds[["amount"]] <- append(listIds[["amount"]], lengthIds)
    
  }
  
  #now construct a dataframe with
  neededLength <- max(listIds$amount)
  
  #loop through the transcription factors, get their EnsIDs, and add "" as necessary for dataframe construction
  for(relevantTFBS in seq_along(finalScoresWholeGenome$informativeTFBS$TFBS)) {
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
  write_csv(dfWGDEnsmblIdsPerTFBS, paste0(dirForRun, "/EnsemblIdsWithPositiveSetEnrichedGeneralisedTFBSPouya.csv"))
  
  #return the df for checking in R
  return (dfWGDEnsmblIdsPerTFBS)
  
  
}

outputPerPositiveSetGeneTFBSAndGeneCharacteristics <- function(dataStructureTFBS      = dataStructureTFBS,
                                                               finalScoresWholeGenome = finalScoresWholeGenome) {
  
  
  Takegenes <- dataStructureTFBS$tFTableMultipleLines %>% filter(EnsemblId %in% dataStructureTFBS$positiveSetTable$EnsemblId) 
  
  #now, for every relevant transcription factor, return the genes that have it
  genelist = list()
  for(relevantTFBS in 1:length(finalScoresWholeGenome$informativeTFBS$TFBS)) {
    
    positiveGenesWithTFBS         <- Takegenes[finalScoresWholeGenome$informativeTFBS$TFBS[relevantTFBS] == as.vector(Takegenes$TFBS),]
    ensemblIds                    <- positiveGenesWithTFBS$EnsemblId
    genelist[[relevantTFBS]]      <- dataStructureTFBS$positiveSetTable %>% filter(EnsemblId  %in%  ensemblIds)
    names(genelist)[relevantTFBS] <- finalScoresWholeGenome$informativeTFBS$TFBS[relevantTFBS]
    genelist[[relevantTFBS]]      <- cbind(TFBS = rep(finalScoresWholeGenome$informativeTFBS$TFBS[relevantTFBS],nrow(genelist[[relevantTFBS]])),
                                           genelist[[relevantTFBS]] )
    if (relevantTFBS == 1) {
      
      write.csv(genelist[[relevantTFBS]], paste0(dirForRun, "/PositiveSetGenesAndAssociatedGeneralisedTFBSPouyaSet.csv"))
      
    } else {
      
      write.table(genelist[[relevantTFBS]], paste0(dirForRun, "/PositiveSetGenesAndAssociatedGeneralisedTFBSPouyaSet.csv"),
                  append = TRUE, sep = ",", col.names = FALSE)
    }
  }
}





######################################################################################################
######################################################################################################
######                                                                                          ######
######                                 Execution of the programme                               ######
######                                                                                          ######
######################################################################################################
######################################################################################################





#data loading
dataStructureTFBS      <- loadLibraryAndData()
#non-crossval run
fisherTests            <- computeFishersTest(dataStructureTFBS$perMotifCountsFullSet,
                                             dataStructureTFBS$perMotifCountsPositiveSet,
                                             runNewTest = TRUE, debug = FALSE)
informativeTFBS        <- fisherTests %>% filter(qValue <= 0.10 & positiveSetCount > 3)


newTFBSData <- data.frame(TFBS = unique(dataStructureTFBS$tFTableMultipleLines$TFBS),
                          Pai = 0, Pbi = 0, avgCountsMHC = 0)

calculatePaiPaibScores <- do.call(rbind, purrr::map(newTFBSData$TFBS, .f = function(nameOfCurrentTFBS, dataStrucTFBS = dataStructureTFBS) {
  
  
  Pai <- nrow(dataStrucTFBS$tFTableMultipleLines %>% filter(EnsemblId %in% dataStrucTFBS$MHCGeneTable$EnsemblId & TFBS == nameOfCurrentTFBS)) / nrow(dataStrucTFBS$tFTableMultipleLines %>% filter(EnsemblId %in% dataStrucTFBS$MHCGeneTable$EnsemblId) %>% distinct(EnsemblId, .keep_all = TRUE))
  Pbi <- nrow(dataStrucTFBS$tFTableMultipleLines %>% filter(EnsemblId %ni% dataStrucTFBS$MHCGeneTable$EnsemblId & TFBS == nameOfCurrentTFBS)) / nrow(dataStrucTFBS$tFTableMultipleLines %>% filter(EnsemblId %ni% dataStrucTFBS$MHCGeneTable$EnsemblId) %>% distinct(EnsemblId, .keep_all = TRUE))
  
  
  
  totalCountsMHC           <- ifelse(nameOfCurrentTFBS %in% (dataStrucTFBS$tFTableMultipleLines %>% filter(EnsemblId %in% dataStrucTFBS$MHCGeneTable$EnsemblId) %>% select(TFBS) %>% `[`(,1)),
                                     dataStrucTFBS$tFTableMultipleLines %>% filter(EnsemblId %in% dataStrucTFBS$MHCGeneTable$EnsemblId & TFBS == nameOfCurrentTFBS) %>% select(Count) %>% sum(), 0)
  avgCountsMHC             <- totalCountsMHC/dataStrucTFBS$MHCGenesInDataset
  
  returnVector <- c(Pai, Pbi, avgCountsMHC)
  names(returnVector) <- c("Pai", "Pbi", 'avgCountsMHC')
  returnVector
  
  
}))
newTFBSData[,2:4] <- calculatePaiPaibScores
newTFBSData$ratioPaiPbi <- newTFBSData$Pai/newTFBSData$Pbi
newTFBSData$absentratioPaiPbi <- (1-newTFBSData$Pai)/(1-newTFBSData$Pbi)
newTFBSData

finalScoresWholeGenomeNonCountsCorrected <- scoresOverRows(dataStrucTFBS = dataStructureTFBS, wholeGenomeDataOneLine = dataStructureTFBS$tFTableOneLine, fisherTests = fisherTests, positiveSetGenesInData = dataStructureTFBS$MHCGenesInDataset, 
                                         runBayesFactor = TRUE, includeCounts = FALSE, scoreStruc = newTFBSData)

finalScoresWholeGenomeCountsCorrected   <- scoresOverRows(dataStrucTFBS = dataStructureTFBS, wholeGenomeDataOneLine = dataStructureTFBS$tFTableOneLine, fisherTests = fisherTests, positiveSetGenesInData = dataStructureTFBS$MHCGenesInDataset, 
                                                          runBayesFactor = TRUE, includeCounts = TRUE)

##do command line
##do 


#make a dataframe that holds per TFBS the Pai and Pbi and the correction


#diagnostics
str(dataStructureTFBS)
dataStructureTFBS$tFTableOneLine
dataStructureTFBS$tFTableMultipleLines
dataStructureTFBS$cvPositiveSetNotTrainedOnList[[1]]
dataStructureTFBS$cvNegativeSetNotTrainedOnList[[2]]
str(fisherTests)
fisherTests
str(informativeTFBS)
informativeTFBS
str(finalScoresWholeGenome)
finalScoresWholeGenome$scoresWholeGenome
finalScoresWholeGenome$informativeTFBS

#plotting
positiveSetGenesTFBSPlot <- plotTFBSPresencePositiveSetGenes(dataStructureTFBS = dataStructureTFBS, informativeTFBS = informativeTFBS)
positiveSetGenesTFBSPlot

barPlotDataTFBS <- createBarPlotData(TFBSdata = dataStructureTFBS, informativeTFBS = informativeTFBS)
str(barPlotDataTFBS)
barPlotDataTFBS
barPlotTFBS <- plotBarPlotData(barPlotDataTFBS)
barPlotTFBS

densityAndHistogramPlot <- plotScoreDistributions(dataStructureTFBS = dataStructureTFBS,
                                                  finalScoresWholeGenome = finalScoresWholeGenome$scoresWholeGenome,
                                                  debug = FALSE)
densityAndHistogramPlot$plottingData
densityAndHistogramPlot$histogramPlot
densityAndHistogramPlot$densityPlot
#crossval run
crossValFisherTestData    <- crossValFisherTests(dataStructureTFBS$cvTotalSetTFBSCountsList, dataStructureTFBS$cvPositiveSetTFBSCountsList)
#gives a strange error that I do not understand??
crossValScoresWholeGenomeData <- crossValScoresOverRows(crossValWholeGenomeDataOneLineList = dataStructureTFBS$cvTotalSetNotTrainedOnWGDOL,
                                                        crossValFisherTestList             = crossValFisherTestData,
                                                        crossValPositiveSets               = dataStructureTFBS$cvPositiveSetTrainedOnList)

crossValScoresOnly <- purrr::map(crossValScoresWholeGenomeData, .f = function(partialScoreDataTable) {
  partialScoreDataTable$scoresWholeGenome
})
totalCrossValScoresAllGenes <- rbindlist(crossValScoresOnly, fill = TRUE)




#diagnostics
str(crossValFisherTestData)
str(crossValFisherTestData[[1]])
crossValFisherTestData[[1]]
crossValFisherTestData[[6]]
crossValFisherTestData[[10]]


str(crossValScoresWholeGenome)
str(crossValScoresWholeGenome[[1]])
str(crossValScoresOnly)
crossValScoresOnly[[1]]
crossValScoresOnly[[5]]
str(totalCrossValScoresAllGenes)
head(totalCrossValScoresAllGenes)


#plotting

plotCrossValHistogramAndDensityPlot <- plotScoreDistributions(dataStructureTFBS = dataStructureTFBS,
                                                              finalScoresWholeGenome = totalCrossValScoresAllGenes
)
plotCrossValHistogramAndDensityPlot$plottingData
plotCrossValHistogramAndDensityPlot$histogramPlotTotalScore
plotCrossValHistogramAndDensityPlot$densityPlotTotalScore
plotCrossValHistogramAndDensityPlot$plotsAndDataPerTFBS$GSC$histogram
plotCrossValHistogramAndDensityPlot$plotsAndDataPerTFBS$NFKB
plotCrossValHistogramAndDensityPlot$plotsAndDataPerTFBS$RAD21
plotCrossValHistogramAndDensityPlot$plotsAndDataPerTFBS$IRF
#close log file connection
close(fileConn)

