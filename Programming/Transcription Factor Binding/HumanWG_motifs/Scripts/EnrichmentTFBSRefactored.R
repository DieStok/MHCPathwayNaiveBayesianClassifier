

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
pacman::p_load(tidyverse, stringr, purrr, readr, LaF, ff, data.table, R.utils, base, utils, reshape2, optparse, cowplot, svglite)
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
MHCIPositiveSetPath   <- "/home/dieter/Documents/Project/Programming/DataForAll/PositiveList_MHCIOnly.csv"
MHCIIPositiveSetPath  <- "/home/dieter/Documents/Project/Programming/DataForAll/PositiveList_MHCIIOnly.csv"
multipleLinesTFBSFile <- "~/Documents/Project/Programming/Transcription Factor Binding/HumanWG_motifs/Data/singleLineFinalDf and Unified table/underscoresRemovedMultipleLinesPouya.csv"
singleLineTFBSFile    <- "~/Documents/Project/Programming/Transcription Factor Binding/HumanWG_motifs/Data/singleLineFinalDf and Unified table/underscoresRemovedSingleLinePouya.csv"

#if you do not wish to load data from source files, input the .rds file locations here
TFBSDataObjectRDS   <- NA
fisherDataObjectRDS <- NA







#Sanity checks

singleLineTFBSFileLoaded <- fread(singleLineTFBSFile)
head(singleLineTFBSFileLoaded)
nrow(singleLineTFBSFileLoaded)
#aha, so not all genes are in there











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


loadLibraryAndDataRevamped     <- function(foldCrossVal = 10, entityName = "TFBS", UI = FALSE, reloadFromSource = TRUE, positiveSetFile = positiveSetFile) {
  
  
  
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
    #saveRDS(returnList, file = paste0(dirForRun, "/TFBSDataObject.rds"))
    
    
    
    
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

computeFishersTestRevamped     <- function(wholeGenomeTFBSCountData, positiveSetTFBSCountData, wholeGenomeMultipleLinesTFBSData, positiveSetEnsemblIDs, runNewTest = TRUE, debug = FALSE) {
  
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
    
    #select from the whole genome data only the positive set genes. Need these later for calculating medium counts of TFBS
    #in Positive Set genes that have them
    wholeGenomeDataMultipleLinesPositiveSetOnly <-  wholeGenomeMultipleLinesTFBSData %>% filter(EnsemblId %in% positiveSetEnsemblIDs)
    
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
    fisherDataFrame <- as.data.table(matrix(ncol = 9, nrow = length(unique(relevantWholeGenome$TFBS))))
    colnames(fisherDataFrame) <- c("TFBS", "pValue", "sampleEstimate", "lowerConfidenceInterval", "upperConfidenceInterval",
                                   "wholeGenomeCount", "positiveSetCount", "medianInPositiveSetGenesThatHaveThisTFBS", 
                                   "sdInPositiveSetGenesThatHaveThisTFBS")
    fisherDataFrame$TFBS <- rep("NA", nrow(fisherDataFrame))
    fisherDataFrame[, colnames(fisherDataFrame)[-1] := lapply(.SD, as.numeric), .SDcols = colnames(fisherDataFrame)[-1]]
    fisherDataFrame
    if(debug == TRUE) {
      print(class(fisherDataFrame$TFBS))
      print(class(fisherDataFrame$pValue))
    }
  

    #cycle over all the rows of the TFBS found in the positive set, and determine the enrichment
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
      #select only those genes in the positive set that have the current TFBS
      relevantPositiveSetGenesForMedian <- wholeGenomeDataMultipleLinesPositiveSetOnly %>% filter(TFBS == positiveSetTFBSCountData$TFBS[bindingSites])
      #calculate the median count and the SD of the count of this TFBS in the positive set
      medianCountInTFBSContainingPositiveSetGenes <- median(relevantPositiveSetGenesForMedian$Count)
      sdCountInTFBSContainingPositiveSetGenes     <- sd(relevantPositiveSetGenesForMedian$Count)
      
      listMatrices[[bindingSites]] <- fisherMatrix
      addVector        <- c(as.character(relevantWholeGenome$TFBS[bindingSites]),
                            pValue,
                            sampleEst,
                            lowerConfInt,
                            upperConfInt,
                            wholeGenomePresence,
                            positiveSetPresence,
                            medianCountInTFBSContainingPositiveSetGenes,
                            sdCountInTFBSContainingPositiveSetGenes
                            )
      names(addVector) <- c("TFBS", "pValue", "sampleEstimate",
                            "lowerConfidenceInterval", "upperConfidenceInterval",
                            "wholeGenomeCount", "positiveSetCount", "medianInPositiveSetGenesThatHaveThisTFBS",
                            "sdInPositiveSetGenesThatHaveThisTFBS")
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

crossValFisherTestsRevamped    <- function(crossValTotalSetTFBSCounts, crossValPositiveTFBSCountList) {
  
  result <- purrr::pmap(list(crossValTotalSetTFBSCounts,crossValPositiveTFBSCountList), .f = computeFishersTestRevamped,
                        wholeGenomeMultipleLinesTFBSData = dataStructureTFBS$tFTableMultipleLines, positiveSetEnsemblIDs = dataStructureTFBS$positiveSetTable$EnsemblId)
}

#scoring over rows functions
scoresOverRows         <- function(wholeGenomeDataOneLine, fisherTests, positiveSetGenesInData, debug = FALSE) {
  
  #find TFBS that are significantly overrepresented in positive set genes
  #assign them a score. In this case, the log2 point estimate from the Fisher tests.
  
  
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
    
    #logical vector that says per informativeTFBS whether it is present in the current gene
    presenceOfInformativeTFBS <- informativeTFBS$TFBS %in% namesTFBS
    #initialse the score. Note the +1 for the totalScore column
    scoreCompendium           <- vector("numeric", length = (length(presenceOfInformativeTFBS) + 1))
    
    for(particularTFBS in 1:length(presenceOfInformativeTFBS)) {
      
    if(presenceOfInformativeTFBS[particularTFBS] == TRUE) {
      
      #if it is there but the score is negative, keep as is. If the score is positive, correct to the average positive set amount
    
      if(informativeTFBS$sampleEstimate[particularTFBS] <= 1) {
        
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
  print("These are their scores:")
  print(informativeTFBS$score)
  
  
  list(
    scoresWholeGenome = wholeGenomeScores,
    informativeTFBS   = informativeTFBS)
  
}

scoresOverRowsRevamped         <- function(wholeGenomeDataOneLine, fisherTests, positiveSetGenesInData, debug = FALSE) {
  
  #find TFBS that are significantly overrepresented in positive set genes
  #assign them a score. In this case, the log2 point estimate from the Fisher tests.
  
  
  informativeTFBS                    <- subset(fisherTests, qValue <= 0.10 & positiveSetCount > 3)
  informativeTFBS                    <- informativeTFBS[order(informativeTFBS$qValue),]
  #rownames(informativeTFBS)         <- seq(1:nrow(informativeTFBS))
  informativeTFBS$score              <- log2(informativeTFBS$sampleEstimate)
  
  
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
                                                                  
                                                                  #logical vector that says per informativeTFBS whether it is present in the current gene
                                                                  presenceOfInformativeTFBS <- informativeTFBS$TFBS %in% namesTFBS
                                                                  #initialse the score. Note the +1 for the totalScore column
                                                                  scoreCompendium           <- vector("numeric", length = (length(presenceOfInformativeTFBS) + 1))
                                                                  
                                                                  for(particularTFBS in 1:length(presenceOfInformativeTFBS)) {
                                                                    
                                                                    if(presenceOfInformativeTFBS[particularTFBS] == TRUE) {
                                                                      
                                                                      #if it is there but the score is negative, keep as is. If the score is positive, correct to the average positive set amount
                                                                      
                                                                      if(informativeTFBS$sampleEstimate[particularTFBS] <= 1) {
                                                                        
                                                                        valueToAdd <- informativeTFBS$score[particularTFBS]
                                                                      } else {
                                                                        
                                                                        particularTFBSCount  <- countsTFBS[informativeTFBS$TFBS[particularTFBS]]
                                                                        #als gemiddeld aantal keer deze TFBS >1 in MHC, corrigeer dan het te scoren gen voor het aantal keer dat het deze TFBS heeft
                                                                        
                                                                          #normal score
                                                                          valueToAdd <- informativeTFBS$score[particularTFBS] *
                                                                            #correction for deviation from median of genes in positive set that have this TFBS
                                                                            (min(informativeTFBS$medianInPositiveSetGenesThatHaveThisTFBS[particularTFBS],particularTFBSCount)/max(informativeTFBS$medianInPositiveSetGenesThatHaveThisTFBS[particularTFBS],particularTFBSCount)) *
                                                                            #correction for the sd (larger sd = less weight should be put on different motif count in gene to be scored)
                                                                            #take the log2 to make more distinction in smaller scores
                                                                            log2(1/informativeTFBS$sdInPositiveSetGenesThatHaveThisTFBS[particularTFBS])
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
  print("These are their scores:")
  print(informativeTFBS$score)
  
  
  list(
    scoresWholeGenome = wholeGenomeScores,
    informativeTFBS   = informativeTFBS)
  
}

crossValScoresOverRows <- function(crossValWholeGenomeDataOneLineList, crossValFisherTestList, crossValPositiveSets) {
  
  positiveSetGenesInData <- lapply(crossValPositiveSets, FUN = function(x) {length(unique(x$EnsemblId))})
  result <- purrr::pmap(list(crossValWholeGenomeDataOneLineList, crossValFisherTestList, positiveSetGenesInData), .f = scoresOverRows)
  result
  
}


crossValScoresOverRowsRevamped <- function(crossValWholeGenomeDataOneLineList, crossValFisherTestList, crossValPositiveSets) {
  
  positiveSetGenesInData <- lapply(crossValPositiveSets, FUN = function(x) {length(unique(x$EnsemblId))})
  result <- purrr::pmap(list(crossValWholeGenomeDataOneLineList, crossValFisherTestList, positiveSetGenesInData), .f = scoresOverRowsRevamped)
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
                                                 qValue = rep(finalScoresWholeGenome$informativeTFBS[finalScoresWholeGenome$informativeTFBS$TFBS == infTFBS[significantTFBS],]$qValue,4))
    
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
                                             debug = FALSE, bw = 2, title = "Distribution of total TFBS score in positive set (MHC genes)\n and negative set (all genes in the human genome)") {
  
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
      geom_histogram(aes(y = ..density.., fill = set), color = "black", position = "dodge", stat = "bin", binwidth = bw) +
      scale_y_continuous(breaks=waiver(), expand=c(0,0)) +
      theme_bw() + ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5))  
    
    densityPlotTotalScore    <- ggplot(data = totalScoreData, aes(x = scores)) +
      geom_density(aes(y = ..density.., fill = set), color = "black",  alpha = 0.3) +
      scale_y_continuous(breaks=waiver()) +
      theme_bw() + ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5))
    
    #now do this also for all the separate TFBS
    TFBSForWhichToScore <- unique(finalPlotDfHist$ScoreFor)
    listPlotsPerTFBS    <- vector("list", length(TFBSForWhichToScore))
    for(TFBS in seq_along(TFBSForWhichToScore)) {
      
      dataToUse     <- finalPlotDfHist %>% dplyr::filter(ScoreFor == TFBSForWhichToScore[TFBS])
      sequenceOfBreaks <- seq(floor(min(dataToUse$scores, na.rm = TRUE)), ceiling(max(dataToUse$scores, na.rm = TRUE)), by = 1)
      #give at least 2 breaks
      if(length(sequenceOfBreaks) <= 1) {sequenceOfBreaks = c(0,1)}
      print(sequenceOfBreaks)
      dataToUseSummarized <- dataToUse %>% group_by(gr = cut(scores, breaks = sequenceOfBreaks), set) %>%
        summarize(Count = n()) %>% group_by(set) %>% mutate(densityPerSet = Count/sum(Count))
      
      histogramPlot <- ggplot(dataToUseSummarized, aes(x = as.factor(gr), y = densityPerSet, fill = set)) +
        geom_bar(stat = "identity", colour = "black", width = 1, position = "dodge") +
        theme_bw() + theme(panel.grid.major.x = element_blank()) +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_x_discrete(expand   = c(0,0)) +
        scale_y_continuous(expand = c(0,0), limits = c(0,1), breaks = seq(0,1,0.1)) +
        ggtitle(label = paste0("Score distribution for TFBS ", TFBSForWhichToScore[TFBS]))
      
      
      # histogramPlot <- ggplot(data = dataToUse, aes(x = scores)) +
      #   geom_histogram(aes(y = ..density.., fill = set), color = "black", position = "dodge", stat = "bin", binwidth = bw) +
      #   scale_y_continuous(breaks=waiver(), expand=c(0,0), limits = c(0,1)) +
      #   theme_bw() + ggtitle(paste0("Score distribution for TFBS ", TFBSForWhichToScore[TFBS])) +
      #   theme(plot.title = element_text(hjust = 0.5))
      
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




#####################################
#       MHC I and II PS combined    #
#####################################

#data loading
dataStructureTFBS      <- loadLibraryAndData()
#non-crossval run
fisherTests            <- computeFishersTestRevamped(dataStructureTFBS$perMotifCountsFullSet,
                                             dataStructureTFBS$perMotifCountsPositiveSet, wholeGenomeMultipleLinesTFBSData = dataStructureTFBS$tFTableMultipleLines,
                                             positiveSetEnsemblIDs = dataStructureTFBS$positiveSetTable$EnsemblId,
                                             runNewTest = TRUE, debug = FALSE)
informativeTFBS        <- fisherTests %>% filter(qValue <= 0.10 & positiveSetCount > 3)
finalScoresWholeGenome <- scoresOverRowsRevamped(dataStructureTFBS$tFTableOneLine, fisherTests, dataStructureTFBS$MHCGenesInDataset)

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
                                                  debug = FALSE, bw = 4.5)
densityAndHistogramPlot$plottingData
densityAndHistogramPlot$histogramPlot
densityAndHistogramPlot$densityPlot

#crossval run
crossValFisherTestData    <- crossValFisherTestsRevamped(dataStructureTFBS$cvTotalSetTFBSCountsList, dataStructureTFBS$cvPositiveSetTFBSCountsList)
crossValScoresWholeGenomeData <- crossValScoresOverRows(crossValWholeGenomeDataOneLineList = dataStructureTFBS$cvTotalSetNotTrainedOnWGDOL,
                                                        crossValFisherTestList             = crossValFisherTestData,
                                                        crossValPositiveSets               = dataStructureTFBS$cvPositiveSetTrainedOnList)

crossValScoresOnly <- purrr::map(crossValScoresWholeGenomeData, .f = function(partialScoreDataTable) {
  partialScoreDataTable$scoresWholeGenome
})
totalCrossValScoresAllGenes <- rbindlist(crossValScoresOnly, fill = TRUE)


#crossValScoresOnly is per test set all the TFBS that contributed to the score in that run, and the totalScore 
#for the EnsemblIDs tested in that run

#totalCrossValScoresAllGenes gives the total scores obtained with cross-valuation and fills in missing ones. That means that
#all TFBS that have a score in at least one crossVal are in there, and the total scores for each gene are too.

#see in how many crossVals a certain TFBS is found as informative (i.e., trend (qvalue<0.1) in Fisher test after multiple testing correction)
#The below displays in how many Crossvals a TFBS is found as informative. more than 5/10 is what I will call a 'core motif'
totalGenesWithScore <- nrow(totalCrossValScoresAllGenes)
howManyCrossValsInformativeTFBS <-  lapply(totalCrossValScoresAllGenes, FUN = function(element) {round(sum(!is.na(element))/totalGenesWithScore*10,0)}) %>% unlist()
howManyCrossValsInformativeTFBS <- sort(howManyCrossValsInformativeTFBS, decreasing = TRUE)
howManyCrossValsInformativeTFBS <- howManyCrossValsInformativeTFBS[-which(names(howManyCrossValsInformativeTFBS) == "totalScore")]
howManyCrossValsInformativeTFBS <- howManyCrossValsInformativeTFBS[-which(names(howManyCrossValsInformativeTFBS) == "EnsemblId")]
howManyCrossValsInformativeTFBS

dataFrameTFBSMHC <- data.frame(Name = names(howManyCrossValsInformativeTFBS), CrossValsInformative = howManyCrossValsInformativeTFBS)
dataFrameTFBSMHC
dataFrameTFBSMHC %>% write_csv("/home/dieter/Documents/Project/presentaties/finalTalk/TFBSDF.csv")



#also get median odds ratio

totalCrossValScoresAllGenes %>% select(-EnsemblId, -totalScore) %>% sapply(median)

#check that not all scores are 0 in a column
colSums(totalCrossValScoresAllGenes[, -c("EnsemblId", "totalScore")], na.rm = TRUE)


#diagnostics
str(crossValFisherTestData)
str(crossValFisherTestData[[1]])
crossValFisherTestData[[1]]
crossValFisherTestData[[6]]
crossValFisherTestData[[10]]


str(crossValScoresWholeGenomeData)
str(crossValScoresWholeGenomeData[[1]])
str(crossValScoresOnly)
crossValScoresOnly[[1]]
crossValScoresOnly[[5]]
str(totalCrossValScoresAllGenes)
head(totalCrossValScoresAllGenes)
nrow(totalCrossValScoresAllGenes)


#plotting

plotCrossValHistogramAndDensityPlot <- plotScoreDistributions(dataStructureTFBS = dataStructureTFBS,
                                                              finalScoresWholeGenome = totalCrossValScoresAllGenes
)
plotCrossValHistogramAndDensityPlot$plottingData
plotCrossValHistogramAndDensityPlot$histogramPlotTotalScore
plotCrossValHistogramAndDensityPlot$histogramPlotTotalScore %>% ggsave(filename = paste0(dirForRun, "/CrossValHistogramPouyaGeneralisedMHCI-IICombinedPS.svg"),
                                                                       dpi = 600)

plotCrossValHistogramAndDensityPlot$histogramPlotTotalScore + scale_fill_manual(breaks = c("positive","negative"),
                                                                                values = c("#686caa", "#eabc4f"),
                                                                                labels = c("PS", "NS"))
plotCrossValHistogramAndDensityPlot$densityPlotTotalScore
plotCrossValHistogramAndDensityPlot$densityPlotTotalScore %>% ggsave(filename = paste0(dirForRun, "/CrossValDensityPlotPouyaGeneralisedMHCI-IICombinedPS.svg"),
                                                                     dpi = 600)


#positiveSetGenesTFBSPlot <- plotTFBSPresencePositiveSetGenes(dataStructureTFBS = dataStructureTFBS, informativeTFBS = informativeTFBS)
#positiveSetGenesTFBSPlot


######################################################################
###                                                                ###
### Output of data file used in ScriptForGeneratingBayesianCSV.R   ###
###                         (MHC I and II)                         ###
######################################################################                                                

#here, output the totalScores of the TFBS for use in the .csv used for the final classifier
TFBSScoreOutput <- totalCrossValScoresAllGenes %>% select(EnsemblId, totalScore, everything())
head(TFBSScoreOutput)
nrow(TFBSScoreOutput)
sum(TFBSScoreOutput$EnsemblId %in% dataStructureTFBS$positiveSetTable$EnsemblId)
write_csv(TFBSScoreOutput, path = paste0("~/Documents/Project/Programming/Transcription Factor Binding/HumanWG_motifs/Data/singleLineFinalDf and Unified table/DataForBayesianCSVMHCIandIIPositiveSet",
                                         gsub(" ", "_", Sys.time(), fixed = TRUE), ".csv"))


##########################
#       MHC I PS only    #
##########################

#data loading for MHCI only positive set

dataStructureTFBSMHCI <- loadLibraryAndDataRevamped(positiveSetFile = MHCIPositiveSetPath)
dataStructureTFBSMHCI$positiveSetTable
dataStructureTFBSMHCI$perMotifCountsPositiveSet


#crossval run
crossValFisherTestDataMHCI    <- crossValFisherTestsRevamped(dataStructureTFBSMHCI$cvTotalSetTFBSCountsList, dataStructureTFBSMHCI$cvPositiveSetTFBSCountsList)
crossValScoresWholeGenomeDataMHCI <- crossValScoresOverRows(crossValWholeGenomeDataOneLineList = dataStructureTFBSMHCI$cvTotalSetNotTrainedOnWGDOL,
                                                        crossValFisherTestList                 = crossValFisherTestDataMHCI,
                                                        crossValPositiveSets                   = dataStructureTFBSMHCI$cvPositiveSetTrainedOnList)

crossValScoresOnlyMHCI <- purrr::map(crossValScoresWholeGenomeDataMHCI, .f = function(partialScoreDataTable) {
  partialScoreDataTable$scoresWholeGenome
})
totalCrossValScoresAllGenesMHCI <- rbindlist(crossValScoresOnlyMHCI, fill = TRUE)


#diagnostics
str(crossValFisherTestDataMHCI)
str(crossValFisherTestDataMHCI[[1]])
crossValFisherTestDataMHCI[[1]]
crossValFisherTestDataMHCI[[6]]
crossValFisherTestDataMHCI[[10]]


str(crossValScoresWholeGenomeDataMHCI)
str(crossValScoresWholeGenomeDataMHCI[[1]])
str(crossValScoresOnlyMHCI)
crossValScoresOnlyMHCI[[1]]
crossValScoresOnlyMHCI[[5]]
str(totalCrossValScoresAllGenesMHCI)
head(totalCrossValScoresAllGenesMHCI)
nrow(totalCrossValScoresAllGenesMHCI)

#plotting
plotCrossValHistogramAndDensityPlotMHCI <- plotScoreDistributions(dataStructureTFBS = dataStructureTFBSMHCI,
                                                              finalScoresWholeGenome = totalCrossValScoresAllGenesMHCI,
                                                              bw = 1,
                                                              title = "Distribution of total TFBS score in positive set (MHC I genes)\n and negative set (all genes in the human genome)"
)
plotCrossValHistogramAndDensityPlotMHCI$plottingData
plotCrossValHistogramAndDensityPlotMHCI$histogramPlotTotalScore
plotCrossValHistogramAndDensityPlotMHCI$histogramPlotTotalScore %>% ggsave(filename = paste0(dirForRun, "/CrossValHistogramPouyaGeneralisedMHCIPS.svg"),
                                                                       dpi = 600)

plotCrossValHistogramAndDensityPlotMHCI$histogramPlotTotalScore + scale_fill_manual(breaks = c("positive","negative"),
                                                                                values = c("#686caa", "#eabc4f"),
                                                                                labels = c("PS", "NS"))
plotCrossValHistogramAndDensityPlotMHCI$densityPlotTotalScore
plotCrossValHistogramAndDensityPlotMHCI$densityPlotTotalScore %>% ggsave(filename = paste0(dirForRun, "/CrossValDensityPlotPouyaGeneralisedMHCIPS.svg"),
                                                                     dpi = 600)

######################################################################
###                                                                ###
### Output of data file used in ScriptForGeneratingBayesianCSV.R   ###
###                            (MHC I)                             ###
######################################################################                                                

#here, output the totalScores of the TFBS for use in the .csv used for the final classifier
TFBSScoreOutputMHCI <- totalCrossValScoresAllGenesMHCI %>% select(EnsemblId, totalScore, everything())
head(TFBSScoreOutputMHCI)
nrow(TFBSScoreOutputMHCI)
sum(TFBSScoreOutputMHCI$EnsemblId %in% dataStructureTFBSMHCI$positiveSetTable$EnsemblId)
write_csv(TFBSScoreOutputMHCI, path = paste0("~/Documents/Project/Programming/Transcription Factor Binding/HumanWG_motifs/Data/singleLineFinalDf and Unified table/DataForBayesianCSVMHCIPositiveSet",
                                         gsub(" ", "_", Sys.time(), fixed = TRUE), ".csv"))



##########################
#       MHC II PS only    #
##########################

#data loading for MHCII only positive set

dataStructureTFBSMHCII <- loadLibraryAndDataRevamped(positiveSetFile = MHCIIPositiveSetPath)
dataStructureTFBSMHCII$positiveSetTable
dataStructureTFBSMHCII$perMotifCountsPositiveSet


#crossval run
crossValFisherTestDataMHCII    <- crossValFisherTestsRevamped(dataStructureTFBSMHCII$cvTotalSetTFBSCountsList, dataStructureTFBSMHCII$cvPositiveSetTFBSCountsList)
crossValScoresWholeGenomeDataMHCII <- crossValScoresOverRows(crossValWholeGenomeDataOneLineList = dataStructureTFBSMHCII$cvTotalSetNotTrainedOnWGDOL,
                                                            crossValFisherTestList                 = crossValFisherTestDataMHCII,
                                                            crossValPositiveSets                   = dataStructureTFBSMHCII$cvPositiveSetTrainedOnList)

crossValScoresOnlyMHCII <- purrr::map(crossValScoresWholeGenomeDataMHCII, .f = function(partialScoreDataTable) {
  partialScoreDataTable$scoresWholeGenome
})
totalCrossValScoresAllGenesMHCII <- rbindlist(crossValScoresOnlyMHCII, fill = TRUE)


#diagnostics
str(crossValFisherTestDataMHCII)
str(crossValFisherTestDataMHCII[[1]])
crossValFisherTestDataMHCII[[1]]
crossValFisherTestDataMHCII[[6]]
crossValFisherTestDataMHCII[[10]]


str(crossValScoresWholeGenomeDataMHCII)
str(crossValScoresWholeGenomeDataMHCII[[1]])
str(crossValScoresOnlyMHCII)
crossValScoresOnlyMHCII[[1]]
crossValScoresOnlyMHCII[[5]]
str(totalCrossValScoresAllGenesMHCII)
head(totalCrossValScoresAllGenesMHCII)
nrow(totalCrossValScoresAllGenesMHCII)

#plotting
plotCrossValHistogramAndDensityPlotMHCII <- plotScoreDistributions(dataStructureTFBS = dataStructureTFBSMHCII,
                                                                  finalScoresWholeGenome = totalCrossValScoresAllGenesMHCII,
                                                                  bw = 1,
                                                                  title = "Distribution of total TFBS score in positive set (MHC II genes)\n and negative set (all genes in the human genome)"
)
plotCrossValHistogramAndDensityPlotMHCII$plottingData
plotCrossValHistogramAndDensityPlotMHCII$histogramPlotTotalScore 
plotCrossValHistogramAndDensityPlotMHCII$histogramPlotTotalScore %>% ggsave(filename = paste0(dirForRun, "/CrossValHistogramPouyaGeneralisedMHCIIPS.svg"),
                                                                           dpi = 600)

plotCrossValHistogramAndDensityPlotMHCII$histogramPlotTotalScore + scale_fill_manual(breaks = c("positive","negative"),
                                                                                    values = c("#686caa", "#eabc4f"),
                                                                                    labels = c("PS", "NS"))
plotCrossValHistogramAndDensityPlotMHCII$densityPlotTotalScore
plotCrossValHistogramAndDensityPlotMHCII$densityPlotTotalScore %>% ggsave(filename = paste0(dirForRun, "/CrossValDensityPlotPouyaGeneralisedMHCIIPS.svg"),
                                                                         dpi = 600)




######################################################################
###                                                                ###
### Output of data file used in ScriptForGeneratingBayesianCSV.R   ###
###                            (MHC II)                             ###
######################################################################                                                

#here, output the totalScores of the TFBS for use in the .csv used for the final classifier
TFBSScoreOutputMHCII <- totalCrossValScoresAllGenesMHCII %>% select(EnsemblId, totalScore, everything())
head(TFBSScoreOutputMHCII)
nrow(TFBSScoreOutputMHCII)
sum(TFBSScoreOutputMHCII$EnsemblId %in% dataStructureTFBSMHCII$positiveSetTable$EnsemblId)
write_csv(TFBSScoreOutputMHCII, path = paste0("~/Documents/Project/Programming/Transcription Factor Binding/HumanWG_motifs/Data/singleLineFinalDf and Unified table/DataForBayesianCSVMHCIIPositiveSet",
                                             gsub(" ", "_", Sys.time(), fixed = TRUE), ".csv"))




######################################################################################################
######################################################################################################
######                                                                                          ######
######                            Manual determination of score bins                            ######
######  (This is a test/relic. Actual calculations are done in ScriptForRunningTheBayesian.R)   ######
######                                                                                          ######
######################################################################################################
######################################################################################################


######################################################################################################
######                                                                                          ######
######    All code from here on was not updated for additional positive sets. This was done     ######
######    because the code is not used in final score calculations, but rather, was the first   ######
######    iteration of work on defining bins and calculating log enrichments, and fitting       ######
######    polynomials to score distributions. For the Bayesian classifier, these functions      ######
######    are handled in ScriptForRunningTheBayesian.R. Therefore, the code below is            ######
######    illustrative and deprecated, but has remained for reference purposes.                 ######
######                                                                                          ######
######################################################################################################


par(mfrow = c(2,2))
breeaks <- c(-14, -4, -2, 2, 4, 6, 8, 22)




getScores <- function(plotBuild) {
  
  logScores <- numeric()
  for(i in seq(2, nrow(plotBuild$data[[1]]), by = 2)) {
    
    logScores = c(logScores, log2(plotBuild$data[[1]][(i-1),]$density/plotBuild$data[[1]][i,]$density))
    
  }
  newPlot <- plotBuild$plot + annotate("text", x = unique(plotBuild$data[[1]]$x), y = 0.25, label = round(logScores, digits = 3), angle = 90, size = 8)
  
  list(plot = newPlot, scores = logScores)
  
}


createScoreHist <- function(brks, data = plotCrossValHistogramAndDensityPlot$plottingData) {
orgPlot <- data %>% filter(ScoreFor == "totalScore") %>% 
  ggplot(aes(x = scores, fill = set, y = ..density..)) + geom_histogram(breaks = brks, colour = "black") +
  scale_x_continuous(breaks = brks) + theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3), limits = c(0, 0.3)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
    theme( legend.justification = c(1.5,1.5), legend.position = c(0.92,0.5), legend.text = element_text(size = 14)) + panel_border(remove = TRUE) +
  #ggtitle("Distribution of total TFBS score in positive set (MHC genes)\n and negative set (all genes in the human genome)") +
  annotate("text", x = 4, y = 0.3, label = "Log2 Enrichment Scores", size = 9)
  
orgPlotBuild <- ggplot_build(orgPlot)
orgPlotBuild
h <- getScores(orgPlotBuild)
return(list(scores = h$scores, plot = ggplot_build(h$plot)))
}


orgBreaks <- c(-14, -4, -2, 0, 2, 4, 6, 8, 22)
betterBreaks <- c(-14, -4, -2, 2, 4, 6, 8, 22)
bestBreaks  <- c(-14, -4, -2, 2, 4, 6, 8, 14, 22)
oneBreaks <- (seq(-14, 22, 1))

k <- createScoreHist(brks = orgBreaks)
k$plot
f <- createScoreHist(brks = betterBreaks)
f$plot
p <- createScoreHist(brks = bestBreaks)
p$plot
cowplot::plot_grid(nrow = 2, ncol = 2, plotlist = list(k$plot$plot, f$plot$plot, p$plot$plot))
orgHist <- createScoreHist(brks = oneBreaks)
orgHist$plot


pietje <- totalCrossValScoresAllGenes %>% group_by(gr = cut(totalScore, breaks = bestBreaks)) %>%
  mutate(finalLog2Score = p$scores[gr])
unique(pietje$finalLog2Score)                                           
range(pietje$finalLog2Score)

saveForCsv <- pietje %>% select(EnsemblId, finalLog2Score) %>% ungroup() %>% arrange(EnsemblId)
saveForCsv
write_csv(saveForCsv, path = paste0("~/Documents/Project/Programming/Transcription Factor Binding/HumanWG_motifs/Data/singleLineFinalDf and Unified table/FinalLog2Scores/binningLog2Scores",gsub(" ", "_", Sys.time(), fixed = TRUE), ".csv"))

plotCrossValHistogramAndDensityPlot$plotsAndDataPerTFBS$GSC$histogram
plotCrossValHistogramAndDensityPlot$plotsAndDataPerTFBS$NFKB
plotCrossValHistogramAndDensityPlot$plotsAndDataPerTFBS$RAD21
plotCrossValHistogramAndDensityPlot$plotsAndDataPerTFBS$IRF












######################################################################################################
######################################################################################################
######                                                                                          ######
######                              Fit a polynomial to scores                                  ######
###### (Note that this is a test/relic, the actual work is done in ScriptForRunningTheBayesian) ######
######                                                                                          ######
######################################################################################################
######################################################################################################


par(mfrow = c(2,2))
positiveSetValues <- totalCrossValScoresAllGenes %>% filter(EnsemblId %in% dataStructureTFBS$positiveSetTable$EnsemblId)
negativeSetValues <- totalCrossValScoresAllGenes %>% filter(EnsemblId %ni% dataStructureTFBS$positiveSetTable$EnsemblId)



positiveSetActualScores <- positiveSetValues$totalScore
negativeSetActualScores <- negativeSetValues$totalScore

positiveSetActualScores
negativeSetActualScores
densityPositive <- density(positiveSetActualScores, bw = 2.5)
densityNegative <- density(negativeSetActualScores, bw = 2.5)
plot(densityPositive)
plot(densityNegative)

dfPositive <- data.frame(x = densityPositive$x, y = densityPositive$y)
polyP <- lm(y ~ poly(x, 20), data = dfPositive)

dfNegative <- data.frame(x=densityNegative$x, y=densityNegative$y)
polyN <- lm(y ~ poly(x, 20), data = dfNegative)


plot(densityPositive$y ~ densityPositive$x)
lines(predict(polyP) ~ densityPositive$x)

plot(densityNegative$y ~ densityNegative$x)
lines(predict(polyN) ~ densityNegative$x)


polyFunction <- function(x, polyModel) {
  
  k <- polyModel$coefficients
  print(k)
  print(k[21])
  print(x^20 * k[21])
  
  func <- k[1] + x * k[2] + x^2 * k[3] + x^3 * k[4] + x^4 * k[5] + x^5 * k[6] +
    x^6 * k[7] + x^7 * k[8] + x^8 * k[9] + x^9 * k[10] + x^10 * k[11] +
    x^11 * k[12] + x^12 * k[13] + x^13 * k[14] + x^14 * k[15] +
    x^15 * k[16] + x^16 * k[17] + x^17 * k[18] + x^18 * k[19] +
    x^19 * k[20] + x^20 * k[21]
  print(func)
  
  func
  
}


positivePolyScoresAllGenes <- predict(polyP, data.frame(x = totalCrossValScoresAllGenes$totalScore))
negativePolyScoresAllGenes <- predict(polyN, data.frame(x = totalCrossValScoresAllGenes$totalScore))
totalLog2ScoresPoly        <- log2(positivePolyScoresAllGenes/negativePolyScoresAllGenes)
totalLog2ScoresPoly
dfTotalLog2Scores <- data.frame(EnsemblId = totalCrossValScoresAllGenes$EnsemblId,
                                totalScore = totalCrossValScoresAllGenes$totalScore,
                                Log2Scores = totalLog2ScoresPoly)

dfTotalLog2Scores <- dfTotalLog2Scores %>% arrange(totalScore)
dfTotalLog2Scores[dfTotalLog2Scores$Log2Scores > 10,]


tail(dfTotalLog2Scores)
#totalPlot
plot(densityNegative$y ~ densityNegative$x, col = "red", cex = 0.2)
lines(predict(polyN) ~ densityNegative$x, col = "darkred")
points(densityPositive$y ~ densityPositive$x, col = "green", cex = 0.2)
lines(predict(polyP) ~ densityPositive$x, col = "darkgreen")
plot(Log2Scores ~ totalScore, data = dfTotalLog2Scores, col = "black")

dfTotalLog2Scores[dfTotalLog2Scores$totalScore >= 3.95 & dfTotalLog2Scores$totalScore <= 4.05,]
dfTotalLog2Scores[dfTotalLog2Scores$totalScore >= 5.95 & dfTotalLog2Scores$totalScore <= 6.05,]
dfTotalLog2Scores[dfTotalLog2Scores$totalScore >= 7.95 & dfTotalLog2Scores$totalScore <= 8.05,]

#find the local maximum that I will use as the max and min values
head(dfTotalLog2Scores, 30)
#minimal value = -3.363383, from ENSG00000167822, totalScore = -10.83280
tail(dfTotalLog2Scores, 50)
#maximal value = 2.406856, from ENSG00000205659, totalScore = 13.22579

#edit the scores
dfTotalLog2ScoresBounded <- dfTotalLog2Scores
minimalBoundGene <- dfTotalLog2ScoresBounded[dfTotalLog2ScoresBounded$EnsemblId == "ENSG00000167822",]
maximalBoundGene <- dfTotalLog2ScoresBounded[dfTotalLog2ScoresBounded$EnsemblId == "ENSG00000205659",]
dfTotalLog2ScoresBounded[dfTotalLog2ScoresBounded$totalScore < minimalBoundGene$totalScore, ]$Log2Scores <- minimalBoundGene$Log2Scores
dfTotalLog2ScoresBounded[dfTotalLog2ScoresBounded$totalScore > maximalBoundGene$totalScore, ]$Log2Scores <- maximalBoundGene$Log2Scores

plot(densityNegative$y ~ densityNegative$x, col = "red", cex = 0.2, main = "score density P and N sets")



scoreSeq <- seq(-20, 20, by = 5)
log2Seq <- numeric(length = length(scoreSeq))
#get the closest value from the list
for(i in seq_len(length(scoreSeq))) {
closestScore <- which(abs(dfTotalLog2ScoresBounded$totalScore-scoreSeq[i])==min(abs(dfTotalLog2ScoresBounded$totalScore-scoreSeq[i])))
print(closestScore)
log2Seq[i]   <- ifelse(length(closestScore) > 1, dfTotalLog2ScoresBounded[sample(closestScore, size = 1),]$Log2Scores, dfTotalLog2ScoresBounded[closestScore,]$Log2Scores)
}
log2Seq
#in the plot you need
#1. the points output from smoothening the actual data with density()
#2. the lines that are generated by the polynomial that is used to calculate the scores
#3. nice area colouration and legend
#4. plot Example scores from -10 to 20
ggPlotDataFrame <- data.frame(y = c(densityNegative$y, densityPositive$y), x = c(densityNegative$x, densityPositive$x),
                              posNeg = c(rep("negative", length(densityNegative$y)), rep("positive", length(densityPositive$y))))
ggPlotDataFramePosOnly <- ggPlotDataFrame %>% filter(posNeg == "positive") %>% mutate(prediction = predict(polyP))
ggPlotDataFrameNegOnly <- ggPlotDataFrame %>% filter(posNeg == "negative") %>% mutate(prediction = predict(polyN))
head(ggPlotDataFrame); tail(ggPlotDataFrame, 30)
polynomialDensityPlot <- ggplot(data = ggPlotDataFrame, aes(y = y, x = x, fill = factor(posNeg, levels = c("negative", "positive")))) +
  geom_point(data = ggPlotDataFramePosOnly, aes(y = y), colour = "darkgreen", size = 0.5) +
  geom_point(data = ggPlotDataFrameNegOnly, aes(y = y), colour = "darkred", size = 0.5) +
  geom_line(data = ggPlotDataFramePosOnly, aes(y = prediction), colour = "darkgreen", size = 0.5) +
  geom_line(data = ggPlotDataFrameNegOnly, aes(y = prediction), colour = "darkred", size = 0.5) +
  scale_fill_manual(breaks = c("negative", "positive"), values = c("positive" = "seagreen", "negative" = "darkred")) +
  geom_area(data = ggPlotDataFramePosOnly, aes(y = prediction, x = x),  alpha = 0.5) +
  geom_area(data = ggPlotDataFrameNegOnly, aes(y = prediction, x = x),  alpha = 0.5) + 
  guides(fill = guide_legend(title = NULL)) +
  scale_x_continuous(name = "total TFBS score", expand = c(0,0)) +
  scale_y_continuous(name = "density of values", expand = c(0,0)) +
  coord_cartesian(ylim = c(0,0.12)) +
  #ggtitle("Smoothened distributions of scores for TFBS")  +
  theme_bw() +  theme(plot.title = element_text(hjust = 0.5), legend.justification = c(1.5,1.5), legend.position = c(0.92,0.5), legend.text = element_text(size = 14)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  annotate("text", scoreSeq, 0.09, label = round(log2Seq, 3), angle = 90, size = 8) + annotate("text", 0, 0.11, label = "Log2 Enrichment Scores", size = 9) + panel_border(remove = TRUE) 
polynomialDensityPlot
ggsave(plot = polynomialDensityPlot,
       filename = paste0(dirForRun, "/PolynomialFittingPlot.pdf"),
       device = "pdf", dpi = 600)


#lines(predict(polyN) ~ densityNegative$x, col = "darkred")
#points(densityPositive$y ~ densityPositive$x, col = "green", cex = 0.2)
#lines(predict(polyP) ~ densityPositive$x, col = "darkgreen")
plot(Log2Scores ~ totalScore, data = dfTotalLog2Scores, col = "black", main = "unbounded Log2Scores", xlab = "total TFBS score", ylim = c(-4,6))
plot(Log2Scores ~ totalScore, data = dfTotalLog2ScoresBounded, col = "black", main = "bounded Log2Scores", xlab = "total TFBS score", ylim = c(-4,6))

orderedLogScores <- dfTotalLog2ScoresBounded %>% arrange(EnsemblId)
head(orderedLogScores)
colnames(orderedLogScores)[3] <- "polynomial-generated bounded log 2 scores"
head(orderedLogScores)
write_csv(orderedLogScores, path = paste0("~/Documents/Project/Programming/Transcription Factor Binding/HumanWG_motifs/Data/singleLineFinalDf and Unified table/FinalLog2Scores/polynomialLog2Scores",
                                          gsub(" ", "_", Sys.time(), fixed = TRUE), ".csv"))










#close log file connection
close(fileConn)

