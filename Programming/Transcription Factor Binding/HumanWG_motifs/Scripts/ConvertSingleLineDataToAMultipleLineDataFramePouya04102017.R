##Dieter Stoker, 4th of October 2017
##Change single line format to multiple line format for Pouya Kheradpour data


options(stringsAsFactors = FALSE)
install.packages("pacman", repos = "http://cran-mirror.cs.uu.nl/")
library("pacman")
pacman::p_load(ggplot2, dplyr, stringr, purrr, readr, LaF, ff, data.table, R.utils, base, utils, devtools, microbenchmark)
#argument: the data frame you want multipleline-nised


convertSingleLineToMultipleLine <- function(fileName, TFBScolumn, Ensemblcolumn) {
  
  

Df <- fread(fileName)
print(head(Df))
nRowsFile <- nrow(Df)
print(nRowsFile)



newDf <- data.frame(EnsemblID = rep("", times = (nRowsFile * 10000)),
                    TFBS = rep("", times = (nRowsFile * 10000)),
                    Count = rep(-1, times = (nRowsFile * 10000))
                    )
nRowsNewDf <- nrow(newDf)
counterPositionInNewDf <- 1
#print(head(newDf))


apply(Df, MARGIN = 1, FUN = function (x){
  
  if(counterPositionInNewDf >= nRowsNewDf)
  {stop("The Data frame assigned for this operation was too small. Please initialise with a larger data frame.")}
  
  
  if(counterPositionInNewDf %% 200 == 0) {
  print(paste0("counterposition", counterPositionInNewDf))
  }
  #print(x[TFBScolumn])
  #Sys.sleep(2)
  
  #might need to change this to str_match!
  presentTFBSAndCounts <- str_match_all(x[TFBScolumn],"([^|]+)-(\\d+)")
  #print(presentTFBSAndCounts[[1]])
  #Sys.sleep(10)
  currentEnsemblID     <- x[Ensemblcolumn]
  #if there are motifs present:
  if (length(head(presentTFBSAndCounts[[1]] > 0))) {
  
  presentTFBS          <- presentTFBSAndCounts[[1]][,2]
  countsTFBS           <- as.numeric(presentTFBSAndCounts[[1]][,3])
  
  
  rowsToFill <- length(presentTFBS)
  #print(rowsToFill)
  #print(rep(currentEnsemblID, rowsToFill))
  #print(newDf[counterPositionInNewDf:(counterPositionInNewDf + rowsToFill-1), "EnsemblID"])
  
  #Sys.sleep(4) 
  newDf[counterPositionInNewDf:(counterPositionInNewDf + rowsToFill-1), "EnsemblID"] <<- rep(currentEnsemblID, rowsToFill)
  #print(head(presentTFBS))
  newDf[counterPositionInNewDf:(counterPositionInNewDf + rowsToFill-1),      "TFBS"] <<- presentTFBS
  #print(head(countsTFBS))
  newDf[counterPositionInNewDf:(counterPositionInNewDf + rowsToFill-1),     "Count"] <<- countsTFBS
  
  #print(head(newDf[counterPositionInNewDf:(counterPositionInNewDf + rowsToFill-1), ]))
  #Sys.sleep(4) 
  
  } else {
    
    rowsToFill <- 1
    newDf[counterPositionInNewDf, "EnsemblID"] <<- currentEnsemblID
    newDf[counterPositionInNewDf, "TFBS"] <<- "NONE"
    newDf[counterPositionInNewDf, "Count"] <<- -1
  }
  
  counterPositionInNewDf                                <<- counterPositionInNewDf + rowsToFill
  

  
})

newDf <- newDf[newDf["EnsemblID"]!= "",]
head(newDf)
tail(newDf)
print("Done!")
newDf
}

#converting various data

PouyaKheradPourENCODETFdataMultipleLines <- convertSingleLineToMultipleLine("/home/dieter/Documents/Project/Programming/Transcription Factor Binding/HumanWG_motifs/Data/singleLineFinalDf and Unified table/combinedMotifFilePouyaKheradpour.csv", 2, 1)
write_csv(PouyaKheradPourENCODETFdataMultipleLines, path = "/home/dieter/Documents/Project/Programming/Transcription Factor Binding/HumanWG_motifs/Data/singleLineFinalDf and Unified table/MultipleLinesFilePouyaKheradPour.csv")
#twentyNineMammalsMyOwnCallingDataMultipleLines <- convertSingleLineToMultipleLine("/home/dieter/Documents/Project/Programming/Transcription Factor Binding/29mammals/CallingMotifsMyself/Data/FinalDfAndUnifiedMotifFile/combinedMotifFile29Mammals.csv", 2, 1)
#write_csv(twentyNineMammalsMyOwnCallingDataMultipleLines, path = "/home/dieter/Documents/Project/Programming/Transcription Factor Binding/29mammals/CallingMotifsMyself/Data/FinalDfAndUnifiedMotifFile//MultipleLines29Mammals.csv")
