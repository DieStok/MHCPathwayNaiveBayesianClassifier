##Dieter Stoker, 4 September 2017
##Change Single Line Domain Info Into multiple lines


options(stringsAsFactors = FALSE)
install.packages("pacman", repos = "http://cran-mirror.cs.uu.nl/")
library("pacman")
pacman::p_load(ggplot2, dplyr, stringr, purrr, readr, LaF, ff, data.table, R.utils, base, utils, devtools, microbenchmark)
#argument: the data frame you want multipleline-nised


convertSingleLineToMultipleLineDomain <- function(fileName, domainColumn, Ensemblcolumn) {
  
  
  
  Df <- fread(fileName)
  print(head(Df))
  nRowsFile <- nrow(Df)
  print(nRowsFile)
  
  
  
  newDf <- data.frame(EnsemblID = rep("", times = (nRowsFile * 100)),
                      Domains = rep("", times = (nRowsFile * 100))
  )
  nRowsNewDf <- nrow(newDf)
  counterPositionInNewDf <- 1
  #print(head(newDf))
  
  
  apply(Df, MARGIN = 1, FUN = function (x){
    
    if(counterPositionInNewDf >= nRowsNewDf)
    {stop("The Data frame assigned for this operation was too small. Please initialise with a larger data frame.")}
    
    #print(paste0("counterposition", counterPositionInNewDf))
    #print(x[TFBScolumn])
    #Sys.sleep(2)
    
    #might need to change this to str_match!
    presentDomains <- str_match_all(x[domainColumn],"([^|]+)")
    #print(presentTFBSAndCounts[[1]])
    #Sys.sleep(10)
    currentEnsemblID     <- x[Ensemblcolumn]
    #if there are motifs present:
    if (length(head(presentDomains[[1]] > 0))) {
      
      presentDomains         <- presentDomains[[1]][,1]
      
      
      
      rowsToFill <- length(presentDomains)
      #print(rowsToFill)
      #print(rep(currentEnsemblID, rowsToFill))
      #print(newDf[counterPositionInNewDf:(counterPositionInNewDf + rowsToFill-1), "EnsemblID"])
      
      #Sys.sleep(4) 
      newDf[counterPositionInNewDf:(counterPositionInNewDf + rowsToFill-1), "EnsemblID"] <<- rep(currentEnsemblID, rowsToFill)
      #print(head(presentTFBS))
      newDf[counterPositionInNewDf:(counterPositionInNewDf + rowsToFill-1), "Domains"] <<- presentDomains

      
      #print(head(newDf[counterPositionInNewDf:(counterPositionInNewDf + rowsToFill-1), ]))
      #Sys.sleep(4) 
      
    } else {
      
      rowsToFill <- 1
      newDf[counterPositionInNewDf, "EnsemblID"] <<- currentEnsemblID
      newDf[counterPositionInNewDf, "Domains"] <<- "NONE"
    }
    
    counterPositionInNewDf                                <<- counterPositionInNewDf + rowsToFill
    
    
    
  })
  
  newDf <- newDf[newDf["EnsemblID"]!= "",]
  head(newDf)
  tail(newDf)
  newDf
}

#various conversions

PfamDomains <- convertSingleLineToMultipleLineDomain(fileName = "Documents/Project/Programming/Domains/Data/PfamDomainsSingleLineTable.csv", 
                                                     domainColumn = 2,
                                                     Ensemblcolumn = 1)
head(PfamDomains, 30)
tail(PfamDomains)
write_csv(PfamDomains, "Documents/Project/Programming/Domains/Data/MultipleLinePfamDomainData.csv")
