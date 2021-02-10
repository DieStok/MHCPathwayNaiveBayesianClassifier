
install.packages("pacman", repos = "http://cran-mirror.cs.uu.nl/")
library("pacman")
pacman::p_load(ggplot2, dplyr, stringr, purrr, readr, LaF, ff, data.table, R.utils, base, utils, devtools, microbenchmark)
options(error=utils::recover)
options(stringsAsFactors = FALSE)


`%ni%` <- Negate(`%in%`)
oldwd <- getwd()
scriptwd <- "/home/dieter/Documents/Project/Programming/Transcription Factor Binding/29mammals/CallingMotifsMyself/Data/FinalDfAndUnifiedMotifFile/"
setwd(scriptwd)


amountOfFiles <- (16-1)
#so you read the first file to create the list,
#then read in the other files, combine them.
listFiles <- list(fread("singleLineFinalDf00.csv"))




for (i in 2:amountOfFiles) {
  
  #load files
  if(i < 10)
    listFiles <- append(listFiles, list(fread(paste0("singleLineFinalDf", "0", i, ".csv")))) 
  else
    listFiles <- append(listFiles, list(fread(paste0("singleLineFinalDf", i, ".csv"))))
  print(i)
  
  ENSIDsNotEmptyFirst   <- listFiles[[1]][which(listFiles[[1]][, "Motifs"] != ""), ][, "Ensembl_ID"][[1]]
  ENSIDsNotEmptySecond  <- listFiles[[i]][which(listFiles[[i]][, "Motifs"] != ""), ][, "Ensembl_ID"][[1]]
  head(ENSIDsNotEmptyToPaste)
  ENSIDsNotEmptyToPaste <- ENSIDsNotEmptySecond[ENSIDsNotEmptySecond %in% ENSIDsNotEmptyFirst]
  
  
  listFiles[[1]][listFiles[[1]]$Ensembl_ID %in% ENSIDsNotEmptyToPaste,][, "Motifs"] <- paste0(listFiles[[1]][listFiles[[1]]$Ensembl_ID %in% ENSIDsNotEmptyToPaste,][, "Motifs"][[1]],
                                                                                                            "|",
                                                                                              listFiles[[i]][listFiles[[i]]$Ensembl_ID %in% ENSIDsNotEmptyToPaste,][, "Motifs"][[1]])
  
  print(head(ENSIDsNotEmptySecond[ENSIDsNotEmptySecond %ni% ENSIDsNotEmptyToPaste]))
  ENSIDsEmptyInFirstToPaste <- ENSIDsNotEmptySecond[ENSIDsNotEmptySecond %ni% ENSIDsNotEmptyToPaste]
  
  listFiles[[1]][listFiles[[1]]$Ensembl_ID %in% ENSIDsEmptyInFirstToPaste, ][, "Motifs"][[1]] <- listFiles[[i]][listFiles[[i]]$Ensembl_ID %in% ENSIDsEmptyInFirstToPaste,][, "Motifs"][[1]]
  
}
  



  barbie <- str_match_all(listFiles[[1]]$Motifs,"([^|]+)-(\\d+)")
  finalSingleLineDf <- data.table(Ensembl_ID = listFiles[[1]]$Ensembl_ID,
                                  Motifs = "")
  
  
  for (i in 1:length(barbie)){
    dt <- as.data.table(barbie[[i]][, c(2,3), drop = FALSE])
    dt$V2 <- as.numeric(dt$V2)
    dt <- dt[, list(motiftotals = sum(V2)), by = V1]
    motifs <- paste0(dt[, "V1"][[1]], "-", dt[, "motiftotals"][[1]], collapse = "|")
    if (motifs == "-") motifs = ""
    finalSingleLineDf[i,][,"Motifs"][[1]] <- motifs
  }
  
  

  head(finalSingleLineDf)  
  write_csv(finalSingleLineDf, path = paste0(scriptwd, "/combinedMotifFile29Mammals.csv"))  
    
 #oldwd
  setwd(oldwd)
  

##############################################################################  
##                              TESTING                                     ##
##############################################################################
  
##TESTING STUFF
##
# testSingleLineFinalDf <- singleLineFinalDf[singleLineFinalDf$Motifs != "",]
# testSingleLineFinalDf[testSingleLineFinalDf$Ensembl_ID == "ENSG00000131584",][,"Motifs"] <- ""
# testSingleLineFinalDf[testSingleLineFinalDf$Ensembl_ID == "ENSG00000078900",][,"Motifs"] <- ""
# testSingleLineFinalDf2 <- singleLineFinalDf[singleLineFinalDf$Motifs != "",]
# testSingleLineFinalDf2[testSingleLineFinalDf2$Ensembl_ID == "ENSG00000131584",][,"Motifs"] <- ""
# testSingleLineFinalDf2[testSingleLineFinalDf2$Ensembl_ID == "ENSG00000008130",][,"Motifs"] <- "Japiekrekelenstein-3"
# testSingleLineFinalDf2[testSingleLineFinalDf2$Ensembl_ID == "ENSG00000157933",][,"Motifs"] <- ""
# ##

##MORE TESTING STUFF

# ENSIDsNotEmptyFirst <- testSingleLineFinalDf[which(testSingleLineFinalDf[, "Motifs"] != ""), ][, "Ensembl_ID"][[1]]
# ENSIDsNotEmptySecond <- testSingleLineFinalDf2[which(testSingleLineFinalDf2[, "Motifs"] != ""), ][, "Ensembl_ID"][[1]]
# ENSIDsNotEmptyToPaste <- ENSIDsNotEmptySecond[ENSIDsNotEmptySecond %in% ENSIDsNotEmptyFirst]
# 
# 
# testSingleLineFinalDf[testSingleLineFinalDf$Ensembl_ID %in% ENSIDsNotEmptyToPaste,][, "Motifs"] <- paste0(testSingleLineFinalDf[testSingleLineFinalDf$Ensembl_ID %in% ENSIDsNotEmptyToPaste,][, "Motifs"][[1]],
#                                                                                                         "|",
#                                                                                                        testSingleLineFinalDf2[testSingleLineFinalDf2$Ensembl_ID %in% ENSIDsNotEmptyToPaste,][, "Motifs"][[1]])
# 
# ENSIDsEmptyInFirstToPaste <- ENSIDsNotEmptySecond[ENSIDsNotEmptySecond %ni% ENSIDsNotEmptyToPaste][[1]]
# 
# testSingleLineFinalDf[testSingleLineFinalDf$Ensembl_ID %in% ENSIDsEmptyInFirstToPaste, ][, "Motifs"][[1]] <- testSingleLineFinalDf2[testSingleLineFinalDf2$Ensembl_ID %in% ENSIDsEmptyInFirstToPaste,][, "Motifs"][[1]]
# 
# barbie <- str_match_all(testSingleLineFinalDf$Motifs,"([^|]+)-(\\d+)")
# 
# 
# for (i in 1:length(barbie)){
#   dt <- as.data.table(barbie[[i]][, c(2,3), drop = FALSE])
#   dt$V2 <- as.numeric(dt$V2)
#   dt <- dt[, list(motiftotals = sum(V2)), by = V1]
#   #testing
#   motifs <- paste0(dt[, "V1"][[1]], "-", dt[, "motiftotals"][[1]], collapse = "|")
#   if (motifs == "-") motifs = ""
#   motifs
#   testFinalSingleLineDf[i,][,"Motifs"][[1]] <- motifs
# }








##old test##
# i = 2 
# dt <- as.data.table(barbie[[2]][, c(2,3), drop = FALSE])
# dt$V2 <- as.numeric(dt$V2)
# dt <- dt[, list(motiftotals = sum(V2)), by = V1]
# #testing
# barpie <- data.table(V1 = "Honderd-3", motiftotals = 3)
# dt<- rbind(dt, barpie)
# testFinalSingleLineDf[i,][,"Motifs"][[1]] <- paste0(dt[, "V1"][[1]], "-", dt[, "motiftotals"][[1]], collapse = "|")

