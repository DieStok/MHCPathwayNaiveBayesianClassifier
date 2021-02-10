##Dieter Stoker
#01-07-2017
#Script for calculating motif enrichments when different motif versions not taken into account
#i.e. STAT_9 and STAT_12 in one gene both become simply STAT in that gene



#packages
#install.packages("pacman", repos = "http://cran-mirror.cs.uu.nl/")
library("pacman")
pacman::p_load(ggplot2, dplyr, stringr, purrr, readr, LaF, ff, data.table, R.utils, base, utils, devtools, microbenchmark)



#create new data: unique (i.e. one gene with GABP_2 and GABP_3 will now contain GABP only once)
oldDataOneLine <- fread("~/Documents/Project/Programming/Transcription Factor Binding/29mammals_unique_tfbs_per_gene_one_line.txt", header = FALSE)
regexSplitList <- lapply(oldDataOneLine[,2], FUN = str_match_all, pattern = "([a-zA-Z0-9-:.]+)(_\\d+){0,1}")[[1]]
head(regexSplitList)
uniqueMotifsOneLine <- unlist(lapply(regexSplitList, FUN = function(x) { uniqueMotifs <- unique(x[,2]);paste0(uniqueMotifs, collapse = "|")}))
ensemblIds <- oldDataOneLine[,1]
uniqueNewMotifDataOneLine <- data.table(EnsemblId = ensemblIds, TranscriptionFactorBindingSites = uniqueMotifsOneLine)
names(uniqueNewMotifDataOneLine) <- c("EnsemblId", "TranscriptionFactorBindingSites")
head(uniqueNewMotifDataOneLine)
head(oldDataOneLine)
write_csv(uniqueNewMotifDataOneLine, "~/Documents/Project/Programming/Transcription Factor Binding/29mammals_unique_tfbs_per_gene_one_line_motifversionsdeleted_unique_DieterStoker_01_07_2017.txt")

#now for the multiple lines

oldDataOneLine <- fread("~/Documents/Project/Programming/Transcription Factor Binding/29mammals_unique_tfbs_per_gene_one_line.txt", header = FALSE)
regexSplitList <- lapply(oldDataOneLine[,2], FUN = str_match_all, pattern = "([a-zA-Z0-9-:.]+)(_\\d+){0,1}")[[1]]
head(regexSplitList)
dataTableCounter = 1
#the below gets the unique motifs per gene from the oldDataOneLine, and then
#repeats the linked ENSG ID the needed number of times, and rbinds the results.
uniqueNewMotifDataMultipleLines <- do.call(rbind, lapply(regexSplitList, FUN = function(x) { uniquemotifs = unique(x[,2])
repsGene <- rep(oldDataOneLine[dataTableCounter,1], length(uniquemotifs))
dataTableCounter <<- dataTableCounter + 1
dataTable <- data.table(repsGene, uniquemotifs)
names(dataTable) <- c("EnsemblId", "TranscriptionFactorBindingSites")
dataTable} ))
head(uniqueNewMotifDataMultipleLines)
uniqueNewMotifDataMultipleLines <- data.frame(sapply(uniqueNewMotifDataMultipleLines, as.character))
head(uniqueNewMotifDataMultipleLines)
tail(uniqueNewMotifDataMultipleLines)
#sapply call because the rbind and lapply above have rendered the columns lists instead of vectors.
write_csv(uniqueNewMotifDataMultipleLines, "~/Documents/Project/Programming/Transcription Factor Binding/29mammals_unique_tfbs_per_gene_multiple_lines_motifversionsdeleted_unique_DieterStoker_01_07_2017.txt")


##########################################################
#Initial options
##########################################################

options(stringsAsFactors = FALSE)
setwd("~/Documents/Project/Programming/Transcription Factor Binding")

##########################################################
#Function to load in the data and required libraries
#reloadFromSource = whether source files like positive list etc. should be selected manually
##########################################################

loadLibraryAndData <- function(reloadFromSource = FALSE) {
  
  #load required libraries
  #if(!require("plyr")) install.packages("plyr")
  if(!require("ggplot2")) install.packages("ggplot2")
  if(!require("reshape2")) install.packages("reshape2")
  if(!require("tcltk")) install.packages("tcltk")
  
  #library("plyr")
  library("ggplot2")
  library("reshape2")
  library("tcltk")
  
  currentDir <- getwd()
  
  options(stringsAsFactors = FALSE)
  
  if(reloadFromSource == TRUE) {
    
    colNamesTFBSTables <- c("EnsemblId","TranscriptionFactorBindingSites")
    
    
    print("Choose the file that contains the conserved TFBS on one line")
    tFTableOneLine <- fread(tk_choose.files(default = "~/Documents/Project/Programming/Transcription Factor Binding/29mammals_unique_tfbs_per_gene_one_line_motifversionsdeleted_unique_DieterStoker_01_07_2017.txt",
                                            caption = "Choose the .txt file that contains the conserved TFBS on one line",
                                            multi = FALSE))
    colnames(tFTableOneLine) <- colNamesTFBSTables
    print(head(tFTableOneLine))
    
    
    tFTableMultipleLines <- fread(tk_choose.files(default = "~/Documents/Project/Programming/Transcription Factor Binding/29mammals_unique_tfbs_per_gene_multiple_lines_motifversionsdeleted_unique_DieterStoker_01_07_2017.txt",
                                                  caption = "Choose the .txt file that contains the conserved TFBS on multiple lines",
                                                  multi = FALSE))
    colnames(tFTableMultipleLines) <- colNamesTFBSTables
    print(head(tFTableMultipleLines))
    
    MHCGeneTable <- read.csv(tk_choose.files(default = "~/Documents/Project/Programming/DataForAll/PositiveList_2.csv",
                                             caption = "Choose the .csv file that contains the positive MHC gene set",
                                             multi = FALSE),header = TRUE)
    
    names(MHCGeneTable)[7] <- "EnsemblId"
    print(head(MHCGeneTable))
    
    #save as .R data files
    
    saveRDS(tFTableOneLine, file = paste(currentDir, "/tFTableOneLineUniqueMotifs.rds", sep=""))
    saveRDS(tFTableMultipleLines, file = paste(currentDir, "/tFTableMultipleLinesUniqueMotifs.rds", sep=""))
    saveRDS(MHCGeneTable, file = paste(currentDir, "/MHCGeneTable.rds", sep=""))
    
  } else if (reloadFromSource == FALSE) {
    
    tFTableOneLine <- readRDS(file = paste(currentDir, "/tFTableOneLineUniqueMotifs.rds", sep=""))
    tFTableMultipleLines <- readRDS(file = paste(currentDir, "/tFTableMultipleLinesUniqueMotifs.rds", sep=""))
    MHCGeneTable <- readRDS(file = paste(currentDir, "/MHCGeneTable.rds", sep=""))
    
  }
  
  #______________________________
  
  amountPositiveSetInList <- sum(MHCGeneTable$EnsemblId %in% tFTableMultipleLines$EnsemblId)
  relevantTFBSDf <- tFTableMultipleLines[tFTableMultipleLines$EnsemblId %in% MHCGeneTable$EnsemblId,]
  print(paste("Amount of positive set genes in 29 mammals data set:", amountPositiveSetInList))
  print(head(relevantTFBSDf))
  #note: relevantTFBSDf thus contains the TFBS linked to MHC pathway related genes
  
  
  #Count data TFBS for MHC gene set
  countDataNames <-c("TFBS","Counts")
  
  positiveSetTFBSCounts <- plyr::count(relevantTFBSDf$TranscriptionFactorBindingSites)
  names(positiveSetTFBSCounts) <- countDataNames
  print("The counts of the TFBS present in MHC genes:")
  print(head(positiveSetTFBSCounts))
  
  #count the totals for the positive set
  totalCountsMHCSet <- sum(positiveSetTFBSCounts$Counts)
  
  #Count data for TFBS in all genes
  
  totalSetTFBSCounts = plyr::count(tFTableMultipleLines$TranscriptionFactorBindingSites)
  names(totalSetTFBSCounts) = countDataNames
  
  totalCountsFullSet = sum(totalSetTFBSCounts$Counts)
  
  #calculate the odds that a random gene in the full set has a certain TFBS
  
  genomeOdds <- apply(totalSetTFBSCounts, FUN = function(x) {
    as.numeric(x[2])/totalCountsFullSet
  }, MARGIN = 1)
  names(genomeOdds) <- totalSetTFBSCounts$TFBS
  
  #calculate the odds that a gene in the MHC set has a certain TFBS
  
  MHCOdds <- apply(positiveSetTFBSCounts, FUN = function(x) {
    as.numeric(x[2])/totalCountsMHCSet
  }, MARGIN = 1)
  names(MHCOdds) <- positiveSetTFBSCounts$TFBS
  
  
  return(list(tFTableOneLine = tFTableOneLine,
              tFTableMultipleLines = tFTableMultipleLines,
              MHCGeneTable = MHCGeneTable,
              MHCGenesInDataset = amountPositiveSetInList,
              genomeOdds = genomeOdds,
              MHCOdds = MHCOdds,
              positiveSetTFBSCounts = positiveSetTFBSCounts,
              totalSetTFBSCounts = totalSetTFBSCounts,
              TFBSinMHCGenes = relevantTFBSDf
  )
  )
  
}



##########################################################
#computerFishersTest
#function takes whole genome count data and positive set count data,
#executes fisher's exact p on all relevant TFSB to find if they are enriched
#in positive set
##########################################################

computeFishersTest = function(wholeGenomeTFBSCountData, positiveSetTFBSCountData, run = FALSE) {
  
  #don't compute again unless explicitly told to. Otherwise, just open the .rds
  if(run == FALSE) {
    
    fisherDataFrame <- readRDS(file = paste(getwd(), "/fisherDataFrameUniqueMotifs.rds", sep = ""))
    return(fisherDataFrame)
    
  }
  else {
    
    #Get only those genes from the whole genome that have TFBS that are also found in the positive MHC set
    relevantWholeGenome <- wholeGenomeTFBSCountData[wholeGenomeTFBSCountData$TFBS %in% positiveSetTFBSCountData$TFBS,] 
    totalCountsWholeGenome <- sum(wholeGenomeTFBSCountData$Counts)
    totalCountsMHCGenes <- sum(positiveSetTFBSCountData$Counts)
    
    #listMatrices is for debug purposes, to check whether correct contingency tables are created
    listMatrices <- list()
    fisherDataFrame <- data.frame(TFBS = character(),
                                  pValue=numeric(),
                                  sampleEstimate=numeric(),
                                  lowerConfidenceInterval=numeric(),
                                  upperConfidenceInterval=numeric())
    
    #cycle over all the rows of 
    for(bindingSites in 1:nrow(positiveSetTFBSCountData)) {
      
      fisherMatrix <- matrix(c(positiveSetTFBSCountData$Counts[bindingSites],
                               totalCountsMHCGenes - positiveSetTFBSCountData$Counts[bindingSites],
                               relevantWholeGenome$Counts[bindingSites],
                               totalCountsWholeGenome - relevantWholeGenome$Counts[bindingSites]),
                             nrow = 2,
                             dimnames = list(TFBS = c(relevantWholeGenome$TFBS[bindingSites], "TotalTFBS"),
                                             Set = c("Positive set", "Whole genome")))
      fisherTestResult <- fisher.test(fisherMatrix)
      pValue = fisherTestResult["p.value"]
      lowerConfInt <- unlist(fisherTestResult["conf.int"])[1]
      upperConfInt <- unlist(fisherTestResult["conf.int"])[2]
      sampleEst <- fisherTestResult["estimate"]
      
      listMatrices <- append(listMatrices,fisherMatrix)
      addVector <- c(as.character(relevantWholeGenome$TFBS[bindingSites]),
                     pValue,
                     sampleEst,
                     lowerConfInt,
                     upperConfInt)
      names(addVector) <- c("TFBS", "pValue", "sampleEstimate",
                            "lowerConfidenceInterval", "upperConfidenceInterval")
      #add to the Df the results of the fisherTest
      fisherDataFrame <- rbind(fisherDataFrame, addVector, make.row.names = FALSE)
      
    }
    #adjust for multiple testing with Benjamini-Hochberg methodology (utilises FDR)
    fisherDataFrame$qValue <- p.adjust(fisherDataFrame$pValue, method="BH")
    saveRDS(fisherDataFrame, file = paste(getwd(), "/fisherDataFrameUniqueMotifs.rds", sep = ""))
    return(fisherDataFrame)
  }
  
  
}

##########################################################
#function scoresOverRows
#takes the dataframe with Fisher tests and sample estimates as well as whole genome data (one line).
#returns dataframe with score per TFBS and total score of all genes in the 29 mammals data.
##########################################################
scoresOverRows <- function(wholeGenomeDataOneLine, fisherTests) {
  
  # find TFBS that are significantly overrepresented in positive set genes
  #assign them a score. log(sample Estimate/qValue). Log to scale,
  #sample estimate to incorporate effect size (i.e. HOW overrepresented are they in MHC genes?)
  #qValue to incorporate probability ( TFBS with p>0.05<0.10 are less sure to be correct than those with p< 0.001)
  
  informativeTFBS <- fisherTests[fisherTests$qValue <= 0.10,]
  informativeTFBS[order(informativeTFBS$qValue),]
  rownames(informativeTFBS) = seq(1:nrow(informativeTFBS))
  informativeTFBS$score <- log(informativeTFBS$sampleEstimate / informativeTFBS$qValue)
  
  
  
  scoreCompendium <- vector()
  dataFrameResults <- data.frame()
  #cycle over rows. Could be refined with an apply function
  for(rows in 1:nrow(wholeGenomeDataOneLine)) {
    
    #check whether any of the informative TFBS are present in a gene's promoter region
    currentRow <- wholeGenomeDataOneLine[rows,]
    currentRowTFs <- currentRow[, 2]
    characterVectorTFs <- unlist(strsplit(as.character(currentRowTFs), split = "|", fixed = TRUE))
    logicalPresence <- informativeTFBS$TFBS %in% characterVectorTFs
    currentRowEnsemblId <- currentRow[, "EnsemblId"]
    
    for(bindingSites in 1:length(logicalPresence)) {
      #for all informative TFBS, if they are present, add their score, otherwise, add 0
      if(logicalPresence[bindingSites] == TRUE) {
        
        valueToAdd <- informativeTFBS$score[bindingSites]
        
        
      } else if (logicalPresence[bindingSites] == FALSE) {
        
        valueToAdd <- 0
        
      }
      
      scoreCompendium <- append(scoreCompendium,valueToAdd)
      names(scoreCompendium)[bindingSites] <-  informativeTFBS$TFBS[bindingSites]
      
    }
    #set rownames as the EnsemblId of the gene from the whole genome
    dataFrameResults <- rbind(dataFrameResults, scoreCompendium)
    rownames(dataFrameResults) <- c(head(rownames(dataFrameResults),-1),currentRow["EnsemblId"])
    scoreCompendium <- vector()
    #print(nrow(dataFrameResults))
    
  }
  colnames(dataFrameResults) <- paste("Score for TFBS", informativeTFBS$TFBS)
  dataFrameResults$totalScore <- rowSums(dataFrameResults)
  dataFrameResults <- dataFrameResults[order(dataFrameResults$totalScore, decreasing = TRUE),]
  print("These are the informative TFBS:")
  print(informativeTFBS$TFBS)
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
dataStructureTFBS <- loadLibraryAndData(reloadFromSource = FALSE)
str(dataStructureTFBS)
dataStructureTFBS$tFTableOneLine
dataStructureTFBS$tFTableMultipleLines
fisherTests <- computeFishersTest(dataStructureTFBS$totalSetTFBSCounts,dataStructureTFBS$positiveSetTFBSCounts, run = TRUE)
str(fisherTests)
informativeTFBS <- subset(fisherTests, qValue <= 0.10)
informativeTFBS

#switch for execution of the selection of all genes in the WGD that have the TFBS sig. overrepresented in positive set genes
#set to TRUE when the positive set has changed and/or the analysis should be run again.

runselection = FALSE

if (runselection == TRUE) {
  finalScoresWholeGenome <- scoresOverRows(as.data.frame(dataStructureTFBS$tFTableOneLine),as.data.frame(fisherTests))
  #Save the final scores of each gene with respect to the amount of MHC-associated TFBS it has, if rerunning
  saveRDS(finalScoresWholeGenome, file = paste(getwd(), "/FinalScoreWholeGenomeUniqueMotifs.rds",sep=""))
}
if (runselection == FALSE) {
  finalScoresWholeGenome <- readRDS(file = paste(getwd(), "/FinalScoreWholeGenomeUniqueMotifs.rds", sep=""))
}
str(finalScoresWholeGenome)
head(finalScoresWholeGenome$scoresWholeGenome)




##################################################################################################
# This piece of code creates histograms of the positive and negative set for the scores of all TFBS
#
##################################################################################################



positiveSetIDs <- dataStructureTFBS$MHCGeneTable$EnsemblId
inPositiveSet <- rownames(finalScoresWholeGenome$scoresWholeGenome) %in% positiveSetIDs
notInPositiveSet <- rownames(finalScoresWholeGenome$scoresWholeGenome) %ni% positiveSetIDs
finalPlotDfHist <- data.frame(scores = numeric(), set = character(), ScoreFor = character())
for (i in 1:ncol(finalScoresWholeGenome$scoresWholeGenome)) {
  finalScoresPositiveSetTotal <- finalScoresWholeGenome$scoresWholeGenome[inPositiveSet, i]
  finalScoresNegativeSetTotal <- finalScoresWholeGenome$scoresWholeGenome[notInPositiveSet, i]
  plotDfHistPositive <- data.frame(scores = finalScoresPositiveSetTotal, set = "positive")
  plotDfHistNegative <- data.frame(scores = finalScoresNegativeSetTotal, set = "negative")
  totalPlotDfHist <- rbind(plotDfHistNegative, plotDfHistPositive)
  totalPlotDfHist$ScoreFor <- colnames(finalScoresWholeGenome$scoresWholeGenome)[i]
  finalPlotDfHist <- rbind(finalPlotDfHist, totalPlotDfHist)
}


ggplot(data = finalPlotDfHist, aes(x = scores)) +
  geom_density(aes(y = ..density.., fill = set), color = "black",  alpha = 0.3) +
  scale_y_continuous(breaks=waiver()) +
  theme_bw() + 
  facet_wrap(~ScoreFor, scales = "free")

bwidth <- 4
chicken <- ggplot(data = finalPlotDfHist, aes(x = scores)) +
  geom_histogram(aes(y = 4 * ..density.., fill = set), color = "black", position = "dodge", stat = "bin", binwidth = 4) +
  scale_y_continuous(breaks=waiver(), expand=c(0,0), trans = "log1p") +
  theme_bw()  + 
  facet_wrap(~ScoreFor, scales = "free") 












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
    proportionInMHC <- sum(TFBSdata$TFBSinMHCGenes$TranscriptionFactorBindingSites %in% infTFBS[significantTFBS])/TFBSdata$MHCGenesInDataset
    proportionNotInMHC <- (TFBSdata$MHCGenesInDataset - sum(TFBSdata$TFBSinMHCGenes$TranscriptionFactorBindingSites %in% infTFBS[significantTFBS]))/TFBSdata$MHCGenesInDataset
    proportionInWGD <- sum(TFBSdata$tFTableMultipleLines$TranscriptionFactorBindingSites %in% infTFBS[significantTFBS])/nrow(TFBSdata$tFTableOneLine)
    proportionNotInWGD <- (nrow(TFBSdata$tFTableOneLine)-sum(TFBSdata$tFTableMultipleLines$TranscriptionFactorBindingSites %in% infTFBS[significantTFBS]))/nrow(TFBSdata$tFTableOneLine)
    tempDf = data.frame(TF = rep(infTFBS[significantTFBS],4),
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
yesNoPlotData <- createGraphData()
#step 2:
#Now plot this data.

plotGraphData <- function(plottingData, run = TRUE) {
  
  if(run == TRUE) {
    green <- rgb(0, 140, 71, maxColorValue = 255)
    red   <- rgb(237, 45, 46, maxColorValue = 255)
    barBorder <- "black"
    initialPlot <-ggplot(data = plottingData,
                         aes(x= factor(yesno, levels = c("yes","no")),
                             y = proportion, fill = factor(genomeMHC, levels = c ("MHC", "genome")))) + 
      facet_wrap(~TF, scales = "free") +
      geom_bar(stat = "identity", position = position_dodge(), colour = barBorder) +
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
    ggsave(plot = annotatedPlot, filename = "TFBSOverrepresentationPlot_binaryUniqueMotifs.pdf", device = "pdf", dpi = 600)
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
    Takegenes <- transcriptionFactorTableMultipleLines[transcriptionFactorTableMultipleLines$EnsemblId %in% positiveMHCGenes$EnsemblId,]
    Takegenes
    
    #now, for every relevant transcription factor, return the genes that have it
    genelist = list()
    for(relevantTFBS in 1:length(finalScoresWholeGenome$informativeTFBS$TFBS)) {
      
      positiveGenesWithTFBS = Takegenes[finalScoresWholeGenome$informativeTFBS$TFBS[relevantTFBS] == as.vector(Takegenes$transcriptionFactorBindingSites),]
      print(positiveGenesWithTFBS)
      EnsemblIds = positiveGenesWithTFBS$EnsemblId
      genelist[[relevantTFBS]] <- positiveMHCGenes[positiveMHCGenes$EnsemblId  %in%  EnsemblIds,]
      names(genelist)[relevantTFBS] <- finalScoresWholeGenome$informativeTFBS$TFBS[relevantTFBS]
      genelist[[relevantTFBS]] <- cbind(TFBS = rep(finalScoresWholeGenome$informativeTFBS$TFBS[relevantTFBS],length(genelist[[relevantTFBS]]$MHC1OR2)),
                                        genelist[[relevantTFBS]] )
      if (relevantTFBS == 1) {
        
        write.csv(genelist[[relevantTFBS]], "PositiveGenesWithTFBSUniqueMotifs.csv")
        
      }
      
      else if (relevantTFBS > 1) {
        
        write.table(genelist[[relevantTFBS]], "PositiveGenesWithTFBSUniqueMotifs.csv", append = TRUE, sep = ",", col.names = FALSE)
        
      }
    }
  }
  
}
#manual switch for running this, set to true if the positive set changes in some form.
outputMHCGeneTFBS(run = TRUE)

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
      Ids <- unique(WGD_relevantgenes$EnsemblId)
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
    write.csv(dfWGDEnsmblIdsPerTFBS, "EnsemblIdsWGDUniqueMotifs.csv")
    
    #return the df for checking in R
    return (dfWGDEnsmblIdsPerTFBS)
  }
  
}

outputEnsIdsDataFrame <- outputEnsIdsWithTFBS(run = TRUE)

###################################################
#                       END                       #
###################################################

