##Finalised on 13_09_2017. Comments last edited 5_02_2018
##Requires a CSV file that contains chromosome, position, gene name, strand etc., with naming x00-x15
## configured to run on mutants with manual input of file to work on


#initialise

commandLineArguments <- commandArgs(trailingOnly = TRUE)
if (!length(commandLineArguments) == 1) {
  stop("Please supply the name of the input match file")
}

#options
options(stringsAsFactors = FALSE)
r <- getOption("repos")
r["CRAN"] <- "http://cran-mirror.cs.uu.nl/"
options(repos <- r)

#packages and associated options
install.packages("pacman", repos = "http://cran-mirror.cs.uu.nl/")
library("pacman")
pacman::p_load(ggplot2, dplyr, stringr, purrr, readr, LaF, ff, data.table, R.utils, base, utils, devtools, microbenchmark)
options(datatable.verbose = FALSE)

#filenames
matchesFile <- paste0("/linuxhome/tmp/dieter/", commandLineArguments[1])
outputFile <- paste0("/linuxhome/tmp/dieter/singleLineFinalDf", str_match(matchesFile, "\\d{2}")[1],".csv")
print(matchesFile)
print(outputFile)

#initialise the function to map TFBS to genes
mapTFBS <- function(x) {
  
  
  #give the row some logical names, convert to correct classes
  #note: this is because I am using apply while instead, I should use something else.
  #at this point, I am sorry to say that that is legacy code. I have improved and incorporated purr and mapply in newer code.
  rowData <- x
  names(rowData) <- c("motifName", "motifChrom", "motifStart", "motifEnd", "motifStrand")
  
  #get the motif name
  
  newMotif <- rowData["motifName"]
  #Just debbuging things below#
  #print(x)
  #print(rowData)
  #print(typeof(rowData[3]))
  #print(paste0("the motif name: ", newMotif))
  
  if(rowData["motifChrom"] %in% names(subsetList)) {
  
  #get only the data on the relevant chromosome and strand
  useDf <- subsetList[[as.character(rowData["motifChrom"])]]
  subsetallGenes <- useDf #[useDf$Strand == as.numeric(rowData["motifStrand"]), ]
  #print(head(subsetallGenes))
  #cycle through every gene on the right chromosome and strand
  
  for(i in 1:nrow(subsetallGenes)) {
    #if the motif start or end are present within that gene 
    
    if( all(as.numeric(rowData["motifStart"]) >= as.numeric(subsetallGenes$Gene.Promoter.Start[i]), as.numeric(rowData["motifEnd"]) <= subsetallGenes$Gene.Promoter.End[i])) {
      
      #select the gene you are at now from all the genes, select the good row in the totalDf
      
      geneToAnnotate <- as.character(subsetallGenes[i,"Ensembl Gene ID"])
      currentGene    <- totalDf[totalDf$Ensembl_ID == geneToAnnotate, ]
      #print(currentGene)
      #print(names(currentGene[,"MotifCounts"][[1]][[1]]))
      #Sys.sleep(2)
      #presentMotifs <- as.character(totalDf[totalDf$Ensembl_ID == geneToAnnotate, "Motifs"])
      
      #if the motif is not yet present in the gene-motif table, or the motif table is as yet empty, add the new motif with count 1
      if(!newMotif %in% names(currentGene[, "MotifCounts"][[1]][[1]])) {
        
        #the line below selects the list item belonging to the gene of interest, and then selects the first non-filled spot in the list, gives it the correct motifname, and sets count to 1 
        names(totalDf[totalDf$Ensembl_ID == geneToAnnotate,][,"MotifCounts"][[1]][[1]])[which(is.na(names(totalDf[totalDf$Ensembl_ID == geneToAnnotate,][,"MotifCounts"][[1]][[1]])) == TRUE)[1]] <<- newMotif
        totalDf[totalDf$Ensembl_ID == geneToAnnotate,][,"MotifCounts"][[1]][[1]][totalDf[totalDf$Ensembl_ID == geneToAnnotate,][, "MotifCounts"][[1]][[1]] == 0][1] <<- 1
        
        
      } else  {
        
        #select the correct count cell and add 1 to it:
        #specifically: select the list element in counts that corresponds to that element in motifs
        countToReplace <- totalDf[totalDf$Ensembl_ID == geneToAnnotate,][, "MotifCounts"][[1]][[1]][which(names(totalDf[totalDf$Ensembl_ID == geneToAnnotate,][,"MotifCounts"][[1]][[1]]) == newMotif)]
        totalDf[totalDf$Ensembl_ID == geneToAnnotate,][, "MotifCounts"][[1]][[1]][which(names(totalDf[totalDf$Ensembl_ID == geneToAnnotate,][,"MotifCounts"][[1]][[1]]) == newMotif)] <<- countToReplace + 1
        #print(countToReplace)
        
      } 
      
    }
    
  }
  }
  #print("it works")
}


allGenes <- fread("/linuxhome/tmp/dieter/Grch37p13_genes.txt", sep = ",", header = T)
#allGenes <- fread("~/Documents/Project/Programming/Transcription Factor Binding/HumanWG_motifs/copytomutant/Grch37p13_genes.txt", sep = ",", header = T)
nrow(allGenes)


#encode the Chromosome Name column as Chromosome.Name
#Take only the 'real' chromosomes (i.e. not the patches etc.)
#Then remove duplicate entries

for (i in 1:length(names(allGenes))) {
  if (names(allGenes)[i] == "Chromosome Name") {
    names(allGenes)[i] <- "Chromosome.Name"}}
allGenes$Chromosome.Name <- as.factor(allGenes$Chromosome.Name)
head(allGenes)

#genes to test on are those lying on the relevant chromosomes
#I will take all the numbered, X, Y and MT chromosomes only, because those are the ones to which genes are mapped
chromosomes <- grep("^\\d+|^MT|^X|^Y", levels(allGenes$Chromosome.Name), value = TRUE)
relevantAllGenes <- allGenes[allGenes$Chromosome.Name %in% chromosomes,]
relevantAllGenes <- relevantAllGenes[!duplicated(relevantAllGenes$"Ensembl Gene ID"), ]
nrow(relevantAllGenes)

#Include promoter regions in the data frame to work with.
#subtract 2000 from the start to include promoter regions. Different on the minus strand.
#Note: check on which strand you are
relevantAllGenes$Gene.Promoter.Start <- rep(0, nrow(relevantAllGenes))
relevantAllGenes$Gene.Promoter.End <- rep(0, nrow(relevantAllGenes))
#positive strand
relevantAllGenes[relevantAllGenes$Strand == 1,]$Gene.Promoter.Start    <- (relevantAllGenes[relevantAllGenes$Strand == 1,]$"Gene Start (bp)" - 2000)
relevantAllGenes[relevantAllGenes$Strand == 1,]$Gene.Promoter.End      <- (relevantAllGenes[relevantAllGenes$Strand == 1,]$"Gene Start (bp)" + 2000)
#negative strand
relevantAllGenes[relevantAllGenes$Strand == -1,]$Gene.Promoter.End     <- (relevantAllGenes[relevantAllGenes$Strand == -1,]$"Gene End (bp)" + 2000)
relevantAllGenes[relevantAllGenes$Strand == -1,]$Gene.Promoter.Start   <- (relevantAllGenes[relevantAllGenes$Strand == -1,]$"Gene End (bp)" - 2000)
print(head(relevantAllGenes))

#pre-subset the data into a list per chromosome for more efficient searching (in theory)
subsetList <- vector("list", length(chromosomes))
for (i in 1:length(chromosomes)) {
  subsetList[[i]] <- relevantAllGenes[relevantAllGenes$Chromosome.Name == chromosomes[i], ]
}
names(subsetList) <- chromosomes

#make the data frame that will hold the final data
maxMotifCount = 250
motifList <- list(rep("", maxMotifCount))
countList <- list(rep(0, maxMotifCount))
names(countList[[1]]) <- NA 
names(countList[[1]])
totalDf <- data.frame(relevantAllGenes[,"Ensembl Gene ID"],
                      I(rep(countList, nrow(relevantAllGenes))))
names(totalDf) <- c("Ensembl_ID", "MotifCounts")
totalDf <- as.data.table(totalDf)
head(totalDf)
names(totalDf[, "MotifCounts"][[1]][[1]] == 0)

#now actual reading begins

#get the file length
motifFileToCount <- file(matchesFile, "rb")
#motifFileToCount <- file("~/Documents/Project/Programming/Transcription Factor Binding/HumanWG_motifs/copytomutant/1000matches.txt", "rb")
totalLines  <- R.utils::countLines(motifFileToCount)
close(motifFileToCount)

numberOfRowsPerTime <- 1000

#open the file for reading
motifFileToRead <- file(matchesFile, "r")
#motifFileToRead <- file("~/Documents/Project/Programming/Transcription Factor Binding/HumanWG_motifs/copytomutant/1000matches.txt", "r")

for (i in 1:ceiling((totalLines[1]/numberOfRowsPerTime))) {
  
  print(paste("Working on slice:", i, "out of", ceiling((totalLines[1]/numberOfRowsPerTime)), "total slices."))

  motifDataSlice            <- read.csv(motifFileToRead, nrow = numberOfRowsPerTime, sep = " ", header = FALSE)
  head(motifDataSlice)
  colnames(motifDataSlice)  <- c("motif","chromosome","start", "end", "strand")
  motifDataSlice$chromosome <- gsub("chr", "", motifDataSlice$chromosome)
  motifDataSlice$strand     <- gsub("+", "1", motifDataSlice$strand, fixed = TRUE)
  motifDataSlice$strand     <- gsub("-", "-1", motifDataSlice$strand, fixed = TRUE)
  head(motifDataSlice)

  apply(motifDataSlice, MARGIN = 1, FUN = mapTFBS)

}

close(motifFileToRead)
head(totalDf)


##concatenate the format to:
#EnsG000345 | "Mfg-2|GGG-3|Horp_6_known-3" --> this is an example of the format. | as separator, -number at the end to indicate amount

makeFinalDfMotifs = function(x,
                             motifSep = "|",
                             countSep = "-",
                             showYourWork = FALSE) {
                               

  currentEnsembl_ID <- x[[1]]
  counts <- x[[2]]
  motifs <- names(counts)
  
  usedCountMotifs <- counts[counts != 0]
  if(length(usedCountMotifs) >= 1 & showYourWork == TRUE) {
    print(usedCountMotifs)
  }
  
  
  #get the final data frame from the outside environment
  #assign the motifs and counts in the correct way
  
  #print(head(df))
  motifCountPairs <- paste0(names(usedCountMotifs), countSep, as.character(usedCountMotifs), collapse = motifSep)
  if (motifCountPairs != "-" & showYourWork == TRUE) {
    print(paste0("motifCountPairs: ", motifCountPairs))
  }
  
  #if there are no motifs in this gene, give it NA. Otherwise, give the motifs.
  if (motifCountPairs == "-") {
    singleLineFinalDf[singleLineFinalDf$Ensembl_ID == currentEnsembl_ID, 2] <<- ""
  } else {
  
  singleLineFinalDf[singleLineFinalDf$Ensembl_ID == currentEnsembl_ID, 2] <<- motifCountPairs
  
  }
  
  #print("Motifs assigned to singleLineFinalDf")
  
}
  



#make the final datafile
singleLineFinalDf <- data.table(Ensembl_ID = relevantAllGenes[,"Ensembl Gene ID"],
                               Motifs = rep("", nrow(relevantAllGenes)))
colnames(singleLineFinalDf) <- c("Ensembl_ID","Motifs")
head(singleLineFinalDf)

apply(totalDf, MARGIN = 1, FUN = makeFinalDfMotifs)

#write.csv(totalDf, file = "/linuxhome/tmp/dieter/HumanWholeGenomeMotifOccurrences_listformat.csv",
#          row.names = FALSE)
write_csv(singleLineFinalDf, path = outputFile)
