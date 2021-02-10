##Human macrophage gene upregulation after activation with LPS and IFN-gamma



###########################################
###########################################
##                                       ##
##      Loading packages and data        ##
##                                       ##
###########################################
###########################################





library("pacman")
pacman::p_load(ggplot2, dplyr, stringr, purrr, readr, LaF, ff, data.table, R.utils, base, utils, devtools, microbenchmark, reshape2, magrittr, tidyr)
source("https://bioconductor.org/biocLite.R")
biocLite("biomaRt")


options(stringsAsFactors = FALSE)
`%ni%` <- Negate(`%in%`)

##step 1: Get EnsemblIDs and normalize all M! repeats to their respective RM repeats at correct time points
##then take log 2 values and save
#output: 3 files, format: SPot ID \t Gene Symbol (ENsemblID) \t timepoint 1 log2 norm. exp value \t tp 2 log2 norm. exp. value \t tp 3 log2 norm. expr. val.
#one file for each repeat

rawMacrophageData <- fread("~/Documents/Project/Programming/Human macrophage activation/Data/GSE57614_series_matrix_editedForLoadingIntoR.csv")
head(rawMacrophageData)
allEnsemblGenes <- fread("~/Documents/Project/Programming/DataForAll/EnsemblBioMartIDs_GRCh38.p10.txt", sep = ",", header = T)


#get only relevant columns
m1RMMacrophageData <- rawMacrophageData %>% select(-dplyr::contains("M2a"), -dplyr::contains("M2c"))
head(m1RMMacrophageData)
tail(m1RMMacrophageData)
ncol(m1RMMacrophageData)
rep1Data <- m1RMMacrophageData %>% select(-dplyr::contains("Rep"))
nrow(m1RMMacrophageData)

##########################################################################################
##########################################################################################
##                                                                                      ##
##        TRANSLATE AGILENT PROBES TO HGNC-SYMBOLS AND GENERATE FORMAT FOR STEM         ##
##                                                                                      ##
##########################################################################################
##########################################################################################


ensembl = useMart("ENSEMBL_MART_ENSEMBL", "hsapiens_gene_ensembl", host = "apr2018.archive.ensembl.org")
listAttributes(ensembl)
attrQuery   <- c("agilent_wholegenome_4x44k_v2", "description", "hgnc_symbol")
filterQuery <- c("agilent_wholegenome_4x44k_v2")
filterValue <- c(m1RMMacrophageData$ID_REF)

ensemblIDsAgilent <- getBM(attrQuery,filterQuery,filterValue,ensembl)
head(ensemblIDsAgilent, 40)
ensemblIDsAgilent <- ensemblIDsAgilent %>% filter(description != "" & hgnc_symbol != "")
head(ensemblIDsAgilent, 100) ; tail(ensemblIDsAgilent, 100)
ensemblIDsAgilent <- ensemblIDsAgilent %>% arrange(agilent_wholegenome_4x44k_v2)
ensemblIDsAgilent
#checked this after outputting files at the bottom to .tsv, to see whether numbering was still correct. Seems to be so! (:)
ensemblIDsAgilent %>% filter(agilent_wholegenome_4x44k_v2 == "A_23_P101434")
#0's are required if no gene is associated with a probe for the program STEM
notInEnsemblDataAgilentIDs <- data.frame(agilent_wholegenome_4x44k_v2 = m1RMMacrophageData$ID_REF[m1RMMacrophageData$ID_REF %ni% ensemblIDsAgilent$agilent_wholegenome_4x44k_v2],
                                         description = "0", hgnc_symbol = "0")
head(notInEnsemblDataAgilentIDs); length(notInEnsemblDataAgilentIDs)

#check this data more thoroughly. Down the line, JPT1 gets an empty row and a row with the JPT1 hgnc symbol. Why?
ensemblIDsAgilent %>% filter(hgnc_symbol == "JPT1")
ensemblIDsAgilent %>% filter(agilent_wholegenome_4x44k_v2 == "A_23_P100632")
#so that just yields two separate probes to JPT1, which seems fine.
#one probe is unique. FIltering step above for description and hgnc_symbol != "" did the trick
#what about the notInEnsemblDataAgilentIDs?
notInEnsemblDataAgilentIDs %>% filter(agilent_wholegenome_4x44k_v2 == "A_23_P100632" | agilent_wholegenome_4x44k_v2 == "A_23_P118435")
#okay so up until now everything is correct and A_23_P100632 is not in there twice


totalProbeEnsemblConversionTable <- ensemblIDsAgilent %>% full_join(notInEnsemblDataAgilentIDs) %>% arrange(agilent_wholegenome_4x44k_v2)
head(totalProbeEnsemblConversionTable, 100)
totalProbeEnsemblConversionTable$hgnc_symbol[duplicated(totalProbeEnsemblConversionTable$hgnc_symbol)]
#the below is okay, that's just two probes to the same gene
totalProbeEnsemblConversionTable[totalProbeEnsemblConversionTable$hgnc_symbol == "VPS72",]

tail(totalProbeEnsemblConversionTable)
#check for duplicates
totalProbeEnsemblConversionTable[duplicated(totalProbeEnsemblConversionTable$agilent_wholegenome_4x44k_v2)|duplicated(totalProbeEnsemblConversionTable$agilent_wholegenome_4x44k_v2, fromLast = TRUE),]
##still not unique, combine double HGNC symbols with ; as sep....possibly with group_by. Also combine descriptions.
totalProbeEnsemblConversionTable <- totalProbeEnsemblConversionTable %>% group_by(agilent_wholegenome_4x44k_v2) %>% 
  mutate(allGenesPerProbe = paste0(hgnc_symbol, collapse = ";"), allDescriptionsPerProbe = paste0(description, collapse = ";"))
head(totalProbeEnsemblConversionTable, 100); totalProbeEnsemblConversionTable %>% filter(agilent_wholegenome_4x44k_v2 == "A_23_P16157")
#now, we have probes, all the genes they map to, and all the descriptions of those genes they map to. Filter on unqiue agilent probes.

uniqueAgilentProbesDescriptionsHGNCSymbols <- totalProbeEnsemblConversionTable[!duplicated(totalProbeEnsemblConversionTable$agilent_wholegenome_4x44k_v2),]
#check that#
uniqueAgilentProbesDescriptionsHGNCSymbols %>% filter(agilent_wholegenome_4x44k_v2 == "A_23_P16157") %>% as.data.frame()
uniqueAgilentProbesDescriptionsHGNCSymbols %>% filter(agilent_wholegenome_4x44k_v2 == "A_23_P104116") %>% as.data.frame()

##########################################################################################
##########################################################################################
##                                                                                      ##
##        DE-LOG10 THE DATA, NORMALIZE TO APPROPRIATE REPEAT CONTROL VALUE, LOG2        ##
##                                                                                      ##
##########################################################################################
##########################################################################################

head(m1RMMacrophageData)
colnames(m1RMMacrophageData)

#problem: log2(negative number) cannot be computed. Remedy by adding 12 to every value (lapply(m1RMMacrophageData, FUN = min)) gives -10.105 as min data

rawDataMacrophageWithoutLog10 <- m1RMMacrophageData %>% select(-ID_REF) %>% as.data.frame()
#quick and dirty loop
for(row in 1:nrow(rawDataMacrophageWithoutLog10)) {
  
  for(col in 1:ncol(rawDataMacrophageWithoutLog10)) {
    
    rawDataMacrophageWithoutLog10[row, col] = 10^(rawDataMacrophageWithoutLog10[row, col])
  }
  
  
}

head(rawDataMacrophageWithoutLog10)
rawDataMacrophageWithoutLog10 <- rawDataMacrophageWithoutLog10 %>% mutate(ID_REF = m1RMMacrophageData$ID_REF) %>% select(ID_REF, everything())
head(rawDataMacrophageWithoutLog10)

rep1Score6h  <- log2(as.numeric(unlist(rawDataMacrophageWithoutLog10[,5]))/as.numeric(unlist(rawDataMacrophageWithoutLog10[,2])))
rep1Score12h <- log2(as.numeric(rawDataMacrophageWithoutLog10[,6])/as.numeric(rawDataMacrophageWithoutLog10[,3]))
rep1Score24h <- log2(as.numeric(rawDataMacrophageWithoutLog10[,7])/as.numeric(rawDataMacrophageWithoutLog10[,4]))

rep2Score6h   <- log2(as.numeric(rawDataMacrophageWithoutLog10[,"GSM1515657_M1_6hRep2_251485022397"])/as.numeric(rawDataMacrophageWithoutLog10[,"GSM1515651_RM_6hRep2_251485022401"]))
rep2Score12h  <- log2(as.numeric(rawDataMacrophageWithoutLog10[,"GSM1515659_M1_12hRep2_25148502239"])/as.numeric(rawDataMacrophageWithoutLog10[,"GSM1515653_RM_12h_Rep2_2514850223"]))
rep2Score24h  <- log2(as.numeric(rawDataMacrophageWithoutLog10[,"GSM1515661_M1_24hRep2_25148505307"])/as.numeric(rawDataMacrophageWithoutLog10[,"GSM1515655_RM_24hRep2_25148505307"]))

rep3Score6h   <- log2(as.numeric(rawDataMacrophageWithoutLog10[,"GSM1515658_M1_6hRep3_251485022399"])/as.numeric(rawDataMacrophageWithoutLog10[,"GSM1515652_RM_6hRep3_251485053074"]))
rep3Score12h  <- log2(as.numeric(rawDataMacrophageWithoutLog10[,"GSM1515660_M1_12hRep3_25148505307"])/as.numeric(rawDataMacrophageWithoutLog10[,"GSM1515654_RM_12hRep3_25148505307"]))
rep3Score24h  <- log2(as.numeric(rawDataMacrophageWithoutLog10[,"GSM1515662_M1_24hRep3_1_4.txt.gz"])/as.numeric(rawDataMacrophageWithoutLog10[,"GSM1515656_RM_24hRep3_25148502238"]))

head(rep3Score6h); head(rawDataMacrophageWithoutLog10[,"GSM1515658_M1_6hRep3_251485022399"]); head(rawDataMacrophageWithoutLog10[,"GSM1515652_RM_6hRep3_251485053074"])



######################################
######################################
##                                  ## 
##            SAVE DATA             ##
##                                  ## 
######################################
######################################





#make the data files
dFRep1 <- data.frame(Spot = rawDataMacrophageWithoutLog10$ID_REF, GeneSymbolHGNC = uniqueAgilentProbesDescriptionsHGNCSymbols$allGenesPerProbe,
                     sixHours = rep1Score6h, twelveHours = rep1Score12h, twentyFourHours = rep1Score24h)

dFRep2 <- data.frame(Spot = rawDataMacrophageWithoutLog10$ID_REF, GeneSymbolHGNC = uniqueAgilentProbesDescriptionsHGNCSymbols$allGenesPerProbe,
                     sixHours = rep2Score6h, twelveHours = rep2Score12h, twentyFourHours = rep2Score24h)

dFRep3 <- data.frame(Spot = rawDataMacrophageWithoutLog10$ID_REF, GeneSymbolHGNC = uniqueAgilentProbesDescriptionsHGNCSymbols$allGenesPerProbe,
                     sixHours = rep3Score6h, twelveHours = rep3Score12h, twentyFourHours = rep3Score24h)

write_tsv(dFRep1, path = "~/Documents/Project/Programming/Human macrophage activation/Data/STEMData1.tsv")
write_tsv(dFRep2, path = "~/Documents/Project/Programming/Human macrophage activation/Data/STEMData2.tsv")
write_tsv(dFRep3, path = "~/Documents/Project/Programming/Human macrophage activation/Data/STEMData3.tsv")


##############################CODE HERE IS MANUAL CHECKING OF THIS SCORE ONLY. SEE BELOW FOR CODE THAT IS USED FOR THE FINAL BAYESIAN###



######################################
######################################
##                                  ## 
##  Checking nrichments of profiles ##
##                                  ## 
######################################
######################################



fileList <- list(higher7  = "~/Documents/Project/Programming/Human macrophage activation/Data/Profile7_STEManalysis1",
                 higher11 = "~/Documents/Project/Programming/Human macrophage activation/Data/Profile11_STEManalysis1",
                 higher15 = "~/Documents/Project/Programming/Human macrophage activation/Data/Profile15_STEManalysis1",
                 lower5   = "~/Documents/Project/Programming/Human macrophage activation/Data/Profile5_STEManalysis1",
                 lower0   = "~/Documents/Project/Programming/Human macrophage activation/Data/Profile0_STEManalysis1",
                 zigzag2  = "~/Documents/Project/Programming/Human macrophage activation/Data/Profile2_STEManalysis1",
                 zigzag9  = "~/Documents/Project/Programming/Human macrophage activation/Data/Profile9_STEManalysis1")
fileList

profileList <-purrr::map( fileList, .f = fread, header = TRUE)




head(profileEleven)

#translate spots to ENSEMBL IDs
attrQueryProf11   <- c("agilent_wholegenome_4x44k_v2", "ensembl_gene_id")
filterQueryProf11 <- c("agilent_wholegenome_4x44k_v2")
filterValueProf11 <- c(profileEleven$Spot)

ensemblIDsAgilentProf11 <- getBM(attrQueryProf11,filterQueryProf11,filterValueProf11,ensembl)

positiveSet <- fread("~/Documents/Project/Programming/DataForAll/PositiveList_2.csv")
positiveSetIdCol <- "Ensembl_accession"
positiveSetAmount <- nrow(positiveSet)
allEnsemblGenes <- fread("~/Documents/Project/Programming/DataForAll/EnsemblBioMartIDs_GRCh38.p10.txt", sep = ",", header = T)
allEnsemblGenesIDCol <- "Gene stable ID"
allEnsemblGenesNonPositive <- allEnsemblGenes %>% filter(`Gene stable ID` %ni% positiveSet$Ensembl_accession)
allEnsemblGenesNonPositiveSetAmount <- nrow(allEnsemblGenesNonPositive)

sum(positiveSet$Ensembl_accession %in% ensemblIDsAgilentProf11$ensembl_gene_id)
sum(ensemblIDsAgilentProf11$ensembl_gene_id %in% allEnsemblGenesNonPositive$`Gene stable ID`)

##make this into a function
#note: uses benefit of the doubt-approach. If multiple Ensembl IDs map to one probe, assumed to all be in there!

rm(data)


listForFisherTest <- purrr::map(profileList, .f = function(data) {
  
  #get corresponding Ensembl IDs
  filterValues <- data$Spot
  attrQueryProf   <- c("agilent_wholegenome_4x44k_v2", "ensembl_gene_id")
  filterQueryProf <- c("agilent_wholegenome_4x44k_v2")
  filterValueProf <- filterValues
  ensemblIdsProf  <- getBM(attrQueryProf,filterQueryProf,filterValueProf,ensembl)
  
  #positive set in profile
  positiveSetInProfile <- sum(positiveSet$Ensembl_accession %in% ensemblIdsProf$ensembl_gene_id)
  #negative set in profile
  negativeSetInProfile <- sum(ensemblIdsProf$ensembl_gene_id %in% allEnsemblGenesNonPositive$`Gene stable ID`)
  #which positive set genes are there?
  positiveSetGenesInProfile <- paste0(positiveSet$Ensembl_accession[positiveSet$Ensembl_accession %in% ensemblIdsProf$ensembl_gene_id], collapse = ";")
  
  fisherMatrix <- matrix(c(positiveSetInProfile, positiveSetAmount - positiveSetInProfile,
                           negativeSetInProfile, allEnsemblGenesNonPositiveSetAmount - negativeSetInProfile),
                         nrow = 2,
                         dimnames = list(Profile = c("In this profile", "in the whole genome"),
                                         Set = c("Positive set", "Whole genome/Negative set")))
  fisherTestResult    <- fisher.test(fisherMatrix)
  pValue              <- fisherTestResult["p.value"]
  lowerConfInt        <- unlist(fisherTestResult["conf.int"])[1]
  upperConfInt        <- unlist(fisherTestResult["conf.int"])[2]
  sampleEst           <- fisherTestResult["estimate"]
  
  
  returndF <- data.frame(positiveSetGenesinProfile = positiveSetInProfile,
                         negativeSetGenesInProfile = negativeSetInProfile,
                         pFisher                   = pValue,
                         lowerCIFisher             = lowerConfInt,
                         upperCIFisher             = upperConfInt,
                         sampleEstFisher           = sampleEst,
                         positiveSetIDS            = positiveSetGenesInProfile)
  
  returndF

  
  
  
  
  
})

finalFisherDf <- do.call(rbind, listForFisherTest)
finalFisherDf$p.value <- p.adjust(finalFisherDf$p.value, method = "BH")
finalFisherDf$log2Scores <- log2(finalFisherDf$estimate)
finalFisherDf <- finalFisherDf %>% mutate(profileName = rownames(finalFisherDf)) %>% arrange(rownames(finalFisherDf)) %>% as.data.frame()
finalFisherDf
sort(rep(finalFisherDf$profileName,2))

 geneAmounts = numeric()
 for(rows in 1:nrow(finalFisherDf)) {
   geneAmounts <- c(geneAmounts, c(as.numeric(finalFisherDf$positiveSetGenesinProfile[rows]), as.numeric(finalFisherDf$negativeSetGenesInProfile[rows])))
 }
 geneAmounts
 geneAmounts[1]; geneAmounts[2]; geneAmounts[3] 
 
 
plottingDFMP <- data.frame(profile = sort(rep(finalFisherDf$profileName,2)),
                           set =   rep(c("PS", "NS"), 7),
                           amountOfGenes = geneAmounts)
plottingDFMP <- plottingDFMP %>% arrange(desc(profile))
plottingDFMP
unique(plottingDFMP$profile)
plotProfileBins <- plottingDFMP %>%
  ggplot(aes(x = factor(profile, levels = unique(plottingDFMP$profile)), fill = set, y = amountOfGenes)) + geom_histogram(stat = "identity", colour = "black") 
plotProfileBins 

ggBuild <- ggplot_build(plotProfileBins) 
ggBuild$data
finalFisherSortedDownwards <- finalFisherDf %>% arrange(desc(profileName))

plotProfileBins + annotate("text", x = unique(ggBuild$data[[1]]$x), y = 900, label = round(finalFisherSortedDownwards$log2Scores, digits = 3), angle = 90)


####################################END OF CODE THAT IS SOLELY FOR CHECKING HERE. NOW COMES WHAT IS USED FOR THE BAYESIAN#################



#################################################################################
#                                                                             ###
#   Redo analysis with command line. Load in those files. Save for Bayesian.  ###
#                                                                             ###
#################################################################################

#files generated by using:
#dieter@genera:~/Documents/Project/Programming/Human macrophage activation:java -mx1024M -jar ./stem/stem.jar -b ./Data/AllProfilesSTEManalysis_2 .
#This uses the programme STEM: Short Time-series Expression Miner, to assign profiles to genes with like gene expression over time.


#read in command-line generated files
locationProfileTable <- "~/Documents/Project/Programming/Human macrophage activation/stemoutput/_profiletable.txt"
profileTable <- fread(locationProfileTable)
profileTable

locationGeneTable <- "~/Documents/Project/Programming/Human macrophage activation/stemoutput/_genetable.txt"
geneTable <- fread(locationGeneTable, header = TRUE)
head(geneTable, 40); tail(geneTable, 40)
colnames(geneTable)
#strange 0 column. Let's remove.
geneTable <- geneTable %>% dplyr::select(-`0`)
head(geneTable)
nrow(geneTable)

#translate agilent to ensemblID
ensembl = useMart("ENSEMBL_MART_ENSEMBL", "hsapiens_gene_ensembl", host = "apr2018.archive.ensembl.org")
attrQueryProf   <- c("agilent_wholegenome_4x44k_v2", "ensembl_gene_id")
filterQueryProf <- c("agilent_wholegenome_4x44k_v2", "biotype")
filterValueProf <- list(geneTable$Spot, "protein_coding")
ensemblIdsTotalGeneTable  <- getBM(attrQueryProf,filterQueryProf,filterValueProf,ensembl)
nrow(ensemblIdsTotalGeneTable)




geneTable2 <- geneTable %>% inner_join(ensemblIdsTotalGeneTable, by = c("Spot" = "agilent_wholegenome_4x44k_v2"))
nrow(geneTable2)
head(geneTable2, 20); tail(geneTable2, 30)

nrow(geneTable2)

#Same EnsemblID can have multiple probes. However, it is always in the same profile. I will be doing the Bayesian on profile level,
#but for other ML approaches, expression data for the different time points might be warranted. Therefore, group by EnsemblID, take
#the median, and use that as new values
geneTable2 %>% arrange(ensembl_gene_id) %>% head(70)
geneTable2 %>% filter(ensembl_gene_id == "ENSG00000005075") #two probes in the original data. Is it correctly summarised?
geneTable2 %>% group_by(ensembl_gene_id) %>% summarise(Spots = paste0(Spot, collapse = ";"),
                                                       Profile = paste0(unique(Profile), collapse = ";"),
                                                       sixHours = signif(median(sixHours, na.rm = TRUE), 3),
                                                       twelveHours = signif(median(twelveHours, na.rm = TRUE),3),
                                                       twentyFourHours = signif(median(twentyFourHours, na.rm = TRUE),3)
                                                       ) %>%
  filter(ensembl_gene_id == "ENSG00000005075")
#yes it is correctly summarised. Use this.
geneTable2 %<>% group_by(ensembl_gene_id) %>% summarise(Spots = paste0(Spot, collapse = ";"),
                                                                    Profile = paste0(unique(Profile), collapse = ";"),
                                                                   sixHours = signif(median(sixHours, na.rm = TRUE), 3),
                                                                   twelveHours = signif(median(twelveHours, na.rm = TRUE),3),
                                                                    twentyFourHours = signif(median(twentyFourHours, na.rm = TRUE),3))
      



nrow(geneTable2)
length(geneTable2$ensembl_gene_id)

ensemblIdsNotInTotalGeneTable <- allEnsemblGenes$`Gene stable ID`[allEnsemblGenes$`Gene stable ID` %ni% geneTable2$ensembl_gene_id] %>% data.frame()
ensemblIdsInTotalGeneTable <- allEnsemblGenes$`Gene stable ID`[allEnsemblGenes$`Gene stable ID` %in% geneTable2$ensembl_gene_id] %>% data.frame()
nrow(ensemblIdsInTotalGeneTable); nrow(ensemblIdsNotInTotalGeneTable)
nrow(ensemblIdsInTotalGeneTable) + nrow(ensemblIdsNotInTotalGeneTable)
head(ensemblIdsNotInTotalGeneTable)

colnames(ensemblIdsNotInTotalGeneTable) <- "EnsemblID"
head(ensemblIdsNotInTotalGeneTable)
nrow(ensemblIdsInTotalGeneTable) + nrow(ensemblIdsNotInTotalGeneTable)
nrow(ensemblIdsNotInTotalGeneTable) + nrow(geneTable2); nrow(allEnsemblGenes)
geneTable3 <- geneTable2 %>% full_join(ensemblIdsNotInTotalGeneTable, by = c( "ensembl_gene_id" = "EnsemblID"))
head(geneTable3,50)
nrow(geneTable3)
length(geneTable3$ensembl_gene_id)
length(unique(geneTable3$ensembl_gene_id))
geneTable3 %>% arrange(ensembl_gene_id) %>% tail(50)
#WHY IS THIS NOT THE SAME?
geneTable3$ensembl_gene_id[geneTable3$ensembl_gene_id %ni% allEnsemblGenes$`Gene stable ID`]
nrow(geneTable3[geneTable3$ensembl_gene_id %in% allEnsemblGenes$`Gene stable ID`,])
#THe answer is that there are two which are not in my total Ensembl data file. I have removed them for now.


finalGeneTableForIncorporationIntoBayesian <- geneTable3[geneTable3$ensembl_gene_id %in% allEnsemblGenes$`Gene stable ID`,] %>%
  dplyr::select(ensembl_gene_id, Profile, sixHours, twelveHours, twentyFourHours) %>% rename("EnsemblID" = "ensembl_gene_id") %>%
  arrange(EnsemblID)
head(finalGeneTableForIncorporationIntoBayesian)
#there is one gene that fits to both profile 1 and profile 7. I have assigned it manually to profile 7. This one gene will not make
#that much difference, I think
finalGeneTableForIncorporationIntoBayesian[4531, "Profile"] <- 7

finalGeneTableForIncorporationIntoBayesian %>% write_csv(path = "~/Documents/Project/Programming/Human macrophage activation/Data/MacrophageActivationForBayesianCSV.csv")
#I will delete those two genes, but I don't know what is happening here...how come they exist?? --> answer: Ensembl seems to have updated

