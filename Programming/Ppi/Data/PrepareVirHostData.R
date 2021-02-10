##make virhostnetdata tractable

source("https://bioconductor.org/biocLite.R")
biocLite("biomaRt")
library(biomaRt)
source("https://bioconductor.org/biocLite.R")
biocLite("org.Hs.eg.db")
library(org.Hs.eg.db)

virHostData <- fread("~/Documents/Project/Programming/Ppi/Data/VirhostNetData", header = F, sep = "\t")
head(virHostData)
colnames(virHostData) <- c("InteractorA", "InteractorB", "AliasA", "AliasB", "OtherAliasA", "OtherAliasB",
                           "MeasurementTechnology", "Authornames", "PublicationIdentifier", "TaxidA",
                           "TaxidB", "InteractionType", "SourceDbForInteractions", "SourceDbIdentifiers", "ConfidenceScore")
head(virHostData)


#filter only human:non-human and non-human:human interactions

viralPPIOnly <- filter(virHostData, TaxidA != "taxid:9606" | TaxidB != "taxid:9606")
head(viralPPIOnly)
viralHumanPPIOnly <- filter(viralPPIOnly, TaxidA != "taxid:9606" & TaxidB == "taxid:9606" | TaxidA == "taxid:9606" & TaxidB != "taxid:9606" )
head(viralHumanPPIOnly)
"taxid:9606" %in% viralHumanPPIOnly$TaxidB

viralHumanPPIOnly$InteractorA <- str_match(viralHumanPPIOnly$InteractorA, "uniprotkb:(.+)" )[,2]
viralHumanPPIOnly$InteractorB <- str_match(viralHumanPPIOnly$InteractorB, "uniprotkb:(.+)" )[,2]
#rownames(viralHumanPPIOnly) <- seq(1, nrow(viralHumanPPIOnly))
head(viralHumanPPIOnly)

#now swap it such that the human ID is always interactor A.

wrongWayRound <- filter(viralHumanPPIOnly, TaxidB == "taxid:9606")
viralInteractor <- wrongWayRound$InteractorA
humanInteractor <- wrongWayRound$InteractorB
viralTaxid      <- wrongWayRound$TaxidA
humanTaxid      <- wrongWayRound$TaxidB
viralAlias     <- wrongWayRound$AliasA
humanAlias     <- wrongWayRound$AliasB
viralOtherAlias     <- wrongWayRound$OtherAliasA
humanOtherAlias     <- wrongWayRound$OtherAliasB

wrongWayRound$InteractorA <- humanInteractor
wrongWayRound$InteractorB <- viralInteractor
wrongWayRound$TaxidA      <- humanTaxid
wrongWayRound$TaxidB      <- viralTaxid
wrongWayRound$AliasA      <- humanAlias
wrongWayRound$AliasB      <- viralAlias
wrongWayRound$OtherAliasA <- humanOtherAlias
wrongWayRound$OtherAliasB <- viralOtherAlias

head(wrongWayRound)

viralHumanPPIOnly[viralHumanPPIOnly$SourceDbIdentifiers %in% wrongWayRound$SourceDbIdentifiers,] <- wrongWayRound
viralHumanPPIOnly[16406,]
"taxid:9606" %in% viralHumanPPIOnly$TaxidB
head(viralHumanPPIOnly)

#so it is now all in the correct orientation.
nrow(viralHumanPPIOnly)
length(unique(viralHumanPPIOnly$InteractorA))

#there are NAs in the data, since some proteins apparently had no uniprot identifiers, but do have
#NP_identifiers (NCBI). See below. 

#are there NA's in the interactorA column?
anyNA(viralHumanPPIOnly$InteractorA)
anyNA(viralHumanPPIOnly$InteractorB)

viralHumanWithUniprot <- filter(viralHumanPPIOnly, is.na(InteractorB) == FALSE)
nrow(viralHumanWithUniprot)
viralHumanWithoutUniprot <- filter(viralHumanPPIOnly, is.na(InteractorB) == TRUE) 
nrow(viralHumanWithoutUniprot)
sum(nrow(viralHumanWithUniprot), nrow(viralHumanWithoutUniprot))
nrow(viralHumanPPIOnly)

#taken together, you see that some have NA, but they have an NP_identifier.

#now to make counts for every gene
viralHumanWithUniprot$InteractionMerge <- paste0(viralHumanWithUniprot$InteractorA, "-", viralHumanWithUniprot$InteractorB)
head(viralHumanWithUniprot)
countViralHumanWithUniprot <- count(viralHumanWithUniprot, InteractionMerge)

#There are instances of things measured more than once
countViralHumanWithUniprot[countViralHumanWithUniprot[,2] >1,]
countViralHumanWithUniprot[countViralHumanWithUniprot[,2] >4,]

#if something is in the data multiple times, what is different?
viralHumanWithUniprot[viralHumanWithUniprot$InteractionMerge == "O60563-P04608",]

#as it turns out, it is measured in different literature sources then, but the miscore
#is consistently the same. This means I can take the unique occurences of interactions
#just fine, and still know something about the trust I should vest in them.

viralHumanWithoutUniprot$InteractionMerge <- paste0(viralHumanWithoutUniprot$InteractorA, "-", viralHumanWithoutUniprot$AliasB)
viralHumanWithoutUniprot$InteractionMerge


#What I want:
#for every gene, whether it interacts with a viral protein or not and the miscore
#for every gene, the viral proteins that it interacts with, and the miscore.

#objective 1:
uniqueInteractionsOnly <- viralHumanPPIOnly[!duplicated(viralHumanPPIOnly$InteractorA),]
head(uniqueInteractionsOnly)
nrow(uniqueInteractionsOnly)

#Get the EnsemblIDs
#not working. Why?
ensemblMart <- useMart("ensembl",dataset="hsapiens_gene_ensembl")



VirHostDataEnsemblDataFrame <- getBM(filters = c("biotype", "uniprot_gn"),
         attributes = c("ensembl_gene_id", "uniprot_gn"),
         values = list("protein_coding", uniqueInteractionsOnly$InteractorA),
         mart = ensemblMart)
VirHostDataEnsemblDataFrame$measuredPPI <- "yes"
head(VirHostDataEnsemblDataFrame)



class(VirHostDataEnsemblDataFrame)
#are there any duplicates caused by uniprot ids mapping to the same ensembl gene id?
length(unique(VirHostDataEnsemblDataFrame$ensembl_gene_id))
#yes, so make unique and save
VirHostDataEnsemblDataFrame <- VirHostDataEnsemblDataFrame[!duplicated(VirHostDataEnsemblDataFrame$ensembl_gene_id),]
nrow(VirHostDataEnsemblDataFrame)
write_csv(VirHostDataEnsemblDataFrame, "~/Documents/Project/Programming/Ppi/Data/virhostnetuniqueForAnalysis.csv", col_names = TRUE)
?write_csv

#Ben hier gebleven. Probleem waar ik nu tegenaan loop is dat uniprot identifiers naar verschillende ensembl IDs gemapt kunnen
#worden. Zoals bijvoorbeeld de uniprotding: 
# http://www.uniprot.org/uniprot/P27824
#met als ENsembl IDs:
#http://www.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000127022;r=5:179678628-179730925 
#http://www.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000283777;r=CHR_HG30_PATCH:179679612-179731798 
#zie ook hier: https://www.biostars.org/p/48361/ 

#Oke we lossen het op door gewoon alle ensemblIDs te accepteren.
#dat is dus een 'benefit of the doubt approach'

#objective 2: not done yet.