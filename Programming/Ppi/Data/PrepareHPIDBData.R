#HPIDB2.0 data analysis

HPIDBInitial <- fread("~/Documents/Project/Programming/Ppi/Data/HPIDB2_completedata.txt", header = T, sep = "\t")
head(HPIDBInitial)
tail(HPIDBInitial)
nrow(HPIDBInitial)
nrow(HPIDBInitial[HPIDBInitial$protein_taxid_1 == "9606",])
HPIDBInitial <- HPIDBInitial[HPIDBInitial$protein_taxid_2_cat == "VIRUS",]

databaseTypes <- str_match(HPIDBInitial$protein_xref_1, "(\\w+):(.+)")

uniqueHPIDB <- HPIDBInitial[!duplicated(HPIDBInitial$protein_xref_1_unique),]
nrow(uniqueHPIDB)
head(uniqueHPIDB)




uniqueDataBaseTypes <- str_match(viralInteractionsonly$protein_xref_1_unique, "(\\w+):(.+)")
uniqueViralProteinXRefs <- str_match(viralInteractionsonly$protein_xref_2_unique, "(\\w+):(.+)")
nrow(uniqueDataBaseTypes)
head(uniqueDataBaseTypes)
unique(uniqueDataBaseTypes[,2])
length(uniqueDataBaseTypes[uniqueDataBaseTypes[,2] == "INTACT",][,2])
nrow(uniqueViralProteinXRefs)
head(uniqueViralProteinXRefs)
unique(uniqueViralProteinXRefs[,2])

humanDataBaseAndIdentifier <- data.frame(humandatabase = rep("", nrow(uniqueDataBaseTypes)),
                                         humanidentifier = rep("", nrow(uniqueDataBaseTypes)))
humanDataBaseAndIdentifier$humandatabase <- uniqueDataBaseTypes[, 2]
humanDataBaseAndIdentifier$humanidentifier <- uniqueDataBaseTypes[, 3]
head(humanDataBaseAndIdentifier)
humanDataBaseAndIdentifier$viralinteractortaxid <- viralInteractionsonly$protein_taxid_2
humanDataBaseAndIdentifier$viralinteractorproteindatabase <- uniqueViralProteinXRefs[,2]
humanDataBaseAndIdentifier$viralinteractorproteinxref <- uniqueViralProteinXRefs[,3]
humanDataBaseAndIdentifier$proteinseq1 <- viralInteractionsonly$protein_seq1
orderedByDataBase <- humanDataBaseAndIdentifier[order(humanDataBaseAndIdentifier$humandatabase),]
head(orderedByDataBase)
tail(orderedByDataBase)

intactInteractions <- orderedByDataBase[orderedByDataBase$humandatabase == "INTACT",]
str(intactInteractions)

#write_csv(intactInteractions, "~/Documents/Project/Programming/Ppi/Data/manuallyAnnotateINTACTIDs.csv")

?write_csv

#Lots of finnicky stuff going on in that database...pseudogenes and what not. See file manuallyAnnotateINTACTIDs

nonINTACTorderedByDataBase <- orderedByDataBase[orderedByDataBase$humandatabase != "INTACT",]
head(nonINTACTorderedByDataBase)
unique(nonINTACTorderedByDataBase$humandatabase)
nonINTACTorderedByDataBase <- nonINTACTorderedByDataBase[order(nonINTACTorderedByDataBase$humanidentifier),]
head(nonINTACTorderedByDataBase)

ensemblMart <- useMart("ensembl",dataset="hsapiens_gene_ensembl")

ensemblIDsHumanGenesHPIDB <- getBM(filters = c("biotype", "uniprot_gn"),
                 attributes = c("ensembl_gene_id", "uniprot_gn"),
                 values = list("protein_coding", nonINTACTorderedByDataBase$humanidentifier),
                 mart = ensemblMart)
head(ensemblIDsHumanGenesHPIDB)
ensemblIDsHumanGenesHPIDB <- ensemblIDsHumanGenesHPIDB[order(ensemblIDsHumanGenesHPIDB$uniprot_gn),]
head(ensemblIDsHumanGenesHPIDB)
nrow(ensemblIDsHumanGenesHPIDB)

#now, I have the ensemblIDs for the human proteins with uniprot identifiers. I add to that the from the INTACT database
INTACTdata <- fread("~/Documents/Project/Programming/Ppi/Data/manuallyAnnotateINTACTIDs.csv")
INTACTdata <- filter(INTACTdata, EnsemblID != "NO")
INTACTdata
INTACTdataClean <- INTACTdata[, c("EnsemblID","UniprotID")]
names(INTACTdataClean) <- c("ensembl_gene_id", "uniprot_gn")


totalEnsemblIDs <- rbind(ensemblIDsHumanGenesHPIDB, INTACTdataClean)
totalEnsemblIDs$measuredPPI <- "yes"
head(totalEnsemblIDs)
str(totalEnsemblIDs)

#check that data is unique and make it so if not
length(unique(totalEnsemblIDs$ensembl_gene_id))
nrow(totalEnsemblIDs)

totalEnsemblIDs <- filter(totalEnsemblIDs, !duplicated(ensembl_gene_id))
length(unique(totalEnsemblIDs$ensembl_gene_id))
nrow(totalEnsemblIDs)
#unique now. save this data.

write_csv(totalEnsemblIDs, "~/Documents/Project/Programming/Ppi/Data/HPIDBDATAforEnrichmentAnalysis.csv")

