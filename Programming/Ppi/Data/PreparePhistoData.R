#Prepare Phistodata for analysis

install.packages("pacman", repos = "http://cran-mirror.cs.uu.nl/")
library("pacman")
pacman::p_load(ggplot2, dplyr, stringr, purrr, readr, LaF, ff, data.table, R.utils, base, utils, devtools, microbenchmark)


Phisto <- fread("~/Documents/Project/Programming/Ppi/Data/phi_data.csv", header = T, sep = ",")
head(Phisto)
tail(Phisto)
?grep
#check if all are virus interactions
sum(grepl("virus", Phisto$Pathogen, fixed = TRUE))
nrow(Phisto)

#no, so let us find out what these stragglers are

nonViral <- Phisto[!grepl("virus", Phisto$Pathogen, fixed = TRUE),]

#phages? The phage lambda one seems ill-supported, as is the C. botulinum one. 
#I will leave them in for now. Discuss with John tomorrow.
#Done, was wrong. Removed these strange entries.
#resolution --> 

nonDuplicatePhisto <- Phisto[ !duplicated(Phisto[,5]), ]
nonDuplicateOnlyViralPhisto <- nonDuplicatePhisto[grepl("virus", nonDuplicatePhisto$Pathogen, fixed = TRUE),]
nrow(nonDuplicateOnlyViralPhisto)
str(nonDuplicateOnlyViralPhisto)
head(nonDuplicateOnlyViralPhisto)
names(nonDuplicateOnlyViralPhisto)[5] <- "humanUniprotID"
head(nonDuplicateOnlyViralPhisto)

ensemblMart <- useMart("ensembl",dataset="hsapiens_gene_ensembl")

finalDfPhisto <-  getBM(filters = c("biotype", "uniprot_gn"),
                                   attributes = c("ensembl_gene_id", "uniprot_gn"),
                                   values = list("protein_coding", nonDuplicateOnlyViralPhisto$humanUniprotID),
                                   mart = ensemblMart)
nrow(nonDuplicateOnlyViralPhisto)
nrow(finalDfPhisto)


print(paste0("Difference between unique uniprots and called ENsemblids: ", nrow(finalDfPhisto)-nrow(nonDuplicateOnlyViralPhisto), " more ensembl IDs" ))
finalDfPhisto$measuredPPI <- "yes"
#remove the duplicate Ensembl IDs and save
finalDfPhisto <- finalDfPhisto[!duplicated(finalDfPhisto$ensembl_gene_id), ]
nrow(finalDfPhisto)

write_csv(finalDfPhisto, "~/Documents/Project/Programming/Ppi/Data/PhistoDataForAnalysis.csv")
