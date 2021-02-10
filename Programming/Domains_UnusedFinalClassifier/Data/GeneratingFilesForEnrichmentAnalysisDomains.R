##Domains analysis
library("pacman")
pacman::p_load(ggplot2, dplyr, stringr, purrr, readr, LaF, ff, data.table, R.utils, base, utils, devtools, microbenchmark, reshape2)


PfamDomains <- fread("~/Documents/Project/Programming/Domains/Data/PfamDomainsBioMart.txt")
colnames(PfamDomains) <- c("Ensembl_Id", "Pfam")
head(PfamDomains)

singleLineTablePfamDomains <- dcast(PfamDomains, Ensembl_Id ~ Pfam)
singleLineTablePfamDomains <- singleLineTablePfamDomains[, -2]
singleLineTablePfamDomains[, 1:30]
singleLineTablePfamDomains[is.na(singleLineTablePfamDomains)] <- ""

headTable <- head(singleLineTablePfamDomains)
paste0(headTable[,], collapse = "|")
headTable[is.na(headTable)] <- ""
head(headTable)[1:6]

lapply(headTable, FUN = function(x) {paste0(x[2:length(x)], collapse = "|") })

FinalSingleLineTablePfamDomains <- data.table(Ensembl_Id = singleLineTablePfamDomains$Ensembl_Id,
                                              Domains = "")
domainsToAdd <- apply(singleLineTablePfamDomains, FUN = function(x) {toPaste <- x[x!= ""]; paste0(toPaste[-1], collapse = "|")}, MARGIN = 1)
tail(domainsToAdd)
length(domainsToAdd)
length(FinalSingleLineTablePfamDomains$Ensembl_Id)
FinalSingleLineTablePfamDomains$Domains <- domainsToAdd
head(FinalSingleLineTablePfamDomains)

write_csv(FinalSingleLineTablePfamDomains, path = "~/Documents/Project/Programming/Domains/Data/PfamDomainsSingleLineTable.csv")
