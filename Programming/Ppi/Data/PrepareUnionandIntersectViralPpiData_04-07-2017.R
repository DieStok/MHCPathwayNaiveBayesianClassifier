##Script for creating the Union of Phisto, HPIDB, and VirHostNet
#then calculating enrichment on that. Note, Will make Union as well as
#just throwing everything in with each other (i.e. the most complete set.)


options(stringsAsFactors = FALSE)
`%ni%` <- Negate(`%in%`)

library("pacman")
pacman::p_load(ggplot2, dplyr, stringr, purrr, readr, LaF, ff, data.table, R.utils, base, utils, devtools, microbenchmark)


PhistoDataFinal <- fread("~/Documents/Project/Programming/Ppi/Data/PhistoDataForAnalysis.csv")
HPIDBFinal      <- fread("~/Documents/Project/Programming/Ppi/Data/HPIDBDATAforEnrichmentAnalysis.csv")
VirHostNetFinal <- fread("~/Documents/Project/Programming/Ppi/Data/virhostnetuniqueForAnalysis.csv")

head(PhistoDataFinal)
head(HPIDBFinal)
head(VirHostNetFinal)

PhistoDataFinal$ensembl_gene_id %in% HPIDBFinal$ensembl_gene_id
PhistoDataFinal$ensembl_gene_id %in% VirHostNetFinal$`Gene stable ID`
HPIDBFinal$ensembl_gene_id %in% VirHostNetFinal$`Gene stable ID`

#they are all present in each other, so presumably, same ensembl gene versions. 
#change name VirHostNetFinal for compatibility, add in source
colnames(VirHostNetFinal) <- colnames(PhistoDataFinal)

#add columns presentInPhisto, presentInHPIDB and presentInVirHostNet to all data
PhistoDataFinal$presentInPhisto <- "yes"
PhistoDataFinal$presentInHPIDB  <- ifelse(PhistoDataFinal$ensembl_gene_id %in% HPIDBFinal$ensembl_gene_id, "yes", "no")
PhistoDataFinal$presentInVirHostNet <- ifelse(PhistoDataFinal$ensembl_gene_id %in% VirHostNetFinal$ensembl_gene_id, "yes", "no")

HPIDBFinal$presentInPhisto <- ifelse(HPIDBFinal$ensembl_gene_id %in% PhistoDataFinal$ensembl_gene_id, "yes", "no")
HPIDBFinal$presentInHPIDB  <- "yes"
HPIDBFinal$presentInVirHostNet <- ifelse(HPIDBFinal$ensembl_gene_id %in% VirHostNetFinal$ensembl_gene_id, "yes", "no")

VirHostNetFinal$presentInPhisto <- ifelse(VirHostNetFinal$ensembl_gene_id %in% PhistoDataFinal$ensembl_gene_id, "yes", "no")
VirHostNetFinal$presentInHPIDB <- ifelse(VirHostNetFinal$ensembl_gene_id %in% HPIDBFinal$ensembl_gene_id, "yes", "no")
VirHostNetFinal$presentInVirHostNet <- "yes"

head(PhistoDataFinal)
head(HPIDBFinal)
head(VirHostNetFinal)

TotalUnionIds <- union(union(PhistoDataFinal$ensembl_gene_id, HPIDBFinal$ensembl_gene_id), VirHostNetFinal$ensembl_gene_id)
TotalUnionDf <- union(union(PhistoDataFinal, HPIDBFinal), VirHostNetFinal)
length(TotalUnionIds)
nrow(TotalUnionDf)
nrow(filter(TotalUnionDf, !duplicated(ensembl_gene_id)))
length(PhistoDataFinal$ensembl_gene_id)
length(HPIDBFinal$ensembl_gene_id)
length(VirHostNetFinal$ensembl_gene_id)

TotalIntersectIds <- intersect(intersect(PhistoDataFinal$ensembl_gene_id, HPIDBFinal$ensembl_gene_id), VirHostNetFinal$ensembl_gene_id)
TotalIntersectDf  <- intersect(intersect(PhistoDataFinal, HPIDBFinal), VirHostNetFinal)

#different lengths, why?
nrow(TotalIntersectDf)
length(unique(TotalIntersectIds))
nrow(filter(TotalIntersectDf, !duplicated(ensembl_gene_id)))
nrow(filter(TotalIntersectDf, !duplicated(uniprot_gn)))
TotalIntersectIds[TotalIntersectIds %ni% TotalIntersectDf$ensembl_gene_id]
TotalIntersectIds[TotalIntersectIds %ni% TotalIntersectDf$ensembl_gene_id] %in% VirHostNetFinal$ensembl_gene_id
remarkableIDs <- TotalIntersectIds[TotalIntersectIds %ni% TotalIntersectDf$ensembl_gene_id]
VirHostNetFinal[ensembl_gene_id %in% remarkableIDs, ]
head(TotalIntersectDf)

##I do not know exactly. I will not spend more time to find out why.
##Anyways, I want the intersect of the IDs, and will manually add the meta-information, giving
##The data.table seen below.

IntersectDfBasedOnIds <- data.table(ensembl_gene_id = TotalIntersectIds)
head(IntersectDfBasedOnIds)
IntersectMerge <- merge(IntersectDfBasedOnIds, PhistoDataFinal)
head(IntersectMerge)
tail(IntersectMerge)
length(unique(IntersectMerge$ensembl_gene_id))
nrow(IntersectMerge)
#save data tables IntersectMerge and TotalUnionDf. Enrichment analysis will be done in a different script.

write_csv(IntersectMerge, "~/Documents/Project/Programming/Ppi/Data/IntersectionPPiData.csv")
write_csv(TotalUnionDf, "~/Documents/Project/Programming/Ppi/Data/UnionPPiData.csv")

