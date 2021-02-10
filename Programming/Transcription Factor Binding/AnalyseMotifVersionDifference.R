######################################################################################################
# Script started 3-08-2017. Author: Dieter Stoker. Function: analyse the difference in motif version
# distribution. I.e. in which (positive set) genes are Nrf2_2 and Nrf2_3 relative to the less specific Nrf?
#
######################################################################################################

EnsemblIdsOriginal29MammalsAnalysis <- fread("~/Documents/Project/Programming/Transcription Factor Binding/EnsemblIdsWGD.csv")

EnsemblIdsOriginalNonUniqueMotifsGeneralised <- fread("~/Documents/Project/Programming/Transcription Factor Binding/EnsemblIdsWGDNonUniqueMotifs.csv")

EnsemblIdsOriginalUniqueMotifsGeneralised <- fread("~/Documents/Project/Programming/Transcription Factor Binding/EnsemblIdsWGDUniqueMotifs.csv")

positiveSet <- fread("~/Documents/Project//Programming/DataForAll/PositiveList_2.csv")
EnsemblIdsPositiveSet <- positiveSet$Ensembl_accession

head(EnsemblIdsOriginal29MammalsAnalysis)
head(EnsemblIdsOriginalUniqueMotifsGeneralised)
head(EnsemblIdsOriginalNonUniqueMotifsGeneralised)

tail(EnsemblIdsOriginal29MammalsAnalysis)
tail(EnsemblIdsOriginalUniqueMotifsGeneralised)
tail(EnsemblIdsOriginalNonUniqueMotifsGeneralised)


Nrf2NonUnique <- EnsemblIdsOriginalNonUniqueMotifsGeneralised$`Nrf-2`[EnsemblIdsOriginalNonUniqueMotifsGeneralised$`Nrf-2` != ""]
Nrf2Unique    <- EnsemblIdsOriginalUniqueMotifsGeneralised$`Nrf-2`[EnsemblIdsOriginalUniqueMotifsGeneralised$`Nrf-2` != ""]
Nrf2_2Original29Mammals <- EnsemblIdsOriginal29MammalsAnalysis$`Nrf-2_2`[EnsemblIdsOriginal29MammalsAnalysis$`Nrf-2_2` != ""]
Nrf2_3Original29Mammals <- EnsemblIdsOriginal29MammalsAnalysis$`Nrf-2_3`[EnsemblIdsOriginal29MammalsAnalysis$`Nrf-2_3` != ""]

length(Nrf2NonUnique)
length(Nrf2Unique)
length(Nrf2_2Original29Mammals) + length(Nrf2_3Original29Mammals)


#welke MHC hebben Nrf-2_3 versus Nrf-2_2?
MHCNrf22 <- EnsemblIdsPositiveSet[EnsemblIdsPositiveSet %in% Nrf2_2Original29Mammals]
positiveSet[positiveSet$Ensembl_accession %in% MHCNrf22,]

MHCNrf23 <- EnsemblIdsPositiveSet[EnsemblIdsPositiveSet %in% Nrf2_3Original29Mammals]
positiveSet[positiveSet$Ensembl_accession %in% MHCNrf23,]
