##Convert specific motifs from the Pouya data to less specific motifs##

library("pacman")
pacman::p_load(ggplot2, dplyr, stringr, purrr, readr, LaF, ff, data.table, R.utils, base, utils, devtools, microbenchmark, reshape2)
options(stringsAsFactors = FALSE)
oldwd <- getwd()
scriptwd <- "/home/dieter/Documents/Project/Programming/Transcription Factor Binding/HumanWG_motifs/Data/singleLineFinalDf and Unified table/"


EncodeMotifsMultiple <- fread(paste0(scriptwd, "MultipleLinesFilePouyaKheradPour.csv"))
EncodeMotifsSingle <- fread(paste0(scriptwd, "combinedMotifFilePouyaKheradpour.csv"))
testEncodeMotifsSingle <- head(EncodeMotifsSingle, 20)


head(EncodeMotifsMultiple)
head(EncodeMotifsSingle)
EncodeMotifsMultiple$TFBS <- str_match(EncodeMotifsMultiple$TFBS,"([^|]+)_([^|]+)")[,2]

head(EncodeMotifsMultiple)
EncodeMotifsMultipleAggregated <- aggregate(Count~EnsemblID+TFBS, EncodeMotifsMultiple, sum)
EncodeMotifsMultipleAggregated <- EncodeMotifsMultipleAggregated[order(EncodeMotifsMultipleAggregated$EnsemblID),]
head(EncodeMotifsMultipleAggregated)

#make a histogram of motif distributions
motifsPerGene <- aggregate(Count~EnsemblID, EncodeMotifsMultipleAggregated, sum)
head(motifsPerGene)
hist(motifsPerGene$Count, breaks = 50)

write_csv(EncodeMotifsMultipleAggregated, paste0(scriptwd, "underscoresRemovedMultipleLinesPouya.csv"))


conversionvectorSingleLine <- character(length = nrow(EncodeMotifsSingle))

for(i in 1:nrow(EncodeMotifsSingle)) {

  print(i)
  if(EncodeMotifsSingle[i,2] == "") {
    next
  }
disentangleSingleLineEncodeMotifs <- str_match_all(as.character(EncodeMotifsSingle[i,2]), '([A-Za-z0-9:_-]+?)(_\\d+|_\\w+?\\d+)-(\\d+)')


  
  dfEncode <- data.frame(disentangleSingleLineEncodeMotifs[[1]])
  dfEncode <- aggregate(as.numeric(unlist(dfEncode["X4"])) , by = dfEncode["X2"], FUN = sum)
  dfEncode["x"] <- paste0("-", as.character(unlist(dfEncode["x"])))
  conversion <- paste0(as.character(unlist(dfEncode["X2"])), as.character(unlist(dfEncode["x"])), collapse = "|")
  conversionvectorSingleLine[i] <- conversion
  
}

EncodeMotifsSingleAltered <- data.frame(Ensembl_ID = EncodeMotifsSingle$Ensembl_ID, Motifs = conversionvectorSingleLine)
head(EncodeMotifsSingleAltered)
EncodeMotifsSingleAltered <- EncodeMotifsSingleAltered[order(EncodeMotifsSingleAltered$Ensembl_ID),]
head(EncodeMotifsSingleAltered)

write_csv(EncodeMotifsSingleAltered, paste0(scriptwd, "underscoresRemovedSingleLinePouya.csv"))

setwd(oldwd)


#to print long motif parts in slices that can be understood:  
#for (i in seq(1, nchar(as.character(EncodeMotifsSingle[112,2])), 40)) 
#{print(substr(as.character(EncodeMotifsSingle[112,2]), i-40, i))}

#test for anomalies via manual inspection with random sampling
sampleEncode <- sample(1:length(conversionvectorSingleLine), 10)
for (i in 1: length(sampleEncode)) {
  print("--------------------------------------------")
  print("")
  print(conversionvectorSingleLine[i])
  
}
