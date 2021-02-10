##Convert specific motifs from the 29mammals data to less specific motifs##
##last updated 29 september 2017

library("pacman")
pacman::p_load(ggplot2, dplyr, stringr, purrr, readr, LaF, ff, data.table, R.utils, base, utils, devtools, microbenchmark, reshape2)
options(stringsAsFactors = FALSE)

twentyNineMammalsMotifsMultiple <- fread("~/Documents/Project/Programming/Transcription Factor Binding/29mammals/CallingMotifsMyself/Data/FinalDfAndUnifiedMotifFile/MultipleLines29Mammals.csv")
twentyNineMammalsMotifsSingle <- fread("~/Documents/Project/Programming/Transcription Factor Binding/29mammals/CallingMotifsMyself/Data/FinalDfAndUnifiedMotifFile/combinedMotifFile29Mammals.csv")
testtwentyNineMammalsMotifsSingle <- head(twentyNineMammalsMotifsSingle, 20)


head(twentyNineMammalsMotifsMultiple, 20)
head(twentyNineMammalsMotifsSingle, 20)
##noot: dit klopt ook niet! zie pou4f3!!!
twentyNineMammalsMotifsMultiple$TFBS <- str_match(twentyNineMammalsMotifsMultiple$TFBS,"([^|]+)_([^|]+)")[,2]

head(twentyNineMammalsMotifsMultiple, 30)
#the NAs that are introduced are removed in the subsequent analysis by filtering on Counts > -1
twentyNineMammalsMotifsMultipleAggregated <- aggregate(Count~EnsemblID+TFBS, twentyNineMammalsMotifsMultiple, sum)
twentyNineMammalsMotifsMultipleAggregated <- twentyNineMammalsMotifsMultipleAggregated[order(twentyNineMammalsMotifsMultipleAggregated$EnsemblID),]
head(twentyNineMammalsMotifsMultipleAggregated, 30)
tail(twentyNineMammalsMotifsMultipleAggregated, 30)


write_csv(twentyNineMammalsMotifsMultipleAggregated, "~/Documents/Project/Programming/Transcription Factor Binding/29mammals/CallingMotifsMyself/Data/FinalDfAndUnifiedMotifFile/GeneralisedMotifs/underscoresRemovedMultipleLines29mammals.csv")


conversionvectorSingleLine <- character(length = nrow(twentyNineMammalsMotifsSingle))

for(i in 1:nrow(twentyNineMammalsMotifsSingle)) {

  print(i)
  if(twentyNineMammalsMotifsSingle[i,2] == "") {
    next
  }
disentangleSingleLinetwentyNineMammalsMotifs <- str_match_all(as.character(twentyNineMammalsMotifsSingle[i,2]), "([^|]+)-([^|]+)")




  print(disentangleSingleLinetwentyNineMammalsMotifs)
  
  dftwentyNineMammals <- data.frame(disentangleSingleLinetwentyNineMammalsMotifs[[1]])
  dftwentyNineMammals$pureMotif <- dftwentyNineMammals$X2
  print(dftwentyNineMammals$X2)
  linesToChange <- grepl("(.+)(_\\d+)$", dftwentyNineMammals$X2 )
  print(linesToChange)
  if(TRUE %in% linesToChange) {
  dftwentyNineMammals$pureMotif[linesToChange] <- str_match(dftwentyNineMammals[linesToChange, "X2"], "(.+)(_\\d+)$")[,2]
  }
  dftwentyNineMammals <- aggregate(as.numeric(unlist(dftwentyNineMammals["X3"])) , by = list(dftwentyNineMammals$pureMotif), FUN = sum)
  print(dftwentyNineMammals)
  dftwentyNineMammals["x"] <- paste0("-", as.character(unlist(dftwentyNineMammals["x"])))
  conversion <- paste0(as.character(unlist(dftwentyNineMammals[,1])), as.character(unlist(dftwentyNineMammals["x"])), collapse = "|")
  print(conversion)
  conversionvectorSingleLine[i] <- conversion
  
}

twentyNineMammalsMotifsSingleAltered <- data.frame(Ensembl_ID = twentyNineMammalsMotifsSingle$Ensembl_ID, Motifs = conversionvectorSingleLine)
head(twentyNineMammalsMotifsSingleAltered, 50)
twentyNineMammalsMotifsSingleAltered <- twentyNineMammalsMotifsSingleAltered[order(twentyNineMammalsMotifsSingleAltered$Ensembl_ID),]
head(twentyNineMammalsMotifsSingleAltered, 50)
tail(twentyNineMammalsMotifsSingleAltered, 50)

write_csv(twentyNineMammalsMotifsSingleAltered, "~/Documents/Project/Programming/Transcription Factor Binding/29mammals/CallingMotifsMyself/Data/FinalDfAndUnifiedMotifFile/GeneralisedMotifs/underscoresRemovedSingleLine29mammals.csv")



#to print long motif parts in slices that can be understood:  
#for (i in seq(1, nchar(as.character(twentyNineMammalsMotifsSingle[112,2])), 40)) 
#{print(substr(as.character(twentyNineMammalsMotifsSingle[112,2]), i-40, i))}

#test for anomalies via manual inspection with random sampling
sampletwentyNineMammals <- sample(1:length(conversionvectorSingleLine), 10)
for (i in 1: length(sampletwentyNineMammals)) {
  print("--------------------------------------------")
  print("")
  print(conversionvectorSingleLine[i])
  
}
#'([A-Za-z0-9_\\-:]+)(_\\d+|_\\w+?\\d+)-(\\d+)'