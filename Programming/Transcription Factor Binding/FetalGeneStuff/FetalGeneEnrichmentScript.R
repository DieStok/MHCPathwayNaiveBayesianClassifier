##script for checking fetal cell stuff

FetalvsAPCZScores <- fread("~/Documents/Project/Programming/Transcription Factor Binding/FetalGeneStuff/nature22795-s10_ZscoresDEGbetweenFetalandAdultAPC.csv")
head(FetalvsAPCZScores)

FetalvsAPCDIffExpressionLevels <- fread("~/Documents/Project/Programming/Transcription Factor Binding/FetalGeneStuff/nature22795-s11_ExpressionLevelsDEGbetweenFetalandAdultAPC.csv")
head(FetalvsAPCDIffExpressionLevels)


positiveSet <- fread("~/Documents/Project/Programming/DataForAll/PositiveList_2.csv")
head(positiveSet)
selectionPositiveDifReg <-positiveSet$Ensembl_name %in% FetalvsAPCZScores$Gene
print(positiveSet$Ensembl_name[positiveSet$Ensembl_name %in% FetalvsAPCZScores$Gene])

FetalGenesOfInterestZscores    <- FetalvsAPCZScores[FetalvsAPCZScores$Gene %in% positiveSet$Ensembl_name,]
FetalGenesOfInterestExprLevels <- FetalvsAPCDIffExpressionLevels[FetalvsAPCDIffExpressionLevels$Gene %in% positiveSet$Ensembl_name,]

FetalGenesOfInterestZscores
FetalGenesOfInterestExprLevels

FetalGenesOfInterestExprLevelsNarrow <- melt(FetalGenesOfInterestExprLevels)
FetalGenesOfInterestExprLevelsNarrow$value
head(FetalGenesOfInterestExprLevelsNarrow)
variable2 <- str_match(FetalGenesOfInterestExprLevelsNarrow$variable, "^\\w+")
variable3 <- str_match(FetalGenesOfInterestExprLevelsNarrow$variable, "^\\w+ (.+)")[,2]
FetalGenesOfInterestExprLevelsNarrow$AdultFetal <- variable2
FetalGenesOfInterestExprLevelsNarrow$CellType <- variable3
FetalGenesOfInterestExprLevelsNarrow
subsetFetalNarrow <- FetalGenesOfInterestExprLevelsNarrow[1:10,]




ggplot(data = FetalGenesOfInterestExprLevelsNarrow, aes(x= FetalGenesOfInterestExprLevelsNarrow$Gene, y= value, fill = AdultFetal)) +
  facet_wrap(~CellType, nrow = 3, scales = "free") + geom_bar(stat = "identity", position = "dodge", colour = "black") + theme_minimal() +
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5))

##questions now:
#how old are these fetuses (and is there earlier data?)
#What is the exact meaning of the z-scores (look into how diff. expressed genes are found, i don't really know)
#



  
