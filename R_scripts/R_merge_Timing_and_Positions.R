#Read Big data table ####
EmbNames <- read.delim("Results/bigTableNamesN.txt", header = F, col.names = c("Asigned","Original"))
EmbNames$Asigned <- sub("-",".",EmbNames$Asigned)
EmbNames <- EmbNames[order(EmbNames$Asigned),]
EmbNames$IDs <- sapply(strsplit(as.character(EmbNames$Original),"_", fixed = T), `[[`, 1)
rownames(EmbNames) <- EmbNames$IDs

AsignedNames <- sub("_Tdeath","",colnames(bigT[grepl("_Tdeath",colnames(bigT[1,]))]))

#match data in bigT dataframe with aligned timing in Data####
if(length(setdiff(EmbNames$IDs,unique(Data$embryo)))==0) print("Embryos in both tables are matching - good")

# parse the big table and transform it to the long format
# for merging with "Data"
for (name in EmbNames$IDs) {
  asign <- EmbNames[name,1]
  subT <- bigT[,grep(paste0(asign,"_"),colnames(bigT))]
  colnames(subT) <- sub(paste0(asign,"_"),"",colnames(subT))
  subT <- subT[apply(subT, 1, function(x) !all(is.na(x))),]
  subT$Cell <- row.names(subT)
  subT$embryo <- name
  print(paste(name," ncols: ",ncol(subT)))
  print(rownames(subT)[1])
  if(name==EmbNames$IDs[1]){LongDF <- subT}
  else {LongDF <- rbind(LongDF, subT)}
}

rownames(LongDF) <- 1:nrow(LongDF)
colnames(LongDF)[colnames(LongDF)=="LifeTime"] <- "LifeTime.Raw"
#override distances and add angles and pOV, aMean in Data from bigT

#Merge positional information with timing in the Data ####
Mdf <- merge(Data,LongDF, by=c("Cell","embryo"),all.x=T) #create merged Data Frame from Timing data and data from positional analysis
Mdf[,c("Outcome","Group","Experiment","AB.rel","Size")]<- Embs[Mdf$embryo,c("Outcome","Group","Experiment","AB_rel","Size")]
Mdf <- Mdf[,c("embryo","Cell","Outcome","Experiment","Group","Divided","AB.rel","Size","LifeTime","StartTime","EndTime","netdis","totdis","pOV","pAP","pDV","pLR","aMean","aAP","aDV","aLR")]
Mdf[is.na(LifeTime),c("netdis", "totdis")] <- NA
Mdf[Divided==0,c("EndTime", "pLR","pAP","pDV")] <- NA

rm(subT, asign, tab, mat, bigT, CellCountList, CellCountNorm, EmbNames)
