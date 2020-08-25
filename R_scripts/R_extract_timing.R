#create new empty list to store data.frames
EmbCycles <- {};

#load raw nuclei data from aggregated starryNite files
for(i in 1:length(Embs$ID)){
  id=as.character(Embs$ID[i])
  
  #look into raw data
  Temp <- Nuclei[[id]] #temporary dataframe
  Temp$IDENTITY <- trimws(Temp$IDENTITY)
  AllCells <- unique(Temp$IDENTITY)
  Temp$frame <- as.numeric(as.factor(Temp$timepoint))
  
  #Retreive startTime from RawData -> first occurence of new Cell Name in an ordered time series 
  StartTimes <- aggregate(Temp$frame,list(Temp$IDENTITY),min)
  EndTimes <- aggregate(Temp$frame,list(Temp$IDENTITY),max)
  Timing <- cbind(StartTimes,EndTimes[2])
  colnames(Timing) <- c("Cell","StartTime","EndTime")
  Timing$LifeTime <- Timing$EndTime - Timing$StartTime
  Timing$frame <- Timing$EndTime
  
  #Create new variable indicating cell division / number of daughters
  #if cell has offspring it will be possible to find at least 2+1 cells at least partialially matching the name of the mother (for cells with regular naming)
  Divided <- sapply(AllCells, function(x) length(grep(paste0(x,".$"), AllCells)))
  
  #Fix cells with irregular naming
  Divided[names(Divided) %in% c("AB","P1","EMS","P2")] <- 2 
  if("Z2" %in% AllCells) Divided[names(Divided)=="P4"]<-2
  if("P4" %in% AllCells) Divided[names(Divided)=="P3"]<-2
  Divided[Divided==2] <- 1
  
  Timing$Divided <- Divided[Timing$Cell]
  
  ##correct for different frame intervals - time is now in minutes
  int <- Embs$interval[Embs$ID==id]/60;
  cat(id, " at interval:", int ," min \n");
  Timing[,c("StartTime","EndTime","LifeTime")]<-Timing[,c("StartTime","EndTime","LifeTime")]*(int);
  
  ## set birth of ABa/pl/r as time zero
  Timing[,2:3] <- Timing[,2:3]-Timing[grep("ABal$",Timing$Cell),"StartTime"]
  Timing$LifeTime[Timing$Divided==0] <- NA
  
  # assign(id, inp);
  EmbCycles[[id]] <- Timing;
}

#### Fix 2 embryos starting at 6C ####
#MM4 is problematic, since movie starts as 7C and thus time is shifted by at least 6 minutes
#P2 really divides at frame 2 that should correspond to start time t = 9 min
MM4 <- EmbCycles$MM4

#add EMS which is not present at the frame 1 
    MM4 <- rbind.all.columns(data.frame(Cell="EMS",StartTime=-2, EndTime=0, LifeTime=0,  Divided=1),MM4)
    
    n <- MM4$Cell
    s <- which(!MM4$Cell %in% c("ABal", "ABar", "ABpl", "ABpr", "ABa", "ABp", "P2", "EMS"))
    
    MM4$StartTime <- MM4$StartTime+9 #ABal/r cycle should be 9 minutes longer >> all cells should start +9 minutes
    MM4$EndTime <- MM4$EndTime+9
    
    MM4$LifeTime[n %in% c("ABal", "ABar", "ABpl", "ABpr")] <- MM4$LifeTime[n %in% c("ABal", "ABar", "ABpl", "ABpr")]+9
    MM4$StartTime[n %in% c("ABal", "ABar", "ABpl", "ABpr")] <- 0
    MM4$StartTime[n %in% c("ABa", "ABp", "EMS", "P2")] <- -10
    
    s <- n %in% c("E","MS")
    MM4$StartTime[s] <- 5
    MM4$LifeTime[s] <- MM4$EndTime[s] - MM4$StartTime[s] +1
    MM4$LifeTime[MM4$Divided==0] <- NA
    
    EmbCycles$MM4 <- MM4

rm(MM4, Temp, Timing, StartTimes, EndTimes)

#### Cell count ####

#calculate number of cells at each timepoint
CellCountList <- {}

for(i in 1:length(EmbCycles)){
  ID <- names(EmbCycles)[i]
  print(ID)
  CellCountList[[ID]] <- CellCount(EmbCycles[[ID]],4)
  #names(CellCountList)[i] <- ID;
}

#plot time course for all embryos to compare them
lim <- max(sapply(CellCountList,length))
      # plot(unlist(CellCountList),type="n",xlim=c(0,lim), xlab ="Time [min]", ylab="Number of cells")
      # 
      # x="PM04"
      # lines(1:length(CellCountList[[x]]),CellCountList[[x]], col="black")
      # 
      # mapply(text,CellCountList, x=mapply(length,CellCountList),y=mapply(max,CellCountList),labels=names(CellCountList),cex=0.5,col=seq_along(CellCountList))
      # 
      # mapply(lines,CellCountList,lty=c(3,1,5)[Embs$Experiment], col=seq_along(CellCountList))
      # legend("topleft",names(CellCountList),lty=c(3,1,5)[Embs$Experiment], cex=0.5, col=seq_along(CellCountList))

#plot time course for all embryos to compare them
leg = c("wild-type 17","lin-5(ev571) 17", "Equalized alive", "Equalized dead")
pdf(paste0(outDir,today,"_Growth_RAW.pdf"), width=7, height=6)

    plot(x=0,y=0,xlim=c(0,lim), ylim=c(4,200), mgp=c(2, .7, 0),
         yaxp=c(0, 175, 7), cex.main=1, cex.lab=1, cex.axis=0.9, tck=-.02,
         xlab ="Time [min]",ylab="Cell count", main="Growth curves of embryos")
    mapply(lines,CellCountList[Embs$Group=="alive"], col=pal[3])
    mapply(lines,CellCountList[Embs$Group=="dead"], col=pal[4])
    mapply(lines,CellCountList[Embs$Experiment=="wt"], col=pal[1])
    mapply(lines,CellCountList[Embs$Experiment=="ctrl"], col=pal[2])
    
    
    
    legend("topleft",c(sapply(leg[1:2], function(x) as.expression(substitute(A~degree~"C",list(A = as.name(x))))), leg[3:4]),
           lty=1, lwd=2, cex=0.8, col=pal[c(1,2,3,4)])
dev.off()
print("Completed")