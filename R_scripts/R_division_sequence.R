df <- Complete
df$Group <- as.character(df$Group)
df$Group[df$embryo %in% Embs$ID[inverted]] <- "inverted"
df <- df[-which(df$embryo %in% Embs$ID[partial]),]
df$Group <- factor(df$Group, levels = c("wt","ctrl","alive","dead","inverted"))

#merge with complete ddftaset
df <- aggregate(df[,c("StartTime", "EndTime","LifeTime", "asynchrony")], by=list(df$Group, df$Cell), function(x)c(mean(x, na.rm=T), sd(x, na.rm=T), length(x[!is.na(x)])))
colnames(df)[1:2] <- c("Group","Cell")

cellstoplot <- c("EMS","ABa","P2","ABal","ABar","MS","E","P3","C","ABala","ABpla","MSa","MSaa","Ea","Eaa","Ca","Caa","D","Da","P4","ABalaa","ABalpa","ABalpaa","ABprppa")

pdf(paste0(outDir,today,"_division_Sequence.pdf"),pointsize=11,width=8,height=4.5)
with(df[(df$Cell %in% cellstoplot),],
  {
  x = EndTime;
  sem = x[,2]/sqrt(x[,3]); #SEM
  #labels = paste(Cell,get.sisterCell(Cell),sep="/")
  labels=Cell;
  
  plot(0,-1,
       yaxt="n", 
       xaxt="n",
       ylim=c(0,6), xlim = c(0,200), 
       main = "Division sequence",
       xlab="Cell Division ± SEM [min]", ylab="Group", 
       cex.lab=0.8
       )
  
  #add gridline and axes
  #abline(v=seq(-30, 200, by=10), lty=3, col="Grey");
  
  axis(2,at=c(1:5),labels=levels(Group), srt=90, cex.axis=0.8, adj=1, las=2)
  axis(1,at=seq(0,200,25), cex.axis=0.8, adj=0, las=1)
  
  #add main time-lines
  segments(x0=by(x[,1],Group,min, na.rm=T),x1=200,y0=1:5, y1=1:5, col=pal[1:5], lwd=2)
  
  # plot sem
  #y <- jitter(as.numeric(Group)-0.1, amount=0.1)
  y <- as.numeric(Group)-0.1
  segments(x0=x[,1]-sem,x1=x[,1]+sem,y0=y, y1=y, col="#231F20", lwd=2)
  
  # plot cell division times
  points(x[,1],as.numeric(Group),
         col=pal[Group],
         cex=1, 
         pch="|"
  )
  
  #label cells
    #for wt
    text(x=x[Group=="wt",1],y=1,labels=Cell[Group=="wt"],cex=0.6, pos=2, offset=c(0,1), srt=90)
    #for inverted
    text(x=x[Group=="inverted",1], y=5, labels=Cell[Group=="inverted"], cex=0.6, pos=4, offset=c(0.1,1), srt=90)

  #connect same cells
  by(x[,1],Cell,function(x)lines(x=x,y=c(1:5),lwd="1",col="#A9ABAE", lty=3))
  }
) 

dev.off()
rm(df)
