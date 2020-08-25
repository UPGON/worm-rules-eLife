# 18.12.2018
# Radek Jankele
#
# Functions for manipulation of C. elegans lineage data
# and for plotting
# Saving and subseting lineages 

# read file with cell fates
CellF <- read.csv("R_Scripts/CellFates.csv")
CellF <- CellF[-grep("l1",CellF$Stage),]

colors <- c("#4d804d", "#56c5c5", "#9B9B9B", "#474747", "#9D3787","#D2D266","#D1833B","#374DA5","#212121", "#212121")
names(colors)<-c("ABal","ABar","ABpl","ABpr","MS","E","C","D","P4","P3")
colors <- colors[sort(names(colors))]

addTrans <- function(color,trans){
  # This function adds transparency to a color.
  # Define transparancy with an integer between 0 and 255
  # 0 being fully transparent and 255 being fully visable
  # Works with either color and trans a vector of equal length,
  # or one of the two of length 1.
  
  if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
  if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
  if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))
  
  num2hex <- function(x)
  {
    hex <- unlist(strsplit("0123456789ABCDEF",split=""))
    return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
  }
  rgb <- rbind(col2rgb(color),trans)
  res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
  return(res)
}

sourcePartial <- function(fn,startTag='#from here',endTag='#to here') {
  lines <- scan(fn, what=character(), sep="\n", quiet=TRUE)
  st<-grep(startTag,lines)
  en<-grep(endTag,lines)
  tc <- textConnection(lines[(st+1):(en-1)])
  source(tc)
  close(tc)
}

starformat <- function(val){
  v<-ifelse(val<0.1,ifelse(val<0.05,ifelse(val<0.01,ifelse(val<0.001,"***","**"),"*"),"."),"n.s.")
  v
}

Tuk <- function(var, groupvar, data) { #TukeyHSD returnig compact letter display in addition for given variable
  forml <- eval(substitute(v ~ gr, list(v = as.name(var), gr = as.name(groupvar))))
  mod <- TukeyHSD(aov(forml, data=data))
  pv <- mod[[groupvar]][,4]
  comp <- multcompLetters(pv)$Letters
  levls <- levels(data[,groupvar])
  comp <- comp[match(levls, names(comp))]
  output <- list(formula = forml, model = mod, comp = comp)
  attributes(output$formula) <- NULL
  return(output)
}

toWideTable <- function(df, anot.cols, var.cols, startCell){
  #restructure data to creat cell.variable combination as columns
  tempDf <- df[grep(paste0(startCell,"$"), df$Cell),]
  colnames(tempDf)[var.cols] <- paste0(startCell,".",colnames(tempDf)[var.cols])
  WideDf <- tempDf[,c(anot.cols, var.cols)]
  
  for (c in sort(unique(df$Cell),decreasing=T)){
    if(c %in% Cellorder[nchar(Cellorder)<=8]){  
      if(c!=startCell){
        tempDf <- df[df$Cell==c,c(which(colnames(df)=="embryo"),var.cols)]
        print(paste(c,":",nrow(tempDf)))
        if(nrow(tempDf)>5){
          colnames(tempDf)[2:ncol(tempDf)] <- paste0(c,".",colnames(tempDf)[2:ncol(tempDf)])
          WideDf <- merge(WideDf,tempDf,by="embryo",all.x=T)
        }
      }
    }
  }
  rownames(WideDf) <- WideDf$embryo
  emptyCols <- which(apply(WideDf,2,function(x) all(is.na(x), na.rm=T)))
  WideDf <- WideDf[,-emptyCols]
  # should be 0 in most embryos
  WideDf <- WideDf[,-which(colnames(WideDf) %in% c("ABal.StartTime","ABar.StartTime","ABpl.StartTime","ABpr.StartTime"))]
  WideDf
}

groupplot <- function(var, test="t.test"){
  outcome <- relevel(outcome,'hatched')
  
  pv.w <- wilcox.test(mat[,var]~outcome)$p.value
  pv.t <- t.test(mat[,var]~outcome)$p.value
  
  pv <- c(pv.w,pv.t)
  pv.print <- ifelse(pv < 0.01, formatC(pv, format = "e", digits = 1), signif(pv,digits=3))
  
  boxplot(mat[,var]~outcome, main=var, ylab=NA,xlab = NA, outpch="", lwd=0.5)
  
  if(test=="wilcox") mtext(paste0("wilcox P=",pv.print[1]), cex=0.6)  else mtext(paste0("t-test P=",pv.print[2]), cex=0.6) 
  if(test=="wilcox") mtext(starformat(pv[1]), cex=1, line=-1.3)  else mtext(starformat(pv[2]), cex=1, line=-1.3) 
  
  
  beeswarm(mat[,var]~outcome, add=T, spacing=1, corral="random", cex=0.8, pwcol=pal[3:4][outcome], pch=16)
  cat(var, "\twilcox p =",pv.print[1],starformat(pv[1]),"\tt.test p =",pv.print[2],starformat(pv[2]),"\n")
}

plot.bxp <- function (filename, data, var1, ylab1, plot.width=1.5, saveit=T,...){
  if(saveit)pdf(paste0(today,filename),width = plot.width, height = 2.1, pointsize = 10, useDingbats=F)
  par(mar=c(2.5,2.5,0.1,0.15), mgp=c(1.25,0.5,0), cex.lab=0.7, cex.axis=0.6, cex.main=0.8, tcl=-.35, las=3)
  
  v1 <- data[,var1]
  g <- data$Group
  boxplot(v1~g, lwd=0.5, outlwd=1, staplewex=0.25, whisklty=1, outpch=NA,ylab=ylab1, xlab=NA,...)
  beeswarm(v1~g, add=T, pwcol=pal[g], pch=16, cex=.6, corral="random", corralWidth=1.15 ,spacing=1)
  if(saveit)dev.off()
}

compare.bxp <- function (filename, data, var1, var2, ylab1, ylab2, plot.width=3, saveit=T){
  if(saveit) pdf(paste0(today,filename),width =plot.width, height = 2.5, pointsize = 10, useDingbats=F)
  par(mar=c(2,2.3,0.2,0.15), mgp=c(1.25,0.5,0), cex.lab=0.7, cex.axis=0.6, cex.main=0.8, tcl=-.35)
  par(mfrow=c(1,2))
  
  v1 <- data[,var1]
  v2 <- data[,var2]
  g <- data$Group
  boxplot(v1~g, lwd=0.5, outlwd=1, staplewex=0.25, whisklty=1,outpch=NA,ylab=ylab1, xlab=NA)
  beeswarm(v1~g, add=T, col=pal[2:4], pch=16, cex=.8, corral="random", spacing=1)
  
  boxplot(v2~g, lwd=0.5, outlwd=1, staplewex=0.25, whisklty=1,outpch=NA,ylab=ylab2, xlab=NA)
  beeswarm(v2~g, add=T, col=pal[2:4], pch=16, cex=.8, corral="random")
  if(saveit) dev.off()
}

plot.corr <- function (filename, data, xvar, yvar, xlab="x", ylab="y", groupvar="Group",title="correlation", plot.width=2.5, save=T, leg=F, xlim=F, ylim=F){
  data <- as.data.frame(data)
  
  if(save==T){ 
    pdf(filename,width = plot.width, height = plot.width*1.05, pointsize = 10, family='Helvetica', useDingbats = FALSE)
    par(mar=c(2.3,2.3,1,0.2), mgp=c(1.25,0.5,0), cex.lab=0.8, cex.axis=0.7, cex.main=0.8, tcl=-.35)
  }
  
  v1 <- data[,xvar]
  v2 <- data[,yvar]
  g <- data[,groupvar]
  
  xlim <- if(all(xlim==F)) xlim <- range(v1, na.rm = T)
  ylim <- if(all(xlim==F)) xlim <- range(v1, na.rm = T)
  
  plot(v2~v1, col=pal[g], pch=20, xlab=xlab, ylab=ylab, main=title, xlim=xlim, ylim=ylim)
  m <- lm(v2~v1)
  ct <- cor.test(v2,v1, method = "pearson")
  
  if(ct$p.value < 0.05){
    abline(m, lwd=1.5)
    
    newx <- seq(min(v1, na.rm = T)*0.9,max(v1, na.rm = T)*1.1, length.out=length(v1))
    prd <- predict.lm(m,newdata=data.frame(v1=newx),interval = c("confidence"),level = 0.90, type="response")
    lines(newx,prd[,2],col="grey",lty=2)
    lines(newx,prd[,3],col="grey",lty=2)
  }
  mtext(paste(c("R=",round(unlist(ct[c("estimate")]), digits=2),", P=", signif(ct$p.value, digits=-1)), collapse = ""), side=3, cex=0.7, line=-1.2) #adj=0.1, padj=0, )
  if(leg==T)legend("bottomright",levels(g),pch=16, cex=0.8, col=pal[factor(levels(g))], inset=0.02)
  
  if(save==T) dev.off()
}

rbind.all.columns <- function(x, y) {
  
  x.diff <- setdiff(colnames(x), colnames(y))
  y.diff <- setdiff(colnames(y), colnames(x))
  
  x[, c(as.character(y.diff))] <- NA
  
  y[, c(as.character(x.diff))] <- NA
  
  return(rbind(x, y))
}

stat.test.by.row <- function(tab, f1, f2, type){
  if(type=="t.test"){
    p.val <- apply(tab[2:(nrow(Embs)+1)], 1, function(x, f1, f2) t.test(x[f1][!is.na(x[f1])],x[f2][!is.na(x[f2])])$p.value, f1=f1, f2=f2)
  } else if(type=="wilcox"){
    p.val <- apply(tab[2:(nrow(Embs)+1)], 1, function(x, f1, f2) wilcox.test(x[f1][!is.na(x[f1])],x[f2][!is.na(x[f2])])$p.value, f1=f1, f2=f2)
  } else print("unknown test")
  names(p.val) <- tab$CellName
  p.val[p.val>0.1] <- NaN
  return(p.val)
}

se <- function(x, na.rm=FALSE){
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}
row.stats <- function(x) c(m=mean(x,na.rm=TRUE), sd=sd(x,na.rm=TRUE), n=length(which(!is.na(x))), sem=se(x, na.rm=T))

getCellFates <- function (CellNames){
  CellFates <- CellF$CellName
  names(CellFates) <- CellF$Generalized
  sapply(CellNames, function(x) paste(unique(names(CellFates)[grep(x, CellFates)]), collapse=", "))
}

dist3D <- function(co) sqrt( (co[4] - co[1])^2 + (co[5] - co[2])^2 + (co[6] - co[3])^2 )

cell.groups <- function(cells){
  group <- rep("Early", length(cells))
  group[grepl("P4|Z",cells)] <- "P4"
  group[grepl("MS",cells)] <- "MS"
  group[grepl("E",cells)] <- "E"
  group[grepl("C",cells)] <- "C"
  group[grepl("D",cells)] <- "D"
  group[grepl("ABar",cells)] <- "ABar"
  group[grepl("ABpr",cells)] <- "ABpr"
  group[grepl("ABpl",cells)] <- "ABpl"
  group[grepl("ABal",cells)] <- "ABal"
  group[cells=="EMS"] <- "Early"
  group <- factor(group)
  group
}

get.motherCell <- function(Daughter){
  #deduce and fix parent names
  irregular <- c("P0","P1","P2","P3","P4","Z1","Z2","EMS","E","MS","C","D","AB")
  Daughter <- as.character(Daughter)
  parent <- mapvalues(Daughter ,c("AB","C" ,"D" ,"E" , "EMS", "MS","P1","P2","P3","P4","Z2","Z3"),
                      c("P0","P2","P3","EMS","P1" ,"EMS","P0","P1","P2","P3","P4","P4"),warn_missing =F)
  s <- ! parent %in% irregular
  parent[s] <- substr(parent[s],1, nchar(parent[s])-1)
  parent
}

get.sisterCell <- function(sibling){
  #deduce and fix parent names
  sibling <- as.character(sibling)
  sister <- mapvalues(sibling ,c("AB","C" ,"D" ,"E" , "EMS", "MS","P1","P2","P3","P4","Z2","Z3","ABal","ABar","ABpl","ABpr","Eal","Ear","Epl","Epr"),
                               c("P1","P3","P4","MS", "P2" , "E", "AB","EMS","C","D","Z3","Z2", "ABar","ABal","ABpr","ABpl","Ear","Eal","Epr","Epl"),warn_missing =F)

  #if sibling doesn't mach any value in the mapvalue, sister==sibling is TRUE
    sister = sub("a$","p",sister)
    b = sister==sibling
    
    sister[b] = sub("p$","a",sister[b])
    
  sister
}

#function expects data frame with column Cell
lineage.subset <- function(dd, root){
  #subset dd for cells in the root sublineage
  if(grepl("(AB)|[CD]", root)){
    sub <- dd[grepl(root,dd$Cell),]
  } else if(grepl("P", root)) {
    if(root=="P0") sub <- dd
    if(root=="P1") sub <- dd[grepl("(P[1234])|[EMZCD]",dd$Cell),]
    if(root=="P2") sub <- dd[grepl("(P[234])|[ZCD]",dd$Cell),]
    if(root=="P3") sub <- dd[grepl("(P[34])|[ZD]",dd$Cell),]
    if(root=="P4") sub <- dd[grepl("(P4)|[Z]",dd$Cell),]
  } else if(root =="EMS"){
    sub <- dd[grepl("[EM]",dd$Cell),]
  } else {
    sub <- dd[grepl(root,dd$Cell)&!dd$Cell=="EMS",]
  }
  sub
} 

CellCount <- function(Tab, trim){
  ##calculate number of cells over time
  startT <- ceiling(min(Tab$StartTime, na.rm = T))
  endT <- ceiling(max(Tab$StartTime, na.rm = T))
  
  ##nCells at the first frame
  nCells=length(which(Tab$StartTime<=startT))
  nCellsVector=nCells;
  
  for(t in (startT+1):endT){
    #index starting from 2:
    i <- t-startT+1;
    ##how many rows has StartTime <= then t?
    nLastT <- length(which(Tab$StartTime<=(t-1)))
    n <- length(which(Tab$StartTime<=t))
    
    ##if nCells is different then at previous timepoint, there was a division. 
    #We take number of new cells and substract half of their number from nCells, because mothers no longer exist
    if(nLastT!=n){
      newDaughters=n-nLastT
      nCells=nCells+ceiling(newDaughters/2) #if there is even number of new cells, we should keep +1
    }
    ##record number of cells at each timepoint
    
    nCellsVector[i] <- nCells
  }
  ##strip the vector so it starts with a single occurence of 4C stage
  FirstDiv=min(which(nCellsVector>trim), na.rm = T)-1
  nCellsVector <- nCellsVector[FirstDiv:length(nCellsVector)]
  return(nCellsVector)
}

rescale <- function (input, len) approx(seq_along(input), input, n = len)$y #fuction that interpolates vector to different lenght 

load.bigTable <-function (filepath){
  
  #filepath <- "/Users/jankele/switchdrive/ev571_embryos/PositionalAnalysis/Lineages2Statisctics/ContrastUnequalEqualized_Results/bigTableN.txt"
  data<-read.table(filepath,sep="\t",header=TRUE,row.names=1)
  
  #Compute lifetime for each cell
  ltSub <-data[,grepl("_Tdeath",colnames(data[1,]))]-data[,grepl("_Tbirth",colnames(data[1,]))] # lifeTimes Subtable
  colnames(ltSub) <- gsub("Tdeath","LifeTime",colnames(ltSub))
  data <- cbind(data,ltSub)
  
  #colnames(data)[grep("dead.1",colnames(data))]
  
  #get all relevant column indices
  Tbirth<-grepl("_Tbirth",colnames(data[1,]))
  Tdeath<-grepl("_Tdeath",colnames(data[1,]))
  
  angleAP<-grepl("_aAP",colnames(data[1,]))
  angleDV<-grepl("_aDV",colnames(data[1,]))
  angleLR<-grepl("_aLR",colnames(data[1,]))
  angleMean<-grepl("_aMean",colnames(data[1,]))
  
  posAP<-grepl("_pAP",colnames(data[1,]))
  posDV<-grepl("_pDV",colnames(data[1,]))
  posLR<-grepl("_pLR",colnames(data[1,]))
  posOV<-grepl("_pOV",colnames(data[1,]))
  
  lifeTime<-grepl("_LifeTime",colnames(data[1,]))
  netdis<-grepl("_netdis",colnames(data[1,]))
  totdis<-grepl("_totdis",colnames(data[1,]))
  
  #set values for undivided cells to NA
  notDivided=is.na(data[,grep("_aAP",colnames(data[1,]))])
  
  sub<-data[,Tdeath]
  sub[notDivided] <- NA
  data[,Tdeath] <- sub
  
  sub<-data[,netdis]
  sub[notDivided]=NA
  data[,netdis]=sub
  
  sub<-data[,totdis]
  sub[notDivided]=NA
  data[,totdis]=sub
  
  sub<-data[,lifeTime]
  sub[notDivided]=NA
  data[,lifeTime]=sub
  
  sub<-data[,posOV]
  sub[notDivided]=NA
  data[,posOV]=sub
  
  data$Cell <- rownames(data)
  data$parent <- get.motherCell(data$Cell)
  data <- data[order(data$Cell),]
  data
}

SaveNuclei <- function (MyEmbryo, Dir){
  NewDir <- paste0(Dir, "/nuclei")
  if(!dir.exists(NewDir)) dir.create(NewDir, recursive = T)
  
  for(i in 1:length(MyEmbryo)){
    write.table(MyEmbryo[[i]], file=paste0(NewDir,"/t",sprintf("%03.f", i),"-nuclei"), quote=F, sep=",", col.names = F, row.names = F)
    
  }
  print(paste0("Nuclei files were saved to: ", getwd(),"/",NewDir))
}

log10Tck <- function(side, type){
  lim <- switch(side, 
                x = par('usr')[1:2],
                y = par('usr')[3:4],
                stop("side argument must be 'x' or 'y'"))
  at <- floor(lim[1]) : ceiling(lim[2])
  return(switch(type, 
                minor = outer(1:9, 10^(min(at):max(at))),
                major = 10^at,
                stop("type argument must be 'major' or 'minor'")
  ))
}
