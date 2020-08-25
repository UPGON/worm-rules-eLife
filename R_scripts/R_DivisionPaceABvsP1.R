#### Functions #### 
split.AB.P.cells <- function(tab){
  AB.indices <- grep("[AN]", tab$Cell) #Will cath all names starting with A or N - thus AB and New or Nucleus will all be picked up
  ABtab <- tab[AB.indices,]
  ABtab <- ABtab[!is.na(ABtab$StartTime),]
  
  Ptab <- tab[-AB.indices,]
  Ptab <- Ptab[!is.na(Ptab$StartTime),]
  
  AB.count <- CellCount(ABtab,2)
  P.count <- CellCount(Ptab,0)
  
  n <- max(length(AB.count), length(P.count))
  length(AB.count) <- n                      
  length(P.count) <- n
  
  mat <- matrix(c(AB.count,P.count),ncol=2, nrow=length(AB.count))
  colnames(mat) <- c("AB", "P")
  return(mat)
}
row.stats <- function(x) c(m=mean(x,na.rm=TRUE), sd=sd(x,na.rm=TRUE), n=length(which(!is.na(x))), sem=sd(x,na.rm=TRUE)/sqrt(length(which(!is.na(x)))) )

#### Transform data to matrices for AB/P1 #### 
Data <- Complete[!(Complete$embryo %in% Embs$ID[partial]),] #remove partial embryos
#Data <- Data[!Data$embryo%in%c("GM37"),]

Data <- Data[!Data$Cell %in% c("P1","AB","ABa","ABp"),]
Data$StartTime[Data$Cell %in% c("ABar","ABal","ABpl","ABpr")] <- 0
Data$StartTime[Data$Cell %in% c("P2","EMS")] <- 0

time.lim=c(5,210)

split <- by(Data[,c("embryo","Cell","StartTime")], factor(Data$embryo), split.AB.P.cells)

ABmat <- sapply(split, "[", TRUE, 1)
ABmat <- lapply(ABmat, na.omit)
seq.max <- max(sapply(ABmat, length))
ABmat <- sapply(ABmat,"[",1:seq.max)

Pmat <- sapply(split, "[", TRUE, 2)
Pmat <- lapply(Pmat, na.omit)
seq.max <- max(sapply(Pmat, length))
Pmat <- sapply(Pmat,"[",1:seq.max)

id=Embs$ID

hatched=id[equalized&alive]
dead=id[equalized&!alive]
equal=id[equalized]
ctrls=id[Embs$Group=="ctrl"]
wt=id[Embs$Group=="wt"]
inv=id[inverted]
inv <- inv[inv!="GM39"]

#### Plots #### 
#calc aggregated stats for each curve in AB by row - mean, sd, number of observations, and SEM
leglabels <- c(paste("Wild-type n =", length(wt)),
               paste("Unequal controls n =", length(ctrls)), 
               paste("Equal. Alive n =", length(hatched)), 
               paste("Equal. Dead n=", length(dead)),
               paste("Inverted n=", length(inv)))

#line widths
lw = c(1.5,1,1)

opar <- par()
pdf(paste0(outDir,today,"_Growth_ABvsP1_raw.pdf"), width=6, height=3, pointsize = 10)
{
  par(mfcol=c(1,2), mar=c(2.5,2.5,2.1,1.1), mgp=c(1.5,0.5,0),cex.lab=0.8, cex.axis=0.7, cex.main=1, lwd=0.5)
      
      # plot AB descendants raw
      
      matplot(1:nrow(ABmat),ABmat[,wt], xlim=time.lim, type="l",col=pal[1],lty=1, 
              main="Division pace in the AB lineage - Raw", 
              xlab ="Time [min]", 
              ylab="Number of cells")
      matplot(1:nrow(ABmat),ABmat[,hatched], type="l", lty=1, col=pal[3], add=T)
      matplot(1:nrow(ABmat),ABmat[,dead], type="l", lty=1, col=pal[4], add=T)
      matplot(1:nrow(ABmat),ABmat[,ctrls], xlim=time.lim, type="l",col=pal[2],add=T) 
      matplot(1:nrow(ABmat),ABmat[,inv], xlim=time.lim, type="l",col=pal[5],add=T) 
      
      
      legend("topleft",leglabels,lty=c(1,1,1), cex=0.7, col=pal)
      
       # plot P1 descendants raw
      matplot(1:nrow(Pmat),Pmat[,wt], 
              xlim=time.lim, type="l",
              col=pal[1],lty=1, 
              main="Division pace in the P1 lineage - Raw", 
              xlab ="Time [min]", 
              ylab="Number of cells")
      matplot(1:nrow(Pmat),Pmat[,dead], type="l", lty=1, col=pal[4], add=T)
      matplot(1:nrow(Pmat),Pmat[,hatched], type="l", lty=3, col=pal[3], add=T)
      matplot(1:nrow(Pmat),Pmat[,ctrls], type="l", lty=1, col=pal[2], add=T)
      matplot(1:nrow(Pmat),Pmat[,inv], type="l", lty=1, col=pal[5], add=T)
dev.off();      
}

pdf(paste0(outDir,today,"_Growth_ABvsP1_avg.pdf"), width=6, height=6, pointsize = 12)
{
  par(mfrow=c(2,2), mar=c(2.5,2.5,2.1,0.2), mgp=c(1.5,0.5,0),cex.lab=0.9, cex.axis=1, cex.main=1.1)
  
  #plot AB averages
    a <- t(apply(ABmat[,dead], 1, row.stats))
    a <- a[a[,"n"]>2,]
    b <- t(apply(ABmat[,hatched], 1, row.stats))
    b <- b[b[,"n"]>2,]
    c <- t(apply(ABmat[,ctrls], 1, row.stats))
    c <- c[c[,"n"]>2,]
    d <- t(apply(ABmat[,wt], 1, row.stats))
    d <- d[d[,"n"]>2,]
    
    matplot(1:nrow(d),d[,"m"] + outer(d[,"sd"], c(0,1,-1)), type="l", 
            xlim=time.lim,ylim=c(1,max(a[,"m"])),  
            col=pal[1],lty=c(1,3,3), 
            lwd=lw,
            main="Division pace in the AB lineage",  
            xlab ="Time [min]", 
            ylab="Number of cells")
    matplot(1:nrow(b),b[,"m"] + outer(b[,"sd"], c(0,1,-1)), type="l", add=T, lwd=lw, col=pal[3],lty=c(1,3,3), yaxt="n", xaxt="n")
    matplot(1:nrow(a),a[,"m"] + outer(a[,"sd"], c(0,1,-1)), type="l", add=T, lwd=lw,  col=pal[4],lty=c(1,3,3), yaxt="n", xaxt="n")
    matplot(1:nrow(c),c[,"m"] + outer(c[,"sd"], c(0,1,-1)), type="l", add=T, lwd=lw,  col=pal[2],lty=c(1,3,3), yaxt="n", xaxt="n")
    
    legend("topleft",leglabels,lty=c(1,1,1), cex=0.7, col=pal)
    
    # zoom
    matplot(1:nrow(d),d[,"m"] + outer(d[,"sd"], c(0,1,-1)), type="l", 
            xlim=c(145,215),ylim=c(60,max(a[,"m"])),  
            col=pal[1],lty=c(1,3,3), 
            lwd=lw,
            main="Division pace in the AB lineage",  
            xlab ="Time [min]", 
            ylab=NA)
    matplot(1:nrow(b),b[,"m"] + outer(b[,"sd"], c(0,1,-1)), type="l", add=T, lwd=lw, col=pal[3],lty=c(1,3,3), yaxt="n", xaxt="n")
    matplot(1:nrow(a),a[,"m"] + outer(a[,"sd"], c(0,1,-1)), type="l", add=T, lwd=lw,  col=pal[4],lty=c(1,3,3), yaxt="n", xaxt="n")
    matplot(1:nrow(c),c[,"m"] + outer(c[,"sd"], c(0,1,-1)), type="l", add=T, lwd=lw,  col=pal[2],lty=c(1,3,3), yaxt="n", xaxt="n")

    #P1 lineage ####
    #calc aggregated stats for each curve
    a <- t(apply(Pmat[,dead], 1, row.stats))
    a <- a[a[,"n"]>2,]
    b <- t(apply(Pmat[,hatched], 1, row.stats))
    b <- b[b[,"n"]>2,]
    c <- t(apply(Pmat[,ctrls], 1, row.stats))
    c <- c[c[,"n"]>2,]
    d <- t(apply(Pmat[,wt], 1, row.stats))
    d <- d[d[,"n"]>2,]
    e <- t(apply(Pmat[,inv], 1, row.stats))
    e <- e[e[,"n"]>2,]
    
    #plot P1 averages
    matplot(1:nrow(d),d[,"m"] + outer(d[,"sd"], c(0,1,-1)), type="l", 
            xlim=time.lim, ylim=c(1,60), 
            col=pal[1],
            lty=c(1,3,3),
            lwd=lw,
            main="Division pace in the P1 lineage",  
            xlab ="Time [min]", 
            ylab="Number of cells")
    matplot(1:nrow(b),b[,"m"] + outer(b[,"sd"], c(0,1,-1)), type="l", add=T, lwd=lw, col=pal[3],lty=c(1,3,3),yaxt="n", xaxt="n")
    matplot(1:nrow(a),a[,"m"] + outer(a[,"sd"], c(0,1,-1)), type="l", add=T, lwd=lw,  col=pal[4],lty=c(1,3,3),yaxt="n", xaxt="n")
    matplot(1:nrow(c),c[,"m"] + outer(c[,"sd"], c(0,1,-1)), type="l", add=T, lwd=lw,  col=pal[2],lty=c(1,3,3),yaxt="n", xaxt="n")
    matplot(1:nrow(e),e[,"m"] + outer(e[,"sd"], c(0,1,-1)), type="l", add=T, lwd=lw,  col=pal[5],lty=c(1,3,3),yaxt="n", xaxt="n")
    
    #zoom
    matplot(1:nrow(d),d[,"m"] + outer(d[,"sd"], c(0,1,-1)), type="l",
            xlim=c(145,215),ylim=c(30,60), 
            col=pal[1],
            lty=c(1,3,3),
            lwd=lw, 
            main="Division pace in the P1 lineage",  
            xlab ="Time [min]", 
            ylab=NA)
    matplot(1:nrow(b),b[,"m"] + outer(b[,"sd"], c(0,1,-1)), type="l", add=T, lwd=lw, col=pal[3],lty=c(1,3,3),yaxt="n", xaxt="n")
    matplot(1:nrow(a),a[,"m"] + outer(a[,"sd"], c(0,1,-1)), type="l", add=T, lwd=lw,  col=pal[4],lty=c(1,3,3),yaxt="n", xaxt="n")
    matplot(1:nrow(c),c[,"m"] + outer(c[,"sd"], c(0,1,-1)), type="l", add=T, lwd=lw,  col=pal[2],lty=c(1,3,3),yaxt="n", xaxt="n")
    matplot(1:nrow(e),e[,"m"] + outer(e[,"sd"], c(0,1,-1)), type="l", add=T, lwd=lw,  col=pal[5],lty=c(1,3,3),yaxt="n", xaxt="n")
    
    dev.off() 
}
print(paste0("Plots of AB and P1 growth saved as: ", paste0(outDir,today,"_Growth_ABvsP1.pdf")))
par <- opar
