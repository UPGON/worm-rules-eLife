
row.stats <- function(x) c(m=mean(x,na.rm=TRUE), sd=sd(x,na.rm=TRUE), n=length(which(!is.na(x))), sem=sd(x,na.rm=TRUE)/sqrt(length(which(!is.na(x)))) )
timestats <- function(mat,subsetrows,k=9){
  m <- t(apply(mat[subsetrows,], 2, row.stats))
  m <- m[m[,"n"]>3,]
  if(k!=0) s=runmed(m[,"m"],k, endrule = "median", algorithm="Turlach") else s=m[,"m"]
  m <- cbind(m, s)
  m
}

#subset source data ####
Complete <- Complete[order(Complete$EndTime),]
#start from 6 cell stage
tempData <- Complete[-grep("(AB|AB[ap]|EMS|P[1-2])$",Complete$Cell),c("pOV","pLR","pDV","pAP","StartTime","EndTime","Cell","embryo","Divided")]
tempData <- tempData[tempData$Divided==1,] #keep only divided cells, since deviation is measured at the end - point
tempData <- tempData[tempData$EndTime<=200,] #end at 180 minutes, after there is too few emryos for most groups

tempData <- tempData[complete.cases(tempData),]  #remove embryos that has only too few timepoints
tempData <- tempData[order(tempData$EndTime),] 

#get cumulative deviation for each embryo ####
deviation <- by(tempData$pOV,tempData$embryo, function(x) cumsum(na.omit(x)))
maxDev <- sapply(deviation, max, na.rm=T) # get max for each embryo

deviation <- deviation[sapply(deviation, length) > 10] # keep only embryos with more than 10 frames of data
endtimes <- by(tempData[,c("pOV","EndTime")], tempData$embryo, function(x) x$EndTime[!is.na(x$pOV)], simplify = T)

# define groups ####
e <- names(deviation)
wt <- which(e %in% Embs$ID[Embs$Group=="wt"])
ctrl <- which(e %in% Embs$ID[Embs$Group=="ctrl"])
dead <- which(e %in% Embs$ID[Embs$Group=="dead"])
equalz.a <- which(e %in% Embs$ID[equalized&alive])
equalz.d <-  which(e %in% Embs$ID[equalized&!alive])
equalz <- which(e %in% Embs$ID[equalized])
invert <- which(e %in% Embs$ID[inverted])
near <-   which(e %in% Embs$ID[partial])
outl <- which(e %in% outliers)
eqd.outl <- setdiff(equalz.d, outl) # equalized dead w.o. outliers

legnd <- function(){
  legend("topleft",c(
  paste("Wildtype n = ", length(wt)), 
  paste("Unequal controls n = ", length(ctrl)), 
  paste("Equalized alive n = ", length(equalz.a)), 
  paste("Equalized dead n = ", length(eqd.outl)),
  paste("Inverted n = ", length(invert))
),
inset=0.01,
lty=1,
lwd=1.5,
cex=0.8, 
col=c(pal[1:4],"#58595B")
)}

# get x axis for ech embryo (time of division)
maxTime <- sapply(endtimes, max, na.rm=T)
minTime <- sapply(endtimes, min, na.rm=T)

#plot raw curves ####
pdf(paste0(outDir,today,"_positional_divergence_raw.pdf"), width=4.5,height = 4.5, pointsize = 10)
{
  #plot deviations
  par(mfcol=c(1,1), mar=c(2.5,3.5,2.1,1.1), mgp=c(1.5,0.5,0), cex.lab=0.8, cex.axis=0.80, cex.main=1)
  plot(x=0,y=0,
    main="Cumulative positional divergence",
    type="n",
    #bty="n",
    #log="y",
    #yaxt="n",
    ylab="Cummulative deviation from \n the mean cell positions in controls [um]",
    xlab="Time from ABa division",
    xlim=c(40,170),
    ylim=c(5,900)# max(maxDev))
  )
  
  # plot individual embryos
  # wt c("GZ01", "GZ02", "GZ03", "GZ04", "GZ05", "GZ06", "GZ07", "GZ08", "GZ09", "GZ10")
  for(id in names(deviation)){
      #if(!id %in% outliers){
    l <- lines(endtimes[[id]],deviation[[id]], col=pal[Embs[id,"Group"]], )
      #}
    }
    
  # overplot outliers
  for(id in outliers){
    #if(!id %in% outliers){
    l <- lines(endtimes[[id]],deviation[[id]], col="black")
    #}
  }
  
  dev.off()
}

#remap values so that they can be coerced to the matrix ####
  #initiate empty matrixes
  devMat <- matrix(nrow=length(deviation), ncol = ceiling(max(maxTime)))
  devMat.norm <- matrix(nrow=length(deviation), ncol = ceiling(max(maxTime)))
  
  rownames(devMat) <- names(deviation)
  rownames(devMat.norm) <- names(deviation)
  
  for(id in names(deviation)){ 
    v <- {}
    v.norm <- {}
    
    for (t in ceiling(minTime[id]):maxTime[[id]]) { # step by 1 minute to get the matrix
      i = max(which(endtimes[[id]] <= t)) #index of value to look for below t 
      v = c(v, deviation[[id]][i])
      v.norm = c(v.norm, deviation[[id]][i]/i) #to normalize by number of cells until this timepoint - N corresponds to i
      names(v)[length(v)] <- t
      names(v.norm)[length(v.norm)] <- t
    }
    devMat[id,as.numeric(names(v))] <- v
    devMat.norm[id,as.numeric(names(v.norm))] <- v.norm
  }
  
  crop.devMatrix <- ceiling(min(minTime))+2

  devMat <- devMat[,crop.devMatrix:ncol(devMat)]
  colnames(devMat) <- crop.devMatrix:(ncol(devMat)+crop.devMatrix-1)
  
  devMat.norm <- devMat.norm[,crop.devMatrix:ncol(devMat.norm)]
  colnames(devMat.norm) <- crop.devMatrix:(ncol(devMat.norm)+crop.devMatrix-1)


  rmFewObservation <- function(mat, lim=4){
    #require more than 4 embryos in each group, otherwise set to NA
    N <- apply(mat, 2, function(x)tapply(x, Embs[rownames(mat),"Group"], function(y) sum(!is.na(y))))
    #expand the logic matrix
    N <- N>lim
    mat[rownames(mat) %in% Embs$ID[Embs$Group=="wt"],!N["wt",]] <- NA
    mat[rownames(mat) %in% Embs$ID[Embs$Group=="ctrl"],!N["ctrl",]] <- NA
    mat[rownames(mat) %in% Embs$ID[Embs$Group=="alive"],!N["alive",]] <- NA
    mat[rownames(mat) %in% Embs$ID[Embs$Group=="dead"],!N["dead",]] <- NA
    mat
  }
  
  devMat <- rmFewObservation(devMat)
  devMat.norm <- rmFewObservation(devMat.norm)
  
#plot averages + SEM #### 
  lw <- c(2.5,1.5,1.5)
  
#calc aggregated stats for each curve

  a <- timestats(devMat, wt)
  b <- timestats(devMat, ctrl)
  c <- timestats(devMat, equalz.a)
  d <- timestats(devMat, equalz.d)
  e <- timestats(devMat, invert)
  f <- timestats(devMat, near)
  g <- timestats(devMat, equalz)
  h <- timestats(devMat, outl)
  j <- timestats(devMat, eqd.outl)# eq. dead and not outliers
  

pdf(paste0(outDir,today,"_cumu_positional_divergence.pdf"), width=4.5,height = 4.5, pointsize = 10)
{
  par(mfcol=c(1,1), mar=c(2.5,2.5,2.1,1.1), mgp=c(1.5,0.5,0), cex.lab=0.8, cex.axis=0.80, cex.main=1)
  #plot averages
  
  #ctrls
  matplot(as.numeric(rownames(b)) , b[,"s"] + outer(b[,"sem"], c(0,1,-1)), type="l", 
          xlim=c(40,165), 
          ylim=c(0,650), 
          col=pal[2],lty=c(1,3,3),
          lwd=lw, 
          # log="y",
          main="Cumulative positional divergence",  
          xlab ="Time from ABa division", 
          ylab="Cumulative divergence from \n mean cell positions in controls [um]")
  #wildtype
  matplot(as.numeric(rownames(a)), a[,"s"] + outer(a[,"sem"], c(0,1,-1)), type="l", add=T, lwd=lw,  col=pal[1],lty=c(1,3,3))

  #equal alive
  matplot(as.numeric(rownames(c)), c[,"s"] + outer(c[,"sem"], c(0,1,-1)), type="l", add=T, lwd=lw, col=pal[3],lty=c(1,3,3))
  #equal dead
  
  matplot(as.numeric(rownames(d)), d[,"s"] + outer(d[,"sem"], c(0,1,-1)), type="l", add=T, lwd=lw, col=pal[4],lty=c(1,3,3))
  
  #inverted
  matplot(as.numeric(rownames(e)), e[,"s"] + outer(e[,"sem"], c(0,1,-1)), type="l", add=T, lwd=lw, col="darkgrey",lty=c(1,3,3))
  
  #Dead without EMS
  #matplot(as.numeric(rownames(j)), j[,"s"] + outer(j[,"sem"], c(0,1,-1)), type="l", add=T, lwd=lw,  col=pal[4],lty=c(1,3,3))
  
  dev.off()
}

#plot normalized divergence by number of cells ####

a <- timestats(devMat.norm, wt)
b <- timestats(devMat.norm, ctrl)
c <- timestats(devMat.norm, equalz.a)
d <- timestats(devMat.norm, equalz.d)
e <- timestats(devMat.norm, invert)
j <- timestats(devMat.norm, setdiff(equalz.d, outl))# eq. dead and not outliers


pdf(paste0(outDir,today,"_positional_divergence_normalized.pdf"), width=2.4,height = 2.4, pointsize = 10)
{
  par(mfcol=c(1,1), mar=c(2.5,2.5,2.1,1.1), mgp=c(1.5,0.5,0), cex.lab=0.8, cex.axis=0.80, cex.main=1)
  
#plot averages

  #ctrls
  matplot(as.numeric(rownames(b)) , b[,"s"] + outer(b[,"sem"], c(0,1,-1)), type="l", 
          xlim=c(35,165), ylim=c(2,7), 
          col=pal[2],lty=c(1,3,3),
          lwd=lw, 
          # log="y",
          main="Normalized positional divergence",  
          xlab ="Time from ABa/p division", 
          ylab="Mean divergence per cell [um]")
  
  #wildtype
  matplot(as.numeric(rownames(a)), a[,"s"] + outer(a[,"sem"], c(0,1,-1)), type="l", add=T, lwd=lw,  col=pal[1],lty=c(1,3,3))
  #equal alive
  matplot(as.numeric(rownames(c)), c[,"s"] + outer(c[,"sem"], c(0,1,-1)), type="l", add=T, lwd=lw, col=pal[3],lty=c(1,3,3))
  #equal dead
  #matplot(as.numeric(rownames(d)), d[,"s"] + outer(d[,"sem"], c(0,1,-1)), type="l", add=T, lwd=lw, col="#BF30A4",lty=c(1,3,3))
  matplot(as.numeric(rownames(d)), d[,"s"] + outer(d[,"sem"], c(0,1,-1)), type="l", add=T, lwd=lw, col=pal[4],lty=c(1,3,3))
  #inverted
  matplot(as.numeric(rownames(e)), e[,"s"] + outer(e[,"sem"], c(0,1,-1)), type="l", add=T, lwd=lw, col="darkgrey",lty=c(1,3,3))
  #dead w.o. outliers
  #matplot(as.numeric(rownames(j)), j[,"s"] + outer(j[,"sem"], c(0,1,-1)), type="l", add=T, lwd=lw,  col=pal[4],lty=c(1,3,3))
  
  dev.off()
}

rm(tempData, maxDev, deviation, endtimes, maxTime, minTime, devMat, devMat.norm, crop.devMatrix, e, a,b,c,d)

