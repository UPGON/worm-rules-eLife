#normalize for the speed differences (due to temp, emb size)
#loop over all curves and find one which gives the hihest overal correlation when alligned with all other control curves

#### Time alignment ####
toAlign <- names(CellCountList) #subset to be aligned against the references

ScaledEmbs <- list() #initialize empty list

##Let's find the best reference first, which gives the best correlation score over all
avgCorr <- {}

for(r in references){
  ref <- CellCountList[[r]]
  
  plot(1:length(ref), ref, col=c("black"), ylim=c(0,100), xlim=c(0,150), type="l", lty=3, lwd=2, main = paste("alignment against", r))
  a <- seq(0.85, 1.2, 0.001)
  scale <- {}; #initialize emty vector to take in scaling factors for each embryos
  maxCorr <- {} #initialize emty vector store max correlation for each embryos
  
  for(n in toAlign){
    id <- n
    
    if(id!=r){
      print(paste("aligning: ",id))
      match <- CellCountList[[n]]
      #trim curve that is longer than reference curve
      if(length(match)>rlength) match <- match[1:rlength]
      
      cc <- 0 ##correlation vector
      
      for(i in 1:length(a)){
        l <- ceiling(a[i]*length(match))
        if(l!=length(match))res <- rescale(match, l)
        trim.ref <- ref[1:length(res)]
        cc[i] <- cor(trim.ref, res)
      }
      
      cc <- na.omit(cc)
      #Retrieve factor a resulting in highest correlation value
      maxCorr[n] <- max(cc)
      amax=a[min(which(cc==maxCorr[n]))]
      
      scale[n] <- amax
      
      best.curve <- rescale(match, ceiling(amax*length(match)))
      lines(1:length(best.curve), best.curve, type="l", lty=c(1,3)[Embs$Outcome[Embs$ID==id]], col="gray")
      #text(length(best.curve),max(best.curve),id,cex=0.6)
      
    } else { # if ID is the reference ID
      scale[n] <- 1
    }
  
    #recompute time for all embryos according to the scale
    tab <- EmbCycles[[id]]
    tab[,2:4] <- tab[,2:4]*scale[n]
    
    #Fill the list with corrected tables
    ScaledEmbs[[n]]<-tab
  }
  avgCorr[r] <- mean(maxCorr, na.rm = T)
} # end for loop references

Embs[names(scale),"timeScaling"] <- scale

#### Normalized cell count and plots ####

CellCountNorm <- list()
for(id in Embs$ID){
  CellCountNorm[[id]] <- CellCount(ScaledEmbs[[id]],4)
}

leg = c("Wild-type 17","Controls 17", "Equalized alive","Equalized dead")

pdf(paste0(outDir,today,"_Growth_aligned_RAW.pdf"), width=7, height=6)

    plot(unlist(CellCountNorm),type="n",xlim=c(0,230), ylim=c(4,180), mgp=c(2, .7, 0),
         yaxp=c(0, 200, 8), cex.main=1, cex.lab=1, cex.axis=0.9, tck=-.02,
         xlab ="Time [min]",ylab="Cell count", main="Growth curves of embryos")
    mapply(lines,CellCountNorm[Embs$ID[Embs$Experiment=="ctrl"]], col=pal[2])
    mapply(lines,CellCountNorm[Embs$ID[Embs$Group=="alive"]], col=pal[3])
    mapply(lines,CellCountNorm[Embs$ID[Embs$Group=="dead"]], col=pal[4])
    mapply(lines,CellCountNorm[Embs$ID[Embs$Experiment=="wt"]], col=pal[1])
    
    # x <- sapply(CellCountNorm, function(x) c(length(x),max(x)))
    # points(x[1,],x[2,],col="red",pch=16)
    # text(x[1,],x[2,],colnames(x),cex=0.5)
    
    legend("topleft",c(sapply(leg[1:2], function(x) as.expression(substitute(A~degree~"C",list(A = as.name(x))))), leg[3:4]),
           lty=1, lwd=2, cex=0.8, col=pal)
    
dev.off()
#print(paste0("Aligned growt curves saved to: ", outDir, today, "_Growth_aligned_RAW.pdf"))

# debug

  # x="GM46"
  # lines(1:length(CellCountNorm[[x]]),CellCountNorm[[x]], col="red")


#aggregate curves to a matrix
maxL <- max(sapply(CellCountNorm, length)) #get length of the longest vector
mat <- sapply(CellCountNorm, '[', seq(maxL))

#plot mean curve for controls
pdf(paste0(outDir,today,"_Growth_aligned_MEAN.pdf"), width=7, height=6)
    plot(unlist(CellCountNorm),type="n",xlim=c(0,230), ylim=c(4,180), mgp=c(2, .7, 0),
         yaxp=c(0, 200, 8), cex.main=1, cex.lab=1, cex.axis=0.9, tck=-.02,
         xlab ="Time [min]",ylab="Number of embryonic cell", main="Growth curves of embryos")
    
    lines(apply(mat[,Embs$Experiment=="wt"],1, mean, na.rm=T), lwd=2, col=pal[1])
    lines(apply(mat[,Embs$Experiment=="ctrl"],1, mean, na.rm=T), lwd=2, col=pal[2])
    lines(apply(mat[,Embs$Group=="alive"],1, mean, na.rm=T), lwd=2, col=pal[3])
    lines(apply(mat[,Embs$Group=="dead"],1, mean, na.rm=T), lwd=2, col=pal[4])
    
    legend("topleft",c(sapply(leg[1:2], function(x) as.expression(substitute(A~degree~"C",list(A = as.name(x))))), leg[3:4]),
           lty=1, lwd=2, cex=0.8, col=pal)
dev.off()
print(paste0("Aligned growt curves saved to: ", outDir, today, "_Growth_aligned_MEAN.pdf"))
