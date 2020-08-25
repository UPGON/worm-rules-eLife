
##hist(Embs$MaxTime)
temp <- Complete[Complete$Divided==1,]
temp <- temp[temp$embryo %in% Embs$ID[Embs$MaxTime>125],]

embs100 <- unique(temp$embryo)
cells100 <- unique(temp$Cell)
cells100 <- cells100[nchar(cells100)<10]

dt100 <- prep.data(TempDf, cells100, outlierRM=T)

(405/(nrow(dt100)*ncol(dt100)-7))*100 #percentage of replaced values - 0.9% 

pdf(paste0(outDir,today,"_PCA_upto100.pdf"), width=8, height=2.8, pointsize = 11)
{ par(mfrow=c(1,3))

    eqlzd <- dt100$embryo %in% Embs$ID[equalized]
    
    p <- embPCA(dt100[,],7:ncol(dt100),F, paste0("PCA at up to ~100-cell stage", ncol(dt100)-7, " variables"),"topright",T) #args: Df, cols, anot, title, legpos
    #add correlation with AB size
    plot.corr(filename="", save=F,
              data=data.frame(PC1=p$x[,1], AB.rel=dt100$AB.rel*100, Group=dt100$Group),
              xvar = "PC1",yvar = "AB.rel", groupvar="Group", 
              xlab="PC1",ylab = "AB size %",
              title="PC1 ~ AB size ")
    
    embPCA(dt100[eqlzd,],7:ncol(dt100),F,paste0("PCA at ~100-cell stage - equalized only"),"topright",F) #args: Df, cols, anot, title, legpos
    
    #with labels
    embPCA(dt100[,],7:ncol(dt100),T,paste0("PCA at ~100-cell stage - all"),"topright",F) #args: Df, cols, anot, title, legpos
    embPCA(dt100[eqlzd,],7:ncol(dt100),T,paste0("PCA at ~100-cell stage - equalized only"),"topright",F) #args: Df, cols, anot, title, legpos
  
  dev.off()
  par(mfrow=c(1,1))
}

#with EMS outliers in -- 
dt100 <- prep.data(TempDf, cells100, outlierRM=F)

(615/(nrow(dt100)*ncol(dt100)-7))*100 #percentage of replaced values - 1.14% 

pdf(paste0(outDir,today,"_PCA_upto100_with_EMSskewed.pdf"), width=8, height=2.8, pointsize = 11)
{ par(mfrow=c(1,3))
  
  eqlzd <- dt100$embryo %in% Embs$ID[equalized]
  
  p <- embPCA(dt100[,],7:ncol(dt100),F, paste0("PCA at up to ~100-cell stage", ncol(dt100)-7, " variables"),"topleft",T) #args: Df, cols, anot, title, legpos
  #add correlation with AB size
  plot.corr(filename="", save=F,
            data=data.frame(PC1=p$x[,1], AB.rel=dt100$AB.rel*100, Group=dt100$Group),
            xvar = "PC1",yvar = "AB.rel", groupvar="Group", 
            xlab="PC1",ylab = "AB size %",
            title="PC1 ~ AB size ")
  
  embPCA(dt100[eqlzd,],7:ncol(dt100),F,paste0("PCA at ~100-cell stage - equalized only"),"topleft",F) #args: Df, cols, anot, title, legpos
  
  embPCA(dt100[,],7:ncol(dt100),T,paste0("PCA at ~100-cell stage - all"),"topleft",F) #args: Df, cols, anot, title, legpos
  embPCA(dt100[eqlzd,],7:ncol(dt100),T,paste0("PCA at ~100-cell stage - equalized only"),"topleft",F) #args: Df, cols, anot, title, legpos
  
    
  dev.off()
  par(mfrow=c(1,1))
}

#with EMS outliers in at 56C ####
dt50 <- prep.data(TempDf, fiftysixC, outlierRM=T)

(170/(nrow(dt50)*ncol(dt50)-7))*100 #percentage of replaced values - 1.5% 

pdf(paste0(outDir,today,"_PCA_56C_witout_EMSskewed.pdf"), width=8, height=2.8, pointsize = 11)
{ par(mfrow=c(1,3))
  
  eqlzd <- dt50$embryo %in% Embs$ID[equalized]
  
  p <- embPCA(dt50[,],7:ncol(dt50),F, paste0("PCA at 56-cell stage", ncol(dt50)-7, " variables"),"topleft",T) #args: Df, cols, anot, title, legpos
  #add correlation with AB size
  plot.corr(filename="", save=F,
            data=data.frame(PC1=p$x[,1], AB.rel=dt50$AB.rel*100, Group=dt50$Group),
            xvar = "PC1",yvar = "AB.rel", groupvar="Group", 
            xlab="PC1",ylab = "AB size %",
            title="PC1 ~ AB size ")
  
  embPCA(dt50[eqlzd,],7:ncol(dt50),F,paste0("PCA at 56-cell stage - equalized only"),"topleft",F) #args: Df, cols, anot, title, legpos
  
  embPCA(dt50[,],7:ncol(dt50),T,paste0("PCA at ~100-cell stage - all"),"topleft",F) #args: Df, cols, anot, title, legpos
  embPCA(dt50[eqlzd,],7:ncol(dt50),T,paste0("PCA at ~100-cell stage - equalized only"),"topleft",F) #args: Df, cols, anot, title, legpos
  
  dev.off()
  par(mfrow=c(1,1))
}
  
