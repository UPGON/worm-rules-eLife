#use all variables up to 15C ####
eqlzd <- pcaDf$embryo %in% Embs$ID[equalized]
cols <- colnames(pcaDf)[7:ncol(pcaDf)]

#PCA with all vars that go into analysis ####

pdf(paste0(outDir,today,"_upto15Cells_PCA_all-vars.pdf"), width=7, height=3.8, pointsize = 10)
{ par(mfrow=c(1,2))
  
  p <- embPCA(pcaDf[,],cols,F, paste("PCA up to 15C - ", length(cols), " variables"),"topleft",T) #args: Df, cols, anot, title, legpos
  #add correlation with AB size
  plot.corr(filename="", save=F,
            data=data.frame(PC1=p$x[,1], AB.rel=pcaDf$AB.rel, Group=pcaDf$Group),
            xvar = "PC1",yvar = "AB.rel", groupvar="Group", 
            xlab="PC1",ylab = "AB size %",
            title="PC1 ~ AB size ")
  
  dev.off()
  par(mfrow=c(1,1))
}

#PCA of equalized embryos ####
  # ids <- pcaDf$embryo[eqlzd]
  # embPCA(pcaDf[ids,],cols,T,"Variable up to 15C",NA,F) #args: Df, cols, anot, title, legpos

#use signif vars for PCA ####  
pdf(paste0(outDir,today,"_PCA_with_signif_vars_upto15C.pdf"), width=7, height = 2.6, pointsize = 11, useDingbats = F)
{
  par(mfrow=c(1,3))
  
  pca <- embPCA(pcaDf[,],signif,F,"Significantly diff. variables H/D up to 15C","bottomleft",T) #args: Df, cols, anot, title, legpos
  
  biplot (pca,cex=0.6,cex.axis=0.7,cex.lab=0.7, pc.biplot=T, 
          expand=0.9, arrow.len=0.07, xlim=c(-2,2))  
  
  plot(pca$x[,1],Embs[row.names(pca$x),"AB_rel"]*100, col=pal[Embs[row.names(pca$x),"Group"]], pch=16, main="PC1 ~ AB relative size", ylab="AB size %", xlab="PC1")      
  
  #plot only equalized
  pca <- embPCA(pcaDf[ids,],signif,F,"Significantly diff. variables H/D up to 15C","topleft",T) #args: Df, cols, anot, title, legpos
  biplot (pca, scale=1 ,cex=0.6,cex.axis=0.7,cex.lab=0.7, pc.biplot=T, 
          expand=0.9, arrow.len=0.06, xlim=c(-2,2))  
  
  plot(pca$x[,1],Embs[row.names(pca$x),"AB_rel"]*100, col=pal[Embs[row.names(pca$x),"Group"]], pch=16, main="PC1 ~ AB relative size", ylab="AB size %", xlab="PC1")      
  
  dev.off()
  par(mfrow=c(1,1))
}
