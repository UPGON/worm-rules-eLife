#prepare data ####
q <- paste("^(",paste0(cells, collapse="|"),")\\..*",sep="")

cols <- grep(q,colnames(TempDf))
#remove EMS outliers
rows <- TempDf$embryo %in% Embs$ID[equalized & !(Embs$ID%in%outliers)]

df <- FilterAndSubset(TempDf, rows, cols, minN=5)

#plot volcano ####
ids <- 1:nrow(df)
cols <- colnames(df)[7:ncol(df)]

  mat <- df[ids,cols]
  outcome <- df[ids,"Outcome"]
  
  pdf(paste0(outDir,today,"_volcano_upTo15C_dead.over.hatched.pdf"), width=3, height = 3.4, pointsize = 11, useDingbats = F)
    
    par(mar=c(2.5,2.5,2,0.15), mgp=c(1.25,0.5,0), cex.lab=0.7, cex.axis=0.6, cex.main=0.8, tcl=-.35)
    resDf1 <- volcano(mat, outcome, "Volcano Dead/Hatched up to 15C")
    
  dev.off()
    
  resDf1 <- resDf1[order(resDf1$p.value),]
  write.csv(resDf1, file=paste0(outDir,today,"_volcano_upTo15C_dead.hatched.csv"))
  signif <- rownames(resDf1)[resDf1$signif==T]
  #"ABara.pDV" "ABara.pOV" "Ca.netdis" "MS.aMean"  "MS.aAP"    "ABpra.aAP" "ABarp.aAP" "ABprp.aAP"
  
  print(paste0("Volcano plot saved as: ", outDir,today,"_volcano_upTo15C_dead.over.hatched.pdf"))
  cat("There is", length(signif),"significant variables: ", signif)