rows <- Complete$Group!="wt" & Complete$Divided==1 & (Complete$Cell %in% Cellorder[c(5:54,87:92)])
df <- Complete[rows,]
df$Cell <- factor(df$Cell)

timeCor <- with(df,
  by(cbind(LifeTime,AB.rel), Cell, function(x){
    s = complete.cases(x);
    if(sum(s)>5) cor.test(x[s,1],x[s,2])[c("estimate","conf.int","p.value")]
  })
)
  #unlist data
  ci <- as.data.frame(t(sapply(timeCor, "[[", 2)))
  colnames(ci) <- c("lower","upper")
  
  r <- sapply(timeCor, "[[", 1)
  pv <- sapply(timeCor, "[[", 3)
  
  c.df <- cbind(ci,r,pv)
  c.df <- c.df[Cellorder[c(5:54,87:92)],]
  c.df$p.adj <- p.adjust(c.df$pv,method="BH")
  c.df$order <- 1:nrow(ci)
  c.df$cell <- factor(rownames(c.df),levels =rownames(c.df))
  

pdf(paste0(outDir,today,"_corr_AB-size_cell-cycle.pdf"),pointsize=10,width=6.5,height=2.2, useDingbats = F)
{
  par(mfcol=c(1,1), mar=c(2.5,2.5,2.1,1.1), mgp=c(1.5,0.6,0),cex.lab=0.9, cex.axis=0.7, cex.main=1)
  
  plot(1,-1.5, 
       main="Correlation of Asymmetry vs Cell cycle duration",
       ylab="Correlation AB size % ~ Cell cycle [min]",
       xlab = "",
       xaxt="n",
       ylim=c(-0.6,1),
       xlim=c(1,nrow(c.df)),
       cex.lab=0.8,
       cex.main=0.9,
       cex.axis=0.7
  )
  
  with(c.df,{
    #c <- colors[cell.groups(cell)]
    c <- ifelse(grepl("AB",c.df$cell),"#2B8E82","#BA2424")
    
    by(c.df,cell, function(x) lines(x=rep(x[,6],2), y=c(x[,1],x[,2]),lwd=0.7, col="grey"))
    points(1:length(r),r,
           cex=0.8, 
           col=c,
           pch=ifelse(p.adj<0.05,16,1),
    )
    abline(h=0,lty=2, lwd=0.5)
    axis(1, at=1:length(r),labels =cell, cex.axis=.6, las=2)
  })
  
dev.off()
}
write.xlsx(c.df, paste0(outDir, today, "_corr_AB-size_cell-cycle.xlsx"))

rm(df, ci, r, pv, timeCor, rows)

