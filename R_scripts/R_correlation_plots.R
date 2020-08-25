#plot correlation heatmap of signif variables ####
df$compression <- Embs[rownames(df),"height"]
vars <- c(signif,"Size","compression")
# no significant correlations with AB.size observed, thus omitted

pdf(paste0(outDir,today,"_Correlation_matrix_signif_vars_in_eqlzd_upto15C_heatmap.pdf"), height=2.5, width=3) 
{
  r <- cor.mtest(df[ids,vars], conf.level=.95)
  M <- cor(df[ids,vars], use="pairwise.complete.obs")
  
  cormat <- corrplot(M, p.mat=r$p,
                     method="color",
                     type="lower",
                     sig.level = c(.001, .01, .05),
                     insig="label_sig",
                     pch.cex=0.8,
                     pch.col="white",
                     order = "hclust",
                     tl.col = "black", 
                     tl.srt = 45, 
                     tl.cex = 0.6,
                     diag = T)
  
  cormat <- corrplot(M, p.mat=r$p, 
                     type="upper", 
                     method="color",
                     order = "hclust", 
                     hclust.method="centroid",
                     sig.level = .01, 
                     insig="blank", 
                     tl.col = "black", 
                     tl.srt = 45, 
                     tl.cex = 0.6)
  
  r <- cor.mtest(df[ids,signif], conf.level=.95)
  M <- cor(df[ids,signif], use="pairwise.complete.obs")
  
  cormat <- corrplot(M, p.mat=r$p, type="lower", 
                     order = "hclust", 
                     hclust.method="centroid",
                     sig.level = .05, insig="blank", 
                     method="color",
                     tl.col = "black", 
                     tl.srt = 45, 
                     tl.cex = 0.6,
                     diag=T)
  
  cormat <- corrplot(M, p.mat=r$p,
                     method="color",
                     type="lower",
                     sig.level = c(.001, .01, .05),
                     insig="label_sig",
                     pch.cex=0.8,
                     pch.col="white",
                     order = "hclust",
                     hclust.method="centroid",
                     tl.col = "black", 
                     tl.srt = 45, 
                     tl.cex = 0.6,
                     diag = T)
  dev.off()
}

#plot signif variables as boxplots ####
mat <- df[ids,signif]
outcome <- df[ids,"Outcome"]

pdf(paste0(outDir,today,"_boxplots_of_signif_vars_w.o.outl.pdf"), width=6, height = 1.8, pointsize = 11, useDingbats = F)
par(mfcol=c(1,8), mar=c(1.5,1.35,2.5,0.35), mgp=c(1.5,0.5,0), cex.lab=0.9, cex.axis=0.8, cex.main=1, tcl=-.35)
{
  for(i in signif){
    groupplot(i, "wilcox")
  }
  dev.off()
  par(mfrow=c(1,1))
}
