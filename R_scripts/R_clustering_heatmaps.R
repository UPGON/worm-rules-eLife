
clusterDf <- lassoDf
levels(clusterDf$Group) <- c("Wildtype","Ctrl","Alive","Dead","Inverted")
# data is scaled lassoDf 

# row annotations
  prediction <- as.numeric(predict(mod1, newx=as.matrix(data[,row.names(mod1$beta)]), type="class"))
  outcome.all <- as.numeric(clusterDf$Outcome)-1
  sum(prediction==outcome.all)/length(prediction)
  #89.8% success on all data including training and test dataset + outlier embryos
  
  Conf(prediction, outcome.all) #Kappa 79.15
  #             Reference
  # Prediction  0  1
  #           0 21 6
  #           1  0 32

  predC <- ifelse(prediction==0,"dead","alive")
  
  Pha4 <- rep("ND",nrow(clusterDf))
  names(Pha4) <- clusterDf$embryo
  Pha4[grep("^P.*", clusterDf$embryo)] <- "normal"
  Pha4[clusterDf$embryo %in% phaDefects] <- "defective"
  
  Extra.P4 <- as.character(Embs[clusterDf$embryo,"Extra.P4"])
  Extra.P4[is.na(Extra.P4)] <- "ND"
  
  r.annot <- cbind(Extra.P4=Extra.P4,
                   Pha4=Pha4,
                   Prediction=predC,
                   clusterDf[,c("Outcome","Group")]
                   )
  #,Partial="#6BB50B"
  a.cols <- list(Prediction=c(dead="black",alive="chartreuse3"),
                 Outcome=c(hatched="#C9F299",died="black"),
                 Extra.P4=c(F="#C9F299",T="#9F1228",ND="grey"),
                 Pha4=c(normal="#C9F299",defective="#AA00AF",ND="grey"),
                 Group=c(Wildtype=pal[1],Ctrl=pal[2],Alive="#E2A807",Dead=pal[4], Inverted="#555555"))

#heatmaps for E-E similarity ####
ids <- clusterDf$embryo
ids <- ids[!ids %in% c("PM23", "GM50")]

mat <- as.matrix(clusterDf[ids,7:ncol(clusterDf)])

#all embryos and all variables
met="ward.D"
{
  sampledist <- dist(scale(mat[ids,]), method = "eucl")

  ph <- pheatmap(sampledist,
                 main="All variables up to 15C - Ward.D",
                 scale = "none",
                 clustering_method = met,
                 annotation_row = r.annot,
                 annotation_colors = a.cols,
                 breaks = seq(floor(min(sampledist)),max(sampledist),length.out = 100),
                 #color=colramp(),
                 border_color="white",
                 fontsize_row=3,
                 fontsize_col=3,
                 fontsize=6,
                 lwd=0.5
  )
}
save_pheatmap_pdf(ph, paste0(outDir,today,"_dist_all_vars.pdf"), 6.5,5.5)

#plot only equalized embryos ####
ids <- clusterDf$embryo %in% Embs$ID[equalized]

sampledist <- as.matrix(dist(scale(as.matrix(clusterDf[ids,signif]))))
ph2 <- pheatmap(as.matrix(sampledist),
                main="Equalized Embs - Sample Distance - Ward.D",
                annotation_row = r.annot,
                annotation_colors = a.cols,
                clustering_method = met,
                border_color="#DCDCDC",
                breaks = seq(floor(min(sampledist)),max(sampledist),length.out = 100),
                fontsize_row=6,
                fontsize=6,
                lwd=0.5
)

save_pheatmap_pdf(ph2, paste0(outDir,today,"_dist_matrix_signif-vars_upTo15C_equalizedE_w.o.outl.pdf"))

ph3 <- pheatmap(clusterDf[clusterDf$Experiment=="meta",signif],
                main="Clustering based on significantly different variables in the volcano plot equalized embryos",
                annotation_row = r.annot,
                annotation_colors = a.cols,
                clustering_method = met,
                border_color="#DCDCDC",
                fontsize_row=6,
                fontsize=8,
                lwd=0.5,
                scale="column"
)

save_pheatmap_pdf(ph3, paste0(outDir,today,"_clustering_signif_vars_from_volcano_upTo15C_all_upshifted.pdf"))
# 
ph4 <- pheatmap(clusterDf[,signif],
                main="Clustering with signif EA/ED vars in equalized embryos",
                annotation_row = r.annot,
                annotation_colors = a.cols,
                clustering_method = met,
                border_color="white",
                fontsize_row=5,
                fontsize=6.5,
                lwd=0.25,
                treeheight_col = 15,
                treeheight_row = 30,
                scale="column"
)

save_pheatmap_pdf(ph4, paste0(outDir,today,"_clustering_signif_vars_from_volcano_upTo15C_all_embs.pdf"), 3.4,5.5)
