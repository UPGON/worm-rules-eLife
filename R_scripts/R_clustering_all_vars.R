#all embryos and all variables
met="ward.D2"

Dt <- prep.data(TempDf, cells=c(eightC,fifteenC,twentyeightC,fifsixC),outlierRM=T)

{ ph <- pheatmap(as.matrix(Dt[7:ncol(Dt)]),
                 scale = "column",
                 main=paste("all variables up to 28C", met),
                 clustering_method = met,
                 annotation_row = r.annot,
                 annotation_colors = a.cols,
                 color=colramp,
                 border_color="white",
                 fontsize_row=4,
                 fontsize_col=4,
                 fontsize=5,
                 treeheight_row=30,
                 lwd=0.5,
                 width=14,
                 height = 6,
                 filename=paste0(outDir,today,"_cluster_all_vars.pdf")
  )
}

ids <- !Dt$embryo %in% c("GM50","EM10", "EM4")
{
  sampledist <- as.matrix(dist(scale(Dt[ids,7:ncol(Dt)]), method = "eucl"))
  sampledist[sampledist==0] <- min(sampledist[sampledist!=0])
  
  ph <- pheatmap(sampledist,
                 scale = "none",
                 main=paste("All variables, excluded GM50, EM10, EM4",met),
                 clustering_method = met,
                 annotation_row = r.annot,
                 annotation_colors = a.cols,
                 color=colramp,
                 border_color="white",
                 fontsize_row=5,
                 fontsize_col=5,
                 fontsize=6,
                 treeheight_row=30,
                 lwd=0.5,
                 width=6,
                 height=5,
                 filename=paste0(outDir,today,"_cluster_all_vars_sample_dist.pdf")
  )
}

