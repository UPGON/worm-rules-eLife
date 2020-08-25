#MS inversion ####
# plot variables related to MS.angle
rows <- TempDf$embryo %in% Embs$ID[!(Embs$ID%in%outliers)] #remove EMS outliers
vars <- c("MS.aMean","MS.aAP","MSa.pOV","MSp.pOV","ABarp.aAP", "ABprp.aAP", "ABpra.aAP")
lab <- c("MS angular deviation from control [°]", 
         "MS AP division angle deviation [°]",
         rep("Positional deviation from control [um]",2),
         "ABarp AP division angle deviation [°]")

pdf(paste0(outDir,today,"_MS_divison_angle.pdf"),pointsize=10,width=4,height=3.5)
{ 
  par(mfcol=c(1,2), mar=c(2.5,2.5,2.1,1.1), mgp=c(1.5,0.5,0), cex.lab=0.9, cex.axis=0.8, cex.main=1, tcl=-.35)
  
  for (i in 1:length(vars)) {
    v <- vars[i]
    pl <- TempDf[rows,c(v,"Group")]
  #frequency in different groups
  boxplot(pl[,v]~pl$Group,
          main=v,
          ylab=lab[i],
          outline=F,
          boxwex=0.75,
          lwd=0.5
  )
  
  beeswarm(pl[,v]~pl$Group,
           cex=0.6, pch=16, 
           pwcol=pal[pl$Group],
           add=T)
  
  #abline(h=35, lty=2)
  cat(file="MS_angle_stats.txt",
    v, "\n", 
    Tuk(v, "Group", pl)$model$Group[], "\n",
    Tuk(v, "Group", pl)$comp,
    append = T
  )
  }
  
  dev.off()
}
