# plot coef. of variance for each variable and group of embryos at different stages
# asumption is that dying embryos accumulate more variance in all directions
library("tidyr")

# which variables to plot
vars <- c("pAP","pLR","pDV",
          "aAP","aLR","aDV",
          "LifeTime","EndTime","pOV",
          "aMean","netdis","totdis")

ylims <- list(c(0,75),c(0,75),c(0,75),
             c(0,120),c(0,120),c(0,120),
             c(0,15),c(0,15),c(0,8),
             c(0,45),c(0,100),c(0,35))

#use with SD
ylims <- list(c(0,22),c(0,22),c(0,22),
              c(0,45),c(0,45),c(0,45),
              c(0,8),c(0,15),c(0,8),
              c(0,45),c(0,4.5),c(0,7))


multiply <- c(T,T,T,
              F,F,F,
              F,F,F,
              F,F,F
              )

ylabb <- c(rep("SD length [%]",3),
           rep("SD angle [degrees]",3),
           rep("SD time [min]",2),
           "SD mispositioning [um]",
           "SD angular deviation [degrees]",
           rep("SD migration [um]",2)
)

names(ylabb) <- names(ylims) <- names(multiply) <- vars


{
# plot stages side by side ####
# plot SD for pDV and positions

timelim=200 #take only cells dividing up to 200 minutes

pdf(paste0(outDir,today,"_variantion_in_embryos_outlineF.pdf"), width=5,height = 8, pointsize = 10, useDingbats = F)
{
  par(mfrow=c(4,3), mar=c(3.5,2.5,2.1,0.5), mgp=c(1.5,0.55,0), cex.lab=0.9, cex.axis=0.8, cex.main=1, tcl=-.4)
  
  for (var in vars) {
    { #pre-process and calculate variation coeficient 
      res <- results[[var]] #get stats for given variable
      timeres <- results[["EndTime"]]
      
      cc <- apply(timeres[,2:7],1,function(x)any(x<timelim, na.rm = T))
      
      #remove too early cells, start at 8C
      res <- res[cc,]
      
      colnms <- colnames(res)
      
      n <- res[,grep("\\.n$", colnms)]
      vc <- res[,grep("\\.sd$", colnms)]
      
      #if(var %in% c("pOV","aMean","pAP","pDV","pLR")) vc <- res[,grep("\\.sd$", colnms)]

      
      #keep only rows with >5 observation
      filter <- which(n$EA.n > 5)
        # inverted have less than 5 measurement for most cells in late time-points 
        #> there is less cells for the late tp.
      
      n <- n[filter,]
      if(multiply[var]) vc <- vc[filter,]*100
    }
    
    cells <- rownames(vc)
    plotgroups <- c("wt","ctrl","alive","dead","inverted")
    colnames(vc) <- c("wt","ctrl","equalized","alive","dead","inverted")
    
    #plot settings
    bw=0.60
    ymax = quantile(vc, 1, na.rm = T)
    ymax = ylims[[var]][2]
    
    outl = F
    
    stages <- list(eightC,fifteenC,twentyeightC, fiftysixC, cells)
    stage_labs <- c("8C","16C", "28C","56C","All")
    
    for (i in 1:5){
      #for a given variable, which combinations are significant? -> feed them into a list that will be queried in the plotting function
      
      td <-vc[stages[[i]],plotgroups] #test data
      td.l <- gather(td, key="Group",value="vc") #to long format
      td.l$Group <- factor(td.l$Group, levels = plotgroups)
      
      #calculate comparisons using Tukey HSD test
      comparisons <- Tuk("vc", "Group", data=td.l) # perform multiple comparisons
      ycoord <- boxplot(vc~Group, data=td.l, plot=F)$stats[5,] #coordinate of the upper whisker
        
      capture.output(cat("---------- ",var, "------------\n\n"),
                     comparisons,
                     cat("\n\n\n"),
                     file = paste0(outDir,today,"_variance_comparisons.txt"),
                     append=T)
        
      if(i==1){
        boxplot(td, cex.axis=0.8,  
                main=paste0(var),
                ylab=ylabb[var], outline= outl, boxwex=bw, 
                xlim=c(0.85,32), 
                ylim=ylims[[var]], 
                outcex=0.4, lwd=0.6, lty=1,
                col = pal, las=3,
                staplewex=0.1,
                boxlwd=0.2)  
        lines(x=c(1,5),y=rep(-.0,2), lwd=1)
        text(x=3,y=ymax*0.035,"8C",cex=0.8)
        
        text(x=1:5, ycoord, comparisons$comp, pos=3, cex=0.8, offset = 0.5)
        
      }
      else {
        xoffset <- (6.5*(i-1))
        boxplot(td, cex.axis=0.8,  
                outline = outl, add=T, 
                boxwex=bw, yaxt="n", 
                at=(1:5)+xoffset, outcex=0.4, lwd=0.6, lty=1,
                col=pal, las=2,
                staplewex=0.1,
                boxlwd=0.2
                )
        lines(x=c(1,5)+xoffset,y=rep(-.0,2), lwd=1)
        text(x=3+xoffset,y=ymax*0.035,stage_labs[i],cex=0.8)
        text(x=(1:5)+xoffset, ycoord, comparisons$comp, pos=3, cex=0.8, offset = 0.5)
      }
    }
  }
dev.off()
}

}
print("Completed and saved")
