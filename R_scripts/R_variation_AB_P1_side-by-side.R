library("tidyr")

# which variables to plot
vars <- c("pAP","pLR","pDV",
          "aAP","aLR","aDV",
          "LifeTime","EndTime","asynchrony",
          "netdis","totdis", "pOV", "aMean")

ylims <- list(c(0,25),c(0,25),c(0,25),
              c(0,120),c(0,120),c(0,120),
              c(0,20),c(0,20),c(0,20),
              c(0,100),c(0,35),c(0,8),c(0,40))

names(ylims) <- vars

toLong <- function(vc, cells, cc){
  td <-vc[cells,plotgroups] #test data
  #td <- td[complete.cases(td),]
  td <- td[apply(td[,1:4], 1, function(x) !any(is.na(x))),]
  td$Cell <- rownames(td)
  td$cellGroup <- as.factor(ifelse(grepl("AB",rownames(td)),"AB","P1"))
  
  td.l <- pivot_longer(td, cols=c(1:5),names_to ="Group", values_to="vc") #to long format
  td.l$Group <- factor(td.l$Group, levels = plotgroups)
  td.l <- as.data.frame(td.l)
  
  #for a given variable, compare groups and AB and P1 within each group
  # -> feed them into a list that will be queried in the plotting function
  #calculate comparisons using Tukey HSD test
  comparisons <- Tuk("vc", "Group", data=td.l) # perform multiple comparisons
  ycoord <- boxplot(vc~Group, data=td.l, plot=F)$stats[5,] #coordinate of the upper whisker
  
  compAB.P1 <- as.numeric(by(td.l,td.l$Group, function(x){t.test(x$vc~x$cellGroup)$p.value}))
  names(compAB.P1) <- plotgroups

  capture.output(
                cat("\r\t\t--------\t",var,"\t---------\n"),
                comparisons,
                cat("\n AB vs P1 t.test comparisons, p-values: \n"),
                compAB.P1,
                file = paste0(outDir,today,"_variance_all_comparisons.txt"), append=T)
  list(td=td.l, y=ycoord, comp=comparisons, AB.P1=compAB.P1)
}

pdf(paste0(outDir,today,"_coefficient_of_variantion_AB_P1.pdf"), width=3.2,height = 8, pointsize = 10, useDingbats = F)
  {
  par(mfrow=c(4,3), mar=c(3.5,2.5,2.1,0.5), mgp=c(1.5,0.55,0), cex.lab=0.9, cex.axis=0.8, cex.main=1, tcl=-.4)
  
  for (var in vars) {
    { #pre-process and calculate variation coeficient 
      res <- results[[var]] #get stats for given variable
      
      #remove too early cells, start at 8C
      res <- res[(3:nrow(res)),] #removes 4C stage
      
      colnms <- colnames(res)
      
      n <- res[,grep("\\.n$", colnms)]
      vc <- res[,grep("\\.vc$", colnms)]
      
      if(var %in% c("pAP","pLR","pDV", "pOV","aMean")) vc <- res[,grep("\\.sd$", colnms)]
      
      #keep only rows with >5 observation
      filter <- which(n$EA.n > 5)
      # inverted have less than 5 measurement for most cells in late time-points 
      # there is less cells for the late tp.
      
      n <- n[filter,]
      if(!var %in% c("pOV","aMean")) vc <- vc[filter,]*100 else vc <-vc[filter,]
      
      colnames(vc) <- c("wt","ctrl","equalized","alive","dead","inverted")
    }
    
    cells <- rownames(vc) %in% Cellorder[3:202]
    plotgroups <- c("wt","ctrl","alive","dead","inverted")

    #plot settings
    bw=0.50
    ymax = quantile(vc, 1, na.rm = T) # where to plot labels
    outl = F
    
    cc <- plotgroups
    td <- toLong(vc, cells, "all") #test data
    ypos <- ylims[[var]][2]
    
        boxplot(vc~cellGroup+Group,
                data=td$td,
                cex.axis=0.6,
                main=paste0(var),
                ylab=ifelse(var %in% c("pAP","pLR","pDV","pOV","aMean"),ifelse(var=="aMean","Standard deviation [degrees]","Standard deviation [% length]"), "Coefficient of variation %"), 
                outline= outl, boxwex=bw, 
                xlab=NA,
                xlim=c(-0.2,10),
                ylim=c(-ypos*0.1,ypos),
                outcex=0.4, lwd=0.6, lty=1,
                col = rep(pal, each=2), las=3,
                at=c(rep(c(-0.4,+0.4),5)+rep(seq(1,9,by=2),each=2)),
                staplewex=0.1,
                boxlwd=0.2,
                names=rep(c("AB","P1"),5)
                )
        # comparisons between AB and P1 in each group
        text(x=rep(seq(1,9,by=2)), -ypos*0.03, starformat(td$AB.P1), pos=1, cex=0.8, offset =0.1) #t.test
        for(i in 1:length(plotgroups))lines(x=c(0.6,0.6,1.4,1.4)+(i-1)*2,y=c(0,-ypos*0.015,-ypos*0.015,0), lwd=.5) 
        
        groupmeans=by(td$td$vc, td$td$Group, mean, na.rm=T)
        segments(x0=-0.5,x1=0, y0=groupmeans,y1=groupmeans, col=pal, lwd=0.7)      
        
        text(x=rep(seq(1,9,by=2)), td$y, td$comp$comp, pos=3, cex=0.8, offset = 0.5) #tukey
  }
  dev.off()
}

