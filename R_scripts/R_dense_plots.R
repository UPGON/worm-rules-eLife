# Author: Radek Jankele 08/07/19

###############################################################################
#source("R_scripts/R_All_stats_export_XLSX.R") this should run before to get access to "results" list with all stats

####  PLOTTING FUNCTIONS #####
tufte.boxplot <- function(x, g, at=0, pal, add=F, xlim=NA, ylim=NA, cex=0.65, lw=0.75, ylab=NA) {
  #x - vector with values for plotting
  #g - grouping factor of the same length as x
  #adapted from https://gist.github.com/even4void/1128764
  k <- nlevels(g)
  crit.val <- tapply(x, g, median, na.rm=T)
  spacing <- seq(1.5,k-0.5, length.out = k)
  
  if(at!=0) at <- at
  if(!any(is.na(xlim))) xlim=xlim else xlim=c(0,k)
  if(!any(is.na(ylim))) ylim=ylim else ylim=ylim=c(min(x)*1.1, max(x)*1.1)  
  
  
  if(!add){ #for drawing initial plot and first grouped variable
    plot(spacing+at, crit.val,
         main="",
         col=pal,
         cex=cex,
         cex.axis=0.85,
         cex.lab=0.85,
         tcl=-0.5,
         xlim=xlim,
         ylim=ylim,
         pch=20,
         xlab=NA,
         ylab=ylab, 
         xaxt="n"
    )
  }
  else points(spacing+at, crit.val, pch=20, col=pal, cex=cex) # for fullowing points
  
  for (i in 1:k) {
    tmp <- boxplot.stats(x[as.numeric(g)==i])
    
    if(at!=0) xi <- (spacing[i]+at) else xi <- spacing[i]
    
    segments(xi, tmp$stats[1], xi, tmp$stats[2], col=pal[i], lwd=lw)
    segments(xi, tmp$stats[2], xi, tmp$stats[4], col=pal[i], lwd=lw*0.5, lty=3)
    segments(xi, tmp$stats[4], xi, tmp$stats[5],col=pal[i], lwd=lw)
    points(rep(xi, length(tmp$out)), tmp$out,col=pal[i], cex=0.8*cex)
  }
}
#horizontal tufte boxplot function
tufte.boxplot.h <- function(x, group, at=0, pal, add=F, xlim=NA, ylim=NA, cex=0.65, lw=0.75, ylab=NA, logscale="") {
  #adapted from https://gist.github.com/even4void/1128764
  k <- nlevels(group)
  crit.val <- tapply(x, group, median, na.rm=T)
  m <- tapply(x, group, mean, na.rm=T)
  
  spacing <- seq(1.5,k-0.5, length.out = k)
  
  if(at!=0) at <- at
  if(!any(is.na(xlim))) xlim=xlim else xlim=c(min(x)*1.1, max(x)*1.1)  
  if(!any(is.na(ylim))) ylim=ylim else ylim=c(k,0)
  
  
  if(!add){ #for drawing initial plot and first grouped variable
    plot(crit.val,spacing+at,
         main="",
         log=logscale,
         col=pal,
         cex=cex,
         cex.axis=0.85,
         cex.lab=0.85,
         tcl=-0.3,
         xlim=xlim,
         ylim=ylim,
         pch=20,
         xlab=ylab,
         ylab=NA, 
         yaxt="n"
    )
    # add meand as | plotting symbol
    points(m, spacing+at, pch="|", col=pal, cex=cex*0.8)
  }
  else {# for adding all following points
    # plot median and mean
    points(crit.val, spacing+at, pch=20, col=pal, cex=cex)
    points(m, spacing+at, pch="|", col=pal, cex=cex*0.8)
  }
  
  for (i in 1:k) {
    tmp <- boxplot.stats(x[as.numeric(group)==i])
    
    if(at!=0) xi <- (spacing[i]+at) else xi <- spacing[i]
    
    segments(tmp$stats[1], xi,  tmp$stats[2], xi, col=pal[i], lwd=lw) # outer bottom bar
    segments(tmp$stats[2], xi, tmp$stats[4], xi,  col=pal[i], lwd=lw*0.8, lty=3) # middle part corresponding to classical boxplot
    segments(tmp$stats[4], xi, tmp$stats[5], xi,  col=pal[i], lwd=lw) # outer top bar 
    points(tmp$out, rep(xi, length(tmp$out)),col=pal[i], cex=0.8*cex, lwd=lw*0.8) #p;lot outliers
  }
}
#plot only selected cells
tufteplot.celllist <- function(data, cells, var){
  
  g = length(levels(data$Group))
  #subset the data and order cells by predefined order
  tempdf <- data[,grep(var, colnames(data[1,]))]
  colnames(tempdf) <- sub(paste0("\\.",var),"",colnames(tempdf[1,]))
  tempdf <- tempdf[,match(cells,colnames(tempdf[1,]), nomatch=0)]
  
  #find min and max of the tempdf (for cell range queried) REPLACE
  minY <- min(tempdf, na.rm = T)
  maxY <- max(tempdf*1.05, na.rm = T)
  
  xlim=c(minY*0.8,maxY)
  
  for (a in 1:ncol(tempdf)){
    cell <- colnames(tempdf[a])
    x <- tempdf[,cell] #subset the data
    
    if(a==1){
      tufte.boxplot.h(x, data$Group, xlim=xlim, ylim=c(ncol(tempdf)*g,1), 
                      at=((a-1)*g),  pal=pal, cex=0.4, 
                      ylab=var.units, lw=lw, logscale=logscale) #lw - line width
      firstpass <- 1
    }
    else {
      tufte.boxplot.h(x, data$Group, at=((a-1)*g), pal=pal, cex=0.4, lw=lw, add=T)
    }
  }
  abline(h=seq(0,g*ncol(tempdf),by=g)+0.5, lty=1, lwd=0.5, col="grey")
  title(main=var.labels, cex=0.8)
  axis(2,at=seq(2.5,g*ncol(tempdf),by=g),labels=colnames(tempdf),cex.axis=0.8, adj=0, srt=90, tick=T, tcl=-0.3, las=2, mgp=c(0,0.5,0))
  
  #plot significance stars
  plot.signif.stars(var,cells,maxY)
}
#produce tufte plot and save it as pdf
densePdf <- function(data ,filename, plotvars, leg, var.labels, var.units){
  
  pdf(file=filename, paper = "a4r", width = 11, height = 8, pointsize=12, useDingbats = F)
  par(mfrow=c(1,length(plotvars)),
      mar=c(3, 0, 1.5, 0.2), # bottom, left, top, right
      oma=c(0, 4, 0, 4), # bottom, left, top, right
      mgp = c(1.75, 0.5, 0)
  ) 
  
  for (i in 1:length(plotvars)){
    v <- plotvars[i]
    g <- nlevels(data$Group)
    
    #subset the data and order cells by predefined order
    tempdf <- data[,grep(v, colnames(data[1,]))]
    colnames(tempdf) <- sub(paste0("\\.",v),"",colnames(tempdf[1,]))
    tempdf <- tempdf[,match(Cellorder,colnames(tempdf[1,]), nomatch=0)]
    
    #find min and max of the tempdf (for cell range queried)
    if (ncol(tempdf)>limC){
      minY <- min(tempdf[,1:limC], na.rm = T)
      maxY <- max(tempdf[,1:limC]*1.05, na.rm = T)
      xlim=c(minY*0.8,maxY)
    } else
    {
      minY <- min(tempdf, na.rm = T)
      maxY <- max(tempdf*1.05, na.rm = T)
      xlim=c(minY*0.8,maxY)
    }
    
    #if cell is not in the data (such as first 4 cells) generate empty plotted column
    firstpass <- 0
    for (a in 1:limC) {
      cell <- Cellorder[a]
      match.col <- which(colnames(tempdf) %in% cell) # if the cell is not in the table this is 0
      
      if(length(match.col)==1){# column with given name exists
        if(!firstpass){# plotting for the first time, new plot needs to be initiates
          x <- tempdf[,cell] #subset the data
          tufte.boxplot.h(x, data$Group, xlim=xlim, ylim=c(limC*g,1), 
                          at=((a-1)*g),  pal=pal, cex=0.5, 
                          ylab=var.units[i], lw=lw, logscale=logscale) #lw - line width
          firstpass <- 1
        }
        else {
          x <- tempdf[,cell]
          tufte.boxplot.h(x, data$Group, 
                          at=((a-1)*g), pal=pal, cex=0.5, lw=lw, add=T)
        }
      }
      
    }
    abline(h=seq(0,g*limC,by=g)+0.5, lty=1, lwd=0.5, col="grey")
    title(main=var.labels[i], cex=0.8)
    
    if(v==plotvars[1]){
      axis(2,at=seq(2.5,g*limC,by=g),labels=Cellorder[1:limC],cex.axis=0.9, adj=0, srt=90, tick=T, tcl=-0.3, las=2, mgp=c(0,0.5,0))
      legend("topleft",leg,cex=0.8, pch=20,col=pal)
    }
  
    plot.signif.stars(var=v,cells=Cellorder[1:limC],maxY)
    #plot significance stars
  }
  
  #add an additional Cell names to the right side of the plot
  axis(4,at=seq(2.5,g*limC,by=g),labels=Cellorder[1:limC],cex.axis=0.9, adj=0, srt=90, tick=T, tcl=-0.3, las=2, mgp=c(0,0.5,0))
  dev.off()
}

plot.signif.stars <- function(var,cells, at){
  #plot significance stars
  pv <- data.frame(results[[var]])
  pv <- pv[rownames(pv)%in%cells,grep("p.val",names(pv[1,]))]
  
  p.adj <- p.adjust(pv$EA.vs.ED.p.val, method='BH')
  pv$signif.stars <- stars.pval(p.adj)
  
  #compare equalized with unequal ev571
  p.adj <- p.adjust(pv$C.vs.EQ.p.val, method='BH')
  pv$signif.stars.ctrl.eq <- stars.pval(p.adj)
  
  #compare wt with unequal ev571
  p.adj <- p.adjust(pv$C.vs.WT.p.val, method='BH')
  pv$signif.stars.ctrl.wt <- stars.pval(p.adj)
  
  #compare inverted with eq. alive ev571
  p.adj <- p.adjust(pv$EA.vs.INV.p.val, method='BH')
  pv$signif.stars.inv.ea <- stars.pval(p.adj)
  
  xi <- match(rownames(pv),cells)-1
  text(x=at, y=xi*g+1.3, pv[,"signif.stars.ctrl.wt"], cex=0.85, col=pal[1])
  text(x=at, y=xi*g+2.4, pv[,"signif.stars.ctrl.eq"], cex=0.85, col=pal[3])
  text(x=at, y=xi*g+3.5, pv[,"signif.stars"], cex=0.85, col=pal[4])
  text(x=at, y=xi*g+4.6, pv[,"signif.stars.inv.ea"], cex=0.85, col=pal[5])
}

#### Output #######
## DataFrame WideDf - contains aligned timing plus all big table data in an expanded format 

######### Generate dense plots  ############
#REMOVE OUTLIERS THAT HAVE EMS DIV. SKEWED 
data <- WideDf[!WideDf$embryo%in%outliers,c(4,7:ncol(WideDf))] # keep just group and all variables
g = length(levels(data$Group))

#to compress display of P4 lifetime substract some constant
data$P4.LifeTime <- data$P4.LifeTime-45

#HORIZONTAL tufte plot ####
#plot dense plot for each cell in an HORIZONTAL tufte plot ####

  #set up variables
    lw <- 0.75 # define line-width #define line-width to be used for the tufte plot
    StartC <- Cellorder[3]
    limC=32
    LastC <- Cellorder[limC]
    
    plotvars <- c("EndTime","LifeTime","netdis","pOV","aMean")
    leg <- c("Wild-type","Unequal controls", "Equalized alive", "Equalized dead")
    var.labels <- c("Time of division","Cell cycle duration", "Displacement", "Positional dev.","Angular dev.")
    var.units <- c("min","min","um", "um","pi")
    
    logscale="" #change to plot in the log scale
    filename <- paste0(outDir,today,"_densePlot_Timing.netdis.pOV.aMean.pdf")
    
densePdf(data ,filename, plotvars, leg, var.labels, var.units)

#plot timing only ####
    plotvars <- c("EndTime","LifeTime","asynchrony")
    var.labels <- c("Time of division","Cell cycle duration","Asynchrony")
    var.units <- c("min","min","sister cell cycle ratio")
    filename <- paste0(outDir,today,"_densePlot_timing.pdf")
    
densePdf(data ,filename, plotvars, leg, var.labels, var.units)

#plot asymmetries####
plotvars <- c("asynchrony","V.rel","sisterRatio","netdis")
var.labels <- c("Sister Asynchrony","Relative volume","Sister Vol. ratio", "Displacement")
var.units <- c("time ratio","cell V/embryo V","volume ratio","um")
filename=paste0(outDir,today,"_densePlot_asymmetries.pdf")

densePdf(data ,filename, plotvars, leg, var.labels, var.units)
aggregate(data$MS.aAP, by=list(data$Group), mean, na.rm=T)
aggregate(data$MS.aAP, by=list(data$Group), function(x) sum(!is.na(x))   )
compare_means(MS.aAP~Group,data = data)

#plot positional vars ####

    plotvars <- c("pOV","pAP","pDV","pLR","netdis")
    var.labels <- c("Deviation (pOV)","AP position","DV position","LR position","Displacement")
    var.units <- c("um","position AP%","position DV%","position LR%","um")
    filename=paste0(outDir,today,"_densePlot_pOV.all_positions.pdf")

densePdf(data ,filename, plotvars, leg, var.labels, var.units)

#plot Angles####
    plotvars <- c("aMean","aAP","aDV","aLR")
    var.labels <- c("angular deviation","AP angle","DV angle","LR angle")
    var.units <- c("deviation °","angle from AP axis °","angle from DV axis °","angle from LR axis °")
    filename=paste0(outDir,today,"_densePlot_angles.pdf")
    
densePdf(data ,filename, plotvars, leg, var.labels, var.units)

#generate plot with AB cells on one side and P1 descendatns on the other ####
var <- c("LifeTime")
var.labels <- c("Cell cycle duration")
var.units <- c("min")


pdf(file=paste0(outDir,today,"_AB-P1_timing_densePlot.pdf"), paper = "a4r", width = 6, height = 6)
  par(mfrow=c(1,2),
      mar=c(2.5, 3.5, 1.5, 0), # bottom, left, top, right
      oma=c(0, 0, 0, 0.25), # bottom, left, top, right
      mgp= c(1,0.30,0),
      cex.main=1
      ) 

  ABcells <- Cellorder[grep("AB",Cellorder)]
  ABcells <- ABcells[4:31]
  
  Pcells <- Cellorder[!grepl("AB",Cellorder)]
  Pcells <- Pcells[1:29]
  
  #plot AB cells
  tufteplot.celllist(data,ABcells,var)
  legend("bottomleft",leg,cex=0.7, pch=20,col=pal)
  
  #plot P1 cells
  tufteplot.celllist(data,Pcells,var)
  
dev.off()

#plot dense plot for each cell in an VERTICAL tufte plot ####

# pdf(file=paste0(outDir,today,"_ABa-Cpp_Timing.netdis.pOV.densePlot.pdf"), paper = "a4", width = 7.5, height = 11)
#     
#     par(mfrow=c(length(plotvars),1),
#     mar=c(3, 3.5, 1.5, 0.1),
#     mgp = c(1.75, 0.5, 0)) # bottom, left, top, right
#     lw=0.75
#     
#   for (v in plotvars) {
#       g <- nlevels(data$Group)
#       
#       #subset the data and order cells by predefined order
#       tempdf <- data[,grep(v, colnames(data[1,]))]
#       colnames(tempdf) <- sub(paste0("\\.",v),"",colnames(tempdf[1,]))
#       tempdf <- tempdf[,match(Cellorder,colnames(tempdf[1,]), nomatch=0)]
#       
#       #find min and max of the tempdf (for cell range queried) REPLACE
#       minY <- min(tempdf[,1:limC], na.rm = T)
#       maxY <- max(tempdf[,1:limC], na.rm = T)
#       ylim=c(minY*0.8,maxY)
#       
#       #if cell is not in the data (such as first 4 cells) generate empty plotted column
#       firstpass <- 0
#       for (a in 1:limC) {
#           cell <- Cellorder[a]
#           match.col <- which(colnames(tempdf) %in% cell) # if the cell is not in the table this is 0
#           
#           if(length(match.col)==1){# column with given name exists
#             if(!firstpass){# plotting for the first time, new plot needs to be initiates
#               x <- tempdf[,cell] #subset the data
#               tufte.boxplot(x, data$Group, xlim=c(1,limC*g), ylim=ylim, 
#                             at=((a-1)*g),  pal=pal, cex=0.5, 
#                             ylab=v, lw=lw)
#               firstpass <- 1
#             }
#             else {
#               x <- tempdf[,cell]
#               tufte.boxplot(x, data$Group, 
#                             at=((a-1)*g), pal=pal, cex=0.5, lw=lw, add=T)
#             }
#           }
#           
#       }
#       axis(1,at=seq(2.5,g*limC,by=g),labels=Cellorder[1:limC],cex.axis=0.8, adj=0, srt=90, tick=T, tcl=-0.3, las=2, mgp=c(0,0.5,0))
#       abline(v=seq(0,g*limC,by=g)+0.5, lty=1, lwd=0.25, col="grey")
#       title(main=v, cex=0.8)
#       
#       #plot significance stars
#       pv <- data.frame(results[[v]])
#       pv <- pv[rownames(pv)%in%Cellorder[1:limC],grep("p.value",names(pv[1,]))]
#           pv$signif.stars <- ""
#           pv$signif.stars[pv$hatched.vs.died.p.value<0.05] <- "*"
#           pv$signif.stars[pv$hatched.vs.died.p.value<0.01] <- "**"
#           
#           xi <- match(rownames(pv),Cellorder[1:limC])-1
#           text(x=xi*g+2.5, y=minY, pv[,"signif.stars"], cex=0.8)
#   }
#   leg = c("Wild-type","Unequal controls", "Equalized alive", "Equalized dead")
#   legend("topleft",leg,cex=0.8, pch=20,col=pal)
# dev.off()

rm(data, limC, plotvars, var.labels, var.units)
