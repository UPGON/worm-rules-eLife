#read data through R_Lineage_PCA.r first
#install.packages(c("rgl", "car", "matlib))

library(car)
library(shape)

EightCellStage=c("ABpl","ABpr","E","MS","ABal","ABar","C","P3")
SixteenC.Stage <- c("ABpla","ABplp","ABpra","ABprp","Ea","Ep","MSa","MSp","ABala","ABalp","ABara","ABarp","Ca","Cp","P4","D")

CellEllipses <- function(tab, cells, axs="XY", plotpoints=F, conf=0.5, sem=F){
  
  tab <- tab[tab$Cell %in% cells,]
  if(axs=="XY")y=tab$pLR else y=tab$pDV
  
  #calculate stats for cell
  ag <- do.call(data.frame, aggregate(tab[,c("pAP","pLR","pDV")], by=list(tab$Cell, tab$Outcome), FUN=row.stats))
  colnames(ag)[1:2] <- c("Cell","Outcome")
  cn <- colnames(ag)
  ag.y <- ag[,if(axs=="XY") grep("pLR",cn) else grep("pDV",cn)]
  
  if(length(cells)==8){
    clr = c(rbind(rep("#D3D3D3",length(cells)), colors[sort(cells)])) 
  } else {
    clr = c(rbind(rep("#D3D3D3",length(cells)), colors[get.motherCell(sort(cells))]))
  }
  
  plotchar=""
  if(plotpoints) plotchar=rep(c(4,1),length(cells))
  
  dataEllipse(tab$pAP, y,
              xlim=c(5,mean(Embs$length)-5),
              ylim=c(2,mean(if(axs=="XY") Embs$width else Embs$height)-2),
              groups=factor(tab$G), 
              levels=c(conf),
              col = clr,
              pch = plotchar,
              cex=1,
              lwd=1, fill=T, fill.alpha=0.05,
              ylab=ifelse(axs=="XY","Y - embryo width [um]","Z - embryo heigth [um]"),
              xlab="X - embryo length [microns]",
              group.labels=c(rbind(rep("",length(cells)),sort(cells))),
              center.cex=0
  )
  if(sem==T){
    #plot sem in x
    arrows(x0=ag$pAP.m+ag$pAP.sem, y0=ag.y[,1], 
           x1=ag$pAP.m-ag$pAP.sem, y1=ag.y[,1],
           length=0.02, angle=90, code=3,
           lwd=1,
           col="black")
    #plot sem in y
    arrows(x0=ag$pAP.m, y0=ag.y[,1]+ag.y[,4], 
           x1=ag$pAP.m, y1=ag.y[,1]-ag.y[,4],
           length=0.02, angle=90, code=3,
           lwd=1,
           col="black")
  }
  Arrows(x0=ag$pAP.m[ag$Outcome=="hatched"], 
         y0=ag.y[ag$Outcome=="hatched",1], 
         x1=ag$pAP.m[ag$Outcome=="died"], 
         y1=ag.y[ag$Outcome=="died",1],
         code=2,
         lwd=2,
         col="#cc0000",
         arr.type="triangle", 
         arr.length=0.1,
         arr.width=0.1
         )
  
  if(plotpoints)legend("topright",c("dead", "alive"), pch=c(4,1), inset=0.01)
}

# 8-CELL STAGE ####
#get Cells for equalized embryos without outliers
tab <- Complete

ids <- Embs$ID[equalized]
ids <- ids[!ids%in%outliers]

tab <- tab[tab$embryo %in% ids,]
tab$pAP <- tab$pAP*mean(Embs$length)
tab$pDV <- tab$pDV*mean(Embs$height)
tab$pLR <- tab$pLR*mean(Embs$width)

tab <- tab[order(tab$Cell, tab$Outcome),]
tab$G <- factor(paste(tab$Cell,tab$Group, sep="."))
tab$Group <- factor(tab$Group)

# Dorsal view 8C ####
pdf(paste0(outDir,today,"_Cell_positions_8C.pdf"), width=6, height=5, pointsize = 10, useDingbats = F)
  par(mar=c(2.5,2.5,0.1,0.1), mgp=c(1.5,0.5,0),cex.lab=0.9, cex.axis=1, cex.main=1)
  layout(c(1,2),widths = c(1,1),heights = c(1,0.65))
  CellEllipses(tab, EightCellStage, conf=0.5)
  CellEllipses(tab, EightCellStage, "XZ", conf=0.5)
dev.off()
  
# Dorsal view 16 ####
pdf(paste0(outDir,today,"_Cell_positions_16C.pdf"), width=6, height=5, pointsize = 10, useDingbats = F)
  par(mar=c(2.5,2.5,0.1,0.1), mgp=c(1.5,0.5,0),cex.lab=0.9, cex.axis=1, cex.main=1)
  layout(c(1,2),widths = c(1,1),heights = c(1,0.65))
  CellEllipses(tab, SixteenC.Stage, "XY", conf=0.25)
  CellEllipses(tab, SixteenC.Stage, "XZ", conf=0.25)
dev.off()


# 3D plot ####
#library(rgl)
#library(matlib)

# m <- ag[,grep("p.*\\.m",colnames(ag))]
# 
# plot3d(m, type="s", size=1,
#        xlab="length",ylab = "width",zlab = "height",
#        xlim = c(5,mean(Embs$length)-5),
#        ylim = c(5,mean(Embs$width)-5),
#        zlim = c(2.5,mean(Embs$height)-2.5),
#        aspect = c(1,mean(Embs$width)/mean(Embs$length),mean(Embs$height)/mean(Embs$length)),
#        col=c(rep("#F93F17",8),colors[sort(EightCellStage)]),
#        box=F
#        )
# for(i in 1:(nrow(m)/2)){
#   arrow3d(p1=m[i,],
#           p0=m[nrow(m)/2+i,],
#           type="rotation",col="black",
#           width=0.5
#           )
# }

# 16-CELL STAGE ####
# 3D plot 16 Cell stage ####
# m <- ag[,grep("p.*\\.m",colnames(ag))]
# plot3d(m, type="s", size=1,
#        xlab="length",ylab = "width",zlab = "height",
#        xlim = c(5,mean(Embs$length)-5),
#        ylim = c(5,mean(Embs$width)-5),
#        zlim = c(2.5,mean(Embs$height)-2.5),
#        aspect = c(1,mean(Embs$width)/mean(Embs$length),mean(Embs$height)/mean(Embs$length)),
#        col=c(rep("#F93F17",16),rep(colors,each=2)[-c(11,17:19)]),
#        box=F
# )
# for(i in 1:(nrow(m)/2)){
#   arrow3d(p1=m[i,],
#           p0=m[nrow(m)/2+i,],
#           type="rotation",col="black",
#           width=0.5
#   )
# }
# texts3d(m[1:16,]+c(2,1,1),texts=sort(SixteenC.Stage))
# writeWebGL()
# 
