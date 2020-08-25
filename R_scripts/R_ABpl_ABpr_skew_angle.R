if(!require(colorRamps)) install.packages("colorRamps")
# TO DO
  # calculate ABal ABar distance and angle and plot
  # plot everything based on normalized X, Y, Z coordinates
  # calculate trajectories and net distance - both normalized and absolute
  # reorient embryos to canonical orientation 
  # aggregate positions of groups of cells to calculate how they move with respect to each other
    # ABar, ABal, MS, E, C, P3

# Notes ####
# ventral side defined by MS is rotating prior to the MS division - 
# MS moves from lateral side to more paralel with the cover slip

# ABpl division happans nearly perpendicular to the c-slip and ABpl ABpr then separate and move
# so that ABpl ends up  on the bottom left (paralel to the c-slip) and ABpr moves laterally to top-right side and more posterior

# measure distance of ABpl and ABpr - normalized to the embryo size
# measure in plane angle between these two cells relative to the AP axis
    # let ABpl be defined by x1 y1 z1 coordinates and ABpr x2 y2 z2 / normalized by length witdh and height of the embryo
    # then calculate euclidian distance of the two spots (long side "c" of the right triangle)
    # adjacent side of the triangle is defined by a = x2 y1 z1 and x2 y2 z - angle betwen the two is arcsin(a/c)
# End Notes ####

angleAndDist <- function(id, sub, cellPair){
  
  sel <- sub$Cell %in% cellPair
  
  #cell cell distance over time
  distance <- rbind(by(sub[sel,1:3], sub$frame[sel], function(x) sqrt( (x$X[1] - x$X[2])^2 + (x$Y[1] - x$Y[2])^2 + (x$Z[1] - x$Z[2])^2 )))
  
  #in plane distance corresponding to the opposite side of the right triangle along the AP axis 
  #(assumes that embryos are aligned correctly along the AP axis)
  
  a <- rbind(by(sub[sel,1:3], sub$frame[sel], function(x) sqrt((x$Y[1] - x$Y[2])^2 + (x$Z[1] - x$Z[2])^2 )))
  #angle in degrees
  angle <- asin(a/distance)*57.3 # get angle in degrees
  
  
  m <- matrix(c(distance, angle), nrow=length(angle), ncol=2)
  m <- cbind(with(sub[sub$Cell %in% cellPair[1],], time), m)
  colnames(m) <- c("time","dist","skew")
  m <- as.data.frame(m)
  m$embryo <- id 
  m
}

subsetData <- function(tab, cells){
  sub <- tab[grep(paste0("^(",paste(cells, collapse = "|"),")$"),tab$Cell), c("X.um","Y.um","Z.um","X.rel","Y.rel","Z.rel","Cell","frame","time")]
  colnames(sub)[1:3] <- c("X","Y","Z")
  sub <- sub[order(Cell),]
  sub$Cell <- factor(sub$Cell)
  sub
}


par(mar=c(3, 3, 1.5, 1), # bottom, left, top, right
    oma=c(0, 0, 0, 0), # bottom, left, top, right
    mgp = c(1.75, 0.5, 0)
) 

# Recalculate nuclei data to microns and relative dimensions ####
time.max <- 61 #restrict up to which time to consider cells for alignment (drift/errors, etc)
newNuclei <- list()
for (id in names(Nuclei)){
  embryoDims <- unlist(Embs[id, c("length","width","height","pixel","Z_slicing","interval","timeScaling")])
  #print(id)
  tab <- Nuclei[[id]]
  
# Rescale and align selected embryo ####
  tab$frame <- as.integer(sub("t","",tab$timepoint))
  tab$Cell <-  trimws(tab$IDENTITY, which="both")
  tab$IDENTITY <- NULL
  
  #calculate real time + set ABa division as time 0
  tab$time <- tab$frame * (embryoDims["timeScaling"]*embryoDims["interval"]/60)
  fr <- with(tab, which(Cell=="ABal" & frame==min(frame[Cell=="ABal"])))
  tab$time <- tab$time - tab$time[fr]
  
  #rescale coordinates to um
  um.coords <- t(apply(tab[,c("X","Y","Z")],1, function(x)x*embryoDims[c("pixel","pixel","Z_slicing")]))
  colnames(um.coords) <- c("X.um","Y.um","Z.um")
  
  #align and center positions 
      # find extreme positions of Nuclei in each direction
      # center nuclei within the physical bounds defined by embryo dims
    bounds <- apply(um.coords[tab$time<=time.max,],2,function(x)range(x))
    
      #dimensions defined by the most extreme Nuclei
    innerDims <- bounds[2,]-bounds[1,]
      
      #how different is this from the physical dimension
    diff <- embryoDims[1:3]-innerDims
    
      #if(any(diff<0)) print("embryo is larger than dimensions in the Embs") break
      if(any(diff<0)) diff[diff<0] <- bounds[1,diff<0]
    
      #diff/2 is how far from left, top and Z should the min coordinates sit in order to center the nuclei inside the shell
      #what is the difference from min bounds? 
    correct <- diff/2 - bounds[1,]
      #apply correct vector like bellow to all indices (in microns) to correct them
    um.coords <- t(apply(um.coords, 1, function(x)x+correct))
    
      #calculate relative coordinates and merge everything
    rel.coords <- t(apply(um.coords, 1, function(x)x/embryoDims[1:3]))
    colnames(rel.coords) <- c("X.rel","Y.rel","Z.rel")
    
    tab <- cbind(tab, um.coords, rel.coords)
    
    #write it in the list
  newNuclei[[id]] <- tab
}

# Plot selected cells with depth cueing  ####

for (id in names(newNuclei)){
  if(id=="MM4") next #MM4 coordinates and timing are wrong, cause it starts at 6C and not at 4C as all other embryos
  plotCells <- c("ABal","ABar","ABpl","ABpr","C","MS","P3")  
  tab <- newNuclei[[id]]
  
  sub <- tab[grep(paste0("^(",paste(plotCells, collapse = "|"),")$"),tab$Cell), c("X.um","Y.um","Z.um","X.rel","Y.rel","Z.rel","Cell","frame","time")]
    colnames(sub)[1:3] <- c("X","Y","Z")
    sub <- sub[order(Cell),]
    sub$Cell <- factor(sub$Cell)
  
  palette <- matlab.like2(30)
  l <- length(palette)

  par(mar=c(3, 3, 1.5, 1), # bottom, left, top, right
      oma=c(0, 0, 0, 0), # bottom, left, top, right
      mgp = c(1.75, 0.5, 0)
  ) 
  
 plot(0,0, type="n", xlim=c(0,1), ylim=c(0,1), xlab="AP axis", ylab="LR axis", main=paste(id, " - ABpl/ABpr separation"), cex.main=1)
    by(sub[,c("X.rel","Y.rel","Z.rel")], sub$Cell, function(x){ 
      lines(x$X.rel,x$Y.rel,col="black", type="l");
      points(x$X.rel,x$Y.rel,col=palette[round(x$Z.rel*l*1.2)], pch=20, cex=0.7)
      }
    )
    label.coords <- by(sub[,c("X.rel","Y.rel","Z.rel", "frame")], sub$Cell, function(x){
      x[which(x$frame==max(frame)),1:3]
    })
    #label cell tracks
    label.coords <- rbindlist(label.coords)
    text(label.coords$X, label.coords$Y, plotCells, cex=.8, pos=1)
    
    #legend scale
    points(seq(from=0.03, to=0.2,length.out=l), rep(0.05,l), col=palette, pch=20)
    text(c(0,0.23,0.11),c(0.05,0.05,0.1), c(0,1,"DV position"), cex=0.7)
}

# Calculate angle between ABpl and ABpr   ####
cellpair=c("ABar","MS")
cellpairname <- paste(cellpair,collapse = "-")
x.lim=c(0,65)
y.lim=c(0,20)
y.lim.angle=c(0,90)

rm(skew, skewtab, skewMat)
distAndSkew <- function(cellpair){
  skew <- list()
  
  pdf(paste0(outDir,today,"_",cellpairname,"_dist_angle.pdf"), width = 3.5, height = 2.5,pointsize = 10)
    par(mar=c(3, 3, 1.5, 3), # bottom, left, top, right
        oma=c(0, 0, 0, 0), # bottom, left, top, right
        mgp = c(1.75, 0.5, 0)
    )
    
  for (id in names(newNuclei)){
    #if(id=="MM4") next;
    sub <- subsetData(newNuclei[[id]], cellpair)
    if(nrow(sub)<2) next;
    m <-  angleAndDist(id, sub, cellpair)
    
    # Plot distance and angle between cells ####
      plot(m[,c("time","dist")], cex=0.8, ylab="cell-cell distance [um]", main=paste(id,"-",cellpairname, "angle and dist."), cex.axis=0.8, cex.lab=1,
           ylim=y.lim,xlim=x.lim)
      loessLine(m[,"time"], m[,"dist"])
      
      par(new = TRUE)
      plot(m[,c("time","skew")], pch=20, type = "l", lwd=2, axes = FALSE, bty = "n", ylab="",xlab = "", col=pal[1],
          ylim=y.lim.angle, xlim=x.lim)
      axis(side=4, col=pal[1], cex.axis=0.8, col.lab=pal[1])
      mtext("skew angle [degrees]", side=4, line=1.7,col=pal[1], cex=1)
        
  skew[[id]] <-m 
  
  }
  dev.off()
    skewtab <- do.call(rbind, skew)
    skewtab$group <- Embs[skewtab$embryo,"Group"]
  
  skewtab  
}

skewtab <- distAndSkew(cellpair)

# Plot all embryos in a single graph ####

plot(skewtab[,c("time","dist")], type="n", cex=0.8, ylab="cell-cell distance [um]", main=paste(cellpairname ,"separation"), cex.axis=0.8, cex.lab=1, 
     xlim=range(skewtab$time, na.rm = T),ylim=range(skewtab$dist, na.rm = T))

by(skewtab, skewtab$embryo, function(df){
 lines(df$time, df$dist, col=pal[as.numeric(df$group)])
})

plot(skewtab[,c("time","skew")], type="n", cex=0.8, ylab="angle [degrees]", main=paste(cellpairname ,"angle"), cex.axis=0.8, cex.lab=1, 
     xlim=range(skewtab$time, na.rm = T),ylim=range(skewtab$skew, na.rm = T))

by(skewtab, skewtab$embryo, function(df){
  lines(df$time, df$skew, col=pal[df$group])
})

# compute average curves and sd + plot ####

skewMat <- matrix(nrow=length(levels(skewtab$group)), ncol = ceiling(max(skewtab$time)/2))
skewMat <- as.data.frame(skewMat)
rownames(skewMat) <- levels(skewtab$group)

#innitiate empty lists
skewMat <- list(mean=skewMat, sd=skewMat, n=skewMat, sem=skewMat)
distMat <- skewMat

for(gr in levels(skewtab$group)){ #go by group
  sel <- skewtab$group==gr
  for (t in 1:ceiling(max(skewtab$time)/2)) {# step by 2 minutes to get the matrix
    # starts at t=0
    lim <- t*2
    ind = with(skewtab, time >= (lim-2) & time < lim)
    v1 = na.omit(skewtab$skew[ind & sel])
    v2 = na.omit(skewtab$dist[ind & sel])
    
    skewMat$mean[gr,t] <- mean(v1)
    skewMat$sd[gr,t] <- sd(v1)
    skewMat$n[gr,t] <- length(v1)
    skewMat$sem[gr,t] <- sd(v1)/sqrt(length(v1))
    

    distMat$mean[gr,t] <- mean(v2)
    distMat$sd[gr,t] <- sd(v2)
    distMat$n[gr,t] <- length(v2)
    distMat$sem[gr,t] <- sd(v2)/sqrt(length(v2))
  }
}

plot(0,0, xlim=c(0,ncol(skewMat$mean)*2),ylim=range(skewMat$mean,na.rm=T)*c(0.9,1.05), type="n", main=paste(cellpairname, "skew angle"), ylab="angle with respect to AP axis [degrees]", xlab="time [min]", cex.main=0.8)
for (row in rownames(skewMat$mean)) {
  lines(0:(ncol(skewMat$mean)-1)*2, skewMat$mean[row,], col=pal[which(rownames(skewMat$mean)==row)], lwd=2)
  lines(0:(ncol(skewMat$mean)-1)*2, skewMat$mean[row,]+skewMat$sem[row,], col=pal[which(rownames(skewMat$mean)==row)], lty=3)
  lines(0:(ncol(skewMat$mean)-1)*2, skewMat$mean[row,]-skewMat$sem[row,], col=pal[which(rownames(skewMat$mean)==row)], lty=3)
}

plot(0,0, xlim=c(8,ncol(distMat$mean)*2),ylim=range(distMat$mean,na.rm=T)*c(0.9,1.05), type="n",  main=paste(cellpairname, "distance"), ylab="distance between cells [um]", xlab="time [min]", cex.main=0.8)
for (row in rownames(distMat$mean)) {
  lines(0:(ncol(distMat$mean)-1)*2, distMat$mean[row,], col=pal[which(rownames(distMat$mean)==row)], lwd=2)
  lines(0:(ncol(distMat$mean)-1)*2, distMat$mean[row,]+distMat$sem[row,], col=pal[which(rownames(distMat$mean)==row)], lty=3)
  lines(0:(ncol(distMat$mean)-1)*2, distMat$mean[row,]-distMat$sem[row,], col=pal[which(rownames(distMat$mean)==row)], lty=3)
}

legend("bottomleft",legend=c("Wild-type","Uequal ctrl","equal. alive", "equal. dead"), col=pal, cex=0.8, lty = 1, lwd=2)

# max distance and mean skew ####

meanSkew <- aggregate(skewtab$skew, by=list(skewtab$embryo), mean, na.rm=T)
names(meanSkew) <- c("embryo","angle")
meanSkew$Group <- Embs[meanSkew$embryo,"Group"]
boxplot(angle~Group, data=meanSkew)
compare_means(angle~Group, meanSkew)
#max angle between cells is different for all ev571 embryos incl. ctrls from wt

maxDist <- aggregate(skewtab$dist, by=list(skewtab$embryo), max, na.rm=T)
names(maxDist) <- c("embryo","maxdist")  
maxDist$Group <- Embs[maxDist$embryo,"Group"]
boxplot(maxdist~Group, data=maxDist)
compare_means(maxdist~Group, maxDist)

# play with angles
      # sin(pi/4) # 0.707 #45 degrees
      # 2/2.83    # 0.707 which corresponds to square of side 2
      # 
      # asin(0.707) #  45 degrees = pi/4 = 0.785
      # # 57.3 degrees per radian
      # asin(0.707)*57.3
      # 
      # sin(pi/2) # 90 degrees
      # sin(pi)
      # 
      # #calculate in degrees: e.g. sin of 90 degrees:
      # sin(90*pi/180)
      # 
      # a <- c(2,2)
      # sqrt(sum(a*a)) # pythagoras theorem - but also norm of the vector (length)

