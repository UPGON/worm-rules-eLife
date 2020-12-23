# color ramp for heatmaps 
colramp <- colorRampPalette(
  c("#67001F", "#B2182B", "#D6604D", "#F4A582",
    "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", 
    "#4393C3", "#2166AC", "#053061"))(200)
colramp<-rev(colramp)

#functions
replace.NAs <- function(vec, group=NA) {
  if(!all(is.na(group))){
    means <- aggregate(vec,by=list(group),mean,na.rm=TRUE)
    means <- means$x
    names(means) <- levels(group)
    na <- is.na(vec)
    vec[na] <- means[group[na]] 
  } else {
    m <- mean(vec, na.rm = TRUE) 
    vec[is.na(vec)] <- m 
  }
  return(vec)
}

rsquare <- function(true, predicted) {
  sse <- sum((predicted - true)^2)
  sst <- sum((true - mean(true))^2)
  rsq <- 1 - sse / sst
  
  # For this post, impose floor...
  if (rsq < 0) rsq <- 0
  
  return (rsq)
}

text.labels <- function (x, y, words, cex = 1, show.lines = TRUE, ...) 
#modified function from textplot from wordcloud package
{
  lay <- wordlayout(x, y, words, cex, ...)
  if (show.lines) {
    for (i in 1:length(x)){
      xl <- lay[i, 1]
      yl <- lay[i, 2]
      w <- lay[i, 3]
      h <- lay[i, 4]
      if (x[i] < xl || x[i] > xl + w || y[i] < yl || y[i] > 
          yl + h) {
        nx <- xl + 0.5 * w
        ny <- yl + 0.5 * h
        lines(c(x[i], nx), c(y[i], ny), col = "grey")
      }
    }
  }
  text(lay[, 1] + 0.5 * lay[, 3], lay[, 2] + 0.5 * lay[, 4], words, cex = cex, ...)
}
embPCA <- function(Df, cols, anot=F, title="PCA", legpos="bottomright", return=F){
  embr.pca <- prcomp(Df[,cols], center = T, scale = T)
  s <- summary(embr.pca)
  
  opar <- par()
  par(mar=c(3,3,3,1),
      cex.main=0.8,
      cex.lab=0.8,
      cex.axis=0.8,
      mgp=c(2,1,0))
  plot(embr.pca$x[,c(1,2)], 
       col = pal[as.numeric(Df$Group)], 
       pch = 16, 
       main=title,
       xlab = paste0("PC1 (",round(s$importance[2], digits = 3)*100," %)"),
       ylab = paste0("PC2 (",round(s$importance[5], digits = 3)*100," %)")
  )
  if(anot) text(embr.pca$x[,1],embr.pca$x[,2], labels = Df$embryo, cex=0.5, pos = 1);
  
  print(title)
  exp.var <- s$importance[6]*100 #cummulative variance explained
  print(paste("cummulative variance explained by first two components: ",exp.var,"%")) 
  mtext(paste("first two principal components explain ", exp.var,"% variance"), cex=0.6,line=0.1)
  
  leg = c("Wild-type","Unequal ctrls", "Equalized alive", "Equalized dead", "Inverted")
  legend(legpos,leg,
         cex=0.7, pch=16, col=pal, inset = 0.01)
  #c(sapply(leg[1], function(x) as.expression(substitute(A~degree~"C",list(A = as.name(x))))), leg[2:4])
  
  # workout variables with highest loading
  m <- embr.pca$rotation[,1]
  names(m) <- rownames(embr.pca$rotation)
  m <- m[order(m, decreasing = T)]
  print("Vars with highest loading for PC1:", head(m,5) ) 
  print(head(m,5))
  print("Vars with highest negative loading for PC1:", head(m[order(m, decreasing = T)],5)) 
  print(tail(m,5))

  par <- opar
  if(return) return(embr.pca)
}

save_pheatmap_pdf <- function(x, filename, width=6.5, height=5) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
  dev.off()
}

volcano <- function(mat, outcome, title="volcano"){
  # compute fold change 
    
    groupmeans <- aggregate(mat, list(outcome), mean, na.rm=T)
    rownames(groupmeans) <- groupmeans$Group.1
    groupmeans$Group.1 <- NULL
    
  #died / hatched
    fc <- unlist(groupmeans[1,]/groupmeans[2,])
    log2fc <- log2(fc)
  
      # compute p-value
    pv <- apply(mat, 2, function(x)t.test(x~outcome)$p.value)  
  
  # Make a basic volcano plot
    xl <- c(min(log2fc, na.rm = T),max(log2fc, na.rm = T));
    yl <- c(0,max(-log10(pv)));
    
    plot(0,-1,
       xlim=xl, 
       ylim=yl, 
       pch=1, 
       main=title,
       ylab="-log10( p-value )",
       xlab="log2 fold-change dead/alive"     
    )
    #set gridlines
    abline(v=seq(min(log2fc,  na.rm = T),max(fc,  na.rm = T),0.25),col="black",lty="dotted", lwd=0.3)
    abline(h=seq(0,max(-log10(pv)),1),col="black",lty="dotted",lwd=0.3)
      
    #set limits 
    abline(v=0,col="#989898",lty="dashed") #mark 0
    abline(h=-log10(plim),col="black",lty="dashed") #mark plim
    abline(v=log2(c(1-fchange,1+fchange)),col="black",lty="dashed")
         
    mtext(paste("p-value cut-off:",plim,"; fold change cut-off:",fchange),1,line=-1, cex=0.6)

  points(log2fc, -log10(pv), pch=16, cex=0.8, col="#000000BF")
  
  # Add colored points: red if padj<p.lim, orange of fold.change > f.change limit
  signif <- pv<plim & (log2fc>log2((1+fchange)) | log2fc<log2((1-fchange)))
  points(log2fc[signif], -log10(pv)[signif], pch=16, col="firebrick2", cex=0.8)
  
#with(subset(res, p.value<plim & log2fc<log2((1-fchange))), points(log2fc, -log10(p.value), pch=16, col="firebrick2", cex=0.8))
  
  #point labels (only significant)
    text.labels(log2fc[signif],-log10(pv[signif]), words=names(log2fc)[signif], cex=0.7, xlim=xl,ylim=yl)
    
  #merge all in the table
    res <- data.frame(t(groupmeans), fold.change=fc, log2fc=log2fc, p.value=pv, p.adj=p.adjust(pv, method="BH"), signif=signif)
    return(res)
} #function to plot and save volcano plot from the variable matrix

getSubset <- function(Df, rows, cols, nalim=3, inpute = T, NAsColumn = 0.25){
  temp <- as.data.frame(Df[rows, cols])
  eqlzd <- Df[rows,"Experiment"]=="meta"
  n.eq <- sum(eqlzd)
  
  # remove variables with 0 variance
  v <- apply(temp, 2, stats::var, na.rm=T)==0
  if(any(v)) temp <- temp[,!v] else cat("no rows with 0 variance\n") 
  
  #remove variables with too few samples especially in EQUALIZED embryos
  nasCol <- apply(temp[eqlzd,],2,function(x)sum(is.na(x))) # how many NAs in each column?
  goodCols <- nasCol/n.eq < NAsColumn #maximum 25%
  cat(sum(!goodCols)," Variables containing more than", NAsColumn*100 ,"% of NAs among equalized embryos removed: \n", ifelse(any(!goodCols), paste(names(goodCols)[!goodCols], collapse=", "),"None"),"\n\n")
  if(any(!goodCols)) temp <- temp[,goodCols]
  
  #how many NAs in each row?
  nas <- apply(temp,1,function(x)sum(is.na(x)))
  goodRows <- nas<ncol(temp)/5 # remove embryos with more than 20% variables being NA
  badRows <- rownames(temp)[!goodRows]
  cat("Samples removed containg more than 20% NAs: ", ifelse(length(badRows)>0,paste(badRows, collapse = ","),"None") ,"\n\n")
  temp <- temp[goodRows,]
  
  # set minimal number of observation in each group of embryos, otherwise exclude the variable
  if(any(is.na(temp))){ 
    # how many NAs in each column by group?
    nasCol <- apply(temp,2,function(x) by(x, factor(Df[row.names(temp),"Group"]),function(y) sum(is.na(y)))) 
    badCols <- apply(nasCol, 2, function(x) any(x>nalim, na.rm=T))
    cat("Variables removed containg more than ", nalim, " NAs in each group: ", ifelse(!any(badCols),"None",names(badCols)[badCols]) ,"\n\n")
    temp <- temp[,!badCols]
  }

  if(inpute & any(is.na(temp))) {
      cat("Imputed ", sum(is.na(temp)), " NAs with group means \n\n" )
      temp <- apply(temp,2,replace.NAs,group=factor(Df[rownames(temp),"Group"])) #fill NAs with group averages
  } else print("No values imputed, table is complete")

  temp <- cbind(Df[rownames(temp),1:6],temp)
  
  temp
}

prep.data <- function(Df, cells, outlierRM=T, inpute=T, variance.th=0.025){
  #removes outliers, 
  #strips variables with <2.5% variance
  #fills NA with group means

  #all variables at a given stage ####
  q <- paste("^(",paste0(cells, collapse="|"),")\\..*",sep="")
  cols <- grep(q,colnames(Df))
  
  if(length(cells)>8) Df <- Df[Df$embryo %in% Embs$ID[Embs$MaxTime>70],] # for later stages remove embryos that were not lineaged/imaged far enough (75 minutes is limit - few guys for volume)
  
  # REMOVE OUTLIERS, 
  #because getSubset fills NA with values from each group, values would get skewed
  if(outlierRM) rows <- !(Df$embryo %in% outliers) else rows <- 1:nrow(Df)
  
  # function from the functions_PCA.R that removes  rows and columns with too many NAs and fills remaining NAs with group means
  Df <- getSubset(Df, rows, cols, 4, inpute=inpute)
  
  # row indices of equalized embyos
  eqlzd <- Df$embryo %in% Embs$ID[equalized]
  
  # keep only variables with more than 5% coeficienn of variance
  dt <- Df[,7:ncol(Df)]
  lowVariance <- apply(dt, 2, function(x) sd(x, na.rm = T)/mean(x,na.rm=T))> variance.th
  cn <- colnames(dt)
  dt <- dt[,lowVariance]
  cat("Removed", sum(!lowVariance) ," variables with variation coeficient < ",variance.th*100," %: \n", paste(cn[!lowVariance], collapse =", " ))
  cat("\n Filtered table contains ", ncol(dt), 'variables')
  
  Df <- cbind(Df[,1:6],dt)
  
  #cols <- colnames(Df)[7:ncol(Df)]
  
  Df
}

#http://www.sthda.com/upload/rquery_cormat.r
rquery.cormat <- function(x, type=c('lower', 'upper', 'full', 'flatten'),
         graph=TRUE, graphType=c("correlogram", "heatmap"),
         col=NULL, ...)
{
  library(corrplot)
  # Helper functions
  #+++++++++++++++++
  # Compute the matrix of correlation p-values
  cor.pmat <- function(x, ...) {
    mat <- as.matrix(x)
    n <- ncol(mat)
    p.mat<- matrix(NA, n, n)
    diag(p.mat) <- 0
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        tmp <- cor.test(mat[, i], mat[, j], ...)
        p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      }
    }
    colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
    p.mat
  }
  # Get lower triangle of the matrix
  getLower.tri<-function(mat){
    upper<-mat
    upper[upper.tri(mat)]<-""
    mat<-as.data.frame(upper)
    mat
  }
  # Get upper triangle of the matrix
  getUpper.tri<-function(mat){
    lt<-mat
    lt[lower.tri(mat)]<-""
    mat<-as.data.frame(lt)
    mat
  }
  # Get flatten matrix
  flattenCorrMatrix <- function(cormat, pmat) {
    ut <- upper.tri(cormat)
    data.frame(
      row = rownames(cormat)[row(cormat)[ut]],
      column = rownames(cormat)[col(cormat)[ut]],
      cor  =(cormat)[ut],
      p = pmat[ut]
    )
  }
  # Define color
  if (is.null(col)) {
    col <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", 
                              "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", 
                              "#4393C3", "#2166AC", "#053061"))(200)
    col<-rev(col)
  }
  
  # Correlation matrix
  cormat<-signif(cor(x, use = "complete.obs", ...),2)
  pmat<-signif(cor.pmat(x, ...),2)
  # Reorder correlation matrix
  ord<-corrMatOrder(cormat, order="hclust")
  cormat<-cormat[ord, ord]
  pmat<-pmat[ord, ord]
  # Replace correlation coeff by symbols
  sym<-symnum(cormat, abbr.colnames=FALSE)
  # Correlogram
  if(graph & graphType[1]=="correlogram"){
    corrplot(cormat, type=ifelse(type[1]=="flatten", "lower", type[1]),
             tl.col="black", tl.srt=45,col=col,...)
  }
  else if(graphType[1]=="heatmap")
    heatmap(cormat, col=col, symm=TRUE)
  # Get lower/upper triangle
  if(type[1]=="lower"){
    cormat<-getLower.tri(cormat)
    pmat<-getLower.tri(pmat)
  }
  else if(type[1]=="upper"){
    cormat<-getUpper.tri(cormat)
    pmat<-getUpper.tri(pmat)
    sym=t(sym)
  }
  else if(type[1]=="flatten"){
    cormat<-flattenCorrMatrix(cormat, pmat)
    pmat=NULL
    sym=NULL
  }
  list(r=cormat, p=pmat, sym=sym)
}

