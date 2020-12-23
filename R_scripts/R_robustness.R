# FUNCTIONS ####
#filter data for testing (keep vars with min number of observations)

FilterAndSubset <- function(Df, rows, cols, minN=5){
  temp <- as.data.frame(Df[rows, cols])
  eqlzd <- Df[rows,"Experiment"]=="meta"
  n.eq <- sum(eqlzd)
  
  # remove variables with minimal variance
  v <- apply(temp, 2, sd, na.rm=T)/apply(temp, 2, mean, na.rm=T)<0.01
  v <- v[!is.na(v)]
  if(any(v)) temp <- temp[,names(v)[!v]] else cat("no rows with 0 variance \n") 
  
  # set minimal number of observation in each group of embryos, otherwise exclude the variable
  obsCol <- apply(temp,2,function(x) by(x, factor(Df[row.names(temp),"Group"]),function(y) sum(!is.na(y)))) # how many observations in each column by group?
  badCols <- apply(obsCol, 2, function(x) any(x<minN, na.rm=T))
  cat(sum(badCols),"variables removed containg less than", minN, "values per group: \n", colnames(temp)[badCols],"\n\n")
  temp <- temp[,!badCols]
  
  temp <- cbind(Df[rownames(temp),1:6],temp)
  temp
}
variable.type <- function(varNames){
  group <- rep("Other", length(varNames))
  group[grepl("Time",varNames)] <- "Timing"
  group[grepl(".\\.p",varNames)] <- "Position"
  group[grepl(".\\.a",varNames)] <- "Angles"
  group[grepl("dis",varNames)] <- "Mvmnt"
  group <- factor(group)
  group
}
plotCats <- function(crosstab, column, heading){
  varNames <- rownames(crosstab)[crosstab[,column]]
  varTyps <- variable.type(varNames)
  varGr <- cell.groups(varNames)
  Signif <- data.frame(vars=varNames,type=varTyps,gr=varGr)
  barplot2(table(Signif[,2:3]), legend.text = T, col=colors[2:5], main=heading)
}

#prepare data ####
Df <- FilterAndSubset(TempDf, rows, cols, minN=5) 
# minN indicates minimal number of observations for each variable in every group
#699 variables removed containg less than 5 values per group

Df <- Df[,c(4,7:ncol(Df))]

cat("filtered table contains: ", ncol(Df)-1, " variables \n") #1608

# Compare equalized against controls for variance and mean change #### 
g <- factor(Df$Group)
g2 <- as.character(Df$Group)
g2[g2 %in% c("dead","alive")] <- "equalized"
g2 <- factor(g2)

test.res <- apply(Df[,2:ncol(Df)], 2, function(x)t.test(x~g2)$p.value) # two sided
var.res <- apply(Df[,2:ncol(Df)], 2, function(x)var.test(x~g2, alternative="less")$p.value) #null control has same variance, alternative: true ratio of variances is less than 1

# better use benjamini-hochberg

test.adj <- p.adjust(test.res, method ="BH")
sum(test.adj<0.05) #163

var.adj <- p.adjust(var.res, method ="BH")
sum(var.adj<0.05) # 37 signif with greater variance, 0 with alternative="greater"

crosstab <- data.frame(mean.test=test.adj<0.05, variance=var.adj<0.05, row.names = colnames(Df[2:ncol(Df)]))
table(crosstab)
#             variance
# mean.test   FALSE TRUE
# FALSE       1422  23
# TRUE        149   14

crosstab$mean.pvalue <- test.adj
  crosstab$mean.crtl <- apply(Df[2:ncol(Df)],2, function(x) by(x, g2, mean, na.rm=T)[[1]])
  crosstab$mean.eq <- apply(Df[2:ncol(Df)],2, function(x) by(x, g2, mean, na.rm=T)[[2]])
  crosstab$mean.ratio <- with(crosstab, mean.eq/mean.crtl)

crosstab$var.pvalue <- var.adj
  crosstab$var.crtl <- apply(Df[2:ncol(Df)],2, function(x) by(x, g2, var, na.rm=T)[[1]])
  crosstab$var.eq <- apply(Df[2:ncol(Df)],2, function(x) by(x, g2, var, na.rm=T)[[2]])
  crosstab$var.ratio <- with(crosstab, var.eq/var.crtl)
  
  crosstab$n.ctrl <- apply(Df[2:ncol(Df)],2, function(x) by(x, g2, function(y) sum(!is.na(y)))[[1]])
  crosstab$n.eq <- apply(Df[2:ncol(Df)],2, function(x) by(x, g2, function(y) sum(!is.na(y)))[[2]])

crosstab$Category <- variable.type(rownames(crosstab))
  
write.xlsx(crosstab, file="Fig.6-Supplement-3-Source-data.xlsx", rowNames=T)
  
par(mar=c(3.5,2.5,3,0.5))

# plot the counts of significant features for variance change and mean change according to the variable type and cell lineage of origin
pdf(paste0(outDir,today,"_mean_and_variance_Change.pdf"),width = 3, height = 4, pointsize = 10, family='Helvetica', useDingbats = FALSE)
  par(mar=c(2.3,2,1,0.2), mgp=c(1.25,0.5,0), cex.lab=0.8, cex.axis=0.7, cex.main=0.8, tcl=-.35, mfrow=c(2,2), las=2)

  plotCats(crosstab,1,paste0("Ctrl/Equal - Mean Change [n=",sum(crosstab[,1]),"/",nrow(crosstab),"]"))
  plotCats(crosstab,2,paste0("Ctrl/Equal - Increased variace [n=",sum(crosstab[,2]),"/",nrow(crosstab),"]"))
  
  #only var increase
  plotCats(crosstab[!crosstab$mean.test & crosstab$variance, ],2,paste0("Ctrl/Equal - Same mean, increased variace [n=23/",nrow(crosstab),"]"))
dev.off()

