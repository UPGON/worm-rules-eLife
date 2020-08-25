
####  FUNCTIONS #####
getTstatistics<-function(data, var, groupcol, ctrl, mutantNames){
  #input data table is in Long format - with columns embryo, Cell, Group, Experiment followed by variables from col 9:22
  out={};
  rows1<-which(data[,groupcol]==ctrl)
  
  for(cond in mutantNames){# for each mutant condition compare it to the control
    #subset table for given mutant and condition
    rows2<-which(data[,groupcol]==cond)
    testdf <- data[c(rows1,rows2),c(1:5,grep(var, colnames(data)))]
    testdf[,groupcol] <- factor(testdf[,groupcol], levels = c(ctrl, cond)) #relevel the factor
    testdf[,"Cell"] <- factor(testdf[,"Cell"]) #relevel the factor
    
    x <- by(testdf, testdf$Cell, test, by=groupcol, var=var)
    x <- do.call("rbind", x)
    
    colnames(x) <- paste(c(rep(ctrl,4),rep(cond,4),rep(paste(ctrl,cond,sep=".vs."),2)),c("mean", "sd","vc","n","mean","sd","vc","n","t.stat","p.val"), sep=".")
    x <- as.data.frame(x)
    x$Cell <- rownames(x)
    
    if(!any(grepl(ctrl,colnames(out)))) out<-x else out<- merge(out, x[,5:11], by="Cell", all=T)
    rownames(out) <- out$Cell
  }
  out;
}
vc.adj <- function(v){
  (1+1/(4*length(v)))*(sd(v)/mean(v))
  }
test <- function(df, by, var){
  l <- levels(df[,by])
  control <- na.omit(df[df[,by]==l[1],var])
  cond <- na.omit(df[df[,by]==l[2],var])
  result <- rep(NA,10)
  if(length(control)>=limitN && length(cond)>=limitN){
    result <- t.test(control,cond)[c(1,3)]; #get absolute value of T statistic
    result <- unlist(result)
    result <- c(mean=mean(control),sd=sd(control),vc=vc.adj(control),n=length(control),mean2=mean(cond),sd2=sd(cond),vc2=vc.adj(cond), n2=length(cond), result)
  }
  result
}

######### Define Groups ############
limitN <- 4 #minimal number of observation in each group to be considered for comparison

# Embryos with EMS tilted were not removed here``
# use the stricter definition of equalized with AB 48-52.5 %
# we can relabel groups with inverted and near.eq 
data <- Complete
data$Group <- as.character(data$Group)
data$Group[data$embryo %in% Embs$ID[inverted]] <- "INV"
data$Group[data$embryo %in% Embs$ID[partial]] <- "PART" #partially equalized
data$Group <- factor(data$Group, levels = c("wt","ctrl","alive","dead","INV","PART"))
levels(data$Group) <- c("WT","C","EA","ED","INV","PART")

data$Experiment <- as.character(data$Experiment)
data$Experiment[data$Group %in% c("EA","ED")] <- "EQ"
data$Experiment[data$Experiment=="ctrl"] <- "C"

#We have clear groups defined in Complete$Group
mutantNames <- c("EA", "ED", "WT", "INV")
ctrl <- "C"

vars <- colnames(Complete)[9:24]

######### Calculate statistics for each variable and each mutant vs control comparison +  ############
results <- list()
for (var in vars){
  rtable <- getTstatistics(data, var, "Group", ctrl, mutantNames) #data, var, groupcol, ctrl, mutantNames
  
    cnames <- colnames(rtable)  
    keep <- grep("WT", cnames)
    keep <- c(1,keep,grep("(t.stat|val)", cnames))
  rtable <- rtable[,unique(keep)]
  
  rtable1 <- getTstatistics(data, var, "Group", "EA", "ED")
  rtable1$ED.vs.EA.effect <- with(rtable1, ED.mean/EA.mean)
  
  rtable2 <- getTstatistics(data, var, "Experiment", "C","EQ") #compare all equalized with ev571 control
  rtable2$EQ.vs.C.effect <- with(rtable2, EQ.mean/C.mean)

  rtable3 <- getTstatistics(data, var, "Group", "EA", "INV")
  rtable3$INV.vs.EA.effect <- with(rtable3, INV.mean/EA.mean)
  rtable3 <- rtable3[,-(1:4)]

  res.table <- Reduce(function(x, y) merge(x, y, by="Cell", all=TRUE), list(rtable, rtable1, rtable2, rtable3))
  row.names(res.table) <- res.table$Cell
  cnames=colnames(res.table)[2:ncol(res.table)]
  cols <- c(grep("WT",cnames),grep("^C",cnames),grep("^EQ",cnames),grep("^EA",cnames),grep("ED",cnames),grep("INV",cnames))
  cols <- unique(cols)
  
  res.table <- res.table[,cols+1]
  cnames=colnames(res.table)
  cols <- c(grep("\\.mean$",cnames),
            grep("\\.effect$",cnames),
            grep("\\.sd$",cnames),
            grep("\\.vc$",cnames),
            grep("\\.n$",cnames),
            grep("\\.stat$",cnames),
            grep("\\.val$",cnames)
  )
  
  res.table <- res.table[,cols]
  
  #order cells as they appear in the embryo
  c <- rownames(res.table)
  res.table <- res.table[c[match(Cellorder,c, nomatch = 0)],]

  #correct p.values for multiple testing with Bonferroni-Hochfer
  cnames=colnames(res.table)
  p.adj <- apply(res.table[,grep("\\.val$",cnames)],2, function(x) p.adjust(x, method="BH"))
  colnames(p.adj) <- sub("val","adj",colnames(p.adj))
  res.table <- cbind(res.table, p.adj)
  
  #convert significance to star notation
  sign <- as.data.frame(starformat(res.table[,grep("adj",colnames(res.table))]))
  colnames(sign) <- sub("p.adj","signif",colnames(sign))
  
  res.table <- cbind(Order=1:nrow(res.table),
                     res.table[,1:6],
                     EQ.over.C=res.table$EQ.vs.C.effect,
                     EQ.C.sign=sign$C.vs.EQ.signif,
                     
                     ED.over.EA=res.table$ED.vs.EA.effect,
                     ED.EA.sign=sign$EA.vs.ED.signif,
                     
                     INV.over.EA=res.table$INV.vs.EA.effect,
                     INV.EA.sign=sign$EA.vs.INV.signif,
                     
                     res.table[,10:ncol(res.table)],
                     sign[,c(1:4)]
                     )
  
  results[[var]]<- res.table
}
hs <- createStyle(textDecoration = "BOLD", fontSize=11, numFmt="0.00")
options("openxlsx.numFmt" = "0.0")

#Save results as xlsx table (items of the list are saved as tabs)
write.xlsx(results,paste0(outDir,today,"_all_Vars.statistics.xlsx"), colNames = TRUE, rowNames=T, borders = "rows",colWidths="auto",firstCol=T,firstRow=T, headerStyle = hs)
write.xlsx(Embs,paste0(outDir,today,"_all_lineaged_embs.xlsx"), colNames = TRUE, rowNames=F, borders = "rows",colWidths="auto",firstCol=T,firstRow=T, headerStyle = hs)

print(paste("Success! xlsx file saved at:",paste0(outDir,today,"_all_Vars.statistics.xlsx")))
