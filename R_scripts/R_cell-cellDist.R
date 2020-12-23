# Let's compute and plot distance between ABar and C ans ABpl and ABpr ####

# function to compute 3D euclidian distance
euc.dist.3 <- function(x1, y1, z1, x2, y2, z2 ) sqrt( (x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2 )
euc.dist.tab <- function(m) sqrt( (m[,1] - m[,4])^2 + (m[,2] - m[,5])^2 + (m[,3] - m[,6])^2 )
pair.dist.plot <- function(cellpair, data=Df){
        
        ind <- grep(paste0("^(",paste(cellpair, collapse = "|"),")\\.p[ADL]"),colnames(Df))
        colnames(data)[ind]        
        
        cell.dist <- euc.dist.tab(Df[,ind])
        tab <- cbind(data[,c("embryo","Group","AB.rel")],cell.dist)
        plot.bxp(filename=paste0(outDir,""), data = tab, var1 = "cell.dist", ylab1 = paste0(paste(cellpair, collapse = " <-> "), " 3D distance [au]"), saveit = F)
        
        if(length(levels(tab$Group)) > 2){
                lab <- Tuk("cell.dist",groupvar ="Group",tab)
                text(x=1:5, y=rep(max(tab$cell.dist,na.rm = T)*0.95, 5), labels=lab$comp)
                pv <- wilcox.test(cell.dist~Group,data=tab[tab$Group%in%c("alive","dead"),])$p.value
        } else pv <- wilcox.test(cell.dist~Group, data=tab)$p.value
        
        text(x=2, y=max(tab$cell.dist,na.rm = T)*0.9, labels=paste0("p =" ,format(pv,digits=3)), cex=0.8)
        
        cat(paste0(cellpair)," wilcox p-value:", format(pv, digits=3), "\n")
        
        colnames(tab)[4] <- paste(cellpair,collapse = ".")
        return(tab)
}

#Is there a difference with relation to the Pha-4 expression pattern

Df <- WideDf[grep("PM",WideDf$embryo),]
Df$Group <- NA
Df$Group <- factor(ifelse(Df$embryo %in% phaDefects2, "Ectopic","Normal")) # embryos with ABala transformed vs the rest

pdf(paste0(outDir,today,"_cell-cell-distances_PHA-4-GFP.pdf"),width =6, height = 4.5, pointsize = 10, useDingbats=F)
par(mar=c(2,2.3,0.2,0.15), mgp=c(1.25,0.5,0), cex.lab=0.9, cex.axis=0.6, cex.main=0.8, tcl=-.35)
par(mfrow=c(2,3))

{
        
        tab <- pair.dist.plot(c("ABal", "MS"),Df) #p=0.62
        tab <- cbind(tab, ABala.MS=pair.dist.plot(c("ABala", "MS"),Df)[,4]) # p=0.96
        tab <- cbind(tab, ABala.MSa=pair.dist.plot(c("ABala", "MSa"),Df)[,4]) # p=0.27 n.s.
        
        tab <- cbind(tab, ABalp.MS=pair.dist.plot(c("ABalp", "MS"),Df)[,4]) # p=0.39 n.s.
        tab <- cbind(tab, ABalp.MSa=pair.dist.plot(c("ABalp", "MSa"),Df)[,4]) # p=0.69
        
        tab <- cbind(tab, ABar.MS=pair.dist.plot(c("ABar", "MS"),Df)[,4]) # p=0.62
        tab <- cbind(tab, ABara.MS=pair.dist.plot(c("ABara", "MS"),Df)[,4]) # p=0.014
        tab <- cbind(tab, ABara.MSa=pair.dist.plot(c("ABara", "MSa"),Df)[,4]) # p=0.005
        
        tab <- cbind(tab, ABar.C=pair.dist.plot(c("ABar", "C"),Df)[,4]) # p=1
        tab <- cbind(tab, ABara.C=pair.dist.plot(c("ABara", "C"),Df)[,4]) # p=0.444
        tab <- cbind(tab, ABarp.C=pair.dist.plot(c("ABarp", "C"),Df)[,4]) # p=0.622
        
        dev.off()        
}

Df <- WideDf[!(WideDf$embryo %in% outliers),]

pdf(paste0(outDir,today,"_cell-cell-distances.pdf"),width =6, height = 4.5, pointsize = 10, useDingbats=F)
        par(mar=c(2,2.3,0.2,0.15), mgp=c(1.25,0.5,0), cex.lab=0.8, cex.axis=0.6, cex.main=0.8, tcl=-.35)
        par(mfrow=c(2,3))
{
        
        tab <- pair.dist.plot(c("ABal", "MS"),Df) #p=0.57
        tab <- cbind(tab, ABala.MS=pair.dist.plot(c("ABala", "MS"),Df)[,4]) # p=0.008 **
        tab <- cbind(tab, ABalp.MS=pair.dist.plot(c("ABalp", "MS"),Df)[,4]) # p=0.23 n.s.
        tab <- cbind(tab, ABala.MSa=pair.dist.plot(c("ABala", "MSa"),Df)[,4]) # p=0.27 n.s.
        tab <- cbind(tab, ABalp.MSa=pair.dist.plot(c("ABalp", "MSa"),Df)[,4]) # p=0.012 *
        
        tab <- cbind(tab, ABar.MS=pair.dist.plot(c("ABar", "MS"),Df)[,4]) # p=0.0019 **
        tab <- cbind(tab, ABala.MS=pair.dist.plot(c("ABara", "MS"),Df)[,4]) # p=0.0026 **
        tab <- cbind(tab, ABara.MSa=pair.dist.plot(c("ABara", "MSa"),Df)[,4]) # p=0.68
        
        tab <- cbind(tab, ABar.C=pair.dist.plot(c("ABar", "C"),Df)[,4]) # p=0.5
        tab <- cbind(tab, ABara.C=pair.dist.plot(c("ABara", "C"),Df)[,4]) # p=0.126
        tab <- cbind(tab, ABarp.C=pair.dist.plot(c("ABarp", "C"),Df)[,4]) # p=0.406

dev.off()        
}

# ggscatterhist(data=data.frame(ABar_MS=ABar.MS, AB.rel=WideDf[ids,"AB_rel"],ID=WideDf[ids,"ID"], Outcome=WideDf[ids,"Outcome"]), 
#               x="AB.rel", y="ABar_MS", group="Outcome", label="ID",
#               margin.plot = "boxplot", margin.plot.size = 0.8, margin.params = list(color="Outcome"),
#               shape=c(4,16)[outcome],palette=c("red","black"), legend="bottom")

