Embs$MSap.flip <- factor(Embs$MSap.flip)
Embs$Extra.P4 <- factor(Embs$Extra.P4)

tab <- Embs[!partial,]

{#generate plots    
  #wihtout partial
  cont <- table(tab$MSap.flip,tab$Group)
  #    wt ctrl alive dead inverted
  # F  0    0    15   10    2
  # T  0    0     4   18    5
  
  # for equalized only: p=0.0067, odds ratio 6.5
  
  #Is the outcome associated with the MSa/p flip phenotype?
  ft <- fisher.test(cont[,3:4]) 
  #definitely MSap flip is more commmon in dead embryos
  
  #Is flip linked to compression?
  meta <- Embs$Experiment=="meta"
  t.test(Embs$height[meta]~Embs$MSap.flip[meta]) #no difference
  t.test(Embs$height[equalized]~Embs$MSap.flip[equalized]) #no difference
  
  #Figure 3E and Figure 4K
  pdf(paste0(outDir,today,"_MSap_and_P4_barplot.pdf"),pointsize=10,width=3,height=3)
  par(mfcol=c(1,2), mar=c(2.5,2.5,2.1,1.1), mgp=c(1.5,0.5,0), cex.lab=0.9, cex.axis=0.8, cex.main=1, tcl=-.35)
  
  barplot((1-prop.table(cont, 2))*100, 
          main ="MSa/p flip", 
          ylab="% embryos with MSa/p flipped", 
          width = 0.8, 
          space=0.2, 
          col=c("#b41a6d","#eeeeee"),
          ylim=c(0,106))
  text(x=0.5:4.5, y=rep(97, 5), labels = colSums(cont), cex=0.7)
  lines(x=c(2.4,2.4,3.6,3.6),y=c(100,101.5,101.5,100), lwd=1.5)
  text(x=3, y=104, labels = paste0("p = ", signif(ft$p.value, digits=1)), cex=0.7)
  
  # Extra divisions in the P4 lineage ####
  #Is the outcome associated with the extra divisions in the P4 lineage?
  cont <- table(tab$Extra.P4,tab$Group)
  
  #   wt ctrl alive dead  inv
  # F 10   11    10   14    0
  # T  0    0     4    8    4
  
  #p = 0.328
  #not significant, but definetely accurs only in equalized embryos
  barplot((1-prop.table(cont, 2))*100, 
          cex.names=0.8,
          main ="Extra P4 divisions", 
          ylab="% embryos with >1 ectopic P4 division", width = 0.8, space=0.2, 
          col=c("#2E5268","#eeeeee"), 
          ylim=c(0,106))
  text(x=0.5:4.5, y=rep(97, 5), labels =colSums(cont), cex=0.7)
  dev.off()    
}

cont <- xtabs(~., data=tab[,c("Extra.P4","Outcome","MSap.flip")], subset = equalized)
ftab <- ftable(cont)

#                  MSap.flip F T
# Extra.P4 Outcome              
# F        died              6 7
#          hatched           7 2
# T        died              2 5
#          hatched           3 1

strucplot(cont, shade=T, type="observed")
fisher.test(ftab) # MSa/p flip and Extra P4 divisions are significantly associated

fisher.test(ftab[c(1,3),1:2])
# so it is not associated in equalized embryos strictly speaking

cont <- xtabs(~., data=Embs[,c("Extra.P4","Outcome","MSap.flip")], subset=!controls)
ftab <- ftable(cont)
#                   MSap.flip F T
# Extra.P4 Outcome              
# F         died              7 8
#           hatched           9 3
# T         died              3 9
#           hatched           3 1

fisher.test(ftab) # for all upshifted embryo there seems to be a trend towards association
#p-value 0.076

