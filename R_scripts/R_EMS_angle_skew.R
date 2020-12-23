# Outlier detection EMS angle ####

#massively misspositioned embryos due to tilted EMS division: GM23, GM24, GM19, PM23, GM6
#these are real outliers with massively messed up cell positioning mostly due to EMS division axis twisted

pdf(paste0(outDir,today,"_EMS_angle_outliers.pdf"),pointsize=10,width=2.1,height=3.5)
{   par(mfcol=c(1,1), mar=c(2.5,2.5,2.1,1.1), mgp=c(1.5,0.5,0), cex.lab=0.9, cex.axis=0.8, cex.main=1, tcl=-.35)
  
  #EMS angle skew ####
  
  r <- Embs$EMS.skew
  angle <- Embs$EMS.angle
  # subset
  s <- !partial #without partial
  
    plot.bxp("",data = Embs[s,],var1 = "EMS.angle", ylab1 = "EMS angular deviation from [°]", saveit = F, ylim=c(0,90))
    abline(h=35, lty=2)
    #points(as.numeric(Embs[r,"Group"]),Embs[r,"EMS.angle"],pch=16, cex=.65);
    
    #frequency in different groups
    ft <- ftable(Embs[s,c("Group","EMS.skew")])
    text(1:5, rep(95), paste(round(prop.table(ft,1)*100)[,2],"%", sep ="" ), cex=0.65, pos=1)
  
  dev.off()
}

Tuk("EMS.angle", "Group", Embs) # no significant difference between groups

#among controls:
#GM31 corrects EMS spindle quickly
#GM54 corrects EMS spindle quickly
#GM14 corrects EMS spindle quickly 
#MM6 corrects EMS spindle, but problem with skewing AB cells left/right 

#PM28 has quite mild phenotype, corrects, but MSa/p is flipped later, kept in the dataset
#GM36 is extremely misspositioned, with ABpl division being skewed and C P3 positions totally off
#PM23 is the most extreme outlier  