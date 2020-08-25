#comparisons variance in timing ####
#Cell cycle duration
timeR <- results$LifeTime
timeR <- timeR[3:nrow(timeR),] #start from the 4-cell stage

vc <- timeR[,grep("vc", colnames(timeR))]
vc$cell=rownames(vc)
apply(vc, 2, mean, na.rm=T)*100
apply(vc, 2, sd, na.rm=T)*100


sum(!is.na(vc$C.vc))

pair=T

p <- t.test(vc$WT.vc,vc$C.vc, paired = pair)$p.value #p = 3.8e-5
p <- c(p,t.test(vc$C.vc,vc$EA.vc, paired = pair)$p.value) # n.s. 0.095
p <- c(p, t.test(vc$ED.vc,vc$C.vc, paired = pair)$p.value)
p <- c(p,t.test(vc$ED.vc,vc$EA.vc, paired = pair)$p.value)
round(p.adjust(p,method = "BH"), digits = 4)

# difference between AB and P1 cells in equalized embryos?
ABc <- grep("AB",rownames(timeR))
t.test(vc$EQ.vc[ABc],vc$EQ.vc[-ABc])
t.test(vc$EA.vc[ABc],vc$EA.vc[-ABc])
t.test(vc$ED.vc[ABc],vc$ED.vc[-ABc])
t.test(vc$C.vc[ABc],vc$C.vc[-ABc])

#Variance in POSITIONS along AP ####
APposR <- results$pAP
LRposR <- results$pLR
DVposR <- results$pDV
OVposR <- results$pOV


APposR <- APposR[7:nrow(APposR),]
LRposR <- LRposR[7:nrow(APposR),]
DVposR <- DVposR[7:nrow(APposR),]
OVposR <- OVposR[7:nrow(APposR),]


vcAP <- APposR[,grep("sd", colnames(APposR))]
vcLR <- LRposR[,grep("sd", colnames(APposR))]
vcDV <- DVposR[,grep("sd", colnames(APposR))]
vcOV <- OVposR[,grep("sd", colnames(OVposR))]

apply(vcAP, 2, mean, na.rm=T)
apply(vcLR, 2, mean, na.rm=T)
apply(vcDV, 2, mean, na.rm=T)
apply(vcOV, 2, mean, na.rm=T)

pair=T

p <- t.test(vcAP$C.sd,vcAP$EA.sd, paired = F)
p <- t.test(vcAP$EA.sd,vcAP$ED.sd, paired = pair)
p <- t.test(vcAP$EA.sd,vcAP$INV.sd, paired = pair)

t.test(vcLR$C.vc,vcLR$EA.vc, paired = pair) #12.5 vs 20.4, 2.9 e-16
t.test(vcLR$EA.vc,vcLR$ED.vc, paired = pair) #20.4 vs 24.9,  0.0001
t.test(vcLR$EA.vc,vcLR$INV.vc, paired = pair) #20.4 vs 19.8, 0.66

t.test(vcDV$C.vc,vcDV$EQ.vc, paired = pair) # 13.7 vs 26%,p=2.2e-16
t.test(vcDV$EA.vc,vcDV$ED.vc, paired = pair) # 25.4 vs 26.8%, p=0.43
t.test(vcDV$EA.vc,vcDV$INV.vc, paired = pair) # 25.4  vs 21%, p=0.055

t.test(vcOV$C.vc,vcOV$WT.vc, paired = pair) # #no difference
t.test(vcOV$C.vc,vcOV$EQ.vc, paired = pair) # more variable in EQ, p-value = 6.008e-13 mean difference 12.6%
t.test(vcOV$EA.vc,vcOV$ED.vc, paired = pair) # 25.4  vs 21%, p=0.055


t.test(vcOV$C.vc[ABc],vcOV$C.vc[-ABc]) # no difference
t.test(vcOV$EQ.vc[ABc],vcOV$EQ.vc[-ABc]) # 6.3% AB, vs 5.5% for P1, significant 

ABcl <- grepl("AB",rownames(vcOV))

boxplot(vcOV$EQ.vc~as.factor(ifelse(ABcl,"AB","P1")))
t.test(vcOV$EQ.vc~as.factor(ifelse(ABcl,"AB","P1")))

boxplot(vcOV$WT.vc~as.factor(ifelse(ABcl,"AB","P1")),ylim=c(0,1.2))
boxplot(vcOV$C.vc~as.factor(ifelse(ABcl,"AB","P1")), ylim=c(0,1.2))
boxplot(vcOV$EA.vc~as.factor(ifelse(ABcl,"AB","P1")),ylim=c(0,1.2))
boxplot(vcOV$ED.vc~as.factor(ifelse(ABcl,"AB","P1")),ylim=c(0,1.2))


t.test(vcOV$ED.vc~as.factor(ifelse(ABcl,"AB","P1"))) #no difference for C, WT, surprisingly also for ED
t.test(vcOV$EA.vc~as.factor(ifelse(ABcl,"AB","P1"))) #signif for EA 0.00015

t.test(vcOV$C.vc,vcOV$EQ.vc) # more variable in EQ, p-value = 6.008e-13 mean difference 12.6%
