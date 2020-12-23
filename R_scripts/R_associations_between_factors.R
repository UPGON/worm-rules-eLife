cont <- table(Embs[Embs$Group %in% c("alive","dead", "inverted"),c("MSap.flip", "Extra.P4")])
fisher.test(cont) #no difference / association between groups p=0.2152

cont <- table(Embs[Embs$Group %in% c("alive","dead", "inverted"),c("MSap.flip", "PhaDefect")])
fisher.test(cont) #no difference / association between groups p=1

cont <- table(Embs[Embs$Group %in% c("alive","dead"),c("Extra.P4", "PhaDefect")])
fisher.test(cont) #n.s. / but looks like there could be an effect p=0.052
#when inverted embryos are included it is significantly associated p=0.0223

#EMS angle deviation over 35º
cont <- table(Embs[Embs$Group %in% c("alive","dead","inverted"),c("MSap.flip", "EMS.skew")]) #p=0.1974
fisher.test(cont) #no difference / association between groups

cont <- table(Embs[Embs$Group %in% c("alive","dead","inverted"),c("Extra.P4", "EMS.skew")]) #p=0.223
fisher.test(cont) #no difference / association between groups

cont <- table(Embs[Embs$Group %in% c("alive","dead","inverted"),c("PhaDefect", "EMS.skew")]) #p=0.1105
cont <- table(Embs[Embs$Group %in% c("alive","dead"),c("PhaDefect", "EMS.skew")]) #p=0.294
fisher.test(cont) #no difference / association between groups

p.adjust(c(0.2152, 1, .0223, .1974, .2233, .1105), method = "BH") 
# not significant anymore

#counts for all different combinations
cont <- table(Embs[Embs$Group %in% c("alive","dead", "inverted"),c("MSap.flip", "Extra.P4", "PhaDefect")])
structable(cont)
mosaic(cont, shade=T)
mosaic(cont, shade=T, type="expected")

