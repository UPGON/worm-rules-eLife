
# how well can we predict MSa/MSp flip? ####
# select variables without filtering

select <- grep(paste0("^(",paste0(c(eightC),collapse ="|"),")\\."),colnames(WideDf))
confounding <- grep("^MS\\.[a][MADL].*",colnames(WideDf)) #we don't want angle variables for MS intself
colnames(WideDf)[confounding]
select <- select[!select %in% confounding] 

na <- is.na(Embs$MSap.flip)


outl <- c("GM6","GM19","GM23","GM24","GM36","PM23")
data <- getSubset(WideDf[!na,],1:nrow(WideDf[!na,]),select) # select all embryos
ids <- rownames(data)[data$Experiment=="meta" & !(data$embryo %in% outl)]


testmat <- as.matrix(data[ids,c(7:ncol(data))])
flip <- Embs[rownames(testmat),"MSap.flip"]

# 1 = T = flipped
# 0 = F = normal

lambdaFlip <- rep.cv.glmnet(testmat, flip, 100, "deviance")

fit <- glmnet(x=testmat,y=flip,family="binomial", alpha=1, dfmax=df.max, lambda = lambdaFlip$min)
cc <- names(sort(fit$beta[fit$beta[,1]!=0,1])); cat(cc)
# best model with lambda min
# MS.pAP P3.pLR C.EndTime ABar.pDV
####      \\\\\  \\\\\\\\

ms.outcome <- as.numeric(flip)-1

testdf <-as.data.frame(cbind(ms.outcome,testmat[,c("MS.pAP","P3.pLR","C.EndTime","ABar.pDV")]))
MSmodel <-  glm(ms.outcome~.,data = testdf, family="binomial")
ms <- predict.glm(MSmodel, type="response") > 0.5
as.numeric(ms) == as.numeric(ms.outcome)

anova(MSmodel)


#performance on the training data
p <- predict(fit, newx=testmat, type="class")
sum(p == flip)/length(p)*100 
# Accuracy 75%

misclass <- rownames(testmat)[!p == flip]; misclass 
length(misclass)
length(flip)
# 13/52 missclassified

r <- roc(predict(fit, newx=testmat, type="response")[,1], response=flip, plot=T, smooth=F, auc=T)
mtext(paste("AUC: ",round(pROC::auc(r)*100, digits=1),"%"),cex=0.8, line=1)
# AUC 85.7

Conf(as.factor(p[,1]),flip)
#it is really quite bad model
# without outliers Kappa 0.44, but 43% sensitivity, yet 97% specificity

#test on unequal controls
p <- predict(fit, newx=as.matrix(data[data$Experiment!="meta",-c(1:6)]), type="class")
1-sum(p=="T")/length(p)
# the model without outliers is able to labels all controls as unflipped, which is really cool

# is there a difference in compression between flipped/not flipped?
# boxplot(height~MSap.flip, data=Embs[ids,])
# beeswarm(height~MSap.flip, data=Embs[ids,], add=T)
# t.test(height~MSap.flip, data=Embs[ids,])
# definitely no difference, compression has no effect on this phenotype

