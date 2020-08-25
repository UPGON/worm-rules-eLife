library(car)
library(glmnet)

# functions
rep.cv.glmnet <- function(mat, outcome, rep=50, type.m="mse", df.max=5, folds=10){
  lambda <- {}
  lambda1se <- {}
  error <- {}
  predictors <- {}
  res <- {}
  
  for (i in 1:rep) {
    cat("iteration:", i,"/",rep,"\n")
    CV <- cv.glmnet(x=mat,y=outcome,family="binomial",type.measure = type.m, alpha=1, dfmax=df.max, nfolds =folds, trace.it = 0)
    co <- coef(CV, lambda=CV$lambda.1se)
    
    predictors <- cbind(predictors, co[,1]!=0)
    
    lambda[i] = CV$lambda.min
    lambda1se[i] = CV$lambda.1se
    error[i] = CV$cvm[CV$lambda==CV$lambda.1se]
  }
  
  #par(mfrow=c(1,3))
  plot(CV)
  plot(lambda1se,error, main="Lambda vs Deviance", ylab=type.m, xlab="Lambda")
  
  mi <- mean(lambda[which(min(error)==error)])
  oneSE <- mean(lambda1se[which(min(error)==error)]) # take 1SE value for a model with the best fit
  ma <- max(lambda[which(max(error)==error)])
  
  res$min <-mi
  res$mean <- mean(lambda)
  res$max <- ma
  res$oneSE <- oneSE
  res$predictors <- predictors
  
  res
}  
df.max = 10
  
#data up to 15C ####
      
      data <- lassoDf[,] # let LASSO choose best variables
      data[,7:ncol(data)] <- scale(data[,7:ncol(data)])
      
      #use only equalized to train the model
      eqlzd <- data$embryo %in% Embs$ID[equalized]
      
      capture.output(
        cat("all embryos after filtering for LASSO at 15C: \n",paste(data$embryo,collapse = ", ")),
        cat("\n Equalized Embryos used for training: \n",
            "Alive:", sum(data$Group=="alive"),"\n",
            "Dead:", sum(data$Group=="dead"),"\n",
            paste(data$embryo[eqlzd],collapse = ", "))
        
        , file=paste0(outDir, today, "_embryos_in_LASSO.txt")
      )
      
      testData <- data[eqlzd,c(2,7:ncol(data))]
      outcome <- as.numeric(testData[,"Outcome"])-1
      testData$Outcome <- outcome
      
      # 1 = hatched
      # 0 = died
      
      testmat <- as.matrix(testData[,2:ncol(testData)])

# perform repeated cross-validation LASSO to get the best lambda estimation ####
  set.seed(2)
  
  nrepeats <- 500
    
  pdf(paste0(outDir,today,"_cross-validation_LASSO_predictors.pdf"), width=5, height=2.5, pointsize = 10, useDingbats = F)
  {
    par(mfcol=c(1,2), mar=c(3,3.5,2,0.5), mgp=c(1.5,0.5,0),cex.lab=0.8, cex.axis=0.6, cex.main=0.8)
  
    #perform 500 x 5-fold cross-validation on equalized embryos
    lambda <- rep.cv.glmnet(testmat, outcome, nrepeats, folds = 5, type.m="mse", df.max = df.max)
    
    pred <- lambda$predictors
    pred <- pred[-1,]
    freq <- apply(pred, 1, sum)/nrepeats
    barplot2(freq[freq!=0], horiz=T, cex.names = 0.7, las=1, space=0.5, plot.grid=T, grid.lwd = 0.5,xlim=c(0,1),
            main="Predictor inclusion frequency", 
            xlab=paste0("Inclussion freq. over ", nrepeats, " iterations"))
    nvars <- apply(pred, 2, sum)
    freqN <- table(nvars)
    freqP <- freqN/nrepeats*100
    
    hist(nvars, freq=T, main="# of predictive variables in the lasso model", 
         xlab="Number of variables", 
         ylim = c(0,nrepeats),
         breaks =c(-0.5:6.5)
         )
    text(as.numeric(names(freqN)),freqN,paste(round(freqP),"%"), adj=c(-0.1,0.5), cex=0.7, srt=90)
    
  dev.off()
  }
  
  # lambda1se = 0.1354668
  mod1 <- glmnet(x=testmat,y=outcome,family="binomial", alpha=1, lambda = lambda$oneSE)
  
  cc <- names(sort(mod1$beta[mod1$beta[,1]!=0,1])); cat(cc)

  # best model contains 5 same variables reproducibly
  # analysis for equalized embryos
    # "ABara.pDV" "Ca.netdis" "ABara.pOV" "ABala.pOV" "ABarp.aLR" "ABpr.aDV" 
  # for lambda1se

  #with EMS skewed embryos kept in:
  # ABara.pDV Ca.netdis ABara.pOV MS.aMean ABpla.pOV ABarp.aLR
  
#performance on the training data
  p <- predict(mod1, newx=testmat, type="class")
  
  sum(p == outcome)/length(p)*100 
  #96.8% accuracy
  #92 % accurate with EMS
misclass <- rownames(testmat)[!p == outcome]; misclass 
# only "PM16"

Conf(as.numeric(p[,1]>0.5),outcome, pos="0")
# is this some kind of overfiting? too good to be true?
# Kappa = 0.937
# Accuracy : 0.97
#Sensitivity : 1
#Specificity : 0.93

#           Reference
# Prediction  0  1
          # 0 17  1
          # 1  0 14

    #test on unequal controls (wildtype + lin-5 controls) ####
      p <- predict(mod1, newx=as.matrix(data[data$Experiment!="meta",-c(1:6)]), type="class")
      sum(p==1)/length(p)
      # 3/23 alive controls labeled as dead
      # 87% accuracy
      
    # all data
      prediction <- predict(mod1, newx=as.matrix(data[,-c(1:6)]), type="class")
      outcome.all <- as.numeric(data[,"Outcome"])-1
      sum(prediction==outcome.all)/length(prediction)
    #93.2% success on all data including training and test dataset + outlier embryos
    #100% sensitive, 89.5% specific
    
    #85% accurate with EMS kept in  
      
    pdf(paste0(outDir,today,"_upto15Cells_LASSO_ROC_all_data.pdf"), width=3, height=3, pointsize = 10)
      r <- roc(predict(mod1, newx=as.matrix(data[,-c(1:6)]), type="response")[,1], 
               response=outcome.all, plot=T, smooth=F, auc=T, ci=T, levels = c(1,0))
      mtext(paste("AUC: ",round(pROC::auc(r)*100, digits=1),"%"),cex=0.8, line=-2)
      # AUC 98
    dev.off()
      
    Conf(as.numeric(prediction[,1]), outcome.all, pos="0")
    # Reference
    # Prediction  0  1
    #          0 21  4
    #          1  0 34
    
    #with EMS kept in
    #Reference
    #Prediction  0  1
    #         0 27  9
    #         1  1 30
      