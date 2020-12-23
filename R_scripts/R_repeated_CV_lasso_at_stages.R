library(car)
library(glmnet)

nrepeats <- 250

# functions
rep.cv.glmnet <- function(mat, outcome, rep=50, type.m="mse", df.max=5, folds=5){
  lambda <- {}
  lambda1se <- {}
  error <- {}
  predictors <- {}
  res <- {}
  conf <- {}
  
  for (i in 1:rep) {
    cat("iteration:", i,"/",rep,"\n")
    CV <- cv.glmnet(x=mat,y=outcome,family="binomial",type.measure = type.m, alpha=1, dfmax=df.max, nfolds =folds, trace.it = 0)
    co <- coef(CV, lambda=CV$lambda.1se)
    
    predictors <- cbind(predictors, co[,1]!=0)
    
    #performance of best model with lambda-1SE from CV on the training data
    mod1SE <- glmnet(x=mat,y=outcome,family="binomial", alpha=1, lambda =CV$lambda.1se)
    
      p <- predict(mod1SE, newx=mat, type="class")
      performance <- Conf(as.numeric(p[,1]>0.5),outcome, pos="0")
      conf <- c(conf, performance$acc)
    
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
  
  points(oneSE,min(error),col="red")
  
  
  res$min <-mi
  res$mean <- mean(lambda)
  res$max <- ma
  res$oneSE <- oneSE
  res$predictors <- predictors
  res$accuracy <- conf
    
  res
}
df.max = 10
  
# data prep

      write.xlsx(lassoDf, file=paste0(outDir,today,"_Table_S7_values.xlsx"))
      lassoDf[,7:ncol(lassoDf)] <- scale(lassoDf[,7:ncol(lassoDf)])
      
      #588 variables, 57 embryos, 236 NA inputed
      stages <- list(four=fourC, eight=eightC, fifteen=fifteenC, twentyeight=twentyeightC)
      acc <- {}
      
pdf(paste0(outDir,today,"_LASSO_CV_report_at_stages.pdf"), width=5, height=8, pointsize = 11, useDingbats = F)
{ par(mfcol=c(6,4), mar=c(2.5,2.5,2,0.5), mgp=c(1.5,0.5,0),cex.lab=0.8, cex.axis=0.6, cex.main=0.8)
  
  wb <- createWorkbook("RJ") #for data recording
  
  for (i in 1:length(stages)){
        stage = names(stages)[i]
        cells <- stages[[i]]
  
        cols <- grep(paste("^(",paste0(cells,collapse ="|"),")\\..*", collapse=""),colnames(lassoDf))
        #colnames(lassoDf)[cols]
        length(cols)
        
        data <- cbind(lassoDf[,1:6],lassoDf[,cols]) # let LASSO choose best variables
        
        #use only equalized to train the model
        eqlzd <- data$embryo %in% Embs$ID[equalized]
        
        capture.output(
          cat("------------------------------------------\n\n
              all embryos after filtering for LASSO at", stage ,"-cell stage: \n",paste(data$embryo,collapse = ", ")),
          cat("\n\n Equalized Embryos used for training: \n",
              "Alive:", sum(data$Group=="alive"),"\n",
              "Dead:", sum(data$Group=="dead"),"\n",
              paste(data$embryo[eqlzd],collapse = ",")
              , "\n\n\n")
          
          , file=paste0(outDir, today, "_LASSO_embryos_included_at_stages.txt"), append = T
        )
        
        addWorksheet(wb, sheetName=stage)
        
        testData <- data[eqlzd,c(2,7:ncol(data))]
        outcome <- as.numeric(testData[,"Outcome"])-1
        testData$Outcome <- outcome
        
        # 1 = hatched
        # 0 = died
        
        testmat <- as.matrix(testData[,2:ncol(testData)])
  
  # perform repeated cross-validation LASSO to get the best lambda estimation ####
    set.seed(2)
      
      #perform 500 x 5-fold cross-validation on equalized embryos
      lambda <- rep.cv.glmnet(testmat, outcome, nrepeats, folds = 5, type.m="mse", df.max = df.max)
      
      pred <- lambda$predictors
      pred <- pred[-1,]
      freq <- apply(pred, 1, sum)/nrepeats
      
      if(all(freq==0)) freq=as.matrix(data.frame(Empty_model=1))
      #par(mar=c(6.1, 2.5, 2.1, 0.5))
      barplot2(freq[freq!=0], horiz=F, cex.names = 0.7, cex.lab=0.8, cex.axis = 0.8, las=2, plot.grid=T, grid.lwd = 0.5,ylim=c(0,1),
              main=paste0("Variable inclusion - ",stage), 
              #xlab=paste0("Inclussion freq. over ", nrepeats, " iterations")
              )
      
      nvars <- apply(pred, 2, sum)
      freqN <- table(nvars)
      freqP <- freqN/nrepeats*100
      
      
      hist(nvars, freq=T, main="# of predictors", 
           xlab="Number of variables", 
           ylim = c(0,max(freqN, na.rm=T)*1.2),
           breaks =c(-1:max(nvars, na.rm=T)+.5)
           )
      
      text(as.numeric(names(freqN)),freqN,paste(round(freqP),"%"), adj=c(-0.1,0.5), cex=0.7, srt=90)

    boxplot(lambda$accuracy, ylim=c(0.5,1), ylab="Accuracy over CV iterations")
      
    acc <- cbind(acc,lambda$accuracy)
      
    mod1 <- glmnet(x=testmat,y=outcome,family="binomial", alpha=1, lambda = lambda$oneSE)
    
    writeData(wb, stage, data.frame(variable=colnames(data[,-(1:6)]),inclusion.frequency=freq, row.names=colnames(data[,-(1:6)]), oneSE.model.coef=mod1$beta[,1]))
    
    cc <- names(sort(mod1$beta[mod1$beta[,1]!=0,1])); cat(cc)
    text(length(cc),max(freqN)*1.1,"*", adj=c(0,3), cex=1, srt=90)
    mtext(paste("BestModel:", paste(cc,collapse =" + ")), cex=0.8)
      
  #performance on the training data
    p <- predict(mod1, newx=testmat, type="class")
    misclass <- rownames(testmat)[!p == outcome]
    
    # all data
    predict.all <- predict(mod1, newx=as.matrix(data[,-c(1:6)]), type="class")
    outcome.all <- as.numeric(data[,"Outcome"])-1
  
  capture.output(
    cat("\n----------------------------", stage,"------------------\n\n"),
    cat("Model performance on TRAINING data at", stage ,"-cell stage: \n\n",
        "Average accuracy on training data: ", 
        round(mean(lambda$accuracy)*100,2),"% ± ",round(sd(lambda$accuracy)*100,2)),
    
    cat("\n\n misclassified embryos:", paste(misclass,collapse =","),"\n\n"),
      Conf(as.numeric(p[,1]>0.5),outcome, pos="0"),
    
    cat("\n\n Performance on all embryos:\n\n"),
      Conf(as.numeric(predict.all[,1]), outcome.all, pos="0"),
    cat("\n\n"),
    file=paste0(outDir, today, "_LASSO_performance_at_stages.txt"),
    append = T
  )
        
        r <- roc(predict(mod1, newx=as.matrix(data[,-c(1:6)]), type="response")[,1], 
                 response=outcome.all, plot=T, smooth=F, auc=T, ci=T, levels = c(1,0))
        mtext(paste("AUC: ",round(pROC::auc(r)*100, digits=1),"%"),cex=0.8, line=-2)
  } # end for loop over stages
  saveWorkbook(wb, paste0(outDir,today,"_Table_S7.xlsx"), overwrite = TRUE)
  dev.off()
} #end pdf

acc <- as.matrix(acc)*100
colnames(acc) <- names(stages)
pdf(paste0(outDir,today,"LASSO_accuracy_at_stages.pdf"),width =2, height = 2.5, pointsize = 10, useDingbats=F)
  par(mar=c(2,2.3,0.2,0.15), mgp=c(1.25,0.5,0), cex.lab=0.7, cex.axis=0.6, cex.main=0.8, tcl=-.35)
  boxplot(acc, ylim=c(50,100), ylab="Accuracy over CV iterations")
  
dev.off()
      