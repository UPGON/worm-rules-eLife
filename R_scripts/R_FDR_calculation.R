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

#prepare data ####
Df <- TempDf

q <- "(AB[aplr]{1,3}|[EMS].{0,1}|(P[2-3])|[C].{0,1})\\..*" #up to 15C
cols <- grep(q,colnames(Df))
rows <- Embs$ID[equalized&!(Embs$ID%in%outliers)] #remove EMS outliers

Df <- FilterAndSubset(Df, rows, cols, minN=5) 
# volume related variables don't have enough data points so they will get removed
# minN indicates minimal number of observations for each variable in every group

Df <- Df[,c(4,7:ncol(Df))]

#df <- prep.data(WideDf,cells=c(fourC,eightC,fifteenC),outlierRM=T,inpute = F)
cat("filtered table contains: ", ncol(Df)-1, " variables \n")

g <- factor(Df$Group)
test.res <- apply(Df[2:ncol(Df)], 2, function(x)t.test(x~g)$p.value)

#to calculate FDR use following process: 
  #shuffle group labels for all variables
  # calculate t.test on scrambled data
  # how many p. values are smaller than x in randomized data? Repeat R times
  # compare with number of positive tests without randomizing
# FDR is number of false positives / all positive results

R <- 100
#Calculate t.test on scrambled data ####
set.seed(3)
sim.res <- matrix(nrow=length(test.res), ncol=R)

  cat("Calculating t.test on randomized data \n")
  pb <- txtProgressBar(min=0,max=R, style=3, label="")
  
  for (i in 1:R) {
    setTxtProgressBar(pb,i)
    scr <- sample(g) # scramble groups
    sim.t <- apply(Df[2:ncol(Df)], 2, function(x)t.test(x~scr)$p.value)
    sim.res[,i] <- sim.t
  }
  close(pb)

# plot experimental and randomized distribution
hist(sim.res, breaks = 20, freq = F, ylim = c(0,3.5))
hist(test.res, breaks = 20, freq = F, add=T, col =addTrans("cadetblue3",80))

#FDR distribution at different alpha levels ####
alpha.range <- seq(0.001, 0.05, by=0.001)

#initiate empty matrix
FDR.m <- matrix(nrow = length(alpha.range), ncol=5, dimnames = list(alpha=alpha.range,vars=c("alpha","n.positives","n.fp","FDR", "fp.fraction")))

for (a in 1:length(alpha.range)) {
  alpha <- alpha.range[a]
  
  n.p <- sum(test.res < alpha) # number of positive results on real data
  false.p <- apply(sim.res, 2, function(x)sum(x<alpha))
  n.fp <- mean(false.p, na.rm=T) # get mean number of false positives in R simulations with p < alpha

  fp.fraction <- n.fp/length(test.res) # given number of test we did
  
  # FDR is number of false positives / all positive results
  FDR <- n.fp/n.p
  FDR.m[a, ] <- c(alpha, n.p, n.fp, FDR, fp.fraction)
}

good <- alpha.range[FDR.m[,"FDR"]<0.1] #which alpha levels result in FDR below 10%
alpha.lim <- max(good)
cat("\n Suitable alpha threshold given the data for FDR < 10% is: ", alpha.lim)

#significant vars given the 10% FDR
s.vars <- names(test.res)[test.res < alpha.lim]
cat("\n There is",length(s.vars),"significant variable at this alpha level are: \n", s.vars)

#fdr tool

q <- p.adjust(test.res, method = "BH")
bh <- names(q)[(q<0.05)]
cat("\n\n There are", length(bh)," significant variable with BH correction:\n", bh)
