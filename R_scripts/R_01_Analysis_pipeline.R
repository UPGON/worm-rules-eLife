#analysis performed in R version 3.6.1

#load libraries ####
  if(!require("plyr")) install.packages("plyr")
  if(!require("reshape2")) install.packages("reshape2")
  
  if(!require("beeswarm")) install.packages("beeswarm")
  if(!require("gridExtra")) install.packages("gridExtra")
  if(!require("ggsignif")) install.packages("ggsignif")
  
  if(!require("readxl")) install.packages("readxl")
  if(!require("openxlsx")) install.packages("openxlsx")
  if(!require("data.table")) install.packages("data.table")
  
  if(!require("gtools")) install.packages("gtools")
  if(!require("gplots")) install.packages("gplots")
  if(!require("pheatmap")) install.packages("pheatmap")
  if(!require("ggpubr")) install.packages("ggpubr")
  
  if(!require("bestglm")) install.packages("bestglm")
  if(!require("glmnet")) install.packages("glmnet")
  
  if(!require("pROC")) install.packages("pROC")
  if(!require("corrplot")) install.packages("corrplot") # for rquery.cormat
  
  if(!require("mplot")) install.packages("mplot")
  
  if(!require("vcd")) install.packages("vcd") #plot categorical variables
  
  if(!require("calibrate")) install.packages("calibrate") 
    # for the textxy function for labeling graph
    # conversion of radians to degrees rad2degree
  if(!require("wordcloud")) install.packages("wordcloud") 
  
  if(!require("DescTools")) install.packages(DescTools) # model performance, Conf function
  # useful to treat outliers by replacing them with 0.05 and 0.95 percentile value
  # also contains function "Conf" for confusion matrix calculation and associated statistics
  
  if(!require(multcompView)){install.packages("multcompView")} # for multicomparison letters

#### set environment ####
par(mar=c(3.1,3.1,2.1,1.1))
today <- format(Sys.Date(), "%y%m%d")

basewd <- "/LocalDocs/switchdrive/ev571_embryos/PROJECTS/worm-rules-eLife/"
outDir <- paste0(basewd,"Results/",today,"/")
if(!dir.exists(outDir)){dir.create(outDir);print("Output directory created");}else print("directory already exists")
setwd(basewd)

#load libraries pre-processed data: ####
source("R_scripts/R_functions.R")

#load data that has been previously aligned and exported in proper format
load(paste0(basewd,"/200723_embryos.RData"))
Embs$AR <- Embs$length/Embs$width

    Complete <- as.data.frame(Complete)
    
    #recalculate radians to degrees
    a.cols <- grep("a[ALDM]",colnames(Complete))
    Complete[,a.cols] <- RadToDeg(Complete[,a.cols]) #from calibrate plugin
    
    #record EMS skew in Embs table
    ems <- Complete[Complete$Cell=="EMS", c("embryo","aMean")]
    names(ems) 
    Embs$EMS.angle <- 0
    Embs[ems$embryo,"EMS.angle"] <- ems$aMean
    Embs$EMS.skew <-  Embs$EMS.angle > 35
    
    #color pallete
    pal<-c("#2E5268", "#64A8A1", "#D19812", "#C93230", "azure4")

#Distribution of AB sizes of embryos in the analysis ####
source(file="R_Scripts/R_plot_lineaged_embs.R", local = T) #contains code for multiple comparisons with compact letter display
  
#Define Groups ####
alive <- Embs$Outcome =="hatched"
inverted <- Embs$AB_rel<0.48 
  sum(inverted) #7
equalized <- Embs$AB_rel<=0.53 & Embs$AB_rel>=0.48 #49
  sum(equalized&alive) #21
  sum(equalized&!alive) #28
  
partial <- Embs$Experiment=="meta"&!(inverted|equalized) #6
controls <- Embs$Experiment %in% c("ctrl","wt")
controls <- Embs$Experiment %in% c("ctrl")

Embs$Group <- as.character(Embs$Group)
Embs$Group[inverted] <- "inverted"
Embs$Group <- factor(Embs$Group, levels=c("wt","ctrl","alive","dead", "inverted"))

# Annotate embryos with Pharynx defects - from Figure 5 - supplement E
phaDefects <- c("PM03","PM04","PM08","PM18","PM21","PM23","PM28")
phaDefects2 <- c("PM08","PM18",	"PM21",	"PM23", "PM28") # with ABala transformation

Embs$PhaDefect <- ifelse(Embs$ID %in% phaDefects, T, F)
Embs$PhaDefect[!grepl("P",Embs$ID)] <- NA 

# Size of AB in equalized alive and dead
  t.test(AB_rel~Group,Embs[equalized,]) #alive and dead equalized have the same asymmetry
  t.test(AR~Group,Embs[equalized,]) #there is a signif., but small diff in AR, alive more elongated, and longer 
  t.test(Size~Group,Embs[equalized,]) #Dead gyus were larger, p < 0.0059 - it might be an artefact of higher compression
    with(Embs,cor.test(height,Size)) #r= -0.57, p=5.75e-09, very stronly correlated
    with(Embs,cor.test(Size,AR)) # r=-0.34, p=0.0012, also aspect ratio is correlated with size
    with(Embs,cor.test(height,AB_rel)) #not correlated, compression does not affect asymmetry

# length of lin-5 embryos and WT is not different #### 
  t.test(length~Lin_5,Embs) #no difference in legth in our data
  aggregate(Embs$length,by=list(Embs$Lin_5),sd)  
    #comparison of aspect ration (AR) with published data   
    # t.test(AR~Lin_5,Embs) #there is a signif difference between N2 and our LIN-5 embryos, but small sample size - controls more elongated
    # ts <- Embs$AR[Embs$Lin_5=="ev571"] # aspect ratio of all lineaged lin-5(embryos)
    # yamamoto <- rnorm(length(ts), 1.6, 0.1) #simulated data according to mean and SD for AR from Yamamoto and Kimura, 2017
    # boxplot(list(yamamoto, ts))
    # t.test(yamamoto, ts) #lin-5(ev571) are more elogated than the N2s from their paper
    # t.test(ts, mu=1.6) # it is different from reported AR
    
# Figure 4E More compressed embryos are more likely to die ####
  t.test(height~Group,Embs[equalized,])
  #Equalzed died p=0.0016, dead 19.78 um, hatched 21.92 um
  
  capture.output(
    with(Embs, {
      cat("Inverted\t",mean(height[inverted]),sd(height[inverted]),sum(inverted),"\n"); #20.75, sd=3.3 n=7
      cat("Eq. alive\t", mean(height[equalized&alive]),sd(height[equalized&alive]),sum(equalized&alive),"\n"); #21.7, sd=1.96, n=19
      cat("Eq. dead\t", mean(height[equalized&!alive]),sd(height[equalized&!alive]),sum(equalized&!alive),"\n"); #19.75, sd=2.09 n=26
      cat("Partial\t",mean(height[partial]),sd(height[partial]),sum(partial),"\n"); #21.5, sd=2.7 n=10
      #controls:
      cat("Ctrl\t", mean(height[Experiment=="ctrl"]),sd(height[Experiment=="ctrl"]),length(height[Experiment=="ctrl"]),"\n"); #21.3 um,  sd=2.6 n=19
      cat("Wt\t",mean(height[Experiment=="wt"]),sd(height[Experiment=="wt"]),length(height[Experiment=="wt"]),"\n") #21.6 um,  sd=1.3 n=10
    }),
  file=paste0(outDir, today, "_Compression_in_groups.txt")
  )
  
# Figure 3A-B: Plot division pace in AB and P1 ####
source("R_Scripts/R_DivisionPaceABvsP1.R", local = T)

#Figure 4H-K MSa / MSp flip and P4 extra divisions####
#3E and Figure 4K barplots
source("R_scripts/R_MSa.p_flip_extra.P4.R")

#(for a reply to reviewers)
source("R_scripts/R_associations_between_factors.R")
  # no significant association between MSa/MSp flip / EMS.skew / P4 extra divisions / Pha-4 expression defects
  
# Table S2 - Compute all statistics for groups ####
source("R_scripts/R_All_stats_export_XLSX_from_long_data.R")
    #save results also as RData
    save(results, file=paste0(outDir,today,"_all_Vars.statistics.RData"))

# Figure 3F: Division sequence ####
source("R_scripts/R_division_sequence.R")

# Figure 3 supplement + Table S3: Correlation of CC duration with AB size % ####
source("R_scripts/R_AB-size_CC_corr.R")
  
# Generate Wide Data for PCA ####
## DataFrame WideDf - contains aligned timing plus all big table data in an expanded format with Cell.variable format 

    df <- Complete
    df$Group <- as.character(df$Group)
    df$Group[df$embryo %in% Embs$ID[inverted]] <- "inverted"
    df <- df[!(df$embryo %in% Embs$ID[partial]),] #remove partial
    df$Group <- factor(df$Group, levels=c("wt","ctrl","alive","dead","inverted"))
    
    df <- df[df$Divided==1,] # remove undivided cells
    df$Divided <- NULL #remove Divided column, not usefull anymore
    
    WideDf <- toWideTable(df, 2:7, 8:ncol(df), "EMS") #df, anot.cols, var.cols, startCell

 # Figure 4F: EMS skew  ####
  source("R_scripts/R_EMS_angle_skew.R")

  #outliers with skewed EMS
  outliers <- c("GM6","GM19","GM23","GM24","GM36","GM39","GM43","MM6","PM23", "PM28","EM7","EM10") # all have tilted EMS division, 
  # EM10 seems to be poorly aligned, far from all others in PCA
  #outliers <- c("EM10") #"PM23", "GM19" seem to be also quite different from all others
 
# Figure 4A-B: Positional divergence over time ####
  source("R_scripts/R_positional_deviation.R")
  source("R_scripts/R_angular_deviation.R")
  
# stages
  #4C - ABp ABa EMS P2 
  #8C - ABar ABal ABpl ABpr MS E C P3
  #15C - ABalp ABpla ABpra ABprp ABarp ABplp ABala ABara MSa MSp Ca Cp Ea Ep P3
  #24C - ABalaaa ABalaap ABalppa ABalppp ABaraaa ABaraap ABarppa ABarppp 
  #      ABplaaa ABplaap ABplppa ABplppp ABpraaa ABpraap ABprppa ABprppp 
  #      Ca Cp Ea Ep MSa MSp D P4
    
  #28C MS anc C cells divide
  
  fourC <- Cellorder[3:6]
  eightC <- c("ABal", "ABar", "ABpl", "ABpr", "MS", "C", "E", "P3")
  fifteenC <- c(paste0(rep(eightC, each=2)[1:14],c("a","p")),"P3")
  twentyfourC <- c(paste0(rep(fifteenC, each=2)[1:16],c("a","p")),fifteenC[9:14],"D","P4")
  twentysixC <- c(paste0(rep(fifteenC, each=2)[1:20],c("a","p")),twentyfourC[19:24])
  twentyeightC <- c(paste0(rep(fifteenC, each=2)[1:24],c("a","p")),twentyfourC[21:24])
  fiftysixC <- c(paste0(rep(twentyeightC, each=2)[1:54],c("a","p")),"Z2","Z3")
  
  
# Figure 4C-D: Variation in groups over time ####
  cells=Cellorder[3:202] #remove AB and P1 - not lineaged - lineaging typically started at 4C
  
  source("R_scripts/R_variance_plots.R") 

# Figure 3G - Variation AB vs P1 ####
  
  source("R_scripts/R_variation_AB_P1_side-by-side.R")
  #source("R_scripts/R_variance_t.tests.R")
  
# Figure 6 - supplementary 2 - Plot 8C and 16C positions - comparison between equalized alive/dead####
  source("R_scripts/R_ellipse_position_plots_xyz.R")
  
# Clean the workspace ####
rm(list=setdiff(ls(), c("Embs","Complete","CellF","WideDf","alive","inverted","equalized",
                        "outliers","phaDefects","phaDefects2","compressed","partial","controls","pal",
                        "Nuclei","EmbCycles", "Cellorder", "today", "outDir","basewd", "outliers", 
                        "results", "volumes", "colors", "fourC","eightC","fifteenC","twentyfourC","twentyeightC","fiftysixC")))

source("R_scripts/R_functions.R")
source("R_scripts/R_functions_PCA.R")
  
# Figure 4 - supplement 1 - A #### 
  
  # calculate FDR #### 
  # imput data (no NAs replaced)
  TempDf <- WideDf[,-grep("StartTime",colnames(WideDf))] # remove StartTime since it is redundant for both sisters + same as end time of mother cell
  TempDf <- TempDf[,-grep("V.rel|sisterRatio|asynchrony",colnames(TempDf))] # remove V.rel since it was measured for very small subset of embryos, asynchrony is not significant for any comparison alive/dead
  TempDf <- TempDf[,-grep("Time",colnames(TempDf))] # remove Time variables

source("R_scripts/R_FDR_calculation.R")
  # for 15 cell stage testing all variables/cells combinations
  # excluding embryos with EMS skew
  
  # returns:
  
  #alpha.lim - p.value threshold for 10% FDR, 18 variables
  #s.vars - signif variable a this alpha
  
  # alpha must be set to 0.008 for FDR below 10% for 252 variables
  # in that case FRD is < 8.8%, and there is 18 variables with p.value below 0.005
  s.vars
  
# Figure 4 - supplement A - VOLCANO variable filtering ####
plim <- alpha.lim #p-value cutoff (uncorrected)
fchange <- 0.15 # minimal fold change 15%

cells <- c(fourC, eightC, fifteenC)

# all time variables removed for volcano
source("R_scripts/R_volcano_upto_15C.R")
  #variables with coeficient of variation < 5% removed
  #plot volcano plots 
  
  #10 signif vars from different from volcano (75 NA values imputed)
    # ABara.pDV ABara.pOV Ca.netdis MS.aMean MS.aAP ABpra.aAP ABarp.aAP 
    # ABprp.aAP ABpl.aMean MSp.pOV

#  Figure 4 - supplement B - Correlation among significantly different variables ####
source("R_scripts/R_correlation_plots.R")

#Figure 4H-K: MS division angle / inversion ####
source("R_scripts/R_MS_division_angle.R")
  
# Predict outcome in equalized embryos ####
# LASSO regression with repeated cross-validation

  # imput data for LASSO (NAs replaced) , no Timing
  cells <- c(fourC, eightC, fifteenC, twentyeightC)
  
  TempDf <- WideDf[,-grep("StartTime",colnames(WideDf))] # remove StartTime since it is redundant for both sisters + same as end time of mother cell
  TempDf <- TempDf[,-grep("ABa.EndTime",colnames(WideDf))] #it is 0, since all time is aligned to it
  TempDf <- TempDf[,-grep("V.rel|sisterRatio|asynchrony",colnames(TempDf))] # remove V.rel since it was measured for very small subset of embryos, #asynchrony is not significant for any comparison alive/dead
  TempDf <- TempDf[TempDf$embryo %in% Embs$ID[Embs$MaxTime>50],] #remove embryos with too few timepoints, they would invoke too many NAs 
  
  lassoDf <- prep.data(TempDf, cells=cells, outlierRM = T)

  # Figure 6 ####
  source("R_scripts/R_repeated_CV_lasso_at_stages.R")
 
# Clusterring Heatmaps (not shown) ####
  #source("R_scripts/R_clustering_heatmaps.R")
  #source("R_scripts/R_clustering_all_vars.R") 

# Figure 6 - supplemenet 3 Robustness ####
 # table of change in variability vs change in mean
 # compare equalized vs controls
 # how many variables have greater variance without change in mean

# imput data:
TempDf <- WideDf[,-grep("StartTime",colnames(WideDf))] # remove StartTime since it is redundant for both sisters + same as end time of mother cell
TempDf <- TempDf[,-grep("V.rel|sisterRatio|asynchrony",colnames(TempDf))] # remove V.rel and sisterRatio since it was measured for very small subset of embryos, asynchrony is derived measure
cells=Cellorder[3:202] # Up to ~100C

q <- paste0("(",paste(cells,collapse ="|"),")\\..*") #up to ~100C
cols <- grep(q,colnames(TempDf)) #2300 variables
rows <- Embs$ID[controls|equalized&!(Embs$ID%in%outliers)] #remove EMS outliers

source("R_scripts/R_robustness.R")

# cell-cell distances in PHA-4::GFP abnormal and normal embryos and among groups (rebuttal)####
source("R_scripts/R_cell-cellDist.R")

# not shown: Plot dense tufte plots with all variables ####

# STATS REPORTED IN DENSE PLOTS ARE TAKEN FROM RESULTS WHICH CONTAIN EMS SKEWED EMBRYOS
source("R_scripts/R_dense_plots.R") #plot dense plots without ouliers defined above

#PCA of all embryos the data (not shown) ####
#up to the 15C
cells <- c(fourC, eightC, fifteenC)

#remove rows with too many NA,
TempDf <- TempDf[TempDf$embryo %in% Embs$ID[Embs$MaxTime>50],] #remove embryos with too few timepoints, they would invoke too many NAs 

# function from the functions_PCA.R
pcaDf <- prep.data(TempDf, cells=cells, outlierRM = T)
#removes rows and columns with too many NAs and fills remaining NAs with group means
#removes variables with <5% variance
#removes outliers based EMS spindle and PCA of 15 cell stage (EM10)

source("R_scripts/R_PCA_upto_15C.R")
source("R_scripts/R_PCA_upto_100C.R")

# predict MS flip (not shown) ####
#source("R_scripts/R_repeated_CV_predict_MS-flip.R")
# using variables from other cells at 8 cell stage prior to the division of MS
# model without outliers is really specific, but has poor sensitivity
# 85% AUC, 77% accuracy (n=53)
# "MS.pAP"    "P3.pLR"    "C.EndTime" "ABar.pDV" 
# it 100% successful on the controls

