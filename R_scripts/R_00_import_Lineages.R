if(!require("beeswarm")) install.packages("beeswarm")
if(!require("ggpubr")) install.packages("ggpubr")
if(!require("readxl")) install.packages("readxl")
if(!require("data.table")) install.packages("data.table")
if(!require("plyr")) install.packages("plyr")

#### colors to use in the plots ####

pal <- c("#2E5268","#64A8A1","#D19812","#C93230")

#### set environment ####
par(mar=c(5.1,4.1,4.1,2.1))
today <- format(Sys.Date(), "%y%m%d")

#change this to your local path
zip_lineages <- "/LocalDocs/switchdrive/ev571_embryos/PROJECTS/worm-rules/zip_lineages"
basewd <- "/LocalDocs/switchdrive/ev571_embryos/PROJECTS/worm-rules/"
outDir <- paste0(basewd, "Results/",today,"/")

setwd(basewd)

if(!dir.exists(outDir)){dir.create(outDir);print("Output directory created");}else print("directory already exists")


### import global function ###
source(file="R_Scripts/R_functions.R", local = T)

##read annotations of embryos ####
  Embs <- as.data.frame(read_excel("200723_all_lineaged_embryos.xlsx"), stringsAsFactors=T)
  Embs$Outcome <- factor(Embs$Outcome)
  Embs$Experiment <- factor(Embs$Experiment)
  Embs$Group <- factor(Embs$Group,levels = c("wt","ctrl","alive","dead"))
  
  Embs <- Embs[order(Embs$ID),]
  rownames(Embs) <- Embs$ID
  
#load zipped_lineages in StarryNite format ####

  #StarryNite file format structure
  #int INDEX = 0, X = 5, Y = 6, Z = 7, IDENTITY = 9, SIZE = 8, WT = 10, STATUS = 1, PRED = 2,SUCC1 = 3, SUCC2 = 4, RWT = 11, RSUM = 12, RCOUNT = 13, ASSIGNEDID = 14;
  #it is absolutely unclear what variables V11-V16 mean

  cols <- c("ID", "ALIVE", "PRED", "SUCC1", "SUCC2", "X","Y","Z","SIZE","IDENTITY","WT","RWT","RSUM","RCOUNT","ASSIGNEDID","V16")

Nuclei <- list(); #create new empty list to store data.frames
zipped_folders <- list.files(zip_lineages, pattern=".zip", full.names=T, recursive = T)

#unzip and process
for(eID in Embs$ID){
  print(eID)
  #create temporary folder to unzip files into
  tempd <- tempdir()
  
  #select correct path for given embryo
  EmbZip=zipped_folders[grep(paste0(eID,"_"), zipped_folders)]
  
  unzip(EmbZip,exdir=tempd)
  nuclei_files <- list.files(tempd, pattern="nuclei$",full.names=T, recursive = T)
  
  #load files
  EmbsTps <- list() #list to store timepoints
  for(i in 1:length(nuclei_files)){
    inp <- read.csv(file=nuclei_files[i], header = F);
    inp <- inp[1:10]
    names(inp) <- cols[1:ncol(inp)]
    cat("Read timepoint:", i, "\trows: ", nrow(inp), " cols: ", ncol(inp), "\n");
    
    EmbsTps[[i]] <- inp;
    names(EmbsTps)[i] <- paste0("t",sprintf("%03.f", i));
  }
  Nuclei[[eID]] <- EmbsTps
  unlink(tempd, recursive=TRUE) #close and delete temporary directory
}

rm(cols,eID,inp,i,EmbZip,EmbsTps,nuclei_files, zipped_folders)

#merge all timepoints into a singe data.frame for each embryo 
Nuclei <- lapply(Nuclei,rbindlist, fill=T, idcol="timepoint") 

#extract timing correct timing ####
source(file="R_Scripts/R_extract_timing.R", local = T)

#Normalise timing ####
references <- c("GZ07") #reference embryo for the time alignmnt is set to GZ07 
rlength <- 130 #max timepoints to use for alignment
source(file="R_Scripts/R_NormaliseTime.R", local = T)    

#merge data into a single table
source(file="R_Scripts/R_aggregate_Data.R", local = T)

#run the cell alignment script from Rob Jelier to get positional statistics ####
system("java -jar Lineages2Statistics.jar iniFile=iniFile.txt")

#load table with positional alignment ####
bigT <- load.bigTable("Results/bigTableN.txt")

#merge timing and positional data####
source("R_scripts/R_merge_Timing_and_Positions.R")

#Write max time of each embryo into a file
  maxTimes <- tapply(Mdf$EndTime, Mdf$embryo, function(x) max(x, na.rm = T), simplify = T)
  Embs$MaxTime <- maxTimes[Embs$ID]
  
# merge Mdf with volumes for each cell - "Complete" df ####
  #load volumes
  VOLUMES <- new.env()
  load('191031_Segmented_lineages.RData', envir=VOLUMES, verbose = T)
  volumes <- VOLUMES$volumes
  
  #calculate relative size of sister cells >> asymmetry of division
  volumes$sisterRatio <- NA
  for (cell in unique(volumes$Cell)) {
    row <- volumes$Cell==cell
    
    sister <- get.sisterCell(cell)
    
    v1 <- volumes$V.rel[row]
    names(v1) <- volumes$embryo[row]
    
    #obtain paired volumes
    v2 <- volumes$V.rel[volumes$Cell==sister]
    names(v2) <- volumes$embryo[volumes$Cell==sister]  
    
    ratio <- v1/v2[names(v1)]
    volumes$sisterRatio[intersect(which(volumes$embryo %in% names(v1)), which(volumes$Cell==cell))] <- ratio
  }
  
#merge the two
  Complete <- merge(Mdf, volumes[,c(1:3,8)], by=c("Cell","embryo"), all=T)
  #set variables to NA for cells that didn't divide
  Complete[Complete$Divided==0, c("EndTime","pAP","pLR","pDV","pOV","aAP","aLR","aDV","aMean")] <- NA
  Complete[is.na(Complete$Divided),c("Divided")] <- 1 # fix 3 cells
  
  #calculate asynchrony of divisions of sister cells ####
  Complete$asynchrony <- NA
  for (cell in unique(Complete$Cell)) {
    row <- Complete$Cell==cell
    sister <- get.sisterCell(cell)
    
    cc1 <- Complete$LifeTime[row]
    names(cc1) <- Complete$embryo[row]
    
    #obtain paired volumes
    cc2 <- Complete$LifeTime[Complete$Cell==sister]
    names(cc2) <- Complete$embryo[Complete$Cell==sister]  
    
    ratio <- cc1/cc2[names(cc1)]
    Complete$asynchrony[intersect(which(Complete$embryo %in% names(cc1)), which(Complete$Cell==cell))] <- ratio
  }

  Complete$V.rel <- Complete$V.rel*100 # to express relative cell volume directly in percentage and not as fraction
  
# Debug - batch effects ####
    # Complete$segmented <- Embs[Complete$embryo,"Segmented"]
    # c="MSa"
    # rows <- Complete$Cell==c
    # boxplot(LifeTime~segmented+Group, data=Complete[rows,],main=c)
    # beeswarm(LifeTime~segmented+Group, data=Complete[rows,], pch=16, pwcol = Complete$segmented[rows]+1, add=T)
    # 
    # Complete$segmented <- NULL
  
#Save data####
  Cellorder <- read.table("CellOrder.txt",sep="\t",header=FALSE)
  Cellorder <- as.character(Cellorder[,2])
  
  filename <- paste0(outDir, today,"_embryos.RData")
  save(Complete, Embs, Nuclei, Data, Cellorder ,CellF, pal, file=filename)
  print(paste0("Data have been aggregated and saved as: ", filename))

rm(list=setdiff(ls(), c("Embs","Complete","Cellorder","Data","pal", "CellF", "Nuclei","EmbCycles")))
