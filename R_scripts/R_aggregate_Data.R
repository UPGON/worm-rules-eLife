##################  Generate a single huge table  ##################
Data <- rbindlist(ScaledEmbs, use.names = T, idcol="embryo")
Data  <- cbind(Data, Embs[Data$embryo,c("Outcome","AB_rel","Group","Size")])

#treat initial cells in each embryo - distances and LifeTimes don't make sence for cells that were not born during the acquisition
ncellAtStart <- unlist(lapply(Nuclei, function(x) sum(x$timepoint=="t001"))) #get number of cells at the first frame
less4 <- names(ncellAtStart)[ncellAtStart==4]

Data[Data$Cell %in% c("AB", "P1"), c("LifeTime","StartTime")] <- NA
Data[Data$Cell %in% c("ABa","ABp", "EMS","P2") & Data$embryo %in% less4, c("LifeTime","StartTime")] <- NA
Data[Data$Cell %in% c("ABa","ABp", "EMS","P2") & Data$embryo=="MM4", c("LifeTime","StartTime")] <- NA

#fix time - add 1 scaled interval to the lifeTime and End Time so that there is no gap in time between EndTime and Birth Time of daughter cells
interval <- Embs$interval*Embs$timeScaling
names(interval) <- Embs$ID
Data[,c("LifeTime","EndTime")] <- Data[,c("LifeTime","EndTime")] + (interval[Data$embryo]/60)
