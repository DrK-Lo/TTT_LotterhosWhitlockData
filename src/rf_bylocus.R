#######################################
#
# rf_bylocus.R
# created by Brenna Forester
# modified by Brett Ford
# Created 20180306
#
# This script replicates analyses in Forester et al. (2018)
# Use repRF.sh to run this script across multiple nodes
# Edit paths to make sure code runs properly
#
# Lotterhos KE, Whitlock MC (2015) Data from: The relative power of genome scans 
# to detect local adaptation depends on sampling design and statistical method. 
# Dryad Digital Repository. http://dx.doi.org/10.5061/dryad.mh67v
#
# Usage: Rscript rf_bylocus.R
#
######################################

############################################################################ 
# Data preparation for all GEA methods used to detect multilocus selection #
############################################################################

# input arguments from the command line
args <- commandArgs(TRUE)

# what simulation
seed <- args[1]

# Load required libraries
library(vegan)
library(tools)
library(randomForest)

# Specify paths
forester_simfiles_path <- "/Users/brettford/Desktop/Northeastern/coding/forester_simulation_code/forester_sim_code/TTT_LotterhosWhitlockData/forester_simfiles/"
forester_results_path <- "/Users/brettford/Desktop/Northeastern/coding/forester_simulation_code/forester_sim_code/TTT_LotterhosWhitlockData/forester_results/"

# Specify working directory
setwd(forester_simfiles_path)

# Coordinates 
# -----------

Coords <- list()
Coords$Pairs <- list()
#BF: Coords for 3 diff landscapes under diff sampling strategies
Coords$Pairs$E453 <- read.table("1351142954_453EnviMatPAIRS_ED.txt")
Coords$Pairs$E988 <- read.table("1351142970_988EnviMatPAIRS_ED.txt")
Coords$Pairs$E950 <- read.table("1351142986_950EnviMatPAIRS_ED.txt")
Coords$Transect <- list()
Coords$Transect$E453 <- read.table("1351142954_453EnviMatTRANSECTS_ED_Design.txt")
Coords$Transect$E988 <- read.table("1351142970_988EnviMatTRANSECTS_ED_Design.txt")
Coords$Transect$E950 <- read.table("1351142986_950EnviMatTRANSECTS_ED_Design.txt")
Coords$Random <- list()
#Set random scheme as its own dataframe, for each landscape, in coords object
#Same random set for each landscape
#Each of the three files (E453, E988, E950) have different coords and different variations of same env
Coords$Random$E453 <- Coords$Random$E988 <- Coords$Random$E950 <- read.table("SchemeRandom1.txt")

#For each landscape:
for(k in 1:length(Coords$Transect))
{
  names(Coords$Transect[[k]])[2:3] <- names(Coords$Pairs[[1]])[2:3]       # use same x/y colnames for Transect as for Pairs
  Coords$Transect[[k]][,7:9][is.na(Coords$Transect[[k]][,7:9])] <- FALSE  # Change NA values to FALSE
  #columns 7-9 are true/false for three different transect types?
}

for(i in 1:length(Coords))
{
  for(k in 1:length(Coords[[i]]))
  {
    b <- order(Coords[[i]][[k]][,3], Coords[[i]][[k]][,2]) # Sort by y, then x
    Coords[[i]][[k]] <- Coords[[i]][[k]][b,]  
  }
}

# Extract env data
# ----------------

Filenames.env <- list.files(path=forester_simfiles_path, pattern="env")
## For all 121 env simfiles, create a vector
Env <- list()
for(i in 1:length(Filenames.env))
{
  tmp <- read.table(paste0(forester_simfiles_path, Filenames.env[[i]]))
  Env[[i]] <- unlist(tmp)
}

# Table of design information for each simulation
# -----------------------------------------------   

Filenames.lfmm <- list.files(path=forester_simfiles_path, pattern="lfmm") 

test <- Reduce(rbind, strsplit(Filenames.lfmm, split="[_=]"))                 
test <- cbind(test, Reduce(rbind, strsplit(test[,ncol(test)], split="[.]")))  
tmp <- strsplit(test[,2], split="[.xs]")
#BF: Split based on transect scheme?
test2 <- matrix(NA, nrow(test), 3)
# made empty matrix with 240 rows and 3 columns
for(i in 1:nrow(test)) test2[i,1:length(tmp[[i]])] <- tmp[[i]]
test2 <- gsub("[A-Z, a-z]","", test2)
test <- cbind(test, test2)
dimnames(test) <- list(NULL, c("Demography", "Design", "ID", "Env", "ID2", 
                               "V6", "NumPops", "V8", "V9", "NumInd", "Design2",
                               "NumPops2", "NumTrans", "NumInd2"))
#BF: This adds number of pops, number of transects, and number of individuals per transect?
Design <- as.data.frame(test[,c(1,2,4,7,10,13,14)])
Design$Env <- ordered(Design$Env, levels=c(453,988,950))
Design$Type <- ordered(substr(Design$Design, 1, 1), levels=c("P", "T", "R"))
Design$Design <- as.character(Design$Design)
Design$Design[Design$Design == "T30.T3x10"] <- "T30.3x10s"
Design$Design[Design$Design == "T30.T6x5s"] <- "T30.6x5s"
Design$NumPops <- as.numeric(as.character(Design$NumPops))


# Selection information
# ---------------------   

Filenames.sum <- list.files(path=forester_results_path, pattern = "Cpval")
#Contains summary statistics for each locus in each dataset from Bayenv2, LFMM, PCAdapt
test <- Reduce(rbind, strsplit(Filenames.sum, split="[_=]"))   
test[,9] <- gsub("Bayenv", "_Bayenv", test[,9] )
test <- cbind(test, Reduce(rbind, strsplit(test[,ncol(test)], split="[.]")))  
tmp <- strsplit(test[,2], split="[.xs]")
test2 <- matrix(NA, nrow(test), 3)
for(i in 1:nrow(test)) test2[i,1:length(tmp[[i]])] <- tmp[[i]]
test2 <- gsub("[A-Z, a-z]","", test2)
test <- cbind(test, test2)
dimnames(test) <- list(NULL, c("Demography", "Design", "ID", "Env", "ID2", 
                               "V6", "NumPops", "V8", "V9", "NumInd", "V11",
                               "NumPops2", "NumTrans", "NumInd2"))
Design.sel <- as.data.frame(test[,c(1,2,4,7,10,13,14)])
Design.sel$Env <- ordered(Design.sel$Env, levels=c(453,988,950))
Design.sel$Type <- ordered(substr(Design.sel$Design, 1, 1), levels=c("P", "T", "R"))
Design.sel$Design <- as.character(Design.sel$Design)
Design.sel$Design[Design.sel$Designl == "T30.T3x10"] <- "T30.3x10s"
Design.sel$Design[Design.sel$Design == "T30.T6x5s"] <- "T30.6x5s"
Design.sel$NumPops <- as.numeric(as.character(Design.sel$NumPops))


# IMPORTANT: Data files match for 90 Pop sims, but not for others!
# ----------------------------------------------------------------

# Use only these sites (the ones with 90 populations):
Sites.90 <- c(1:nrow(Design))[Design$NumPops==90]
Sel.90   <- c(1:nrow(Design.sel))[Design.sel$NumPops==90]


# Code to replicate analyses in Forester et al. (submitted) 

##########################################
# Random Forest 
##########################################

# Data prep
# ---------

j = as.numeric(seed)  # Simulation data set #2 (see notes above)
verify <- "Do the files match?"

i=Sites.90[j]
Site <- rep(1:Design$NumPops[i],each=as.numeric(as.vector(Design$NumInd)[i])) 
tmp <- read.table(paste0(forester_simfiles_path, Filenames.lfmm[[i]])) # the genetic data
   # number of loci

z=Sel.90[j]
doubleCheck <- identical(unlist(Design[i,]), unlist(Design.sel[z,]))# sanity check that the files match up
verify <- rbind(verify, doubleCheck)
file1 <- as.data.frame(read.table(paste0(forester_results_path, Filenames.sum[[z]]), header=T))
colnames(tmp) <- file1$SNPnames[file1$SNPIncluded==TRUE]
nloci <- ncol(tmp)

tmp.Site <- split(tmp, Site)

# SNPs are assigned sequentially in the sims, so they need to be randomized within populations.
# Use same seed for each population (1:90)
for (rs in 1:length(tmp.Site)) {               
  set.seed(rs)
  tmp.Site[[rs]] <- apply(as.data.frame(tmp.Site[[rs]]),2,sample)
}

snps <- as.data.frame(do.call("rbind", tmp.Site))
snps[snps==2] <- 1    # set 2s to 1 for snmf, LEA
colnames(snps) <- c(1:ncol(snps))

coord.full <- data.matrix(Coords[[as.numeric(Design$Type[i])]][[as.numeric(Design$Env[i])]])
des <- coord.full[ , grepl( Design$Design[i] , colnames(coord.full) ) ]
coord.full.des <- cbind(coord.full[,2:3],des)
coord.des <- coord.full.des[coord.full.des[,3]==1,]
x <- unsplit(coord.des[,1], Site); y <- unsplit(coord.des[,2], Site)
coord <- cbind(x,y); colnames(coord) <- c("X_Pops","Y_Pops")

Hab <- unname(Env[[i]])
predictors <- as.data.frame(scale(Hab, center=T, scale=T)); colnames(predictors) <- c("Hab")

# Random Forest
# -------------

# NOTE: run time is 3-10 days depending on size of data set and computational capacity

# Full RF runs - 3 replicates for habitat and x
# mtry default for *regression* is p/3 = ~3300 loci at each split
# NOTE: mtry default would need to be dramatically increased for classification!

rf.h1 <- randomForest(x=snps, y=predictors[,1], importance=TRUE , proximity=TRUE, ntree=10000) 
rf.h2 <- randomForest(x=snps, y=predictors[,1], importance=TRUE , proximity=TRUE, ntree=10000) 
rf.h3 <- randomForest(x=snps, y=predictors[,1], importance=TRUE , proximity=TRUE, ntree=10000) 

# Get importance (mean decrease in accuracy) for each full run 
imp.h1 <- data.frame(importance(rf.h1, type=1))
imp.h2 <- data.frame(importance(rf.h2, type=1))
imp.h3 <- data.frame(importance(rf.h3, type=1))

message("Pull best % of loci from...")
# Pull best % of loci from 0.5%, 1%, 1.5%, 2% 

best.loci <- c(0.995,0.990,0.985,0.980)  
PVEs <- data.frame(HabPVE = rep(NA, length(best.loci))); rownames(PVEs) <- best.loci

for(b in 1:length(best.loci)) {
  names.best1 <- rownames(imp.h1)[which(imp.h1$X.IncMSE > quantile(imp.h1$X.IncMSE, probs=best.loci[b]))]
  names.best2 <- rownames(imp.h2)[which(imp.h2$X.IncMSE > quantile(imp.h2$X.IncMSE, probs=best.loci[b]))]
  names.best3 <- rownames(imp.h3)[which(imp.h3$X.IncMSE > quantile(imp.h3$X.IncMSE, probs=best.loci[b]))]
  names.best.unique <- unique(c(names.best1, names.best2, names.best3))
  
  best.dat <- snps[,colnames(snps) %in% names.best.unique]
  
  best.rf1 <- randomForest(x=best.dat, y=predictors[,1], importance=TRUE , proximity=TRUE, ntree=10000)
  best.rf2 <- randomForest(x=best.dat, y=predictors[,1], importance=TRUE , proximity=TRUE, ntree=10000)
  best.rf3 <- randomForest(x=best.dat, y=predictors[,1], importance=TRUE , proximity=TRUE, ntree=10000)
  
  PVEs[b,1] <- mean(c(tail(best.rf1$rsq, n=1), tail(best.rf2$rsq, n=1), tail(best.rf3$rsq, n=1)))
}


h.start <- as.numeric(rownames(PVEs[which.max(PVEs[,1]),,drop=F]))  # drop=F keeps row names in single column data frame
h.start <- (1-(2*(1-h.start)))                                      # double the number of loci detected as "important"


### Purging for Habitat:

names.best1 <- rownames(imp.h1)[which(imp.h1$X.IncMSE > quantile(imp.h1$X.IncMSE, probs=h.start))]
names.best2 <- rownames(imp.h2)[which(imp.h2$X.IncMSE > quantile(imp.h2$X.IncMSE, probs=h.start))]
names.best3 <- rownames(imp.h3)[which(imp.h3$X.IncMSE > quantile(imp.h3$X.IncMSE, probs=h.start))]
names.best.unique <- unique(c(names.best1, names.best2, names.best3))
n <- length(names.best.unique)
dat.purge <- snps[,colnames(snps) %in% names.best.unique]

# Run 3 iterations of random forest with these loci
rf.purge1 = randomForest(x=dat.purge, y=predictors[,1], importance=TRUE , proximity=TRUE, ntree=10000)
rf.purge2 = randomForest(x=dat.purge, y=predictors[,1], importance=TRUE , proximity=TRUE, ntree=10000)
rf.purge3 = randomForest(x=dat.purge, y=predictors[,1], importance=TRUE , proximity=TRUE, ntree=10000)
message("Initiate backward purging")

# Initiate backward purging 
names_hab <- list(); names_hab[[n]] <- names.best.unique
imp_hab <- list()

pve_hab <- data.frame(V1=rep(NA,n), V2=rep(NA,n), V3=rep(NA,n))
pve_hab[n,] <- c(tail(rf.purge1$rsq, n=1), tail(rf.purge2$rsq, n=1), tail(rf.purge3$rsq, n=1))

for(p in 1:(n-1)) {
  imp.purge1 <- data.frame(importance(rf.purge1, type=1))
  imp.purge2 <- data.frame(importance(rf.purge2, type=1))
  imp.purge3 <- data.frame(importance(rf.purge3, type=1))
  all_imp <- cbind(imp.purge1,imp.purge2,imp.purge3)
  
  all_imp$average <- apply(all_imp[,1:3], 1, mean)
  imp_hab[[n-p+1]] <- all_imp$average
  dont_keep <- which(all_imp[,'average']==min(all_imp[,'average']))    # row no. of locus with lowest average importance
  
  if (length(dont_keep)==1) {
    table_keep <- all_imp[-dont_keep,]                              
  } else {
    table_keep <- all_imp[-dont_keep[sample(x=dont_keep, size=1)],]    # if 2 loci are as unimportant, pick one
  }
  
  names_keep <- rownames(table_keep)
  names_hab[[n-p]] <- names_keep
  
  dat.purge <- snps[,colnames(snps) %in% names_keep]
  rf.purge1 = randomForest(x=as.data.frame(dat.purge), y=predictors[,1], importance=TRUE , proximity=TRUE, ntree=10000) # add as.data.frame for final single locus
  rf.purge2 = randomForest(x=as.data.frame(dat.purge), y=predictors[,1], importance=TRUE , proximity=TRUE, ntree=10000)
  rf.purge3 = randomForest(x=as.data.frame(dat.purge), y=predictors[,1], importance=TRUE , proximity=TRUE, ntree=10000)
  
  pve_hab[n-p,] <- c(tail(rf.purge1$rsq, n=1), tail(rf.purge2$rsq, n=1), tail(rf.purge3$rsq, n=1))
}
message("Grab last importance value")
# Grab last importance value:
imp.purge1 <- data.frame(importance(rf.purge1, type=1))
imp.purge2 <- data.frame(importance(rf.purge2, type=1))
imp.purge3 <- data.frame(importance(rf.purge3, type=1))
all_imp <- cbind(imp.purge1,imp.purge2,imp.purge3)

all_imp$average <- apply(all_imp[,1:3], 1, mean)
imp_hab[[n-p]] <- all_imp$average #FINAL MEAN IMPORTANCE VALUE FOR THE MOST IMPORTANT LOCI

# Pull set of loci with maximum PVE
pve_hab$Average<-apply(pve_hab,1,mean)
pve_hab_max <- which(pve_hab$Average==max(pve_hab$Average, na.rm=T)) # Maximum percent variance explained
pve_hab_maximum <- max(pve_hab$Average, na.rm=T)

cand.RFh <- as.numeric(names_hab[[pve_hab_max]])
cand.imp.RFh <- imp_hab[[pve_hab_max]]

# Variable importance:
hout <- as.data.frame(cbind(cand.RFh, cand.imp.RFh)); names(hout) <- c("RFh.loci","RFh.mse")

# ----------------------------------
# Detection levels using full output
# ----------------------------------
# Because RF only produces a set of SNPs (without scores), we don't used 
# cutoffs or calculate empirical p-values as we do for the other methods.
# Instead we pull the SNPs that are TP detections and use those to calculate metrics:
message("TPR and FPR")
# TPR and FPR
# -----------
# NOTE: number of neutral loci is always 9900 in these simulations

cand.RF <- as.numeric(cand.RFh)

TPR <- sum(cand.RF > 9900) / (nloci-9900)  # also the "empirical p-value" output
FPR <- sum(cand.RF <= 9900) / 9900


# Selection levels detected
# -------------------------

#tp <- cand.RF[cand.RF > 9900]

#z=Sel.90[j]
#tmp <- read.table(paste0(base_path,"/SummaryFiles/", Filenames.sum[[z]]), header=T) # the summary file
tmp <- tmp[(tmp$IsNeut=="Sel" & tmp$SNPIncluded==TRUE), 5]
#n.sel <- length(tmp)
#tmp <- as.data.frame(cbind(c(9901:(9900+n.sel)), tmp))
#names(tmp) <- c("SNP","Sel")

#sel001 <- tmp[tmp$Sel==0.001, 1]
#sel005 <- tmp[tmp$Sel==0.005, 1]
#sel01  <- tmp[tmp$Sel==0.01, 1]
#sel1   <- tmp[tmp$Sel==0.1, 1]

#TPR001 <- sum(tp %in% sel001) / length(sel001)
#TPR005 <- sum(tp %in% sel005) / length(sel005)
#TPR01  <- sum(tp %in% sel01) / length(sel01)
#TPR1   <- sum(tp %in% sel1) / length(sel1)


# -----------------
# Write out results
# -----------------
message("Writing results")
rf_total <- as.data.frame(cbind(id=c(1:nloci), SNPnames=file1$SNPnames[file1$SNPIncluded==TRUE], rf_log=rep(NA, nloci)))

for (q in 1:nloci) {
  if(rf_total[q,1] %in% cand.RF) {
    rf_total$rf_log[q] <- "TRUE"
  }
  else {
    rf_total$rf_log[q] <- "FALSE"
  }
}

imp_total <- as.data.frame(cbind(id=c(1:length(imp.h1)), imp_1=imp.h1$X.IncMSE, imp_2=imp.h2$X.IncMSE, imp_3=imp.h3$X.IncMSE))

i=Sites.90[j]

simulation <- file_path_sans_ext(Filenames.lfmm[[i]])

write.table(rf_total, file=paste0(forester_results_path, "rf_results/RF_locus_stats_", simulation, ".txt"), sep = " ", row.names = F)
write.table(imp_total, file=paste0(forester_results_path, "rf_results/RF_importance_", simulation, ".txt"), sep = " ", row.names = F)
write.table(hout, file=paste0(forester_results_path, "rf_results/RF_bestloci_imp_", simulation, ".txt"), sep = " ", row.names = F)
write.table(verify, file= paste0(forester_results_path, "rf_results/", j, "_verifying_match_RF.txt"), sep = " ", row.names =F)
