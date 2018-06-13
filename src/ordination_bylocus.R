#######################################
#
# ordination_by_locus.R
# created by Brenna Forester
# modified by Brett Ford
# Created 20180308
#
# This script replicates analyses in Forester et al. (2018)
# Use repORD.sh to run this script across multiple nodes
# Edit paths to make sure code runs properly
#
# Lotterhos KE, Whitlock MC (2015) Data from: The relative power of genome scans 
# to detect local adaptation depends on sampling design and statistical method. 
# Dryad Digital Repository. http://dx.doi.org/10.5061/dryad.mh67v
#
# Usage: Rscript ordination_by_locus.R
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

Filenames.env <- list.files(path=forester_simfiles_path, pattern=".env")
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
#BF: Splits based on transect scheme

test2 <- matrix(NA, nrow(test), 3)

for(i in 1:nrow(test)) test2[i,1:length(tmp[[i]])] <- tmp[[i]]
test2 <- gsub("[A-Z, a-z]","", test2)
test <- cbind(test, test2)

dimnames(test) <- list(NULL, c("Demography", "Design", "ID", "Env", "ID2", 
                               "V6", "NumPops", "V8", "V9", "NumInd", "Design2",
                               "NumPops2", "NumTrans", "NumInd2"))

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

##########################################
# Constrained Ordinations: RDA and dbRDA #
##########################################

# Function to calculate empirical p-values
# ----------------------------------------
# Code provided by K. Lotterhos
# NOTE: does not replicate empPvals in qvalue package

getEmpP <- function(Obs, sort.Null){
  
  # Obs is a single observed value
  # sort.Null is a list of the null values in ascending order
  # To apply to multiple SNPS:
  # sapply(AllObservedValues, getEmpP, sort(NonCodingValues))
  
  if(is.na(Obs)){
    return(NA)
  }else{
    options(warn=-1)
    out = max(which(sort.Null<=Obs))/length(sort.Null)
    if(is.infinite(out)==TRUE){out=0} ## above code is undefined when XXT < min(sort.Null)
    return(out)
  }
}


# Function to detect outliers from RDA, dbRDA, set cutoff as z
# ------------------------------------------------------------
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)                   # find +/- z sd from mean loading     
  as.numeric(names(out <- x[x < lims[1] | x > lims[2]]))   # locus names in these tails
}

# Data prep
# ---------

j = as.numeric(seed)

verify <- "Do the files match?"

i=Sites.90[j] #the number of file (x out of 72) in all 242 simulation files
Site <- rep(1:Design$NumPops[i],each=as.numeric(as.vector(Design$NumInd)[i])) 
tmp <- read.table(paste0(forester_simfiles_path, Filenames.lfmm[[i]])) # the genetic data

z=Sel.90[j]
doubleCheck <- identical(unlist(Design[i,]), unlist(Design.sel[z,])) # sanity check that the files match up
verify <- rbind(verify, doubleCheck)
file1 <- read.table(paste0(forester_results_path, Filenames.sum[[z]]), header=T)
colnames(tmp) <- file1$SNPnames[file1$SNPIncluded==TRUE]
nloci <- ncol(tmp)   # number of loci; not exactly 10,000 because those with low He were removed

tmp.Site <- split(tmp, Site)

# SNPs are assigned sequentially in the sims, so they need to be randomized within populations.
for (rs in 1:length(tmp.Site)) {               
  set.seed(rs)
  tmp.Site[[rs]] <- apply(as.data.frame(tmp.Site[[rs]]),2,sample)
}

snps <- as.data.frame(do.call("rbind", tmp.Site))
snps[snps==2] <- 1    # set 2s to 1 for snmf, LEA
colnames(snps) <- c(1:ncol(snps))

coord.full <- data.matrix(Coords[[as.numeric(Design$Type[i])]][[as.numeric(Design$Env[i])]])
des <- coord.full[ , grepl( Design$Design[i] , colnames(coord.full) ) ] #search for P90 in coord.full columns and set equal to des
coord.full.des <- cbind(coord.full[,2:3],des) #bind corrdinates with des
coord.des <- coord.full.des[coord.full.des[,3]==1,] #making sure there are only 3 columns?
x <- unsplit(coord.des[,1], Site); y <- unsplit(coord.des[,2], Site) #sorts by y then x coords
coord <- cbind(x,y); colnames(coord) <- c("X_Pops","Y_Pops")

Hab <- unname(Env[[i]])
predictors <- scale(Hab, center=T, scale=T)

# RDA and dbRDA
# -------------

# NOTE: code at end of file includes correction for population structure & would be substituted here.

# RDA
snps.scale <- scale(snps, center=T, scale=T)              # scale and center for RDA
snp.rda <- rda(snps.scale, predictors, scale=F)           # could also use scale=T instead of scaling beforehand                  
load.rda <- summary(snp.rda)$species[,1]                  # RDA loadings
load.rda_df <- as.data.frame(load.rda)

# dbRDA
snps.bray <- vegdist(snps, method="bray")                 # bray-curtis for dbRDA
snp.dbrda <- capscale(snps.bray ~ predictors, comm=snps)  
load.dbrda <- summary(snp.dbrda)$species[,1]              # dbRDA loadings
load.dbrda_df <- as.data.frame(load.dbrda)

emp_total <- cbind(id=c(1:nloci), SNPnames=file1$SNPnames[file1$SNPIncluded==TRUE], load.rda_df)
emp_total_db <- cbind(id=c(1:nloci), SNPnames=file1$SNPnames[file1$SNPIncluded==TRUE], load.dbrda_df)
load.out <- cbind(load.rda, load.dbrda)  # loadings for RDA and dbRDA

# -----------------------------
# Power from empirical p-values
# -----------------------------

# RDA
neu.ep <- abs(load.rda[1:9900])

sel.ep <- abs(load.rda[9901:nloci])

emp.p_r <- sapply(c(neu.ep,sel.ep), getEmpP, sort(neu.ep))  # p-value from cumulative distribution

emp.p_r <- 1-emp.p_r      # p-value for qvalue
emp.p_r_df <- as.data.frame(emp.p_r)
emp_total <- cbind(emp_total, emp.p_r_df)

library(qvalue)
#source("https://bioconductor.org/biocLite.R")
#biocLite("qvalue")

emp.q_r <- qvalue(emp.p_r)$qvalues
emp.q_r_df <- as.data.frame(emp.q_r)
emp_total <- cbind(emp_total, emp.q_r_df)

#rT.empirical <- as.numeric(names(which(emp.q_r[9901:nloci] < 0.01)))
#rTP.empirical <- sum(emp.q_r[9901:nloci] < 0.01, na.rm=TRUE)
#rTPR.empirical <- mean(emp.q_r[9901:nloci] < 0.01, na.rm=TRUE)
emp_total$emp.p_r_log <- "NA"

for (i in 1:nloci) {
  if(emp_total[i,5] < 0.01) {
    emp_total$emp.p_r_log[i] <- "TRUE"
  }
  else {
    emp_total$emp.p_r_log[i] <- "FALSE"
  }
}

# dbRDA
neu.ep <- abs(load.dbrda[1:9900])
sel.ep <- abs(load.dbrda[9901:nloci])

emp.p_db <- sapply(c(neu.ep,sel.ep), getEmpP, sort(neu.ep))  # p-value from cumulative distribution
emp.p_db <- 1-emp.p_db         # p-value for qvalue
emp.p_db_df <- as.data.frame(emp.p_db)
emp_total_db <- cbind(emp_total_db, emp.p_db_df)

emp.q_db <- qvalue(emp.p_db)$qvalues
emp.q_db_df <- as.data.frame(emp.q_db)
emp_total_db <- cbind(emp_total_db, emp.q_db_df)
#head(emp_total_db)
#dT.empirical <- as.numeric(names(which(emp.q_db[9901:nloci] < 0.01)))
#dTP.empirical <- sum(emp.q_db[9901:nloci] < 0.01, na.rm=TRUE)
#dTPR.empirical <- mean(emp.q_db[9901:nloci] < 0.01, na.rm=TRUE)

emp_total_db$emp.p_db_log <- "NA"

for (i in 1:nloci) {
  if(emp_total_db[i,5] < 0.01) {
    emp_total_db$emp.p_db_log[i] <- "TRUE"
  }
  else {
    emp_total_db$emp.p_db_log[i] <- "FALSE"
  }
}

# Selection info
# --------------
z=Sel.90[j]
#doubleCheck <- identical(unlist(Design[i,]), unlist(Design.sel[z,]))  # sanity check that the files match up
#tmp <- read.table(paste0(forester_results_path, Filenames.sum[[z]]), header=T) # the summary file
#tmp <- tmp[(tmp$IsNeut=="Sel" & tmp$SNPIncluded==TRUE), 5]
#The above line of code states to subset the data to only include loci under selection (IsNeut=="Sel") and those that are used (IncludeSNP=="TRUE)
#n.sel <- length(tmp)
#tmp <- as.data.frame(cbind(c(9901:(9900+n.sel)), tmp))
#names(tmp) <- c("SNP","Sel")

#sel001 <- tmp[tmp$Sel==0.001, 1]
#sel005 <- tmp[tmp$Sel==0.005, 1]
#sel01  <- tmp[tmp$Sel==0.01, 1]
#sel1   <- tmp[tmp$Sel==0.1, 1]


# Levels of selection detected from emp pvalue detections
# -------------------------------------------------------
#eTPR001r <- sum(rT.empirical %in% sel001) / length(sel001)
#eTPR005r <- sum(rT.empirical %in% sel005) / length(sel005)
#eTPR01r  <- sum(rT.empirical %in% sel01) / length(sel01)
#eTPR1r   <- sum(rT.empirical %in% sel1) / length(sel1)

#eTPR001d <- sum(dT.empirical %in% sel001) / length(sel001)
#eTPR005d <- sum(dT.empirical %in% sel005) / length(sel005)
#eTPR01d  <- sum(dT.empirical %in% sel01) / length(sel01)
#eTPR1d   <- sum(dT.empirical %in% sel1) / length(sel1)


#View(emp_total)
# ------------------------------
# Detection levels using cutoffs
# ------------------------------

# TPR and FPR
# -----------
# NOTE: number of neutral loci is always 9900 in these simulations

cand.rda25 <- outliers(load.rda, 2.5)   #2.5
cand.rda3  <- outliers(load.rda, 3)     #3

cand.dbrda25 <- outliers(load.dbrda,2.5) #2.5
cand.dbrda3  <- outliers(load.dbrda,3)   #3

emp_total$cand.rda25 <- "NA"
emp_total$cand.rda3 <- "NA"
emp_total_db$cand.dbrda25 <- "NA"
emp_total_db$cand.dbrda3 <- "NA"

for (m in 1:nloci) {
  if(emp_total[m,1] %in% cand.rda25) {
    emp_total$cand.rda25[m] <- "TRUE"
  }
  else {
    emp_total$cand.rda25[m] <- "FALSE"
  }
}

for (n in 1:nloci) {
  if(emp_total[n,1] %in% cand.rda3) {
    emp_total$cand.rda3[n] <- "TRUE"
  }
  else {
    emp_total$cand.rda3[n] <- "FALSE"
  }
}

for (o in 1:nloci) {
  if(emp_total_db[o,1] %in% cand.dbrda25) {
    emp_total_db$cand.dbrda25[o] <- "TRUE"
  }
  else {
    emp_total_db$cand.dbrda25[o] <- "FALSE"
  }
}

for (p in 1:nloci) {
  if(emp_total_db[p,1] %in% cand.dbrda3) {
    emp_total_db$cand.dbrda3[p] <- "TRUE"
  }
  else {
    emp_total_db$cand.dbrda3[p] <- "FALSE"
  }
}

# Selection levels detected
# -------------------------

#tp25r <- cand.rda25[cand.rda25 > 9900]
#tp3r  <- cand.rda3[cand.rda3 > 9900]
#tp25d <- cand.dbrda25[cand.dbrda25 > 9900]
#tp3d  <- cand.dbrda3[cand.dbrda3 > 9900]
#head(emp_total)
# 2.5 SD results
#TPR001r25 <- sum(tp25r %in% sel001) / length(sel001)
#TPR005r25 <- sum(tp25r %in% sel005) / length(sel005)
#TPR01r25  <- sum(tp25r %in% sel01) / length(sel01)
#TPR1r25   <- sum(tp25r %in% sel1) / length(sel1)

#TPR001d25 <- sum(tp25d %in% sel001) / length(sel001)
#TPR005d25 <- sum(tp25d %in% sel005) / length(sel005)
#TPR01d25  <- sum(tp25d %in% sel01) / length(sel01)
#TPR1d25   <- sum(tp25d %in% sel1) / length(sel1)

# 3 SD results
#TPR001r3 <- sum(tp3r %in% sel001) / length(sel001)
#TPR005r3 <- sum(tp3r %in% sel005) / length(sel005)
#TPR01r3  <- sum(tp3r %in% sel01) / length(sel01)
#TPR1r3   <- sum(tp3r %in% sel1) / length(sel1)

#TPR001d3 <- sum(tp3d %in% sel001) / length(sel001)
#TPR005d3 <- sum(tp3d %in% sel005) / length(sel005)
#TPR01d3  <- sum(tp3d %in% sel01) / length(sel01)
#TPR1d3   <- sum(tp3d %in% sel1) / length(sel1)

# -----------------
# Write out results
# -----------------

emp_total_final <- cbind(emp_total, emp_total_db[, 3:8], corresponding_file=Filenames.sum[[z]])
i=Sites.90[j]
simulation <- file_path_sans_ext(Filenames.lfmm[[i]])
write.table(emp_total_final, file= paste0(forester_results_path,"ordination_results", "/ORD_locus_stats_", simulation, ".txt"), sep = " ", row.names = F)
j= as.numeric(seed)
write.table(verify, file= paste0(forester_results_path,"ordination_results/", j, "verifying_match.txt"), sep = " ", row.names =F)
