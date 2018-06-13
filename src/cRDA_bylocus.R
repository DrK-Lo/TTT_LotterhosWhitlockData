#######################################
#
# cRDA_by_locus.R
# created by Brenna Forester
# modified by Brett Ford
# Created 20180309
#
# This script replicates analyses in Forester et al. (2018)
# Use repORD.sh to run this script across multiple nodes
# Edit paths to make sure code runs properly
#
# Lotterhos KE, Whitlock MC (2015) Data from: The relative power of genome scans 
# to detect local adaptation depends on sampling design and statistical method. 
# Dryad Digital Repository. http://dx.doi.org/10.5061/dryad.mh67v
#
# Usage: Rscript cRDA_by_locus.R
#
######################################

# input arguments from the command line
args <- commandArgs(TRUE)

# what simulation
seed <- args[1]

# Load required libraries
library(vegan)
library(tools)
library(hornpa)
library(pracma)
library(qvalue)

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

Filenames.env <- list.files(path= forester_simfiles_path, pattern="env")
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

###########################################
# Redundancy analysis of components: cRDA #
###########################################

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


# Parallel Analysis: run once for 540 & 1800 individuals
# ------------------------------------------------------

 #PA540  <- hornpa(k=540, size=540,   p=0.99, reps=1000, seed=100, method="pca")
 #PA1800 <- hornpa(k=1800, size=1800, p=0.99, reps=1000, seed=100, method="pca")
 
 #write.csv(PA540[,3],  paste0(forester_results_path,"PA_540.csv"),  row.names=F)
 #write.csv(PA1800[,3], paste0(forester_results_path,"PA_1800.csv"), row.names=F)

# Data prep
# ---------

j = as.numeric(seed)

verify <- "Do the files match?" 

i=Sites.90[j]
Site <- rep(1:Design$NumPops[i],each=as.numeric(as.vector(Design$NumInd)[i])) 
tmp <- read.table(paste0(forester_simfiles_path, Filenames.lfmm[[i]])) # the genetic data
z=Sel.90[j]
doubleCheck <- identical(unlist(Design[i,]), unlist(Design.sel[z,]))# sanity check that the files match up
verify <- rbind(verify, doubleCheck)
file1 <- read.table(paste0(forester_results_path, Filenames.sum[[z]]), header=T)
colnames(tmp) <- file1$SNPnames[file1$SNPIncluded==TRUE]
nloci <- ncol(tmp)   # number of loci

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
des <- coord.full[ , grepl( Design$Design[i] , colnames(coord.full) ) ]
coord.full.des <- cbind(coord.full[,2:3],des)
coord.des <- coord.full.des[coord.full.des[,3]==1,]
x <- unsplit(coord.des[,1], Site); y <- unsplit(coord.des[,2], Site)
coord <- cbind(x,y); colnames(coord) <- c("X_Pops","Y_Pops")

Hab <- unname(Env[[i]])
predictors <- scale(Hab, center=T, scale=T)


# cRDA
# ----

# Step 1: PCA
# -----------
snps.scale <- scale(snps, center=T, scale=T)
ncomp <- dim(snps.scale)[1]                   
p <- prcomp(snps.scale, center=F, scale=F) 
p.eigval <- (p$sdev)^2                        # eigenvalues for 540/1800 PCs depending on sample size

# Step 2: Select PCs to retain
# ----------------------------

# Parallel analysis:
PA <- read.csv(paste0(forester_results_path, "PA_", ncomp, ".csv"))
retainPA <- sum(apply(as.data.frame(p.eigval), 2, ">" , PA), na.rm=TRUE)

# Step 3: Varimax rotation
# ------------------------
rawLoadings <- p$rotation[,1:ncomp] %*% diag(p$sdev, ncomp, ncomp) # loadings * eigenvalues
rotatedLoadings <- varimax(rawLoadings)$loadings
invLoadings <- t(pracma::pinv(rotatedLoadings)) #This takes a while to run
scores <- scale(snps.scale) %*% invLoadings   # scores computed via rotated loadings
rcPA <- scores[,1:retainPA]

# Step 4: RDA with rotated components
# -----------------------------------
rda.PA <- rda(rcPA, predictors, scale=F)
#rda.PA2 <- rda(rcPA, predictors, scale=F)
#load.rda.PA <- summary(rda.PA)$species[,1]                  # RDA loadings
#load.rda.PA_df <- as.data.frame(load.PA.rda)

rda.PA.ca <- rda.PA$CCA$u    # constrained axes

# Step 5: Correlation between retained comp & constrained axes
# ------------------------------------------------------------
n <- ncomp                        # sample size of individuals

# Calculate correlation cutoff
a5.PA     <- 0.05/retainPA        # Bourret starts with 0.05
qcrit5.PA <- qnorm(1-a5.PA/2)
r5.PA     <- tanh(qcrit5.PA/sqrt(n-3))  

rda.PA.corr <- apply(rcPA, 2, function(x) cor.test(x, rda.PA.ca[,1]))
signif.PA <- which(lapply(rda.PA.corr, function(x) abs(x$estimate) > r5.PA)==TRUE)

# Step 6: Correlate SNPS with significant retained components
# -----------------------------------------------------------

alpha <- c(0.05, 0.01, 0.001)    # uncorrected cutoffs

out <- list(a0.05=list(snpPA=NA),a0.01=list(snpPA=NA),a0.001=list(snpPA=NA))
loadPA <- matrix(data=NA, nrow=nloci, ncol=length(signif.PA))

for (a in 1:length(alpha)) {
  ac    <- alpha[a]/retainPA
  qcrit <- qnorm(1-ac/2)
  rPA   <- tanh(qcrit/sqrt(n-3))   # correlation cutoff for PA retained components
  
  for (comp in 1:length(signif.PA)) {
    snp.srcPA.corr <- apply(snps.scale, 2, function(x) cor.test(x, rcPA[,signif.PA[comp]]))
    all.snp.srcPA <- unlist(lapply(snp.srcPA.corr, function(x) unname(x$estimate)))
    loadPA[,comp] <- scale(all.snp.srcPA, center=T, scale=T)
    
    # cutoffs
    signif.snp.srcPA <- which(lapply(snp.srcPA.corr, function(x) abs(x$estimate) > rPA)==TRUE) # SHORTER because alpha is more restricted.
    signif.snp.srcPA <- unname(signif.snp.srcPA)
    out[[a]]$snpPA <- c(out[[a]]$snpPA, signif.snp.srcPA)
  }
}  

# -----------------------------
# Power from empirical p-values
# -----------------------------

nloadPA <- dim(loadPA)[2]  # these loadings are correlations between SNPs and retained components...
names(loadPA) <- c(1:nloci)


emp_total_cRDA <- as.data.frame(cbind(id=c(1:nloci), SNPnames=file1$SNPnames[file1$SNPIncluded==TRUE], loadings=loadPA[,1]))
emp_total_cRDA$SNPnames <- file1$SNPnames[file1$SNPIncluded==TRUE]

neu.ep <- abs(loadPA[1:9900])
sel.ep <- abs(loadPA[9901:nloci])

emp.p <- sapply(c(neu.ep,sel.ep), getEmpP, sort(neu.ep))  # p-value from cumulative distribution
emp.p <- 1-emp.p     # p-value for qvalue
emp.p_crda <- as.data.frame(emp.p)

emp_total_cRDA <- cbind(emp_total_cRDA, emp.p_crda)

emp.q <- qvalue(emp.p)$qvalues
emp.q_crda <- as.data.frame(emp.q)
emp_total_cRDA <- cbind(emp_total_cRDA, emp.q_crda)

#paT.empirical <- as.numeric(names(which(emp.q[9901:nloci] < 0.01)))
#paTP.empirical <- sum(emp.q[9901:nloci] < 0.01, na.rm=TRUE)
#paTPR.empirical <- mean(emp.q[9901:nloci] < 0.01, na.rm=TRUE)
emp_total_cRDA$emp.q_cRDA_log <- "NA"

for (i in 1:nloci) {
  if(emp_total_cRDA[i,5] < 0.01) {
    emp_total_cRDA$emp.q_cRDA_log[i] <- "TRUE"
  }
  else {
    emp_total_cRDA$emp.q_cRDA_log[i] <- "FALSE"
  }
}

# Selection info
# --------------
#z=Sel.90[j]
#doubleCheck <- identical(unlist(Design[i,]), unlist(Design.sel[z,]))  # sanity check that the files match up
#tmp <- read.table(paste0(forester_simfiles_path, Filenames.sum[[z]]), header=T) # the summary file
#tmp <- tmp[(tmp$IsNeut=="Sel" & tmp$SNPIncluded==TRUE), 5]
#n.sel <- length(tmp)
#tmp <- as.data.frame(cbind(c(9901:(9900+n.sel)), tmp))
#names(tmp) <- c("SNP","Sel")

#sel001 <- tmp[tmp$Sel==0.001, 1]
#sel005 <- tmp[tmp$Sel==0.005, 1]
#sel01  <- tmp[tmp$Sel==0.01, 1]
#sel1   <- tmp[tmp$Sel==0.1, 1]


# Levels of selection detected from emp pvalue detections
# -------------------------------------------------------
#eTPR001pa <- sum(paT.empirical %in% sel001) / length(sel001)
#eTPR005pa <- sum(paT.empirical %in% sel005) / length(sel005)
#eTPR01pa  <- sum(paT.empirical %in% sel01) / length(sel01)
#eTPR1pa   <- sum(paT.empirical %in% sel1) / length(sel1)


# ------------------------------
# Detection levels using cutoffs
# ------------------------------

# TPR and FPR
# -----------
# NOTE: number of neutral loci is always 9900 in these simulations

#out.NArm <- rapply(out, function(x) x[!is.na(x)], how="list")     # remove NAs
#out.final <- rapply(out.NArm, function(x) unique(x), how="list")  # remove duplicates

#TPR <- rapply(out.final, function(x) sum(x > 9900)/(nloci-9900))
#FPR <- rapply(out.final, function(x) sum(x <= 9900) / 9900)


# Selection levels detected
# -------------------------
#TP <- rapply(out.final, function(x) x[x > 9900], how="list")

#TPR001 <- rapply(TP, function(x) sum(x %in% sel001) / length(sel001))
#TPR005 <- rapply(TP, function(x) sum(x %in% sel005) / length(sel005))
#TPR01  <- rapply(TP, function(x) sum(x %in% sel01) / length(sel01))
#TPR1   <- rapply(TP, function(x) sum(x %in% sel1) / length(sel1))


# -----------------
# Write out results
# -----------------

i=Sites.90[j]
emp_cRDA_final <- cbind(emp_total_cRDA, corresponding_file=Filenames.sum[[z]])

library(tools)
i=Sites.90[j]
simulation <- file_path_sans_ext(Filenames.lfmm[[i]])
write.table(emp_cRDA_final, file=paste0(forester_results_path,"constrainedordination_resultsordination_results/cRDA_locus_stats_", simulation, ".txt"), sep = " ", row.names = F)
write.table(verify, file= paste0(forester_results_path, "constrainedordination_results/cRDA_verifying_match.txt"), sep = " ", row.names =F)
