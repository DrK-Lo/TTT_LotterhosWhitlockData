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
root_path <- "/Users/brettford/Desktop/Northeastern/coding/forester_simulation_code/forester_sim_code/TTT_LotterhosWhitlockData/"
forester_simfiles_path <- paste0(root_path, "forester_simfiles/")
forester_results_path <- paste0(root_path, "forester_results/")

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
print(verify)
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
load.rda <- summary(snp.rda)$species[,1:5]                  # RDA loadings
load.rda_df <- as.data.frame(load.rda)

# dbRDA
snps.bray <- vegdist(snps, method="bray")                 # bray-curtis for dbRDA
snp.dbrda <- capscale(snps.bray ~ predictors, comm=snps)  
load.dbrda <- summary(snp.dbrda)$species[,1:5]              # dbRDA loadings
load.dbrda_df <- as.data.frame(load.dbrda)

# Create dataframes for rda and dbrda with SNPnames
rda_df <- cbind(id=c(1:nloci), SNPnames=file1$SNPnames[file1$SNPIncluded==TRUE], load.rda_df)
dbrda_df <- cbind(id=c(1:nloci), SNPnames=file1$SNPnames[file1$SNPIncluded==TRUE], load.dbrda_df)

# ------------------
# Write out results
# ------------------

rda_final <- cbind(rda_df, dbrda_df[, 3:7], corresponding_file=Filenames.sum[[z]])
i=Sites.90[j]
simulation <- file_path_sans_ext(Filenames.lfmm[[i]])
write.table(rda_final, file= paste0(forester_results_path,"ordination_results/", "ORD_locus_stats_", simulation, ".txt"), sep = " ", row.names = F)

