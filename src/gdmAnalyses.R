
# Load libraries ---------------------------------------------------------------
library(adegenet)
library(gdm)
library(gradientForest)
library(foreach)
library(doParallel)
library(OutFLANK)
library(pbapply)
library(gdata)
library(data.table)
library(PresenceAbsence)
library(ROCR)
library(modEvA)
library(ggplot2)
library(grid)
library(gridExtra)
library(gtools)
################################################################################
options(scipen=999)

# FUNCTIONS & SOURCED SCRIPTS --------------------------------------------------
gfOutObj <- setClass("gfOutObj", slots = c(alFreq="data.frame", imp="list"))

# Helper functions for Generalized Diss. Models
# Work in progress ignore for now
source("/Users/mfitzpatrick/code/plantGenome/FstByRowtoGDMmatrices.R")


#####
# function to calculate pop-level allele counts / frequencies
# alleleTab = input allele data
# popSize = number of individuals sampled
# numPops = number of locations sampled = nrow(geo)
alleleDat <- function(alleleTab, popSize, numPops){
  mat <- matrix(NA, nrow=numPops, ncol=ncol(alleleTab))
  ind <- data.frame(i1=seq(1, popSize*numPops, by=popSize),
                    i2=seq(0, popSize*numPops, by=popSize)[-1])
  for(i in 1:numPops){
    mat[i,] <- apply(alleleTab[ind[i,1]:ind[i,2],], 2, sum)
  }
  colnames(mat) <- colnames(alleleTab)
  return(list(ifelse(1-(mat/popSize)<0.5, 1-(mat/popSize), mat/popSize), mat))}
#####


#####
# function to create a matrix for each locus of the counts of each 
# allele (columns) in each population (rows)
buildMats <- function(counts, popSize){
  bMat <- cbind(counts, popSize-counts)
  return(bMat)}
#####


#####
# function to calculate population pairwise Fst for a single locus
pwFst <- function(tab, numPops, ind1, ind2){
  mat <- matrix(0, numPops, numPops)
  for(i in 1:length(ind1)){
    mat[ind1[i], ind2[i]] <- WC_FST_FiniteSample_Haploids_2AllelesB_MCW(tab[c(ind1[i], ind2[i]),])[3]
    mat[ind1[i], ind2[i]] <- ifelse(mat[ind1[i], ind2[i]]<0, 0, mat[ind1[i], ind2[i]])
  }
  return(upperTriangle(mat, diag=F))}
#####


##### GDM function - ignore for now
# Function to add genetic distance to site-pair, remove NAs, and scale if desired.
# Scaling can improve model fitting in some instances by increasing the range 
# of Fst values.
finalPrepFst <-  function(x, sitePair, scale=F){
  fst <- sitePair
  fst$distance <- x
  fst <- na.omit(fst)
  if(scale==T){
    return(scaleDist(fst)) # function sourced from above
  } else {return(fst)}}
#####
################################################################################


####### START PREP DATA AND RUN GDM #############################################
# Load and prep data (env, allelic) ---------------------------------------
# number of cores to use for parallel processing
cores <- 12

# sims
# Updated to new forester sim folder 11/28/18
simFiles <- list.files(path="/Volumes/localDrobo/Projects/activeProjects/testingTheTests/fitzLab-AL_TTT_LotterhosWhitlockData/forester_simfiles",
                       full.names=T)
simIDs <- unique(sapply(strsplit(sapply(strsplit(simFiles, "_NumPops"), function(x){
  x[1]}), "/forester_simfiles/", fixed=T), function(x){
    x[2]}))

simIDs <- simIDs[-grep(".txt", simIDs)]
simIDs <- simIDs[-grep("README.md", simIDs)]

# for testing lapply function below
# simID <- simIDs[3]
# numInds <- 20

# lapply through each simulation
# Steps are:
# 1. Prep data for GF
# 2. Fit GF models to each simulation
# 3/ Write results to file
lapply(simIDs, function(simID, numInds=c(20, 6)){
  # x-y and environment
  sim <- simFiles[grep(simID, simFiles)]
  
  # simulations with either 20 or 6 individuals
  # select using numInds argument
  sim <- sim[grep(paste0("NumInd=", numInds), sim)]
  
  # Updated to new forester sim folder 11/28/18
  # read in cpVal table
  cpValFile <- list.files(path="/Volumes/localDrobo/Projects/activeProjects/testingTheTests/fitzLab-AL_TTT_LotterhosWhitlockData/forester_results", 
                          pattern=simID, full.names=T)
  cpValFile <- cpValFile[grep(paste0("NumInd=", numInds), cpValFile)]
  #cpValFile <- cpValFile[grep(".Cpval", cpValFile)]
  if(length(cpValFile)>1){cpValFile <- cpValFile[-grep("gradientforests", cpValFile)]}
  cpVal <- read.table(cpValFile, header=T) #should always have 10000 loci
  # select only those included (I think becuase some go to fixation...)
  cpVal.use <- subset(cpVal, SNPIncluded==TRUE)
  
  # env gradient(s)
  envSelect <- read.table(sim[grep("env", sim)])
  names(envSelect) <- "envSelect"
  
  geo <- read.table("/Volumes/localDrobo/Projects/activeProjects/testingTheTests/fitzLab-AL_TTT_LotterhosWhitlockData/forester_simfiles/SchemeRandom1.txt")
  geo <- geo[which(geo$R90==TRUE),]
  
  # allele data
  # columns = loci
  # rows = total # of individuals (#populations x #inds sampled)
  allelic <- fread(sim[grep("lfmm", sim)], header=F, data.table=F)
  
  # create character name for SNPs
  snpID <- paste0("S", cpVal.use$SNPnames) #paste("S", row.names(cpVal.use), sep="")
  names(allelic) <- snpID
  
  # data stats - used for indexing, etc
  popSize <- nrow(allelic)/nrow(geo)
  numPops <- nrow(geo)
  
  # build env table w/ x-y coords
  xC <- geo$X_Pops[sort(rep(1:numPops, popSize))]
  yC <- geo$Y_Pops[sort(rep(1:numPops, popSize))]
  env.ind <- data.frame(x=xC, y=yC, envSelect=envSelect$envSelect)
  envPop <- unique(env.ind)
  envPop <- data.frame(popID=geo$PopID, envPop)
  
  # Minor allele frequencies
  alFreq.x <- alleleDat(allelic, popSize, numPops)
  alFreq <- alFreq.x[[1]]
  
  # Allele counts
  alCount <- alFreq.x[[2]]
  rm(alFreq.x)
  
  # returns a list n loci long, each element is a 2-column matrix with 
  # major & minor allele counts
  locusMats <- lapply(1:ncol(alCount), function(x, tab){
    buildMats(tab[,x], popSize)
  }, tab=alCount)
  
  # Now, calculate pairwise Fst values between all populations
  # Inefficient to do all possible pairs, so reduce to only unique
  allPairs <- expand.grid(1:numPops, 1:numPops)
  for(i in 1:nrow(allPairs)){
    if(allPairs[i,1] > allPairs[i,2]){allPairs[i,] <- allPairs[i,c(2,1)]}
  }
  
  # find unique combinations of population indices
  uniqPairs <- unique(allPairs)
  ind1 <- uniqPairs[,1]
  ind2 <- uniqPairs[,2]
  
  # calculate population pairwise Fst for all loci
  # calculate population pairwise Fst for all loci
  Fst <- mclapply(locusMats, function(tab){pwFst(tab, numPops, ind1, ind2)},
                  mc.cores=cores, mc.cleanup = T)
  
  # build site-pair tables needed for GDM
  # easier to start with fake data then add Fst values
  fstMat <- matrix(0, numPops, numPops)
  fstMat <- data.frame(popID=envPop$popID, fstMat)
  sitePair <- formatsitepair(bioData=fstMat, bioFormat=3, siteColumn="popID", 
                             predData=envPop, XColumn="x", 
                             YColumn="y")
  
  # loop to add genetic distance to site-pair, remove NAs, and scale if desired.
  # Scaling can improve model fitting in some instances by increasing the range 
  # of Fst values.
  finalPrepFst <-  function(x, sitePair, scale=F){
    fst <- sitePair 
    fst$distance <- x
    fst <- na.omit(fst)
    if(scale==T){
      return(scaleDist(fst)) # function sourced from above
    } else {return(fst)}
  }
  
  # create complete site-pair tables for all Fst matrices
  fstGDM <- pblapply(Fst, finalPrepFst, sitePair, scale=T)
  
  # Fit GDM to each locus
  gdmLoci.geoFalse <- mclapply(fstGDM, gdm, geo=F, mc.cores=cores, mc.cleanup = T)
  gdmLoci.geoTrue <- mclapply(fstGDM, gdm, geo=T, mc.cores=cores, mc.cleanup = T)
  
  #modTest <- gdm.varImp(fstGDM[[4]], geo=F, nPerm=50, parallel=F, cores=10)
  
  devExp.geoFalse <- do.call(rbind, pblapply(gdmLoci.geoFalse, function(x){
    if(!is.null(x)){dev <- x$explained} else{ dev <- NA}
    return(dev)}))
  devExp.geoFalse[devExp.geoFalse<0] <- NA
  
  devExp.geoTrue <- do.call(rbind, pblapply(gdmLoci.geoTrue, function(x){
    if(!is.null(x)){dev <- x$explained} else{ dev <- NA}
    return(dev)}))
  devExp.geoTrue[devExp.geoTrue<0] <- NA
  
  
  devExplained <- data.frame(SNPnames=snpID, 
                             sitepairs=sapply(fstGDM, function(x){nrow(x)}),
                             deviance.geoF=devExp.geoFalse,
                             deviance.geoT=devExp.geoTrue)
  
  devExplained$SNPnames <- gsub("S", "", devExplained$SNPnames)
  
  notUsed <- cpVal$SNPnames[which((cpVal$SNPnames %in% devExplained$SNPnames)==FALSE)]
  binder <- devExplained[1:length(notUsed),]
  binder[] <- NA
  binder[,1] <- notUsed
  devExplained <- rbind(binder, devExplained)
  devExplained <- devExplained[match(cpVal$SNPnames, devExplained$SNPnames),]
  
  
  #   # Fit GDM to each locus
  #   gdmLoci <- mclapply(fstGDM, function(x){
  #     gdmLocus <- gdm(x, geo=T)
  #     if(!is.null(gdmLocus)){
  #       iSpline <- isplineExtract(gdmLocus)
  #       #maxHeight <- max(iSpline$y)
  #       data.frame(deviance=gdmLocus$explained)
  #     }}, mc.cores=cores, mc.cleanup = T)
  #   
  #   # get result is the model is not = NULL (no fit)
  #   if(!is.null(gfLocus)){
  #     # cumulative importance values (to plot turnover functions)
  #     cImp <- cumimp(gfLocus, "envSelect", type="Species")
  #     cImp <- data.frame(rbindlist(cImp, idcol="allele"))
  #     # return data frame with allele ID, cumulative importance results,
  #     # and R2 value for GF model
  #     data.frame(cImp, r2=gfLocus$result)
  #   }
  # }
  
  
  write.table(devExplained, 
              paste0(getwd(), "/", strsplit(strsplit(sim[1], ".env")[[1]][1], "forester_simfiles/")[[1]][2],
                     "_GDM.Cpval"), row.names=F)
  
  # save(gfAllele.freq, 
  #             file=paste0(getwd(), "/", strsplit(strsplit(sim[1], ".env")[[1]][1], "forester_simfiles/")[[1]][2],
  #                    "_gradientforests_cImp_alleleFreq.Rdata"))
}, numInds=6)
  
  
####### END PREP DATA AND RUN GF ###############################################


  # Prep table of importance values for each SNP (rows) 
  # and each var (columns)
  gfAF <- lapply(gfAllele.freq, function(x){
    if(!is.null(x)){
      ttt <- x[1,c("r2", "allele")]
      return(data.frame(imps=ttt$r2, SNPnames=ttt$allele, row.names = "envSelect",
                        stringsAsFactors=F))}})
  
  
########## START MODEL EVALUATION ##############################################
# sims
simFiles <- list.files(path=paste(getwd(), "/simfiles", sep=""), full.names=T)
simIDs <- unique(sapply(strsplit(sapply(strsplit(simFiles, "_NumPops"), function(x){
  x[1]}), "/simfiles/", fixed=T), function(x){
    x[2]}))
simID <- simIDs[1]

gfAUC <- function(simID){
  
  gfResults <- list.files(path=paste(getwd(), "/gradientForestResults", sep=""), 
                          pattern=paste(simID, "csv", sep="."), full.names = T)

  gfResult.freq <- read.csv(gfResults[grep("Freq", gfResults)])
  gfResult.pa <- read.csv(gfResults[grep("PresAbs", gfResults)])
  
  cpVal <- read.table(list.files(path=paste(getwd(), "/results", sep=""), 
                                 pattern=simID, full.names=T), header=T)
  cpVal <- subset(cpVal, UseSNP==TRUE)
  
  Observed.freq <- c(rep(0, max(which(cpVal$s_high==0))), rep(1, (nrow(gfResult.freq)-max(which(cpVal$s_high==0)))))
  Predicted.freq <- gfResult.freq$envR2
  
  Observed.pa <- c(rep(0, max(which(cpVal$s_high==0))), rep(1, (nrow(gfResult.pa)-max(which(cpVal$s_high==0)))))
  Predicted.pa <- gfResult.pa$envR2
  
  # Run the AUC calculations
  ROC_perf.pa <- performance(prediction(Predicted.pa, Observed.pa),"tpr","fpr")
  ROC_sens.pa <- performance(prediction(Predicted.pa, Observed.pa),"sens","spec")
  ROC_auc.pa <- performance(prediction(Predicted.pa, Observed.pa),"auc")
  
  ROC_perf.freq <- performance(prediction(Predicted.freq, Observed.freq),"tpr","fpr")
  ROC_sens.freq <- performance(prediction(Predicted.freq, Observed.freq),"sens","spec")
  ROC_auc.freq <- performance(prediction(Predicted.freq, Observed.freq),"auc")
  
  # Make plot data
  plotDat.pa <- data.frame(FP=ROC_perf.pa@x.values[[1]], TP=ROC_perf.pa@y.values[[1]],
                        CUT=ROC_perf.pa@alpha.values[[1]], POINT=NA, 
                        AUC=ROC_auc.pa@y.values[[1]],
                        Sens=ROC_sens.pa@y.values[[1]],
                        Spec=ROC_sens.pa@x.values[[1]],
                        type="Allele pres-abs")
  plotDat.pa[unlist(lapply(c(1, 0.99, 0.95, 0.90),function(x){which.min(abs(plotDat.pa$CUT-x))})),"POINT"] <- c(1, 0.99, 0.95, 0.90)
  
  
  plotDat.freq <- data.frame(FP=ROC_perf.freq@x.values[[1]], TP=ROC_perf.freq@y.values[[1]],
                           CUT=ROC_perf.freq@alpha.values[[1]], POINT=NA, 
                           AUC=ROC_auc.freq@y.values[[1]],
                           Sens=ROC_sens.freq@y.values[[1]],
                           Spec=ROC_sens.freq@x.values[[1]],
                           type="Minor allele freq")
  plotDat.freq[unlist(lapply(c(1, 0.99, 0.95, 0.90),function(x){which.min(abs(plotDat.freq$CUT-x))})),"POINT"] <- c(1, 0.99, 0.95, 0.90)
  
  plotDat <- rbind(plotDat.pa, plotDat.freq)
  
  # Plot the curve
  ggplot(plotDat, aes(x=FP,y=TP,col=TP)) + 
    scale_colour_gradientn("",colours=rainbow(14)[1:11]) + facet_grid(~type) +
    geom_abline(intercept=0, slope=1) +
    geom_line(lwd=1) + 
    #geom_point(data=plotDat[!is.na(plotDat$POINT),], aes(x=FP,y=TP,fill=POINT), pch=21, size=3, col="black") +
    #geom_text(data=plotDat[!is.na(plotDat$POINT),], aes(x=FP,y=TP,fill=POINT), label=plotDat$POINT[!is.na(plotDat$POINT)], hjust=1, vjust=0, col="black") +
    scale_x_continuous("False Positive Rate", limits=c(0,1)) +
    scale_y_continuous("True Positive Rate", limits=c(0,1)) +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 18),
          legend.text = element_text(size = 15),
          legend.title = element_text(size=16),
          strip.text = element_text(size=20)) +
    scale_fill_gradientn("Threhsold Cutoff",colours=rainbow(14)[1:11]) +
    geom_text(data=data.frame(x=c(0.8, 0.8), y=c(0.05, 0.05), 
                              label=paste("AUC=", round(unique(plotDat$AUC), 3), sep=""), 
                                   type=c("Allele pres-abs", "Minor allele freq")), 
                   aes(x, y, label=label), inherit.aes=FALSE, size=8) +
    ggtitle(simID) + 
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste(getwd(), "/gradientForestResults/AUC_", simID, ".pdf", sep=""), 
         width = 16, height = 10, units = "in", dpi=300)
}
########## END MODEL EVALUATION ####################################################





