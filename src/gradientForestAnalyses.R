
# Load libraries ---------------------------------------------------------------
#library(devtools)
#install_github("whitlock/OutFLANK")

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


# FUNCTIONS & SOURCED SCRIPTS --------------------------------------------------

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
# build output GF data frames
gfR2tab <- function(gfMods.list, alFreqs){
  i=1
  while(is.null(gfMods.list[[i]])){i=i+1}
  tab <- do.call(rbind, gfMods.list)
  vrNm <- rep(row.names(tab)[1:nrow(gfMods.list[[i]])], 
              nrow(tab)/nrow(gfMods.list[[i]]))
  tab <- data.frame(variable=vrNm, tab)
  tab <- dcast(tab, SNPnames~variable, value.var="imps")
  envR2 <- rowSums(data.frame(tab[,-1]))
  R2Tab <- data.frame(tab, envR2=envR2)
  
  # get name of SNP if it has a positive R2
  posR2 <- unlist(lapply(gfMods.list, function(x){
    return(as.character(unique(x[,2])))}))
  
  # Find which loci have R2 < 0 (no GF model for those) & assign R2=0
  negR2 <- !(colnames(alFreqs) %in% posR2)
  negR2 <- colnames(alFreqs)[negR2] 
  
  noGF <- data.frame(matrix(0, nrow=length(negR2), ncol=ncol(R2Tab)))
  colnames(noGF) <- colnames(R2Tab)           
  noGF$SNPnames <- negR2
  
  R2Tab <- rbind(R2Tab, noGF)
  return(R2Tab[order(R2Tab$SNPnames),])}
#####
################################################################################

zoot <- TRUE
if(zoot){
  path <- "/pool/home/mfitzpatrick/Projects/activeProjects/TTT_LotterhosWhitlockData/"
} else {
  path <- "/Volumes/localDrobo/Projects/activeProjects/testingTheTests/fitzLab-AL_TTT_LotterhosWhitlockData/"
}

####### START PREP DATA AND RUN GF #############################################
# Load and prep data (env, allelic) ---------------------------------------
# number of cores to use for parallel processing
cores <- 11

# "background" environment
# this file is currently in the old repo
# need to check if it is to be used with new forester sims?
#bgEnv <- read.table("/Volumes/localDrobo/Projects/activeProjects/testingTheTests/fitzLab-AL_TTT_LotterhosWhitlockData/results_AdaptreeEnviFor_R90.txt") 

# sims
# Updated to new forester sim folder 11/28/18
simFiles <- list.files(path=paste0(path, "forester_simfiles"), full.names=T)
simIDs <- unique(sapply(strsplit(sapply(strsplit(simFiles, "_NumPops"), function(x){
  x[1]}), "/forester_simfiles/", fixed=T), function(x){
    x[2]}))

simIDs <- simIDs[-grep(".txt", simIDs)]
simIDs <- simIDs[-grep("README.md", simIDs)]

# for testing lapply function below
#simID <- simIDs[1]

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
  cpValFile <- list.files(path=paste0(path, "forester_results"), pattern=simID, full.names=T)
  cpValFile <- cpValFile[grep(numInds, cpValFile)]
  cpValFile <- cpValFile[grep(".Cpval", cpValFile)]
  if(length(cpValFile)>1){cpValFile <- cpValFile[-grep("gradientforests", cpValFile)]}
  if(length(cpValFile)>1){cpValFile <- cpValFile[-grep("GDM", cpValFile)]}
  cpVal <- read.table(cpValFile, header=T) #should always have 10000 loci
  # select only those included (I think because some go to fixation...)
  cpVal.use <- subset(cpVal, SNPIncluded==TRUE)
  
  # env gradient(s)
  envSelect <- read.table(sim[grep("env", sim)])
  names(envSelect) <- "envSelect"
  
  # allele data
  # columns = loci
  # rows = total # of individuals (#populations x #inds sampled)
  allelic <- fread(sim[grep("lfmm", sim)], header=F, data.table=F)
  #allelic <- allelic[,cpVal.use$SNPIncluded]
  
  # create character name for SNPs
  snpID <- paste0("S", cpVal.use$SNPnames) #paste("S", row.names(cpVal.use), sep="")
  names(allelic) <- snpID
  
  # data stats - used for indexing, etc
  popSize <- nrow(allelic)/nrow(unique(envSelect))#nrow(bgEnv)
  numPops <- nrow(unique(envSelect))
  
  popID <- sort(rep(1:numPops, popSize))
  
  # build data tables for individuals & populations
  datInd <- data.frame(popID=popID, envSelect=envSelect, allelic)
  
  # aggregate to population-level allele freqs
  datPop <- aggregate(. ~ popID+envSelect, data=datInd, FUN=function(x, popSize){
    sum(x)/popSize}, popSize=popSize)
  alFreq <- datPop[,-c(1, 2)] #"popID", "envSelect"

  envPop <- data.frame(envSelect = datPop$envSelect)
  
  ##############################################
  # Chunk to fit GF models to minor allele frequencies at the level of
  # populations
  # GF is fit to each SNP individually to 
  # ease computational / memory burden
  print(paste("Minor allele freq", simID, sep="::"))
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  # use for each to run in parallel - each SNP modeled by GF on different
  # processor
  gfAllele.freq <- foreach(k=1:ncol(alFreq), .verbose=F, .packages=c("gradientForest", "data.table")) %dopar% {
    # get locus k & name with SNP ID
    locus <- data.frame(alFreq[,k])
    names(locus) <- colnames(alFreq)[k]
    
    # here's the GF modeling function
    gfLocus <- gradientForest(data=data.frame(envPop, locus),
                              predictor.vars=colnames(envPop), 
                              response.vars=colnames(alFreq)[k], 
                              corr.threshold=0.5, 
                              ntree=500, 
                              trace=F)
    
    # get result is the model is not = NULL (no fit)
    if(!is.null(gfLocus)){
      # cumulative importance values (to plot turnover functions)
      cImp <- cumimp(gfLocus, "envSelect", type="Species")
      cImp <- data.frame(rbindlist(cImp, idcol="allele"))
      # return data frame with allele ID, cumulative importance results,
      # and R2 value for GF model
      data.frame(cImp, r2=gfLocus$result)
    }
  }
  stopCluster(cl)
  ##############################################
  
  # Prep table of importance values for each SNP (rows) 
  # and each var (columns)
  gfAF <- lapply(gfAllele.freq, function(x){
    if(!is.null(x)){
      ttt <- x[1,c("r2", "allele")]
      return(data.frame(imps=ttt$r2, SNPnames=ttt$allele, row.names = "envSelect",
                        stringsAsFactors=F))}})
  
  # run function to convert to table # then write to file
  gfAllele.R2.af <- gfR2tab(gfAF, alFreq)
  gfAllele.R2.af <- gfAllele.R2.af[match(mixedsort(gfAllele.R2.af$SNPnames), 
                                         gfAllele.R2.af$SNPnames),]
  gfAllele.R2.af$SNPnames <- gsub("S", "", gfAllele.R2.af$SNPnames)
  
  notUsed <- cpVal$SNPnames[which((cpVal$SNPnames %in% gfAllele.R2.af$SNPnames)==FALSE)]
  binder <- gfAllele.R2.af[1:length(notUsed),]
  binder[] <- NA
  binder[,1] <- notUsed
  gfAllele.R2.af <- rbind(binder, gfAllele.R2.af)
  gfAllele.R2.af <- gfAllele.R2.af[match(cpVal$SNPnames, gfAllele.R2.af$SNPnames),]

  ##############################################
  # Chunk to fit GF models to allele presence / absence at the level of
  # individuals
  print(paste("Presence-Absence", simID, sep="::"))
  
  # Allele presence-absence using GF
  # Need to change 0/1 to factor for classification 
  # (otherwise defaults to regression) 
  allelic <- apply(allelic, 2, factor)
  # create dummy variable to avoid error in GF fitting 
  dummyEnv <- runif(length(envSelect[,1]), -0.001, 0.001)
  envInd <- data.frame(envSelect=envSelect, dummyEnv=dummyEnv)
  #envInd <- data.frame(envSelect=envSelect)
  
  # same as for MAF above, but now for allele pa/
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  gfAllele.pa <- foreach(k=1:ncol(allelic), .verbose=F, .packages=c("gradientForest", "data.table")) %dopar%{
    locus <- data.frame(allelic[,k])
    names(locus) <- colnames(allelic)[k]
    
    gfLocus <- gradientForest(data=data.frame(envInd, locus),
                              predictor.vars=colnames(envInd), 
                              response.vars=colnames(allelic)[k], 
                              corr.threshold=0.5, 
                              ntree=500, 
                              trace=F)
    
    if(!is.null(gfLocus)){
      cImp <- cumimp(gfLocus, "envSelect", type="Species")
      cImp <- data.frame(rbindlist(cImp, idcol="allele"))
      data.frame(cImp, r2=gfLocus$result)
    }
  }
  stopCluster(cl)
  ##############################################
  
  # table of importance values for each SNP (rows) and each var (columns)
  # run function to convert to table # then write to file
  gfPA <- lapply(gfAllele.pa, function(x){
    if(!is.null(x)){
      ttt <- x[1,c("r2", "allele")]
      return(data.frame(imps=ttt$r2, SNPnames=ttt$allele, row.names = "envSelect",
                        stringsAsFactors=F))}})
  
  # run function to convert to table # then write to file
  gfAllele.R2.pa <- gfR2tab(gfPA, alFreq)
  gfAllele.R2.pa <- gfAllele.R2.pa[match(mixedsort(gfAllele.R2.pa$SNPnames), 
                                         gfAllele.R2.pa$SNPnames),]
  gfAllele.R2.pa$SNPnames <- gsub("S", "", gfAllele.R2.pa$SNPnames)
  
  gfAllele.R2.pa <- rbind(binder, gfAllele.R2.pa)
  gfAllele.R2.pa <- gfAllele.R2.pa[match(cpVal$SNPnames, gfAllele.R2.pa$SNPnames),]
  
  finalOut <- data.frame(SNPnames=gfAllele.R2.af$SNPnames, 
                         #cImp.envSelect.alleleFreq=gfAllele.R2.af$envSelect,
                         r_squared.alleleFreq=gfAllele.R2.af$envR2,
                         #cImp.envSelect.allelePres_Abs=gfAllele.R2.pa$envSelect,
                         r_squared.allelePres_Abs=gfAllele.R2.pa$envR2)
  
  write.table(finalOut, 
            paste0(getwd(), "/", strsplit(strsplit(sim[1], ".env")[[1]][1], "forester_simfiles/")[[1]][2],
                   "_gradientforests.Cpval"), 
            row.names=F)
  
  save(gfAllele.freq, 
              file=paste0(getwd(), "/", strsplit(strsplit(sim[1], ".env")[[1]][1], "forester_simfiles/")[[1]][2],
                     "_gradientforests_cImp_alleleFreq.Rdata"))
  
  #save(gfAllele.pa, 
  #     file=paste0(getwd(), "/forester_results/", strsplit(strsplit(sim[1], ".env")[[1]][1], "forester_simfiles/")[[1]][2],
  #            "_gradientforests_cImp_allelePres_Abs.Rdata"))
}, numInds=6)
################################################################################

# GF on all alleles simultaneously ---------------------------------------
# sims
# Updated to new forester sim folder 11/28/18
simFiles <- list.files(path=paste0(path, "forester_simfiles"), full.names=T)
simIDs <- unique(sapply(strsplit(sapply(strsplit(simFiles, "_NumPops"), function(x){
  x[1]}), "/forester_simfiles/", fixed=T), function(x){
    x[2]}))

simIDs <- simIDs[-grep(".txt", simIDs)]
simIDs <- simIDs[-grep("README.md", simIDs)]

# for testing lapply function below
#simID <- simIDs[1]

# mclapply through each simulation
mclapply(simIDs, function(simID, numInds=c(20, 6)){
  # x-y and environment
  sim <- simFiles[grep(simID, simFiles)]
  
  # simulations with either 20 or 6 individuals
  # select using numInds argument
  sim <- sim[grep(paste0("NumInd=", numInds), sim)]
  
  # Updated to new forester sim folder 11/28/18
  # read in cpVal table
  cpValFile <- list.files(path=paste0(path, "forester_results"), pattern=simID, full.names=T)
  cpValFile <- cpValFile[grep(paste0("NumInd=", numInds), cpValFile)]
  cpValFile <- cpValFile[grep(".Cpval", cpValFile)]
  if(length(cpValFile)>1){cpValFile <- cpValFile[-grep("gradientforests", cpValFile)]}
  if(length(cpValFile)>1){cpValFile <- cpValFile[-grep("GDM", cpValFile)]}
  cpVal <- read.table(cpValFile, header=T) #should always have 10000 loci
  # select only those included (I think because some go to fixation...)
  cpVal.use <- subset(cpVal, SNPIncluded==TRUE)
  
  # env gradient(s)
  envSelect <- read.table(sim[grep("env", sim)])
  names(envSelect) <- "envSelect"
  
  # allele data
  # columns = loci
  # rows = total # of individuals (#populations x #inds sampled)
  allelic <- fread(sim[grep("lfmm", sim)], header=F, data.table=F)
  #allelic <- allelic[,cpVal.use$SNPIncluded]
  
  # create character name for SNPs
  snpID <- paste0("S", cpVal.use$SNPnames) #paste("S", row.names(cpVal.use), sep="")
  names(allelic) <- snpID
  
  # data stats - used for indexing, etc
  popSize <- nrow(allelic)/nrow(unique(envSelect))#nrow(bgEnv)
  numPops <- nrow(unique(envSelect))
  
  popID <- sort(rep(1:numPops, popSize))
  
  # build data tables for individuals & populations
  datInd <- data.frame(popID=popID, envSelect=envSelect, allelic)
  
  # aggregate to population-level allele freqs
  datPop <- aggregate(. ~ popID+envSelect, data=datInd, FUN=function(x, popSize){
    sum(x)/popSize}, popSize=popSize)
  alFreq <- datPop[,-c(1, 2)] #"popID", "envSelect"
  
  envPop <- data.frame(envSelect = datPop$envSelect)
  
  ##############################################
  # Chunk to fit GF models to minor allele frequencies 
  # all loci
  gfMAF <- gradientForest(data=data.frame(envPop, alFreq),
                              predictor.vars=colnames(envPop), 
                              response.vars=colnames(alFreq), 
                              corr.threshold=0.5, 
                              ntree=500, 
                              trace=T)
  cImpMAF <- cumimp(gfMAF, "envSelect", type="Overall")
  cImpMAF <- data.frame(x=cImpMAF$x, y=cImpMAF$y)
  save(cImpMAF, 
       file=paste0(getwd(), "/", strsplit(strsplit(sim[1], ".env")[[1]][1], "forester_simfiles/")[[1]][2],
                   "_gradientforests_cImp_alleleFreq_overall_All_Loci.Rdata"))
  
  # selected loci
  seld <- which(names(alFreq) %in% paste0("S", cpVal$SNPnames[which(cpVal$s_high>0)]))
  gfMAF.sel <- gradientForest(data=data.frame(envPop, alFreq[,seld]),
                          predictor.vars=colnames(envPop), 
                          response.vars=colnames(alFreq)[seld], 
                          corr.threshold=0.5, 
                          ntree=500, 
                          trace=T)
  cImpMAF.sel <- cumimp(gfMAF.sel, "envSelect", type="Overall")
  cImpMAF.sel <- data.frame(x=cImpMAF.sel$x, y=cImpMAF.sel$y)
  save(cImpMAF.sel, 
       file=paste0(getwd(), "/", strsplit(strsplit(sim[1], ".env")[[1]][1], "forester_simfiles/")[[1]][2],
                   "_gradientforests_cImp_alleleFreq_overallSelected.Rdata"))
  
  # neutral loci
  neut <- which(names(alFreq) %in% paste0("S", cpVal$SNPnames[which(cpVal$s_high==0)]))
  gfMAF.neut <- gradientForest(data=data.frame(envPop, alFreq[,neut]),
                              predictor.vars=colnames(envPop), 
                              response.vars=colnames(alFreq)[neut], 
                              corr.threshold=0.5, 
                              ntree=500, 
                              trace=T)
  cImpMAF.neut <- cumimp(gfMAF.neut, "envSelect", type="Overall")
  cImpMAF.neut <- data.frame(x=cImpMAF.neut$x, y=cImpMAF.neut$y)
  save(cImpMAF.neut, 
       file=paste0(getwd(), "/", strsplit(strsplit(sim[1], ".env")[[1]][1], "forester_simfiles/")[[1]][2],
                   "_gradientforests_cImp_alleleFreq_overallNeutral.Rdata"))
  
  #gfMAF.comb <- combinedGradientForest(selected=gfMAF.neut, neutral=gfMAF.neut, method=2)
  ##############################################
  
}, numInds=6, mc.cores = 18, mc.cleanup = T)
################################################################################


####### PLOT GF TURNOVER / CUMULATIVE IMPORTANCE FUNCTIONS #####################
simFiles <- list.files(path=paste0(path, "forester_simfiles"), full.names=T)

impData <- list.files(pattern=".Rdata",
                      path=getwd(),
                      full.names=T,
                      recursive=T)

for(i in 1:length(impData)){
  load(impData[i])
  plotTitle <- strsplit(strsplit(impData[i], "_grad")[[1]][1], "results/")[[1]][2]
  cat(plotTitle, fill=T)
  
  sim <- simFiles[grep(plotTitle, simFiles)]

  # read in cpVal table
  cpValFile <- list.files(path=paste0(path, "forester_results"), pattern=plotTitle, full.names=T)
  cpValFile <- cpValFile[grep(".Cpval", cpValFile)]
  if(length(cpValFile)>1){cpValFile <- cpValFile[-grep("gradientforests", cpValFile)]}
  if(length(cpValFile)>1){cpValFile <- cpValFile[-grep("GDM", cpValFile)]}
  cpVal <- read.table(cpValFile, header=T) #should always have 10000 loci
  # select only those included (I think because some go to fixation...)
  cpVal.use <- subset(cpVal, SNPIncluded==TRUE)
  
  # env gradient(s)
  envSelect <- read.table(sim[grep("env", sim)])
  names(envSelect) <- "envSelect"
  
  impDatList <- gfAllele.freq[unlist(lapply(gfAllele.freq, function(x){!is.null(x)}))]
  impDat <- do.call(rbind, impDatList)
  
  x <- sort(unique(envSelect)[,1])
  
  ggCand <- impDat
  strSel <- factor(cpVal.use$s_high[match(ggCand$allele, paste0("S", cpVal.use$SNPnames))])
  
  ggCand <- data.frame(ggCand, strSel)
  
  #snpID[cpVal.use$IsNeut=="Sel"]
  #ggCand <- data.frame(ggCand, isNeut=as.character(cpVal.use$IsNeut[match(ggCand$allele, paste("X", cpVal.use$SNPnames, sep=""))]))
  
  maxs <- NULL
  maxs[1] <- max(impDatList[[1]]$y)
  for(j in 2:length(impDatList)){
    maxs[j] <- max(impDatList[[j]]$y)
  }
  
  p.imp <- ggplot() + geom_line(data=ggCand, aes(x=x, y=y, group=allele),
                                colour=rgb(0,0,0,0.4), lwd=0.5) +
    facet_grid(. ~ strSel) +
    labs(y="Cumulative Importance", x="Environment") +
    geom_line(data=cImpMAF, aes(x=x, y=y),
              colour=rgb(1,0,0, 0.75), lwd=1) +
    theme(plot.margin = unit(c(1.25,1.25,1.25,1.25), "in")) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 18, colour = "grey60"),
          axis.title.x = element_text(size=24)) +
    theme(axis.text.y = element_text(size = 16, colour = "grey60"),
          axis.title.y = element_text(size=24, vjust=1)) +
    theme(strip.text = element_text(size=16)) +
    ggtitle(plotTitle) +
    theme(plot.title = element_text(size=14, face="bold.italic"))

  pngName <- paste0(getwd(), "/", plotTitle, "_gradientforests_cImp_alleleFreq2.png")
  
  ggsave(pngName, device="png", width = 16, height = 10, units = "in", dpi=300, p.imp)
}

# ##### plot cImp for PAF models #####
# impDatList <- gfAllele.pa[unlist(lapply(gfAllele.pa, function(x){!is.null(x)}))]
# impDat <- do.call(rbind, impDatList)
# 
# x <- sort(unique(envSelect)[,1])
# 
# ggCand <- impDat #rbind(ttt, nnn)
# strSel <- factor(cpVal.use$s_high[match(ggCand$allele, paste("S", row.names(cpVal.use), sep=""))])
# 
# ggCand <- data.frame(ggCand, strSel)
# 
# snpID[cpVal.use$IsNeut=="Sel"]
# ggCand <- data.frame(ggCand, isNeut=as.character(cpVal.use$IsNeut[match(ggCand$allele, paste("X", cpVal.use$SNPnames, sep=""))]))
# 
# #ggCand$colorSNP <- as.character(ggCand$colorSNP)
# #ggCand$colorSNP[ggCand$allele %in% paste("V", 9900:10000, sep="")] <- "red"
# 
# maxs <- NULL
# maxs[1] <- max(impDatList[[1]]$y)
# for(j in 2:length(impDatList)){
#   maxs[j] <- max(impDatList[[j]]$y)
# }
# 
# p.imp <- ggplot() + geom_line(data=ggCand, aes(x=x, y=y, group=allele),
#                               colour=rgb(0,0,0,0.4), lwd=0.5) +
#   facet_grid(. ~ strSel) +
#   labs(y="Cumulative Importance", x="Environment") +
#   
#   #ylim(0, max(maxs)*1.2) +
#   theme(plot.margin = unit(c(1.25,1.25,1.25,1.25), "in")) +
#   theme_bw() +
#   theme(axis.text.x = element_text(size = 18, colour = "grey60"),
#         axis.title.x = element_text(size=24)) +
#   theme(axis.text.y = element_text(size = 16, colour = "grey60"),
#         axis.title.y = element_text(size=24, vjust=1)) +
#   theme(strip.text = element_text(size=16))
# 
# ggsave(paste(getwd(), "/gradientForestResults/facetPA_cImp_", simID, ".png", sep=""),
#        device="png", width = 16, height = 10, units = "in", dpi=300, p.imp)
########################### END PLOTTING #######################################


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


sapply(simIDs, gfAUC)

quants <- quantile(modResults$envSelect, probs=probs)

fnr <- sapply(quants, function(x){
  inds <- 9900:length(modResults$envSelect)
  negs <- sum(modResults$envSelect[inds]<x)
  return(negs/length(inds))
})

fpr <- sapply(quants, function(x){
  inds <- 1:9899
  fpos <- sum(modResults$envSelect[inds]>=x)
  return(fpos/length(inds))
})

plot(fpr, 1-fnr, type="l", col=)
  
  1-(test SNP rank/[number intergenic SNPs+1])

# loop through each simulation and plot gf and gdm models
for(k in 1:length(sims)){
  modResults <- read.csv(list.files(pattern=simIDs[1]))
  
  # alpha to calculate false postives / false negatives
  alpha <- seq(0, 0.1, 0.001)
  quants <- quantile(modResults$envSelect, probs=probs)
  
  # false neg/pos rates at alpha = 0.01
  fnrPlot <- sum(modResults$envSelect[9900:length(modResults$envSelect)] < quantile(gfAlfreq$R2, probs=0.99))/sum(gfAlfreq$ind>9900)
  fprPlot <- sum(gfAlfreq$R2[gfAlfreq$ind<=9900] >= quantile(gfAlfreq$R2, probs=0.99))/sum(gfAlfreq$ind<=9900)
  
  fnr <- NULL
  fpr <- NULL
  for(w in 1:length(quants)){
    fnr[w] <- sum(gfAlfreq$R2[gfAlfreq$ind>9900] < quants[w])/sum(gfAlfreq$ind>9900)
    fpr[w] <- sum(gfAlfreq$R2[gfAlfreq$ind<=9900] >= quants[w])/sum(gfAlfreq$ind<=9900)
  }
  
  errorRates[[k]][[2]] <- data.frame(fnr, fpr)
  
  # colors for plotting
  bg <- (gfAlfreq$R2 >= quantile(gfAlfreq$R2, probs=0.99))+1
  bg <- ifelse(bg==1, rgb(0,0,0,0.15), rgb(1,0,0,0.5))
  bg[gfAlfreq$ind>9900] <- rgb(0,0,1,0.5)
  
  # plot & save 
  filename <- paste("gf.allele.freq_", sims[k], ".tif", sep="")
  tiff(file=filename, width=18, height=12, units="in", compression="lzw", res=150,
       type="cairo")
  
  par(mai=c(1.25,1.25,1.25,1.25))
  
  col <- rgb(0,0,0,0.5)
  
  plot(gfAlfreq$R2, bg=bg, pch=21, xlab="locus", ylab=expression(R^2), cex.lab=3.5, 
       cex.axis=2, col=col, ylim=c(0, 1), main="GF on minor allele frequencies", 
       col.axis="grey50", cex.main=2)
  
  text(5000, 1, sims[k], cex=1.5)
  text(2000, 0.8, paste("false negative rate =", round(fnrPlot, 3), sep=" "), cex=3) 
  text(2000, 0.74, paste("false positive rate =", round(fprPlot, 4), sep=" "), cex=3)  
  dev.off()
  ####################